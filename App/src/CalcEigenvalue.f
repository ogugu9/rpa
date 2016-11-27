c###################################################################
c##### calcEigen: calc Eigen vector and value in all k-points
c###################################################################

      subroutine calcEigenvalue(nx, ny, ispin)

      use common, only : Nkx, Nky, Nqx, Nqy, Eall, Zpsiall, Pi
      implicit none

      integer, intent(in) :: nx, ny, ispin
      integer :: kx, ky

      real*8 :: akx, aky

      do kx = 0, nx-1 ; do ky = 0, ny-1
         ! ** akxとakyはそのままで正しいか？ → 1/2にしてkx=-nx:nxとするべきか？
         akx = DBLE(kx) / DBLE(Nkx)
         aky = DBLE(ky) / DBLE(Nky)
         call  diagHamiltonian(Eall(kx,ky,:,ispin),Zpsiall(kx,ky,:,:,ispin)
     &        ,akx,aky,ispin)
      end do ; end do

c     ** Not Para Magnetic Mode
      if ((nx /= Nkx).or.(ny /= Nky).or.(Nqx /= 1)) then
c         call expandEigen(Nkx,Nky,nx,ny,
c     &        Eall(:,:,:,ispin),Zpsiall(:,:,:,:,ispin))
      end if

      call boundary()

      return
      end

c###################################################################

c###################################################################
c##### expandEigenk: expand calcEigenk for not PM order
c###################################################################

      subroutine expandEigenvalue(mx,my,mmx,mmy,eig,zpsi)

      use common, only : Nqx, Nqy, Nband
      implicit none

      integer, intent(in) :: mx, my, mmx, mmy
      real*8, intent(inout) :: eig(0:mx,0:my,Nband*Nqx)
      complex*16, intent(inout) :: zpsi(0:mx,0:my,Nband*Nqx,Nband*Nqx)

      integer :: kkx, kky, kx, ky
      integer :: nn0, nn1, iQ, iqq

      do iQ = 1, Nqx-1
         do kx = 0, mmx - 1 ; do ky = 0, mmy - 1
            kkx = MOD(kx + iQ * mx / Nqx, mx)
            kky = MOD(ky + iQ * my / Nqy, my)

            do iqq = 0, Nqx-1
               nn0 = Nband * iqq
               nn1 = Nband * MOD(iqq+iQ,Nqx)
               eig(       kkx, kky, :)
     &              = eig( kx,  ky, :)

               zpsi(    kkx, kky, 1+nn0:Nband+nn0, :)=
     &              zpsi(kx,  ky, 1+nn1:Nband+nn1, :)
            end do
         end do ; end do
      end do

      return
      end

c###################################################################

c###################################################################
c##### boundary:
c###################################################################

      subroutine boundary()


      use common, only : Nkx, Nky, Eall, Zpsiall
      implicit none

      Eall(  Nkx, 0:Nky, :,:) = Eall(0    , 0:Nky, :,:)
      Eall(0:Nkx,   Nky, :,:) = Eall(0:Nkx, 0    , :,:)
      Eall(Nkx,Nky,:,:) = Eall(0,0,:,:)

      Zpsiall(  Nkx, 0:Nky,:,:,:) = Zpsiall(0    , 0:Nky,:,:,:)
      Zpsiall(0:Nkx,   Nky,:,:,:) = Zpsiall(0:Nkx, 0    ,:,:,:)
      Zpsiall(  Nkx,   Nky,:,:,:) = Zpsiall(0    , 0    ,:,:,:)

      return
      end

c###################################################################

c###################################################################
c##### calcEigenBand: calc Eigen vector and value for bandplot
c###################################################################

      subroutine calcBandplot(nx, ny, ispin)

      use common, only : Nk, Nky, Nqx, Nqy, Eband, Zpsiband, Pi
     &      Nkpath, Nkmesh
      implicit none

      integer, intent(in) :: nx, ny, ispin
      integer :: nkpt

      real*8 :: dkx, dky

      do nkpt = 0, Nkpath*Nkmesh
         dkx = klist(nkpt)%kx
         dky = klist(nkpt)%ky
         call  diagHamiltonian(Eband(nkpt,:,ispin),Zpsiband(nkpt,:,:,ispin)
     &        ,dkx,dky,ispin)
      end do ; end do

c     ** Not Para Magnetic Mode
      if ((nx /= Nkx).or.(ny /= Nky).or.(Nqx /= 1)) then
c         call expandEigen(Nkx,Nky,nx,ny,
c     &        Eall(:,:,:,ispin),Zpsiall(:,:,:,:,ispin))
      end if

c      call boundary()

      return
      end

c###################################################################
