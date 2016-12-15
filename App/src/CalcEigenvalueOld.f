c###################################################################
c##### calcEigen: calc Eigen vector and value in all k-points
c###################################################################

      subroutine calcEigenvalue(nx, ny, nz, ispin)

      use common, only : Nkx, Nky, Nkz, Nband, Nqx, Nqy, Nqz,
     &      Eall, Zpsiall, Pi
      implicit none

      integer, intent(in) :: nx, ny, nz, ispin
      integer :: kx, ky, kz

      real*8 :: dkx, dky, dkz

      do kx = 0, nx-1 ; do ky = 0, ny-1 ; do kz = 0, nz-1
         ! ** akxとakyはそのままで正しい?→1/2にしてkx=-nx:nxとするべき?
         dkx = DBLE(kx) / DBLE(Nkx)
         dky = DBLE(ky) / DBLE(Nky)
         dkz = DBLE(kz) / DBLE(Nkz)
         call  diagH(Eall(kx,ky,kz,:,ispin),Zpsiall(kx,ky,kz,:,:,ispin),
     &   dkx,dky,dkz,ispin)
         !** debug write(*,'(5F8.3)'), Eall(kx,ky,:,ispin)
      end do ; end do ; end do

c     ** Not Para Magnetic Mode
      if ((nx /= Nkx).or.(ny /= Nky).or.(nz /= Nkz).or.(Nqx /= 1)) then
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

      subroutine expandEigenvalue(mx,my,mz,mmx,mmy,mmz,eig,zpsi)

      use common, only : Nqx, Nqy, Nqz, Nband
      implicit none

      integer, intent(in) :: mx, my, mz, mmx, mmy, mmz
      real*8, intent(inout) :: eig(0:mx,0:my,0:mz,Nband*Nqx)
      complex*16, intent(inout) ::
     &      zpsi(0:mx,0:my,0:mz,Nband*Nqx,Nband*Nqx)

      integer :: kkx, kky, kkz, kx, ky, kz
      integer :: nn0, nn1, iQ, iqq

      do iQ = 1, Nqx-1
         do kx = 0, mmx - 1 ; do ky = 0, mmy - 1 ; do kz = 0, mmz - 1
            kkx = MOD(kx + iQ * mx / Nqx, mx)
            kky = MOD(ky + iQ * my / Nqy, my)
            kkz = MOD(kz + iQ * mz / Nqz, mz)

            do iqq = 0, Nqx-1
               nn0 = Nband * iqq
               nn1 = Nband * MOD(iqq+iQ,Nqx)
               eig(       kkx, kky, kkz, :)
     &              = eig( kx,  ky,  kz, :)

               zpsi(    kkx, kky, kkz, 1+nn0:Nband+nn0, :)=
     &              zpsi(kx,  ky,  kz, 1+nn1:Nband+nn1, :)
            end do
         end do ; end do ; end do
      end do

      return
      end

c###################################################################

c###################################################################
c##### boundary: ** ブリルアンゾーンをどう考慮するか
c###################################################################

      subroutine boundary()


      use common, only : Nkx, Nky, Nkz, Eall, Zpsiall
      implicit none

      Eall(Nkx, 0:Nky, 0:Nkz, :,:) = Eall(0, 0:Nky, 0:Nkz, :,:)
      Eall(0:Nkx, Nky, 0:Nkz, :,:) = Eall(0:Nkx, 0, 0:Nkz, :,:)
      Eall(0:Nkx, 0:Nky, Nkz, :,:) = Eall(0:Nkx, 0:Nky, 0, :,:)
      Eall(Nkx,Nky,Nkz,:,:) = Eall(0,0,0,:,:)

      Zpsiall(Nkx, 0:Nky, 0:Nkz,:,:,:) = Zpsiall(0, 0:Nky, 0:Nkz,:,:,:)
      Zpsiall(0:Nkx, Nky, 0:Nkz,:,:,:) = Zpsiall(0:Nkx, 0, 0:Nkz,:,:,:)
      Zpsiall(0:Nkx, 0:Nky, Nkz,:,:,:) = Zpsiall(0:Nkx, 0:Nky, 0,:,:,:)
      Zpsiall(Nkx, Nky, Nkz,:,:,:) = Zpsiall(0, 0, 0,:,:,:)

      return
      end

c###################################################################

c###################################################################
c##### calcEigenBand: calc Eigen vector and value for bandplot
c###################################################################

      subroutine calcBandplot(nx, ny, nz, ispin)

      use common, only : Nkx, Nky, Nkz, Nqx, Nqy, Nqz,
     &      Eband, Zpsiband, Pi, Nkpath, Nkmesh, dkfrac
      implicit none

      integer, intent(in) :: nx, ny, nz, ispin
      integer :: nkpt

      real*8 :: dkx, dky, dkz

      do nkpt = 0, Nkpath*Nkmesh
         dkx = dkfrac(nkpt,1)
         dky = dkfrac(nkpt,2)
         dkz = dkfrac(nkpt,3)
         call  diagH(Eband(nkpt,:,ispin),Zpsiband(nkpt,:,:,ispin),
     &   dkx,dky,dkz,ispin)
         !write(*,'(8F8.3)'), dkx, dky, dkz, Eband(nkpt,:,ispin)
      end do

c     ** Not Para Magnetic Mode
      if ((nx /= Nkx).or.(ny /= Nky).or.(ny /= Nkz).or.(Nqx /= 1)) then
c         call expandEigen(Nkx,Nky,nx,ny,
c     &        Eall(:,:,:,ispin),Zpsiall(:,:,:,:,ispin))
      end if

c      call boundary()

      return
      end

c###################################################################
