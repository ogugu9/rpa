c################ 2011 Kaneshita ################

c#########################################################
c#### hamiltoniank:
c#########################################################

      subroutine diagHamiltonian(eig,zpsi,dkx,dky,ispin)

      use common, only : Nkx, Nky, Nband, Nqx, Nqy, Check
      implicit none

      integer, intent(in) :: ispin
      real*8, intent(in) :: dkx, dky
      real*8, intent(out) :: eig(Nband*Nqx)
      complex*16, intent(out) :: zpsi(Nband*Nqx,Nband*Nqx)

      integer :: i, j, nm
      integer :: info, lwork
      real*8 :: ene(Nband*Nqx)
      real*8 :: rwork(3*Nband*Nqx-2)
      complex*16 :: zdummy
      complex*16 :: work((Nband*Nqx+1)*Nband*Nqx)

      complex*16 :: zhmat(Nband*Nqx,Nband*Nqx)

      !** ispin = ispnとしなかったことで intent(in)であるはずのispinがmakeH_int
      !の中で変更されてしまっているためにコンパイルエラーが起きる可能性がある。
      if ((ispin /= 1).and.(ispin /=2)) then
         write(*,*) 'Error in hamiltonian  -- invalid spin'
         stop
      end if

      nm = Nband * Nqx
      lwork = (nm + 1) * nm

      call makeH(zhmat,dkx,dky,ispin)

c## HERMITIAN CHECK[

      if (Check == 'y') then
         call checkHermitian(nm,zhmat,info)
         if(info < 0) stop
         if (info == 1) write(*,*)'Hamiltonian is real',dkx,dky
      end if

c## HERMITIAN CHECK]


c## DIAGONALIZATION[

      Info = 0
      call zheev('V','U',nm,zhmat(1:nm,1:nm),nm,ene(1:nm),
     &     work,lwork,rwork,info)
      if(Info /= 0) then
         write(*,*) 'Lapack ZHEEV: Info=',Info
         stop
      end if

      eig(1:nm) = ene(1:nm)
      zpsi(1:nm,1:nm) = zhmat(1:nm,1:nm)

c## DIAGONALIZATION]

      return
      End

c################################################################


c#########################################################
c#### makeH:
c#########################################################

      subroutine makeH(zhmat,dkx,dky,ispin)

      use common, only : Nkx, Nky, Nband, Nsite, Nqx, Nqy,
     &     Pi, Zi, t, Eorbit, Isitex, Isitey
      implicit none

      integer, intent(in), optional :: ispin
      real*8, intent(in) :: dkx, dky
      complex*16, intent(out) :: zhmat(Nband*Nqx,Nband*Nqx)

      integer :: ix, iy, is
      integer :: mu, nu
      integer :: iQ, lQ1
      real*8 :: akx, aky
      complex*16 :: zbloch(0:Nqx,0:Nsite-1)

      if ((ispin /= 1).and.(ispin /=2)) then
         write(*,*) 'Error in makeH  -- invalid spin'
         stop
      end if

c######################################################################
c     interaction.pdf: Eq. (12)
c       h_0 + epsilon_mu * delta_(mu,nu)
c       dkx,dky(=akx,aky) <-> k
c       ix,iy   <-> delta_x, delta_y
c######################################################################

      zhmat(:,:) = 0.0d0

      do iQ = 0, Nqx-1
         akx = dkx + DBLE(iQ) / DBLE(Nqx)
         aky = dky + DBLE(iQ) / DBLE(Nqy)
         call calcBlochFactor(zbloch(iQ,:),akx,aky)
      end do

      do is = 0, Nsite - 1
            do mu = 1,Nband ; do nu = 1,Nband
               do iQ = 0, Nqx-1
                  lQ1 = Nband * iQ
                  zhmat(mu+lQ1,nu+lQ1) =  zhmat(mu+lQ1,nu+lQ1)
     &                 + t(is,mu,nu) * zb(iQ,is)
               end do
         end do ; end do
      end do

c      call makeH_int(zhmat,ispin,Nband,Nqx)

      return

      End

c################################################################

c#########################################################
c#### bloch: calculate bloch factor at specific k-point
c#########################################################

      subroutine calcBlochFactor(zbloch,dkx,dky)

      use common, only: Nkx, Nky, Nsite, Isitex, Isitey, Pi, Zi
      implicit none

      real*8, intent(in) :: dkx, dky
      complex*16, intent(out) :: zbloch(0:Nsite-1)

      zbloch(:) = 0.0d0

      do is = 0, Nsite-1
         ix = isitex(is)
         iy = isitey(is)
         !** orthorhombicもこれでいいのか→よさそう
c         zbloch(is) = EXP( 2.0d0*Pi*Zi * (ix*dkx + iy*dky) )
         zbloch(is) = EXP( Pi*Zi*(ix*dkx+iy*dky) )
      end do

      return
      end

c################################################################

c#########################################################
c#### checkHermitian:
c#########################################################

      subroutine checkHermitian(nm,zhmat,info)

      implicit none

      integer, intent(out) :: info
      integer, intent(in) :: nm
      complex*16, intent(in) :: zhmat(nm,nm)

      integer :: i, j, idummy
      complex*16 :: zdummy

      info = 0
c## HERMITIAN CHECK[
      idummy = 0
      do i = 1, nm
         do j = i, nm
            if(AIMAG(zhmat(i,j)) /= 0) then
               idummy = 1
            end if
            zdummy = zhmat(i,j) - CONJG(zhmat(j,i))
            if(ABS(zdummy) > 0.1d-6)then
               write(*,*) 'ERROR -- H IS NOT HERMITIAN'
               write(*,*) i, j, zdummy
               write(*,*) zhmat(i,j), zhmat(j,i)
               write(*,*) '***'
               info = -1        ! Hamiltonian is not hermit
            end if
         end do
      end do
      if (idummy == 0) then
         info = 1               ! Hamiltonian is real
      end if
c## HERMITIAN CHECK]

      return
      End

c################################################################
