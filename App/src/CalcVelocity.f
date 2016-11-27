
c#########################################################
c#### velocity:
c#########################################################

      subroutine velocity(vel,zpsi,dkx,dky,xory)

      use common, only : llx, lly, lln,
     &     Nkx, Nky, Ns, Ns2, Nqx, Nqy,
     &     Pi, Zi, U, t, DJ, Eall, Zpsiall,
     &     Isitex, Isitey, Check, Zop
      implicit none

      character, intent(in) :: xory*1
      real*8, intent(out) :: vel(Ns*Nqx)
      complex*16, intent(in) :: zpsi(Ns*Nqx,Ns*Nqx)
      real*8, intent(in) :: dkx, dky

      integer :: i, j, k, ix, iy, idummy, is
      integer :: mu, nu, nm
      integer :: istest(-llx:llx,-lly:lly)
      real*8 :: dummat1(Ns*Nqx,Ns*Nqx), dummat2(Ns*Nqx,Ns*Nqx)
      real*8 :: akx, aky, dummy
      complex*16 :: zb1(-llx:llx,-lly:lly)
      complex*16 :: zb2(-llx:llx,-lly:lly)
      complex*16 :: zb3(-llx:llx,-lly:lly)
      complex*16 :: zb4(-llx:llx,-lly:lly)
      complex*16 :: zhmu(Ns*Nqx,Ns*Nqx)
      complex*16 :: zblochx, zblochy, zdummy
      complex*16 :: work((Ns*Nqx+1)*Ns*Nqx)

      nm = Ns * Nqx

      zhmu(:,:) = 0.0d0

      call blochk_bibun(zb1,dkx,dky,xory)
      if (Nqx >= 2) then
         akx = dkx + 2.0d0 * Pi / DBLE(Nqx)
         aky = dky + 2.0d0 * Pi / DBLE(Nqy)
         call blochk_bibun(zb2,akx,aky,xory)
      end if
      if (Nqx >= 3) then
         akx = dkx + 2.0d0 * Pi / DBLE(Nqx) * 2.0d0
         aky = dky + 2.0d0 * Pi / DBLE(Nqy) * 2.0d0
         call blochk_bibun(zb3,akx,aky,xory)
      end if
      if (Nqx >= 4) then
         akx = dkx + 2.0d0 * Pi / DBLE(Nqx) * 3.0d0
         aky = dky + 2.0d0 * Pi / DBLE(Nqy) * 3.0d0
         call blochk_bibun(zb4,akx,aky,xory)
      end if

      istest(:,:) = 0
      istest(0,0) = 1
      do is = 1, lln - 1
         ix = isitex(is)
         iy = isitey(is)
         istest(ix,iy) = istest(ix,iy) + 1
         do mu = 1, Ns ; do nu = 1, Ns
            zhmu(mu,nu) =  zhmu(mu,nu)
     &           + t(ix,iy,mu,nu) * zb1(ix,iy)
            if (Nqx >= 2) then
               zhmu(mu+Ns,nu+Ns) =  zhmu(mu+Ns,nu+Ns)
     &              + t(ix,iy,mu,nu) * zb2(ix,iy)
            end if
            if (Nqx >= 3) then
               zhmu(mu+Ns*2,nu+Ns*2) =  zhmu(mu+Ns*2,nu+Ns*2)
     &              + t(ix,iy,mu,nu) * zb3(ix,iy)
            end if
            if (Nqx >= 4) then
               zhmu(mu+Ns*3,nu+Ns*3) =  zhmu(mu+Ns*3,nu+Ns*3)
     &              + t(ix,iy,mu,nu) * zb4(ix,iy)
            end if
         end do ; end do
      end do


      vel(:) =0.0d0
      do i = 1, nm ; do j = 1, nm ; do k = 1, nm
         vel(i) = vel(i) + CONJG(zpsi(j,i)) * zhmu(j,k) * zpsi(k,i)
      end do ; end do ; end do

      return

      dummat1(:,:) =  MATMUL(zhmu(:,:), zpsi(:,:))

      dummat2(:,:) =0.0d0
      do i = 1, nm ; do j = 1, nm ; do k = 1, nm
         dummat2(i,j) = dummat2(i,j) + CONJG(zpsi(k,i)) * dummat1(k,j)
      end do ; end do ; end do

      do i = 1, nm
         vel(i) =  dummat2(i,i)
      end do

c## HERMITIAN CHECK[
      if (Check == 'y') then
         dummy = SUM(ABS(dummat2(:,:))) - SUM(ABS(vel(:)))
         if (dummy > 1.0d-10) then
            write(*,*) 'error in velocity', dummy
            stop
         end if
      end if
c## HERMITIAN CHECK]


      return
      End

c################################################################

c#########################################################
c#### bloch_bibun: calculate bloch factor
c#########################################################

      subroutine blochk_bibun(zb,dkx,dky,xory)

      use common, only: llx, lly, Nkx, Nky, Pi, Zi, Dkshift
      implicit none

      character, intent(in) :: xory*1
      complex*16, intent(out) :: zb(-llx:llx,-lly:lly)
      real*8, intent(in) :: dkx, dky

      integer :: kx, ky, i
      real*8 :: akx, aky
      complex*16 :: zblochx, zblochy


c# bloch factor:  e^(i * k * a) for a=r-r'
      akx = dkx
      aky = dky

c## need a (pi,pi) shift when using the hopping integral data of Imada group
c      akx = akx + Dkshift(1)
c      aky = aky + Dkshift(2)

      zblochx = EXP(Zi * akx)
      zblochy = EXP(Zi * aky)

      zb( 0, 0) = 1.0d0

      zb( 1, 0) = zblochx
      zb(-1, 0) = CONJG(zblochx)
      zb( 0, 1) = zblochy
      zb( 0,-1) = CONJG(zblochy)

      zb( 2, 0) = zb( 1, 0) * zb( 1, 0)
      zb(-2, 0) = zb(-1, 0) * zb(-1, 0)
      zb( 0, 2) = zb( 0, 1) * zb( 0, 1)
      zb( 0,-2) = zb( 0,-1) * zb( 0,-1)

      zb( 1, 1) = zb( 1, 0) * zb( 0, 1)
      zb( 1,-1) = zb( 1, 0) * zb( 0,-1)
      zb(-1, 1) = zb(-1, 0) * zb( 0, 1)
      zb(-1,-1) = zb(-1, 0) * zb( 0,-1)

      zb( 1, 2) = zb( 1, 0) * zb( 0, 2)
      zb( 1,-2) = zb( 1, 0) * zb( 0,-2)
      zb(-1, 2) = zb(-1, 0) * zb( 0, 2)
      zb(-1,-2) = zb(-1, 0) * zb( 0,-2)

      zb( 2, 1) = zb( 2, 0) * zb( 0, 1)
      zb( 2,-1) = zb( 2, 0) * zb( 0,-1)
      zb(-2, 1) = zb(-2, 0) * zb( 0, 1)
      zb(-2,-1) = zb(-2, 0) * zb( 0,-1)

      zb( 2, 2) = zb( 2, 0) * zb( 0, 2)
      zb( 2,-2) = zb( 2, 0) * zb( 0,-2)
      zb(-2, 2) = zb(-2, 0) * zb( 0, 2)
      zb(-2,-2) = zb(-2, 0) * zb( 0,-2)

c==== bibun[
      if (xory == 'x') then
         do i = -2, 2
            zb(i,:) = Zi * DBLE(i) * zb( i, :)
         end do
      else if (xory == 'y') then
         do i = -2, 2
            zb(:,i) = Zi * DBLE(i) * zb( :, i)
         end do
      end if
c==== bibun]

      return
      end

c################################################################

c###################################################################
c##### Svelo:
c###################################################################

      subroutine Svelo(akx, aky, vx, vy, ispn)

      use common, only : Nqx, Nqy, Ns, Ns2
      implicit none

      integer, intent(in) :: ispn
      real*8, intent(in) :: akx, aky
      real*8, intent(out) :: vx(Ns*Nqx), vy(Ns*Nqx)

      integer :: kkx, kky, kx, ky
      real*8 :: eig(Ns*Nqx,2)
      complex*16 :: zpsi(Ns*Nqx,Ns*Nqx,2)

      character :: chara*15
      integer :: i, j, ispin

       ispin = ispn
      call diag(eig(:,ispin),zpsi(:,:,ispin),akx,aky,ispin)
      call velocity(vx,zpsi(:,:,ispin),akx,aky,'x')
      call velocity(vy,zpsi(:,:,ispin),akx,aky,'y')


      return
      end

c###################################################################
