c################ 2011 Kaneshita ################

c#########################################################
c#### hamiltoniank:
c#########################################################

      subroutine diag(eig,zpsi,dkx,dky,ispn)
      
      use common, only : Nkx, Nky, Ns, Nqx, Nqy, Check
      implicit none

      integer, intent(in) :: ispn
      real*8, intent(out) :: eig(Ns*Nqx)
      complex*16, intent(out) :: zpsi(Ns*Nqx,Ns*Nqx)
      real*8, intent(in) :: dkx, dky


      integer :: i, j, ispin, nm
      integer :: info, lwork
      real*8 :: ene(Ns*Nqx)
      real*8 :: rwork(3*Ns*Nqx-2)
      complex*16 :: zdummy
      complex*16 :: work((Ns*Nqx+1)*Ns*Nqx)

      complex*16 :: zhmat(Ns*Nqx,Ns*Nqx)


      ispin = ispn

      if ((ispin /= 1).and.(ispin /=2)) then
         write(*,*) 'Error in hamiltonian  -- invalid spin'
         stop
      end if


      nm = Ns * Nqx

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
c#### hamiltoniank:
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


c#########################################################
c#### makeH:
c#########################################################

      subroutine makeH(zhmat,dkx,dky,ispin)
      
      use common, only : llx, lly, lln, 
     &     Nkx, Nky, Ns, Ns2, Nqx, Nqy,
     &     Pi, Zi, U, t, DJ, Eall, Zpsiall, Eneorb,
     &     Isitex, Isitey, Check, Zop, Zop0, Zop2
      implicit none

      integer, intent(in), optional :: ispin
c      real*8, intent(out) :: eig(Ns*Nqx,2)
c      complex*16, intent(out) :: zpsi(Ns*Nqx,Ns*Nqx,2)
      real*8, intent(in) :: dkx, dky
      complex*16, intent(out) :: zhmat(Ns*Nqx,Ns*Nqx)


      integer :: i, j, ix, iy, idummy, is, ispin2
      integer :: iQ1, iQ2
      integer :: mu, nu, nm
      integer :: iQ, lQ1, lQ2
      integer :: info, lwork
      integer :: istest(-llx:llx,-lly:lly)
      real*8 :: akx, aky
      real*8 :: ep(Ns)
      complex*16 :: zb(0:Nqx,-llx:llx,-lly:lly)
      complex*16 :: zb1(-llx:llx,-lly:lly)
      complex*16 :: zb2(-llx:llx,-lly:lly)
      complex*16 :: zb3(-llx:llx,-lly:lly)
      complex*16 :: zb4(-llx:llx,-lly:lly)

      complex*16 :: zdummy

      complex*16 :: znop(Ns,Ns,2), znmu(Ns)
      complex*16 :: zn, znup


      if ((ispin /= 1).and.(ispin /=2)) then
         write(*,*) 'Error in makeH  -- invalid spin'
         stop
      end if

      if (ispin == 1) ispin2 = 2
      if (ispin == 2) ispin2 = 1

      nm = Ns * Nqx
      
      zhmat(:,:) = 0.0d0

c## on-site energy: !xyz
c      ep(1) = 10.75d0
c      ep(2) = 10.96d0
c      ep(3) = 10.96d0
c      ep(4) = 10.62d0
c      ep(5) = 11.12d0

c## on-site energy: !xyz
      do mu = 1, Ns
         ep(mu) = Eneorb(mu)
      end do
      

      do mu = 1, Ns
         do iQ = 0, Nqx-1
            zhmat(mu+Ns*iQ,mu+Ns*iQ) = zhmat(mu+Ns*iQ,mu+Ns*iQ) + ep(mu)
         end do
      end do


      do iQ = 0, Nqx-1
         akx = dkx + 2.0d0 * Pi / DBLE(Nqx) * DBLE(iQ)
         aky = dky + 2.0d0 * Pi / DBLE(Nqy) * DBLE(iQ)
         call blochk_any(zb(iQ,:,:),akx,aky)
      end do

      istest(:,:) = 0
      istest(0,0) = 1
      do is = 1, lln - 1
         ix = isitex(is)
         iy = isitey(is)
         istest(ix,iy) = istest(ix,iy) + 1
            do mu = 1, Ns ; do nu = 1, Ns
               do iQ = 0, Nqx-1
                  lQ1 = Ns * iQ
                  zhmat(mu+lQ1,nu+lQ1) =  zhmat(mu+lQ1,nu+lQ1)
     &                 + t(ix,iy,mu,nu) * zb(iQ,ix,iy)
               end do
         end do ; end do
      end do


      if (Check == 'y') then
         is = 1
         do ix = -llx, lly ; do iy = -lly, lly
            is = istest(ix,iy) * is
         end do ; end do
         if (is /= 1) then
            write(*,*) 'ERROR in hamiltonian'
            write(*,*) 'hopping term is not correct', is
         end if
      end if


      call makeH_int(zhmat,ispin,Ns,Nqx)
      return

      End

c################################################################

c#########################################################
c#### makeH_int:
c#########################################################

      subroutine makeH_int_old(zhmat,is,norb,nq)
      
      use common, only : Zi, U, DJ, Check, Zop, Zop0, Zop2
      implicit none

      integer, intent(in) :: is, norb, nq
      complex*16, intent(inout) :: zhmat(norb*nq,norb*nq)


      integer :: i, j, ix, iy, idummy, is1, is2
      integer :: iQ1, iQ2
      integer :: mu, nu, nm
      integer :: iQ, lQ1, lQ2

      complex*16 :: zdummy

      complex*16 :: znop(norb,norb,2), znmu(norb)
      complex*16 :: zn, znup
      complex*16, external ::  zham_U

      if ((is /= 1).and.(is /=2)) then
         write(*,*) 'Error in makeH  -- invalid spin'
         stop
      end if

      is1 = is
      if (is1 == 1) is2 = 2
      if (is1 == 2) is2 = 1

      nm = norb * nq
      

c==== interaction U&J[

      do iQ2 = 0, nq-1 ; do iQ1 = iQ2, nq-1
         lQ1 = norb * iQ1
         lQ2 = norb * iQ2

         if (iQ1-iQ2 == 0)then
            znop(:,:,:) = Zop0(:,:,:)
         else if (iQ1-iQ2 == 1)then
            znop(:,:,:) = Zop(:,:,:)
         else if (iQ1-iQ2 == 2)then
            znop(:,:,:) = Zop2(:,:,:)
         else if (iQ1-iQ2 == 3)then
            do mu = 1, norb ; do nu = 1, norb
               znop(nu,mu,:)= CONJG(Zop(mu,nu,:))
            end do ; end do
         end if

         znup = 0.0d0
         zn = 0.0d0
         do mu = 1, norb
            znup = znup + znop(mu,mu,is1)
            znmu(mu) = znop(mu,mu,is1) + znop(mu,mu,is2)
            zn = zn + znmu(mu)
         end do

         do mu = 1, norb ; do nu = 1, norb            

            if (mu == nu) then
c#### U
               zdummy =  - U * znop(mu,nu,is1)
               zdummy = CONJG(zdummy)
               call makeH_int_set(nm,zhmat,mu,nu,lQ1,lQ2,zdummy)
               
               zdummy = - DJ * (znup - znop(mu,mu,is1))
               zdummy = CONJG(zdummy)
               call makeH_int_set(nm,zhmat,mu,nu,lQ1,lQ2,zdummy)
            else

c#### inter-orbital
               zdummy = - (U - 2.0d0 * DJ) * znop(mu,nu,is1)
               zdummy = CONJG(zdummy)
               call makeH_int_set(nm,zhmat,mu,nu,lQ1,lQ2,zdummy)

c#### pair hopping
               zdummy = - DJ * DCONJG( Zop(nu,mu,is1))
               zdummy = CONJG(zdummy)
               call makeH_int_set(nm,zhmat,mu,nu,lQ1,lQ2,zdummy)

            end if

         end do ; end do
      end do ; end do
c==== interaction U&J]



      return
      End

c################################################################


c#########################################################
c#### makeH_int:
c#########################################################

      subroutine makeH_int(zhmat,is,norb,nq)
      
      use common, only : Zi, U, DJ, Check, Zop, Zop0, Zop2
      implicit none

      integer, intent(in) :: is, norb, nq
      complex*16, intent(inout) :: zhmat(norb*nq,norb*nq)


      integer :: i, j, ix, iy, idummy, is1, is2
      integer :: iQ1, iQ2
      integer :: mu, nu, nm
      integer :: iQ, lQ1, lQ2

      complex*16 :: zdummy

      complex*16 :: znop(norb,norb,2), znmu(norb)
      complex*16 :: zn, znup
      complex*16, external ::  zham_U

      if ((is /= 1).and.(is /=2)) then
         write(*,*) 'Error in makeH  -- invalid spin'
         stop
      end if

      is1 = is
      if (is1 == 1) is2 = 2
      if (is1 == 2) is2 = 1

      nm = norb * nq
      

c==== interaction U&J[

      do iQ2 = 0, nq-1 ; do iQ1 = iQ2, nq-1
         lQ1 = norb * iQ1
         lQ2 = norb * iQ2

         if (iQ1-iQ2 == 0)then
            znop(:,:,:) = Zop0(:,:,:)
         else if (iQ1-iQ2 == 1)then
            znop(:,:,:) = Zop(:,:,:)
         else if (iQ1-iQ2 == 2)then
            znop(:,:,:) = Zop2(:,:,:)
         else if (iQ1-iQ2 == 3)then
            do mu = 1, norb ; do nu = 1, norb
               znop(nu,mu,:)= CONJG(Zop(mu,nu,:))
            end do ; end do
         end if

         znup = 0.0d0
         zn = 0.0d0
         do mu = 1, norb
            znup = znup + znop(mu,mu,is1)
            znmu(mu) = znop(mu,mu,is1) + znop(mu,mu,is2)
            zn = zn + znmu(mu)
         end do

         do mu = 1, norb ; do nu = 1, norb            

            zdummy = zham_U(mu,nu,zn,znop(mu,nu,is1))
            call makeH_int_set(nm,zhmat,mu,nu,lQ1,lQ2,zdummy)

            if (mu == nu) then
c#### intra-orbital J
               zdummy = - DJ * ( 2.0d0 * zn + znup
     &              - 2.0d0 * znmu(mu) - znop(mu,mu,is1))
               zdummy = CONJG(zdummy)
               call makeH_int_set(nm,zhmat,mu,nu,lQ1,lQ2,zdummy)

            else                ! inter-orbital (mu /= nu)
c#### inter-orbital J
               zdummy = DJ * (3.0d0 * znop(mu,nu,is1)
     &              + znop(mu,nu,is2) + znop(nu,mu,is2))
               zdummy = CONJG(zdummy)
               call makeH_int_set(nm,zhmat,mu,nu,lQ1,lQ2,zdummy)

            end if

         end do ; end do
      end do ; end do
c==== interaction U&J]



      return
      End

c################################################################

c#########################################################
c#### makeH_int_set:
c#########################################################

      subroutine makeH_int_set(nsize,zhmat,mu,nu,lQ1,lQ2,zvalue)
      
      use common, only : Zi, U, DJ, Check, Zop, Zop0, Zop2
      implicit none

      integer, intent(in) :: nsize,mu, nu, lQ1, lQ2
      complex*16, intent(inout) :: zhmat(nsize,nsize)
      complex*16, intent(in) :: zvalue


      zhmat(mu+lQ1,nu+lQ2) = zhmat(mu+lQ1,nu+lQ2) + zvalue
      
      if (lQ1 /= lQ2) then
         zhmat(nu+lQ2,mu+lQ1) = zhmat(nu+lQ2,mu+lQ1) + CONJG(zvalue)
      end if


      return
      End

c################################################################

c#########################################################
c#### makeH_int_set:
c#########################################################

      complex*16 function zham_U(mu,nu,zn,znop)
      
      use common, only : Zi, U, DJ, Check, Zop, Zop0, Zop2
      implicit none

      integer, intent(in) :: mu, nu
      complex*16, intent(in) :: zn,znop


      if (mu == nu) then        !intra-orbital
         zham_U = U * (zn - znop) 
         
      else                      ! inter-orbital (mu /= nu)
         zham_U = - U * znop
      end if
      
      zham_U = CONJG(zham_U)

      return
      End

c################################################################

c#########################################################
c#### Hcheck:
c#########################################################

      subroutine Hcheck(dkx,dky,ispn,chara)
      
      use common, only : llx, lly, lln, 
     &     Nkx, Nky, Ns, Ns2, Nqx, Nqy,
     &     Pi, Zi, U, t, DJ, Eall, Zpsiall,
     &     Isitex, Isitey, Check, Zop, Zop0, Zop2
      implicit none

      character(len=*), intent(in) :: chara
      integer, intent(in) :: ispn
      real*8, intent(in) :: dkx, dky

      integer :: nm
      real*8 :: akx, aky
      complex*16 :: zhmat(Ns*Nqx,Ns*Nqx)


      nm = Ns * Nqx

      if ((ispn /= 1).and.(ispn /=2)) then
         write(*,*) 'Error in Hcheck  -- invalid spin'
         stop
      end if

      call makeH(zhmat,dkx,dky,ispn)

c## OUTPUT[
      open(10, file='Hmat.dat.'//chara)
      call writemtrx(zhmat,nm,nm,10,'(100g20.6)')
      close(10)
c## OUTPUT]

      return
      End

c################################################################


c#########################################################
c#### Hcheckall:
c#########################################################

      subroutine Hcheckall(nx,ny,ispn,chara)
      
      use common, only : llx, lly, lln, 
     &     Nkx, Nky, Ns, Ns2, Nqx, Nqy,
     &     Pi, Zi, U, t, DJ, Eall, Zpsiall,
     &     Isitex, Isitey, Check, Zop, Zop0, Zop2
      implicit none

      character(len=*), intent(in) :: chara
      integer, intent(in) :: nx, ny
      integer, intent(in) :: ispn

      character :: chx*4, chy*4, ch*40
      integer :: kx, ky, kkx, kky
      real*8 :: dkx, dky

      character(LEN=40), external :: wordcomb

      if ((ispn /= 1).and.(ispn /=2)) then
         write(*,*) 'Error in Hcheckall -- invalid spin'
         stop
      end if

      do kkx = 0, Nkx-1 ; do kky = 0, Nky-1
         kx = kkx
         ky = kky
         dkx = 2.0d0 * Pi * DBLE(kx) / DBLE(Nkx)
         dky = 2.0d0 * Pi * DBLE(ky) / DBLE(Nky)
         write(chx ,'(i4)') kx
         write(chy ,'(i4)') ky
         ch = wordcomb(chx,chy)
         ch = wordcomb(chara,ch)
         call Hcheck(dkx,dky,ispn,trim(ch))
      end do ; end do

      return
      End

c################################################################


c#########################################################
c#### eigenfncheck:
c#########################################################

      subroutine eigenfncheck(kx,ky,ispn,chara)
      
      use common, only : Ns, Nqx, Eall, Zpsiall, Zop
      implicit none

      character(len=*), intent(in) :: chara
      integer, intent(in), optional :: ispn, kx, ky

      character :: chsrfx*40, frmt*40
      integer :: idfile, nm, ispin
      complex*16 :: zmtrx(Ns*Nqx,Ns*Nqx)


      chsrfx = chara
      nm = Ns * Nqx

      if (present(ispn)) then
         ispin = ispn
      else
         ispin = 1
      end if

      if ((ispin /= 1).and.(ispin /=2)) then
         write(*,*) 'Error in hamiltonian  -- invalid spin'
         stop
      end if

      
      frmt = '(100f11.6)'
      idfile = 10
      open(idfile, file='psi.dat.'//chsrfx)
      call writemtrx( Zpsiall(kx,ky,1:nm,1:nm,ispin),nm,nm,idfile,frmt)
      close(idfile)

      open(idfile, file='ene.dat.'//chsrfx)
      zmtrx(1:nm,1) = CMPLX(Eall(kx,ky,1:nm,ispin),0.0d0)
      call writemtrx( zmtrx(1:nm,1),nm,1,idfile,frmt)
      close(idfile)

      open(idfile, file='op.dat.'//chsrfx)
      call writemtrx( Zop(1:Ns,1:Ns,1),Ns,Ns,idfile,frmt)
      close(idfile)


      return
      end

c################################################################

c#########################################################
c#### writemtrx:
c#########################################################

      subroutine writemtrx(zmtrx,nsize1,nsize2,fid,frmt)

      character, intent(in) :: frmt*40
      integer, intent(in) :: nsize1, nsize2, fid
      complex*16, intent(in) :: zmtrx(nsize1,nsize2)

      character :: chform*40
      integer :: i, j


c## OUTPUT[
      write(fid,'(a)') '#real_part'
      do i = 1, nsize1
         write(fid,frmt) (DBLE(zmtrx(i,j)), j = 1,nsize2)
      end do
      write(fid,'(a)') '#imaginary_part'
      write(fid,*) 
      do i = 1, nsize1
         write(fid,frmt) (AIMAG(zmtrx(i,j)), j = 1,nsize2)
      end do
c## OUTPUT]


      return
      End

c################################################################


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
c#### bloch: calculate bloch factor
c#########################################################

      subroutine bloch(zb,kkx,kky,mkx,mky)
      
      use common, only: llx, lly, Nkx, Nky, Pi, Zi, Dkshift
      implicit none

      integer, intent(in) :: kkx, kky, mkx, mky
      complex*16, intent(out) :: zb(-llx:llx,-lly:lly)

      integer :: kx, ky
      real*8 :: akx, aky
      complex*16 :: zblochx, zblochy

      kx = MOD(kkx,mkx)
      ky = MOD(kky,mky)

c# bloch factor:  e^(i * k * a) for a=r-r' 
      akx = 2.0d0 * Pi / DBLE(mkx) * DBLE(kx)
      aky = 2.0d0 * Pi / DBLE(mky) * DBLE(ky)
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

      return
      end

c################################################################


c#########################################################
c#### bloch: calculate bloch factor
c#########################################################

      subroutine blochk_mesh(zb,kkx,kky,mkx,mky,dkx,dky)
      
      use common, only: llx, lly, Nkx, Nky, Pi, Zi, Dkshift
      implicit none

      integer, intent(in) :: kkx, kky, mkx, mky
      complex*16, intent(out) :: zb(-llx:llx,-lly:lly)
      real*8, intent(in) :: dkx, dky

      integer :: kx, ky
      real*8 :: akx, aky
      complex*16 :: zblochx, zblochy

      kx = MOD(kkx,mkx)
      ky = MOD(kky,mky)

c# bloch factor:  e^(i * k * a) for a=r-r' 
      akx = 2.0d0 * Pi / DBLE(mkx) * DBLE(kx) + dkx
      aky = 2.0d0 * Pi / DBLE(mky) * DBLE(ky) + dky
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

      return
      end

c################################################################

c#########################################################
c#### bloch: calculate bloch factor
c#########################################################

      subroutine blochk_any(zb,dkx,dky)
      
      use common, only: llx, lly, Nkx, Nky, Pi, Zi, Dkshift
      implicit none

      complex*16, intent(out) :: zb(-llx:llx,-lly:lly)
      real*8, intent(in) :: dkx, dky

      integer :: kx, ky
      real*8 :: akx, aky
      complex*16 :: zblochx, zblochy


      zb(:,:) = 0.0d0

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

      return
      end

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


c#########################################################
c#### hamiltoniank:
c#########################################################

      subroutine h_test(eig,zpsi,dkx,dky,ispn)
      
      use common, only : llx, lly, lln, 
     &     Nkx, Nky, Ns, Ns2, Nqx, Nqy,
     &     Pi, Zi, U, t, DJ, Eall, Zpsiall,Eneorb,
     &     Isitex, Isitey, Check, Zop, Zop0, Zop2
      implicit none

      integer, intent(in), optional :: ispn
      real*8, intent(out) :: eig(Ns*Nqx,2)
      complex*16, intent(out) :: zpsi(Ns*Nqx,Ns*Nqx,2)
      real*8, intent(in) :: dkx, dky


      integer :: i, j, ix, iy, idummy, is, ispin, jspin, ispin2
      integer :: iQ1, iQ2
      integer :: mu, nu, nm
      integer :: iQ, lQ1, lQ2
      integer :: info, lwork
      integer :: istest(-llx:llx,-lly:lly)
      real*8 :: rwork(3*Ns*Nqx-2)
      real*8 :: dnup, dndn
      real*8 :: Eu(Ns*Nqx), Ed(Ns*Nqx)
      real*8 :: akx, aky
      real*8 :: ep(Ns)
      complex*16 :: zb(0:Nqx,-llx:llx,-lly:lly)
      complex*16 :: zb1(-llx:llx,-lly:lly)
      complex*16 :: zb2(-llx:llx,-lly:lly)
      complex*16 :: zb3(-llx:llx,-lly:lly)
      complex*16 :: zb4(-llx:llx,-lly:lly)
      complex*16 :: zhmu(Ns*Nqx,Ns*Nqx)
c      complex*16 :: zhmd(Ns2,Ns2)
      complex*16 :: zblochx, zblochy, zdummy
      complex*16 :: work((Ns*Nqx+1)*Ns*Nqx)

      complex*16 :: znop(Ns,Ns,2), znmu(Ns)
      complex*16 :: zn, znup

c      pair = 2.0 + 0.8

      if (present(ispn)) then
         ispin = ispn
      else
         ispin = 1
      end if

      if ((ispin /= 1).and.(ispin /=2)) then
         write(*,*) 'Error in hamiltonian  -- invalid spin'
         stop
      end if

      if (ispin == 1) ispin2 = 2
      if (ispin == 2) ispin2 = 1

      nm = Ns * Nqx

      lwork = (nm + 1) * nm
      
      zhmu(:,:) = 0.0d0

c## on-site energy:
c      ep(1) = 10.75d0
c      ep(2) = 10.96d0
c      ep(3) = 10.96d0
c      ep(4) = 11.12d0
c      ep(5) = 10.62d0

c## on-site energy: !xyz
c      ep(1) = 10.75d0
c      ep(2) = 10.96d0
c      ep(3) = 10.96d0
c      ep(4) = 10.62d0
c      ep(5) = 11.12d0
      
      do mu = 1, Ns
         ep(mu) = Eneorb(mu)
      end do

      do mu = 1, Ns
         do iQ = 0, Nqx-1
            zhmu(mu+Ns*iQ,mu+Ns*iQ) = zhmu(mu+Ns*iQ,mu+Ns*iQ) + ep(mu)
         end do
      end do

      do iQ = 0, Nqx-1
         akx = dkx + 2.0d0 * Pi / DBLE(Nqx) * DBLE(iQ)
         aky = dky + 2.0d0 * Pi / DBLE(Nqy) * DBLE(iQ)
         call blochk_any(zb(iQ,:,:),akx,aky)
      end do

      istest(:,:) = 0
      istest(0,0) = 1
      do is = 1, lln - 1
         ix = isitex(is)
         iy = isitey(is)
         istest(ix,iy) = istest(ix,iy) + 1
            do mu = 1, Ns ; do nu = 1, Ns
               do iQ = 0, Nqx-1
                  lQ1 = Ns * iQ
                  zhmu(mu+lQ1,nu+lQ1) =  zhmu(mu+lQ1,nu+lQ1)
     &                 + t(ix,iy,mu,nu) * zb(iQ,ix,iy)
               end do
         end do ; end do
      end do

      if (Check == 'y') then
         is = 1
         do ix = -llx, lly ; do iy = -lly, lly
            is = istest(ix,iy) * is
         end do ; end do
         if (is /= 1) then
            write(*,*) 'ERROR in hamiltonian'
            write(*,*) 'hopping term is not correct', is
         end if
      end if


c==== interaction U&J[

      do iQ2 = 0, Nqx-1 ; do iQ1 = iQ2, Nqx-1
         lQ1 = Ns * iQ1
         lQ2 = Ns * iQ2

         if (iQ1-iQ2 == 0)then
            znop(:,:,:) = Zop0(:,:,:)
         else if (iQ1-iQ2 == 1)then
            znop(:,:,:) = Zop(:,:,:)
         else if (iQ1-iQ2 == 2)then
            znop(:,:,:) = Zop2(:,:,:)
         else if (iQ1-iQ2 == 3)then
            do mu = 1, Ns ; do nu = 1, Ns
               znop(nu,mu,:)= DCONJG(Zop(mu,nu,:))
            end do ; end do
         end if

         znup = 0.0d0
         zn = 0.0d0
         do mu = 1, Ns
            znup = znup + znop(mu,mu,ispin)
            znmu(mu) = znop(mu,mu,ispin) + znop(mu,mu,ispin2)
            zn = zn + znmu(mu)
         end do

         do mu = 1, Ns ; do nu = 1, Ns            

            if (mu == nu) then

c#### intra-orbital U
               zdummy = U * (zn - znop(mu,mu,ispin))

               zhmu(mu+lQ1,mu+lQ2) =  zhmu(mu+lQ1,mu+lQ2)
     &              + DCONJG(zdummy)

               if (lQ1 /= lQ2) then
                  zhmu(mu+lQ2,mu+lQ1) = zhmu(mu+lQ2,mu+lQ1) + zdummy
               end if
               
c#### intra-orbital J
               zdummy = - DJ * ( 2.0d0 * zn + znup
     &              - 2.0d0 * znmu(mu) - znop(mu,mu,ispin))

               zhmu(mu+lQ1,mu+lQ2) =  zhmu(mu+lQ1,mu+lQ2)
     &              + DCONJG(zdummy)

               if (lQ1 /= lQ2) then
                  zhmu(mu+lQ2,mu+lQ1) = zhmu(mu+lQ2,mu+lQ1) + zdummy
               end if
            else

c#### inter-orbital U
               zdummy = - U * znop(mu,nu,ispin)

               zhmu(mu+lQ1,nu+lQ2) =  zhmu(mu+lQ1,nu+lQ2)
     &              + DCONJG(zdummy)

               if (lQ1 /= lQ2) then
                  zhmu(nu+lQ2,mu+lQ1) = zhmu(nu+lQ2,mu+lQ1) + zdummy
               end if

c#### inter-orbital J
               zdummy = DJ * (3.0d0 * znop(mu,nu,ispin)
     &              + znop(mu,nu,ispin2) + znop(nu,mu,ispin2))

               zhmu(mu+lQ1,nu+lQ2) =  zhmu(mu+lQ1,nu+lQ2)
     &              + DCONJG(zdummy)

               if (lQ1 /= lQ2) then
                  zhmu(nu+lQ2,mu+lQ1) = zhmu(nu+lQ2,mu+lQ1) + zdummy
               end if

            end if

         end do ; end do
      end do ; end do
c==== interaction U&J]



c## HERMITIAN CHECK[
      if (Check == 'y') then
         idummy = 0
         do i = 1, Ns*Nqx
            do j = i, Ns*Nqx
               if(DIMAG(zhmu(i,j)) /= 0) then
                  idummy = 1
               end if
               zdummy = zhmu(i,j) - DCONJG(zhmu(j,i))
               if(ABS(zdummy) > 0.1d-6)then
                  write(*,*) 'ERROR -- H IS NOT HERMITIAN'
                  write(*,*) i, j, zdummy
                  write(*,*) zhmu(i,j), zhmu(j,i)
                  write(*,*) '***'
                  stop
               end if
c     write(*,*)i,j,DBLE(zhm(i,j))
            end do
         end do
         if (idummy == 0) then
            write(*,*)'Hamiltonian is real',dkx,dky
         end if
      end if
c## HERMITIAN CHECK]


c## DIAGONALIZATION[
      Info = 0
      call  zheev('V','U',nm,zhmu(1:nm,1:nm),nm,Eu(1:nm),
     &     work,lwork,rwork,info)
      if(Info /= 0) then
         write(*,*) 'Lapack ZHEEV: Info=',Info
         stop
      end if

      eig(1:nm,ispn) = Eu(1:nm)
      zpsi(1:nm,1:nm,ispin) = zhmu(1:nm,1:nm)

c## DIAGONALIZATION]


      return
      End

c################################################################


c#########################################################
c#### energy_check               EK 2015
c#########################################################

      subroutine energy_check()

      use common, only : Nkx, Nky, Ns, Pi, Eall, Zpsiall, 
     &     Dmu, Nqx, Nqy
      implicit none

      integer :: ikx, iky, kx, ky
      real*8 :: akx, aky, pi2
      real*8 :: eig(Ns*Nqx,2)
      complex*16 :: zpsi(Ns*Nqx,Ns*Nqx,2)

      pi2 = 2.0d0 * Pi
      do ikx = 0, Nkx ; do iky = 0, Nky
         kx = MOD(ikx + Nkx, Nkx)
         ky = MOD(iky + Nky, Nky)
         
         akx = DBLE(kx) / DBLE(Nkx) * pi2
         aky = DBLE(ky) / DBLE(Nky) * pi2
         call  diag(eig(:,1),zpsi(:,:,1),akx,aky,1)
         
         write(*,'(20f10.6)') eig(:,1)-Dmu
         write(*,'(20f10.6)') Eall(kx,ky,:,1)-Dmu
         write(*,'(20f10.6)') eig(:,1)-Eall(kx,ky,:,1)
         write(*,'(20f10.6)') SUM(ABS(eig(:,1)-Eall(kx,ky,:,1)))
         write(*,'(a)') ' '
      end do ; end do

      return
      end 
c###################################################################

c%%%%
c%%%% END_OF_FILE:hamiltonian2.f
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
