
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
