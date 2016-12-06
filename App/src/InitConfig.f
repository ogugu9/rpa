c#########################################################
c#### initConfig: initial configuration
c#########################################################
      subroutine initConfig(ierr)

      use common
      implicit none

      integer, intent(out) :: ierr

      integer :: magneticPattern, i, mu, nu
      real*8 :: dummy, dphase
      complex*16 :: zdummy

c######
c#    nuu: density of electrons with up spin
c#    ndd: density of electrons with dn spin
c######

      ierr = 0

      write(*,*)'  0: uniform   1: stripe   2: checkerboard'
      write(*,*)'  3: double stripe'
      write(*,*)' -1: continue'
      read(5,*,err=999) magneticPattern

      if(magneticPattern == -1) then
         ierr = -1
         return
      end if

c===========================================[
c== uniform ================================
c==
      if(magneticPattern == 0) then

         call reset
         Nqx = 1
         Nqy = 1
         Nqz = 1

         Nredx = 4 / Nqx
         Nredy = 4 / Nqy
         Nredz = 4 / Nqz

         Zop(:,:,:) =  0.0d0
         Zop2(:,:,:) =  0.0d0
c===========================================[
c== stripe =================================
c==
      else if(magneticPattern == 1) then

         call reset
         Nqx = 2
         Nqy = 1
         Nqz = 1

         !** Nredzはどうしたらいいか
         Nredx = 2 / Nqx
         Nredy = 2 / Nqy
         Nredz = 2 / Nqz

         call deallocation()
         call allocation()

c## INITIAL CONFIGURATION[
         dummy = 0.1d0

         Dnuu(:) = dummy
         Dndd(:) = -dummy

         ! ** ここは変更しなくて良いのか
         do mu = 1, Nband
            Zop(mu,mu,1) =  dummy
            Zop(mu,mu,2) = -dummy
         end do

         do mu = 1, Nband ; do nu = mu, Nband
            if (mu /= nu) then
               Zop(mu,nu,1) = cmplx(dummy,dummy)
               Zop(nu,mu,1) = cmplx(dummy,-dummy)
               Zop(mu,nu,2) = cmplx(-dummy,-dummy)
               Zop(nu,mu,2) = cmplx(-dummy,dummy)
            end if
         end do ; end do

         Zop(3,3,:) = Zop(2,2,:)

         Zop(1,2,:) = 0.0d0
         Zop(1,3,:) = 0.0d0
         Zop(1,4,:) = 0.0d0
         Zop(2,1,:) = 0.0d0
         Zop(2,4,:) = 0.0d0
         Zop(2,5,:) = 0.0d0
         Zop(3,1,:) = 0.0d0
         Zop(3,4,:) = 0.0d0
         Zop(3,5,:) = 0.0d0
         Zop(4,1,:) = 0.0d0
         Zop(4,2,:) = 0.0d0
         Zop(4,3,:) = 0.0d0
         Zop(4,5,:) = 0.0d0
         Zop(5,2,:) = 0.0d0
         Zop(5,3,:) = 0.0d0
         Zop(5,4,:) = 0.0d0
c         do mu = 1, Ns ; do nu = mu, Ns
c            if (mu /= nu) then
c               Zop(mu,nu,1) = dummy
c               Zop(nu,mu,1) = dummy
c               Zop(mu,nu,2) = -dummy
c               Zop(nu,mu,2) = -dummy
c            end if
c         end do ; end do

c## INITIAL CONFIGURATION]
c===========================================]

      else if(magneticPattern == 2) then
c===========================================[
c== checkerboard =================================
c==
         call reset

         Nqx = 2
         Nqy = 2
         Nqz = 2
         call deallocation()
         call allocation()

         Nredx = 2 / Nqx
         Nredy = 2 / Nqy
         Nredz = 2 / Nqz
         dummy = 0.1d0

c## INITIAL CONFIGURATION[
         Dnuu(:) = dummy
         Dndd(:) = -dummy

         do mu = 1, Nband
            Zop(mu,mu,1) =  dummy
            Zop(mu,mu,2) = -dummy
         end do

         do mu = 1, Nband ; do nu = mu, Nband
            if (mu /= nu) then
               Zop(mu,nu,1) = cmplx(dummy,dummy)
               Zop(nu,mu,1) = cmplx(dummy,-dummy)
               Zop(mu,nu,2) = cmplx(-dummy,-dummy)
               Zop(nu,mu,2) = cmplx(-dummy,dummy)
            end if
         end do ; end do

         Zop(3,3,1) = Zop(2,2,1)
c## INITIAL CONFIGURATION]
c===========================================]

      else if(magneticPattern == 3) then
c===========================================[
c== double stripe ==========================
c==
         call reset
         Nqx = 4
         Nqy = 4
         Nqz = 4
         call deallocation()
         call allocation()

         Nredx = 4 / Nqx
         Nredy = 4 / Nqy
         Nredz = 4 / Nqz


c## INITIAL CONFIGURATION[
         dummy = 0.1d0

         Zop2(:,:,:) = 0.0d0
         Zop0(:,:,:) = 0.0d0
         Zop(:,:,:) = 0.0d0

         Dnuu(:) = dummy
         Dndd(:) = -dummy

         dphase = 0.25d0 * Pi
         zdummy = EXP(Zi * dphase)
         do mu = 1, Nband
            Zop(mu,mu,1) =  dummy * zdummy
            Zop(mu,mu,2) = -Zop(mu,mu,1)

c            Zop2(mu,mu,:) = cmplx(dummy,dummy)
c            Zop2(mu,mu,:) = dummy
            Zop0(mu,mu,:) = Dne / DBLE(Nband) / 2.0d0
         end do

         do mu = 1, Nband ; do nu = mu, Nband
            if (mu /= nu) then
c               Zop(mu,nu,1) = cmplx(dummy,dummy)
c               Zop(nu,mu,1) = cmplx(dummy,-dummy)
c               Zop(mu,nu,2) = cmplx(-dummy,-dummy)
c               Zop(nu,mu,2) = cmplx(-dummy,dummy)

c               Zop2(mu,nu,:) = cmplx(dummy,dummy)
c               Zop2(nu,mu,:) = cmplx(dummy,-dummy)
c               Zop0(mu,mu,:) = cmplx(dummy,dummy)
c               Zop0(nu,mu,:) = cmplx(dummy,-dummy)
            end if
         end do ; end do



c## INITIAL CONFIGURATION]
c===========================================]

      else
         ierr = 1
      end if


      return


 999  ierr = 1
      return
      end
c########################################################

c#########################################################
c#### reset:
c#########################################################
      Subroutine reset

      use common, only: Zop
      implicit none

      Zop(:,:,:) = 0.0d0

      return
      End
c########################################################
