c#########################################################
c#### inicon: Initial Configuration
c#########################################################
      subroutine inicon(ierr)

      use common
      implicit none

      integer, intent(out) :: ierr
 
      integer :: idummy, i, mu, nu
      real*8 :: dummy, dphase
      complex*16 :: zdummy


c######
c#    Dnuu: density of electrons with up spin
c#    Dndd: density of electrons with dn spin
c######

      ierr = 0

      write(*,*)'  0: uniform   1: stripe   2: checkerboard'
      write(*,*)'  3: double stripe'
      write(*,*)' -1: continue'
      read(5,*,err=999) idummy

      if(idummy == -1) then
         ierr = -1
         return
      end if

      if(idummy == 0) then

c===========================================[
c== uniform ================================
c==
         call reset
         Nqx = 1
         Nqy = 1
         
         Nredx = 4 / Nqx
         Nredy = 4 / Nqy

         Zop(:,:,:) =  0.0d0
         Zop2(:,:,:) =  0.0d0
         
         
c===========================================]

      else if(idummy == 1) then

c===========================================[
c== stripe =================================
c==
         call reset
         Nqx = 2
         Nqy = 1
         
         Nredx = 2 / Nqx
         Nredy = 2 / Nqy
         
         call deallocation()
         call allocation()

c## INITIAL CONFIGURATION[
         dummy = 0.1d0
         
         Dnuu(:) = dummy
         Dndd(:) = -dummy

         do mu = 1, Ns
            Zop(mu,mu,1) =  dummy
            Zop(mu,mu,2) = -dummy
         end do

         do mu = 1, Ns ; do nu = mu, Ns
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

      else if(idummy == 2) then
c===========================================[
c== checkerboard =================================
c==
         call reset

         Nqx = 2
         Nqy = 2
         call deallocation()
         call allocation()

         Nredx = 2 / Nqx
         Nredy = 2 / Nqy
         dummy = 0.1d0
         
c## INITIAL CONFIGURATION[
         Dnuu(:) = dummy
         Dndd(:) = -dummy

         do mu = 1, Ns
            Zop(mu,mu,1) =  dummy
            Zop(mu,mu,2) = -dummy
         end do

         do mu = 1, Ns ; do nu = mu, Ns
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
         
      else if(idummy == 3) then
c===========================================[
c== double stripe ==========================
c==
         call reset
         Nqx = 4
         Nqy = 4
         call deallocation()
         call allocation()
         
         Nredx = 4 / Nqx
         Nredy = 4 / Nqy
         

c## INITIAL CONFIGURATION[
         dummy = 0.1d0

         Zop2(:,:,:) = 0.0d0
         Zop0(:,:,:) = 0.0d0
         Zop(:,:,:) = 0.0d0
         
         Dnuu(:) = dummy
         Dndd(:) = -dummy

         dphase = 0.25d0 * Pi
         zdummy = EXP(Zi * dphase)
         do mu = 1, Ns
            Zop(mu,mu,1) =  dummy * zdummy
            Zop(mu,mu,2) = -Zop(mu,mu,1)

c            Zop2(mu,mu,:) = cmplx(dummy,dummy)
c            Zop2(mu,mu,:) = dummy
            Zop0(mu,mu,:) = Dne / DBLE(Ns) / 2.0d0
         end do

         do mu = 1, Ns ; do nu = mu, Ns
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

      use common, only: Dnuu, Dndd, Zop
      implicit none
      
      Zop(:,:,:) = 0.0d0
                  
      return
      End
c########################################################
