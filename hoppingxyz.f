c################ 2009 Kaneshita ################


c#########################################################
c#### hopping_integral:
c#########################################################

      subroutine hopping_integralxyz(finput)
      
      use common, only : Ns, Nrs, Nb, Eneorb,
     &     t, isitex, isitey, lln, llx, lly, Dkshift, Pi
      implicit none

      real*8, parameter :: tinit = 1.0d10

      character(len=*), intent(in) :: finput
      character :: chara*20
      integer :: i, mu, nu, ix, iy, is
      real*8 :: dummy(1:Ns)
      integer, external :: inverse, isigma_y

      open(10,file = finput, status = 'old')

      read(10,*) chara
      if(chara(1:3) == '##X') then
         call hopping_integral(finput)
         close(10)
         return
      end if
      backspace(10)


      t(:,:,:,:) = tinit
      t(0,0,:,:) = 0.0d0
      i = 0
      read(10,*) chara
      backspace(10)
      do while (chara(1:3) /= 'END')
         do while (chara(1:1) == '#')
            read(10,*) chara
         end do
         backspace(10)
         if (chara(1:3) /= 'END') then
            dummy(:) = 0.0d0
c            read(10,*) mu, nu, t(1,0,mu,nu), t(1,1,mu,nu),
c     &           t(2,0,mu,nu), t(2,1,mu,nu), t(2,2,mu,nu)
            read(10,*) mu, nu, (dummy(is), is= 1,Nb)
            t(1,0,mu,nu) = dummy(1)
            t(1,1,mu,nu) = dummy(2)
            t(2,0,mu,nu) = dummy(3)
            t(2,1,mu,nu) = dummy(4)
            t(2,2,mu,nu) = dummy(5)
            i = i + 1
         end if
         read(10,*) chara
      end do


      if (i /= (Ns + 1) * Ns / 2) then
         write(*,*) 'ERROR in hpping_integral'
         write(*,*) 'Wrong input data'
         stop
      end if

c================
      do while (chara(1:7) /= 'on-site')
         read(10,*) chara
      end do

      i = 0
c      read(10,*) chara
c      backspace(10)
      chara = '###'
      do while (chara(1:3) /= 'END')
         do while (chara(1:1) == '#')
            read(10,*) chara
         end do
c         write(*,*) chara
         backspace(10)
         if (chara(1:3) /= 'END') then
            read(10,*) nu, dummy(1)
            Eneorb(nu) = dummy(1)
         end if
         read(10,*) chara
      end do
c================

c================
      rewind(10)
      Dkshift(:) = 0.0d0
      read(10,*) chara
      do while (chara(1:11) /= 'END_OF_DATA')
         read(10,*) chara
         if (chara(1:7) == 'Dkshift') then
            backspace(10)
            read(10,*) chara, Dkshift(1), Dkshift(2)
            exit
         end if
      end do
c================


      do is = 1, Nrs
         ix = isitex(is)
         iy = isitey(is)
c#### ++ 00
         t(ix,iy,:,:) = t(ix,iy,:,:) * 0.1d0
c#### +- 00  = sigma_y(++00)
         call sigma_d( t( ix, iy,:,:), t( iy, ix,:,:))
         do mu = 1, Ns ; do nu = mu, Ns
c#### -- 00  = inverse(++00)
            t(-ix,-iy,mu,nu) = inverse(mu,nu) * t( ix, iy,mu,nu)
c#### -+ 00  = inverse(-+00)
            t(-iy,-ix,mu,nu) = inverse(mu,nu) * t( iy, ix,mu,nu)
c#### ++ 10  = sigma_y(++00)
            t( ix,-iy,mu,nu) = isigma_y(mu,nu) * t( ix, iy,mu,nu)
c#### -+ 10  = sigma_y(-+00)
            t( iy,-ix,mu,nu) = isigma_y(mu,nu) * t( iy, ix,mu,nu)
c#### +- 10  = sigma_y(+-00)
            t(-iy, ix,mu,nu) = isigma_y(mu,nu) * t(-iy,-ix,mu,nu)
c#### -- 10  = sigma_y(--00)
            t(-ix, iy,mu,nu) = isigma_y(mu,nu) * t(-ix,-iy,mu,nu)


c#### -- 01  = ++00
            t(-ix,-iy,nu,mu) = t( ix, iy,mu,nu)
c#### +- 01  = -+00
            t( ix,-iy,nu,mu) = t(-ix, iy,mu,nu)
c#### -+ 01  = +-00
            t(-ix, iy,nu,mu) = t( ix,-iy,mu,nu)
c#### ++ 01  = --00
            t( ix, iy,nu,mu) = t(-ix,-iy,mu,nu)
c#### -- 11  = ++10
            t(-iy,-ix,nu,mu) = t( iy, ix,mu,nu)
c#### +- 11  = -+10
            t( iy,-ix,nu,mu) = t(-iy, ix,mu,nu)
c#### -+ 11  = +-10
            t(-iy, ix,nu,mu) = t( iy,-ix,mu,nu)
c#### ++ 11  = --10
            t( iy, ix,nu,mu) = t(-iy,-ix,mu,nu)
         end do ; end do
      end do

c      write(*,*) SUM(t)
      do is = 0, lln - 1
         ix = isitex(is)
         iy = isitey(is)
         do mu = 1, Ns ; do nu = 1, Ns
            if (t(ix,iy,mu,nu) == tinit) then
               write(*,*) 'ERROR in hopping integrals'
               write(*,*) 't is not assigned'
               write(*,*) 'ix, iy, mu, nu =', ix, iy, mu, nu
            end if
         end do ; end do
      end do

c      write(*,*) SUM(t)
      do is = 0, lln - 1
         ix = isitex(is)
         iy = isitey(is)
         dummy(1) = (-1.0d0) ** (ix * Dkshift(1) + iy * Dkshift(2)) 
         t(ix,iy,:,:) = t(ix,iy,:,:) * dummy(1)
      end do
      Dkshift(:) = 0.0d0

      Dkshift(:) = Dkshift(:) * Pi

caaaaaaaaaaaaaaaaaaaaa
c      open(10, file='out.hop.dat')
c      do ix = -llx,llx ; do iy = -lly, lly 
c         write(10,'(a,2i3)') 'x,y=',ix,iy
c         write(10,'(25g8.2)') 
c     &        ( (REAL(t(ix,iy,mu,nu)*10.0d0),nu = mu, Ns),mu=1,Ns)
c      end do ; end do
c      close(10)

      close (10)
      return
      end


c#########################################################
c#### isigma_y: symmetry operation (x,y) --> (x,-y) for xyz coordinate
c#########################################################
      
      integer function isigma_y(mu,nu)
      
      use common, only : Ns
      implicit none
      
      integer, intent(in) :: mu, nu
      
      if (mu > nu) then
         write(*,*) 'Error in isigma_y'
         stop
      end if
      
      isigma_y = 1
      
      if (mu == nu) return
      if (mu == 1) then
         if ((nu == 2).or.(nu == 5)) then
            isigma_y = -isigma_y
         end if
         return
      end if
      if (mu == 2) then
         if ((nu == 3).or.(nu == 4)) then
            isigma_y = -isigma_y
         end if
         return
      end if
      if (mu == 3) then
         if (nu == 5) then
            isigma_y = -isigma_y
         end if
         return
      end if
      if (mu == 4) then
         if (nu == 5) then
            isigma_y = -isigma_y
         end if
         return
      end if
      
      return
      end
      
c#########################################################
c#### sigma_d: symmetry operation (x,y) --> (y,x) for xyz coordinate
c#########################################################

      subroutine sigma_d(tin,tout)
      
      use common, only: Ns
      implicit none

      real*8, intent(in) :: tin(Ns,Ns)
      real*8, intent(out) :: tout(Ns,Ns)
      integer :: mu,nu

      tout(1,1) =  tin(1,1)
      tout(1,2) =  tin(1,3)
      tout(1,3) =  tin(1,2)
      tout(1,4) = -tin(1,4)
      tout(1,5) =  tin(1,5)

      tout(2,2) =  tin(3,3)
      tout(2,3) =  tin(2,3)
      tout(2,4) = -tin(3,4)
      tout(2,5) =  tin(3,5)

      tout(3,3) =  tin(2,2)
      tout(3,4) = -tin(2,4)
      tout(3,5) =  tin(2,5)

      tout(4,4) =  tin(4,4)
      tout(4,5) = -tin(4,5)
      tout(5,5) =  tin(5,5)

      return
      end

c#########################################################

