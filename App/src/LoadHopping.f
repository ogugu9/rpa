c#########################################################
c#### loadHopping
c#########################################################

      subroutine loadHopping(finput)

c     ** Nbっていらなくない？
      use common, only : Nband, Nsite, Eorbit, t,
     &     isitex, isitey
      implicit none

      real*8, parameter :: tinit = 1.0d10
      real*8 :: val

      character(len=*), intent(in) :: finput
      character :: chara*20
      integer :: i, mu, nu, ix, iy, is

      open(10,file = finput, status = 'old')

c     // Number Of Wannier
c      read(10,*) chara
c     // Number of Sites
c      read(10,*) chara
c     // Number of Degeneracy
c     read(10,('15I5')) (ndegen(i),i=1,Nir)

c     // Init hopping
      t(:,:,:) = tinit

      write(*,*) 'loading hopping data...'

      is = 0
c     // Read hopping
c      write(*,*) '   is   ix   iy   mu   nu   t'
c      write(*,*) '-----------------------------------------'
      do
         read (10,'(4I5,F12.6)',end=180) ix, iy, mu, nu, val
c           // Siteに対する添え字のラベリングをinputから行う
c           if (ix <= 2 && iy <= 2) then
         if((mu.eq.1).and.(nu.eq.1)) then
            is = is + 1
            isitex(is) = ix
            isitey(is) = iy
         end if
         t(is,mu,nu) = val
c         write(*,'(5I5,F12.6)') is, ix, iy, mu, nu, val
c           end if
      end do

c     // バンドに対する対角成分を軌道エネルギーEorbitとして保持
180   do mu = 1,Nband
         Eorbit(mu) = t(1,mu,mu)
      end do

      if (is .ne. Nsite) then
         write(*,*)
         write(*,*) '!! ERROR in LoadHopping !!'
         print "('# sites = ', I3)", is
      else
         write(*,*) 'loading hopping data is done.'
      end if

      return
      end
c#########################################################
