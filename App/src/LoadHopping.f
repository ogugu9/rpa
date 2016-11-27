c#########################################################
c#### loadHopping
c#########################################################

      subroutine loadHopping(finput)

c     ** Nbっていらなくない？
      use common, only : Nband, Nsite, Eorbit, t,
     &     isitex, isitey
      implicit none

      real*8, parameter :: tinit = 1.0d10
      real*8, parameter :: val

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
      t(:,:,:,:) = tinit
      t(0,0,:,:) = 0.0d0
      Nsite = 0

      write(*,*)
      write(*,*) 'loading hopping data...'
      write(*,*) '---------------------------------'

c     // Read hopping
      do
         do mu = 1,Nband ; do nu = 1,Nband
            read (10,('4I5,2F12.6'),end=180) ix, iy, mu, nu, val
c           // Siteに対する添え字のラベリングをinputから行う
c           if (ix <= 2 && iy <= 2) then
            if (mu == 1 && nu == 1) then
               isitex(is) = ix
               isitey(is) = iy
               is = is + 1
            end if
            t(ix,iy,mu,nu) = val
c           end if
         end do
      end do

180   write(*,*) 'loading hopping data is done.'
      write(*,*) '---------------------------------'

c     // バンドに対する対角成分を軌道エネルギーEorbitとして保持
      do mu = 1,Nband
         Eorbit(mu) = t(0,0,mu,mu)
      end do

      if (is /= Nsite) then
         write(*,*) 'ERROR in labelsite'
         write(*,*) '# of neighbors=', is
      end if

c#########################################################
