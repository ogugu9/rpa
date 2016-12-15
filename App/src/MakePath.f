c###################################################################
c##   makePath: make k-point path for bandplot or suscept
c###################################################################

      subroutine makePath()

      use common, only : req, Nkmesh, Nkpath,
     &      dkfrac, kpath
      implicit none

      integer :: i, j, k, ik
      real*8 :: kini(1:3), kfin(1:3)
      character*128 :: chara

      write(*,*) 'loading kpath for bandplot ...'

      read(5,*,err=110) Nkmesh
      print "(' => Nkmesh  = ', I5)", Nkmesh
      read(5,*,err=110) Nkpath
      print "(' => Nkpath  = ', I5)", Nkpath

      call deallocPath()
      call allocPath()

110   write(*,*)
      write(*,*) '#k-point path for bandplot'
      read(5,*) chara
      do i = 1,Nkpath
         read(5, '(A4,3F6.2,A4,3F6.2)'),
     &      kpath(i)%iniName, kpath(i)%iniPosition,
     &      kpath(i)%finName, kpath(i)%finPosition
         write(*, '(A4,3F6.2,3X,A4,3F6.2)'),
     &      kpath(i)%iniName, kpath(i)%iniPosition,
     &      kpath(i)%finName, kpath(i)%finPosition
      end do
      read(5,*) chara ! chara = ## KPOINT PATH ##
      write(*,*)

200   write(*,*) '#k-point list for bandplot'
      ik = 0
      do i = 1, Nkpath
         kini(:) = kpath(i)%iniPosition
         kfin(:) = kpath(i)%finPosition
         do j = 0,Nkmesh-1
            dkfrac(ik,:) = kini(:)
     &      + (kfin(:)-kini(:)) * DBLE(j)/DBLE(Nkmesh)
            write(*, '(I3,3F6.2)') ik, (dkfrac(ik,k),k=1,3)
            ik = ik + 1
            if (ik == Nkpath*Nkmesh) then
               dkfrac(Nkpath*Nkmesh,:) = dkfrac(0,:)
               write(*, '(I3,3F6.2)') ik, (dkfrac(ik,k),k=1,3)
            end if
         end do
      end do

      write(*,*)
      write(*,*) 'loading kpath for bandplot is done ...'

      return
      end

c###################################################################
