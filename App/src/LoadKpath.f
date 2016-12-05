c###################################################################
c##   loadKpath:
c###################################################################

      subroutine loadKpath()

      use common, only : request, Nkmesh, Nkpath,
     &      kpointList, kpathList
      implicit none

      integer :: i, j, k, nkpt
      real*8 :: kvecIni(1:3), kvecFin(1:3)
      character*128 :: chara

      write(*,*) 'loading kpath for bandplot ...'

      read(5,*,err=110) Nkmesh
      print "(' => Nkmesh  = ', I5)", Nkmesh
      read(5,*,err=110) Nkpath
      print "(' => Nkpath  = ', I5)", Nkpath

      call deallocKpath()
      call allocKpath()

110   write(*,*)
      write(*,*) '#k-point path for bandplot'
      read(5,*) chara
      do i = 1,Nkpath
         read(5, '(A4,3F6.2,A4,3F6.2)'),
     &      kpathList(i)%iniName, kpathList(i)%iniPosition,
     &      kpathList(i)%finName, kpathList(i)%finPosition
         write(*, '(A4,3F6.2,3X,A4,3F6.2)'),
     &      kpathList(i)%iniName, kpathList(i)%iniPosition,
     &      kpathList(i)%finName, kpathList(i)%finPosition
      end do
      read(5,*) chara ! chara = ## KPOINT PATH ##
      write(*,*)

200   write(*,*) '#k-point list for bandplot'
      nkpt = 0
      do i = 1, Nkpath
         kvecIni(:) = kpathList(i)%iniPosition
         kvecFin(:) = kpathList(i)%finPosition
         do j = 0,Nkmesh-1
            kpointList(nkpt,:) = kvecIni(:)
     &      + (kvecFin(:)-kvecIni(:)) * DBLE(j)/DBLE(Nkmesh)
            write(*, '(I3,3F6.2)') nkpt, (kpointList(nkpt,k),k=1,3)
            nkpt = nkpt + 1
            if (nkpt == Nkpath * Nkmesh) then
               kpointList(Nkpath*Nkmesh,:) = kpointList(0,:)
               write(*, '(I3,3F6.2)') nkpt, (kpointList(nkpt,k),k=1,3)
            end if
         end do
      end do

      write(*,*)
      write(*,*) 'loading kpath for bandplot is done ...'

      return
      end

c###################################################################
