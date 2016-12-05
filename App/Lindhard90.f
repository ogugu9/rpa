c###################################################################
c##   Lindhard90
c###################################################################

      program Lindhard90

      use common
      implicit none

      Xmode = 'auto'    ! Manual / Auto

      fnameHop = 'hop.in'
      fnameInit = 'init.out'
      fnameOut = 'lin90.out'
      Zi = DCMPLX(0.0,1.0)
      Pi = ACOS(-1.0d0)

c      call set_hopping()

      call showMenu()                  ! Show Start Menu
      do while (request /= 'q')
 10      if ((request.ne.' ').and.(request.ne.'!')) then
            write(*,*)
         endif
         read (5,'(A1)',err=10) request   ! Read Input File
         call toLowerCase()        ! Convert 'A' to 'a'
         call Router()             ! Run Each Menu
      end do

 999  continue
      end

c###################################################################

c###################################################################
c##   Router:
c###################################################################

      subroutine Router()

      use common, only : request, fnameHop
      implicit none

      integer :: ierr

      if (request == 'p') then
         call setParameter()
      else if (request == 'c') then
         call setConfig()
      else if (request == 'b') then
         call loadKpath()
      else if (request == '1') then
         !call initConfig()
      else if (request == 'd') then
         call display()
      else if (request == 'e') then
         call setEnergyRange()
      else if (request == 'i') then
         call setIteration()
      else if (request == 'h') then
         write(*,*) 'Reading hopping from:',fnameHop
         write(*,*) 'Do you want to change?'
         read(5,*,err=999) request
         if (request == 'y') then
            call getFilename('filename? ',fnameHop)
            call loadHopping(fnameHop)
         end if
      else if (request == 'l') then
         !call loadData()
      else if (request == 'w') then
         !call saveData()
      else if (request == 'm') then
         call startMeanField(ierr)
      else if (request == 'o') then
         !** メソッド名変えたい
         !call orbital()
      else if (request == 'f') then
         !call drawFermiSurface()
      else if (request == 'k') then
         !Optical Conductivity
         !call startOptics()
      else if (request == 'b') then
         !Interband Excitation
         !call startInterband()
      else if (request == 's') then
         !Susceptibility
         !call startSuscepts()
      else if (request == 'r') then
         !call calcResistivity()
      else if (request == 'v') then
         !call calcFermiVelocity()
      else if (request == 'x') then
         !call startRIXS()
      end if

 999  return
      end

c###################################################################

c###################################################################
c##   setParameter:
c###################################################################

      subroutine setParameter()

      use common, only : U, J, kT
      implicit none

      write(*,*) 'setting parameter ...'

      read(5,*,err=110) U
      print "(' => U  = ', F5.3)", U

      J = U * 0.2d0
110   read(5,*,err=120) J
      print "(' => J  = ', F5.3)", J

120   read(5,*,err=130) kT
130   print "(' => kT = ', F5.3)", kT

      write(*,*) 'setting parameter is done ...'

      return
      end

c###################################################################

c###################################################################
c##   setConfig: **Nredx関連がよくわからない　＋　Nkx=Nkyでいいのか？
c###################################################################

      subroutine setConfig()

      use common, only : Nkx, Nky, Nqx, Nqy, Nsite,
     &                  Nredx, Nredy, Nband, Dne, Dnuu, Dndd
      implicit none

      write(*,*) 'setting config ...'

      read(5,*,err=999) Nkx
      if (MOD(Nkx,4) /= 0) then
         Nkx = Nkx - MOD(Nkx,4)
         write(*,*) 'Nkx is set to', Nkx
      end if
      print "(' => Nkx   = ', I5)", Nkx

      Nky = Nkx

      read(5,*,err=999) Nqx
      print "(' => Nqx   = ', I5)", Nqx

      if (Nqx <= 3) Nqy = 1
      if (Nqx == 4) Nqy = 4

      Nredx = 4 / Nqx
      Nredy = 4 / Nqy

      read(5,*,err=999) Nsite
      print "(' => Nsite = ', I5)", Nsite

      read(5,*,err=999) Nband
      print "(' => Nband = ', I5)", Nband

      call deallocation()
      call allocation()
      Dnuu(:) = 0.0d0
      Dndd(:) = 0.0d0

      print "(' #electrons/sites <=',I3)", 2 * Nband
      read(5,*,err=999) Dne
      do while ((Dne < 0.0d0).or.(Dne > DBLE(Nband)*2.0d0))
         write(*,*) 'Enter 0 < n < Nband*2'
         read(5,*,err=999) Dne
      end do
      print "(' => Dne   = ', F5.2)", Dne

      write(*,*) 'setting config done ...'

999   return
      end

c###################################################################

c###################################################################
c##   display:
c###################################################################

      subroutine display()

      use common, only : Nkx, Nky, Nqx, Nqy, Nband,
     &     Dne, U, J, kT, Erange, maxOmega, Dnuu, Dens
      implicit none

      write(*,*) '======== Current Setting ========'
      print "('  Nkx =',I4,'  Nky =',I4)", Nkx, Nky
      print "('  Nqx =',I4,'  Nqy =',I4)", Nqx, Nqy
      print "('  #electrons = ',F5.3)", Dne
      print "('  U =',F5.3,'  J =',F5.3,'  kT =',F5.3)", U, J, kT
      !write(*,'(a,x,5f10.3)') 'op n(XYZ)=',Dnuu(1:Nband)
      !write(*,'(a,x,5f10.3)') 'dens1(XYZ)=',Dens(1:Nbnad,1)
      !write(*,'(a,x,5f10.3)') 'dens2(XYZ)=',Dens(1:Nband,2)
      !write(*,'(a,x,f12.8)') 'Fen=', Fen
      print "('  Erange = ',F5.3,'  MaxOmega =',I5)", Erange, MaxOmega

      return
      end

c###################################################################

c###################################################################
c##   setIteration:
c###################################################################

      subroutine setIteration()

      use common, only : maxIter, conv
      implicit none

c      write(*,*) 'Enter #iterations (>1)'
      read (5,*,err=300) maxIter
      maxIter = MAX(maxIter,1)
c      write(*,*) 'Enter convergence level (>1.0e-15)'
      read (5,*,err=300) conv
      conv = MAX(conv,1.0e-15)

      write(*,*) 'setting iterations ...'
300   print "(' => maxIter = ', I5)", maxIter
      print "(' => conv    = ', ES8.1)", conv
      write(*,*) 'setting iterations is done ...'

      return
      end

c###################################################################

c###################################################################
c##   setEnergyRange:
c###################################################################

      subroutine setEnergyRange()

      use common, only: Erange, maxOmega, minOmega, initOmega,eta

      write(*,*) 'setting energy-range ...'
c      write(*,*) ' Enter energy range? ', REAL(Erange)
      read (5,*,err=400) Erange
c      write(*,*) ' Enter #slices? ', maxOmega
      read (5,*,err=400) maxOmega

      eta = erng / DBLE(maxOmega) * 2.0d0
c      write(*,*) ' Broadening? ', REAL(eta)
      read (5,*,err=400) eta
c      write(*,*) ' initial energy? ', REAL(initOmega)
      read (5,*,err=400) initOmega
c      write(*,*) ' initial step of omega? ', minOmega
      read (5,*,err=400) minOmega

400   print "(' => MaxEnergy  = ', F5.3)", Erange
      print "(' => initOmega  = ', F5.3)", initOmega
      print "(' => minOmega   = ', I5)", minOmega
      print "(' => maxOmega   = ', I5)", maxOmega
      print "(' => eta (blur) = ', F5.3)", eta

      write(*,*) 'setting energy-range is done ...'

      return
      end

c###################################################################

c###################################################################
c##   startMeanField: **saveData()
c###################################################################

      subroutine startMeanField(ierr)

      use common, only: request, fnameInit, fnameOut, Dne, Dne1,
     &         Nkx, Nky
      implicit none

      integer, intent(out) :: ierr
      real*8 :: diff
      integer :: nk0

      call initConfig(ierr)

      if(ierr /= 1) then        !ierr=1 : Input Error

         write(*,*)
         if(ierr == 0) then
            write(*,*)' New config is set.'
         else if(ierr == -1) then
            write(*,*)' Old config is used.'
         end if
         write(*,*)
         !call saveData(fnameInit)
         write(*,*)

         if (request /= '1') then
            call calcMeanField()
            write(*,*) ''
            call display
            write(*,*) '---------------------------------'
            !call saveData(fnameOut)
            write(*,*) '---------------------------------'
            write(*,*) 'saving the eigenenergies...'
            call calcBandplot(Nkx,Nky,1)
            call calcBandplot(Nkx,Nky,2)
            write(*,*) '---------------------------------'
            write(*,*) 'calculating eigenvalue for bandplot is done...'
            call outEigenvalue()
            write(*,*) '---------------------------------'
            write(*,*) 'outputing eigenvalue for bandplot is done...'
            nk0 = 100
            write(*,*) '---------------------------------'
            write(*,*) 'calculating DOS...'
            !call calcDOS3(nk0*2)! 状態密度の計算 > optcond3.f
         end if
      end if

c## check[
      diff = ABS(Dne1 - Dne)
      if(diff >= 1.0d-8)then
         write(*,*) '---------------------------------'
         write(*,*)'### calculation may be wrong ##'
         write(*,*)'Dne=',Dne
         write(*,*)'Dne1=',Dne1
         write(*,*)'|Dne1-Dne|=',diff
      end if
c## check]

 999  return
      end
c####################################################################
