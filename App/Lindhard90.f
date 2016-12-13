c###################################################################
c##   Lindhard90
c###################################################################

      program Lindhard90

      use common
      implicit none

      Xmode = 'auto'    ! Manual / Auto

      fnameHop = 'ato.hop'
      fnameInit = 'init.out'
      fnameOut = 'lin90.out'
      Zi = DCMPLX(0.0,1.0)
      Pi = ACOS(-1.0d0)

c      call set_hopping()

      call showMenu()                  ! Show Start Menu
      do while (req /= 'q')
 10      if ((req.ne.' ').and.(req.ne.'!')) then
            write(*,*)
         endif
         read (5,'(A1)',err=10) req   ! Read Input File
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

      use common, only : req, fnameHop
      implicit none

      integer :: ierr

      if (req == 'p') then
         call setParameter()
      else if (req == 'c') then
         call setConfig()
      else if (req == 'b') then
         call loadKpath()
      else if (req == '1') then
         !call initConfig()
      else if (req == 'd') then
         call display()
      else if (req == 'e') then
         call setEnergyRange()
      else if (req == 'i') then
         call setIteration()
      else if (req == 'h') then
         write(*,*) 'Reading hopping from:',fnameHop
         write(*,*) 'Do you want to change?'
         read(5,*,err=999) req
         if (req == 'y') then
            call getFilename('filename? ',fnameHop)
            call loadHopping(fnameHop)
         end if
      else if (req == 'l') then
         !call loadData()
      else if (req == 'w') then
         !call saveData()
      else if (req == 'm') then
         call startMeanField(ierr)
      else if (req == 'o') then
         !** メソッド名変えたい
         !call orbital()
      else if (req == 'f') then
         !call drawFermiSurface()
      else if (req == 'k') then
         !Optical Conductivity
         !call startOptics()
      else if (req == 'b') then
         !Interband Excitation
         !call startInterband()
      else if (req == 's') then
         !Susceptibility
         !call startSuscepts()
      else if (req == 'r') then
         !call calcResistivity()
      else if (req == 'v') then
         !call calcFermiVelocity()
      else if (req == 'x') then
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

      use common, only : Nkx, Nky, Nkz, Nqx, Nqy, Nqz, Nsite,
     &                  Nredx, Nredy, Nredz, Nband, Dne, Dnuu, Dndd,
     &                  EF, la, lb, lc, recipLat
      implicit none
      integer :: i, j
      character*1 :: chara

      write(*,*) 'setting config ...'

      read(5,*,err=999) Nkx
      if (MOD(Nkx,4) /= 0) then
         Nkx = Nkx - MOD(Nkx,4)
         write(*,*) 'Nkx is set to', Nkx
      end if
      print "(' => Nkx   = ', I5)", Nkx

      read(5,*,err=999) Nky
      if (MOD(Nky,4) /= 0) then
         Nky = Nky - MOD(Nky,4)
         write(*,*) 'Nky is set to', Nky
      end if
      print "(' => Nky   = ', I5)", Nky

      read(5,*,err=999) Nkz
      if (MOD(Nkz,4) /= 0 .and. Nkz /= 1) then
         Nkz = Nkz - MOD(Nkz,4)
         write(*,*) 'Nkz is set to', Nkz
      end if
      print "(' => Nkz   = ', I5)", Nkz

      ! ** z成分のオーダリングベクトルはどうするか
      read(5,*,err=999) Nqx
      print "(' => Nqx   = ', I5)", Nqx

      if (Nqx <= 3) then
         Nqy = 1
         Nqz = 1
      elseif (Nqx == 4) then
         Nqy = 4
         Nqz = 4
      endif

      Nredx = 4 / Nqx
      Nredy = 4 / Nqy
      Nredz = 4 / Nqz

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

      read(5,*,err=999) EF
      print "(' => E_F = ', F6.3)", EF

      write(*,*)
      write(*,*) '#lattice constant'
      read(5,*,err=999) la
      read(5,*,err=999) lb
      read(5,*,err=999) lc
      print "(' => la, lb, lc = ', 3F6.3)", la, lb, lc

      read(5,*,err=999) chara
      write(*,*)
      write(*,*) '#reciprocal lattice'
      do i = 1,3
         read(5,'(A5,3F11.6)',err=999) chara, recipLat(i,:)
         print "(' => b_', I1, ' = ', 3F7.3)", i, recipLat(i,:)
      end do
      read(5,*,err=999) chara

      write(*,*) 'setting config done ...'

999   return
      end

c###################################################################

c###################################################################
c##   display:
c###################################################################

      subroutine display()

      use common, only : Nkx, Nky, Nkz, Nqx, Nqy, Nqz, Nband,
     &     Dne, U, J, kT, Erange, maxOmega, Dnuu, Dens
      implicit none

      write(*,*) '======== Current Setting ========'
      print "('  Nkx =',I4,'  Nky =',I4,'  Nkz =',I4)", Nkx, Nky, Nkz
      print "('  Nqx =',I4,'  Nqy =',I4,'  Nqz =',I4)", Nqx, Nqy, Nqz
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

      use common, only: req, fnameInit, fnameOut, Dne, Dne1,
     &         Nkx, Nky, Nkz
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

         if (req /= '1') then
            call calcMeanField()
            write(*,*) ''
            call display
            write(*,*) '---------------------------------'
            !call saveData(fnameOut)
            write(*,*) '---------------------------------'
            write(*,*) 'saving the eigenenergies...'
            call calcBandplot(Nkx,Nky,Nkz,1)
            call calcBandplot(Nkx,Nky,Nkz,2)
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
