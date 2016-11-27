c###########################################################
c##   Lindhard90
c###########################################################

      program Lindhard90

      use common
      implicit none

      Xmode = 'auto'    ! Manual / Auto

c      call set_hopping()
      call deallocation()
      call allocation()

      call Router()
      do while (request /= 'q')
 10      write(*,*)
         call showMenu()                  ! Show Start Menu

         read (5,'(A1)',err=10) request   ! Read Input File
         call toLowerCase()        ! Convert 'A' to 'a'
         call Router()             ! Run Each Menu
      end do

 999  continue
      end

c###########################################################

c###########################################################
c##   Router:
c###########################################################

      subroutine Router()

      use common, only : request, fnameHop
      implicit none

      integer :: ierr

      if (request == 'p') then
         call setParameter()
      else if (request == 'c') then
         call setConfig()
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
         !call startMeanField(ierr)
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

c###########################################################

c###########################################################
c##   setParameter:
c###########################################################

      subroutine setParameter()

      use common, only : U, J, kT
      implicit none

      write(*,*) '   U = ',REAL(U)
      write(*,*) 'Enter : U'
      read(5,*,err=110) U
      write(*,*) '=> U =', U

      J = U * 0.2d0
110   write(*,*) '   J = ',REAL(J)
      write(*,*) 'Enter : J'
      read(5,*,err=120) J
      write(*,*) '=> J =', J

120   write(*,*) '   kT =',REAL(kT)
      write(*,*) 'Enter : kT'
      read(5,*,err=130) kT
130   write(*,*) '=> kT =', kT

      return
      end

c###########################################################

c###########################################################
c##   setConfig: **Nredx関連がよくわからない　＋　Nkx=Nkyでいいのか？
c###########################################################

      subroutine setConfig()

      use common, only : Nkx, Nky, Nqx, Nqy, Nsite, Nredx, Nredy, Nband,
     &                  Dne
      implicit none

      write(*,*) 'Enter : Nkx', Nkx
      read(5,*,err=999) Nkx
      if (MOD(Nkx,4) /= 0) then
         Nkx = Nkx - MOD(Nkx,4)
         write(*,*) 'Nkx is set to', Nkx
      end if
      write(*,*) '=> Nkx =', Nkx

      Nky = Nkx

      write(*,*) 'Enter : Nqx ', Nqx
      read(5,*,err=999) Nqx
      write(*,*) '=> Nqx =', Nqx

      if (Nqx <= 3) Nqy = 1
      if (Nqx == 4) Nqy = 4

      Nredx = 4 / Nqx
      Nredy = 4 / Nqy

      write(*,*) 'Enter : Nsite ', Nsite
      read(5,*,err=999) Nsite
      write(*,*) '=> Nsite =', Nsite

      write(*,*) 'Enter #electrons/sites (<=', 2 * Nband, ')'
      write(*,*) 'Current setting: Dne=', Dne
      read(5,*,err=999) Dne
      do while ((Dne < 0.0d0).or.(Dne > DBLE(Nband)*2.0d0))
         write(*,*) 'Enter 0 < n < Nband*2'
         read(5,*,err=999) Dne
      end do
      write(*,*) '=> #electrons Dne =', Dne

      end

c###########################################################

c###########################################################
c##   display:
c###########################################################

      subroutine display()

      use common, only : Nkx, Nky, Nband, Dne, U, J, kT,
     &     Nqx, Erange, maxOmega, Dnuu, Dens
      implicit none

      write(*,*) '       ========= Current Setting ========'
      write(*,*) '  #kx=', Nkx, '  #ky=', Nky, '  #electrons=',REAL(Dne)
      write(*,*) '  U=',REAL(U), '  J=',REAL(J)
      write(*,*) '  kT = ', REAL(kT), '  Nqx =', Nqx
      !write(*,'(a,x,5f10.3)') 'op n(XYZ)=',Dnuu(1:Nband)
      !write(*,'(a,x,5f10.3)') 'dens1(XYZ)=',Dens(1:Nbnad,1)
      !write(*,'(a,x,5f10.3)') 'dens2(XYZ)=',Dens(1:Nband,2)
      !write(*,'(a,x,f12.8)') 'Fen=', Fen
      write(*,'(a,x,f12.6,a,x,i5)')'  Erange=', Erange, '  MaxOmega=', maxOmega

      return
      end

c###########################################################

c###########################################################
c##   setIteration:
c###########################################################

      subroutine setIteration()

      use common, only : maxIter, conv
      implicit none

      write(*,*) 'Enter #iterations (>1)'
      read (5,*,err=300) maxIter
      maxIter = MAX(maxIter,1)
      write(*,*) 'Enter convergence level (>1.0e-15)'
      read (5,*,err=300) conv
      conv = MAX(conv,1.0e-15)

300   write(*,*) 'max. #iterations=',maxIter,'  conv. level=',conv

      return
      end

c###########################################################

c###########################################################
c##   setEnergyRange:
c###########################################################

      subroutine setEnergyRange()

      use common, only: Erange, maxOmega, minOmega, initOmega,eta

      write(*,*) ' Enter energy range? ', REAL(Erange)
      read (5,*,err=400) Erange
      write(*,*) ' Enter #slices? ', maxOmega
      read (5,*,err=400) maxOmega

      eta = erng / DBLE(maxOmega) * 2.0d0
      write(*,*) ' Broadening? ', REAL(eta)
      read (5,*,err=400) eta
      write(*,*) ' initial energy? ', REAL(initOmega)
      read (5,*,err=400) initOmega
      write(*,*) ' initial step of omega? ', minOmega
      read (5,*,err=400) minOmega

c      write(*,*) ' Enter energy range? ', REAL(Erng)
c      read (5,*,err=400) Erng
c      write(*,*) ' Enter #slices? ', Maxom
c      read (5,*,err=400) Maxom
c      write(*,*) ' Deta? ', REAL(Deta)
c      read (5,*,err=400) Deta
c      write(*,*) ' initom? ', REAL(Ominit)
c      read (5,*,err=400) Ominit
c      write(*,*) ' initial step of omega? ', Minom
c      read (5,*,err=400) Minom

400   write(*,*) 'MaxEnergy = ', REAL(Erange)
      write(*,*) 'Ominit = ',REAL(initOmega)
      write(*,*) 'MaxOmega = ', maxOmega
      write(*,*) 'eta = ', REAL(eta)

      return
      end

c###########################################################

c###########################################################
c##   startMeanField: **fnameInitとfnameOutが未定義, saveData()つくってない
c###########################################################

      subroutine startMeanField(ierr)

      use common, only: request, fnameInit, fnameOut, Dne, Dne1
      implicit none

      integer, intent(out) :: ierr
      real*8 :: diff
      integer :: nk0

      call initConfig(ierr)

      if(ierr /= 1) then        !ierr=1 : Input Error

         write(*,*) '---------------------------------'
         if(ierr == 0) then
            write(*,*)'New config is set.'
         else if(ierr == -1) then
            write(*,*)'Old config is used.'
         end if
         write(*,*) '---------------------------------'
         call saveData(fnameInit)
         write(*,*) '---------------------------------'

         if (request /= '1') then
            call calcMeanField()
            write(*,*) ''
            call display
            write(*,*) '---------------------------------'
            call saveData(fnameOut)
            write(*,*) '---------------------------------'
            write(*,*) 'saving the eigenenergies...'
            call outEigenvalue()
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
         write(*,*)'|ne1-ne|=',diff
      end if
c## check]

 999  return
      end
c####################################################################
