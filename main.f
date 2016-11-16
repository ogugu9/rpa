c###################################################
c##     Fe-based superconductor
c##
c##                        E. Kaneshita 2009
c###################################################

      program Fe_based_superconductor

      use common
      implicit none

      character :: what*1
      integer :: idummy
      real*8 :: tpd, D, alpha, beta

c## file names[
      fname(2) = 'out.dat        ' !ouput parameters
      fname(3) = 'init.dat       ' !ouput initial parameters
      fname(4) = 'in.hop1111xyz.dat' !input hopping integral
      fname(5) = 'sti.dat        ' !ouput initial spin config
      fname(6) = 'out.green.dat  '
      fname(7) = 'out.greenk.dat '
      fname(8) = 'out.chi0.dat   '
c## file nemes]

c## constants[
      Zi = DCMPLX(0.0,1.0)
      Pi = ACOS(-1.0d0)
c## constants]

      Xmode = 'manual'
      Xmode = 'auto'

      Check = 'y'
      Erng = 4.0d0              !0.4d0
      Maxom = 100               !20
      Minom = 0
      Deta =  0.04d0                  !0.001d0
      Ominit = 0.0d0

      Nb = 5

      Maxit = 10000               ! number of iterations
      Conv = 1.0d-14

c==== Dtemp = k_B T [eV]
      Dtemp = 0.010d0           ! temperature ( [ k_B / t ] )
      Dtemp = 0.0d0
c      Deltae = 1.0d-14

      Ns = llorb                ! number of states (orbitals)
      Ns2 = Ns * 2

      Nkx = 40
      Nky = 40

      Nqx = 4
      Nqy = 4

      Nredx = 4 / Nqx
      Nredy = 4 / Nqy

      U = 1.6d0
      DJ = 0.32d0
      DJ = U * 0.2d0

c      U = 1.2d0 * 2.0d0
c      DJ = 0.25d0 * 2.0d0

c      Dne = 6.0d0 + 0.11d0
      Dne = 6.0d0
      
      what = 'd'
      if (Xmode(1:1) /= 'm') what = '0'

         

      call set_hopping()!ホッピングパラメータの読み込み > mfcalc2.f

      call deallocation()
      call allocation()
      allocate( Dnuu(Ns) )
      allocate( Dndd(Ns) )
      allocate( Dens(Ns,2) )
      allocate( Dnuuxyz(Ns) )
      allocate( Densxyz(Ns,2) )
      Dnuu(:) = 0.0d0
      Dndd(:) = 0.0d0

      call menu(what)
      do while (what /= 'q')
c## main menu[
 10      write(*,*)

         if (Xmode(1:1) == 'm') then
            write(*,*) '             -- MAIN MENU -- '
            write(*,*) '           ENTER THE CHARACTER'
            write(*,*)
            write(*,*) 'SETUP: (P)aram  (C)onfig  (I)teration'
            write(*,*) '     : (1)st-config  (D)isplay'
            write(*,*) '     : (E)nergy-range (#)system size'
            write(*,*) '     : (2)phase factor'
            write(*,*)
            write(*,*) ' FILE: (L)oad   (W)rite'
            write(*,*)
            write(*,*) ' CALC: (M)F  (S)uscept  inter(B)and'
            write(*,*) '     : (F)ermi (V)elocity'
            write(*,*) '     : (O)bital'
            write(*,*) '     : (Q)uit'
         else
            write(*,*) '---------------------------------'
         end if
         
         read (5,'(A1)',err=10) what
         
         if (what == 'Q') then
            what='q'            !quit
         else if (what == 'B') then
            what='b'            !interband excitation
         else if (what == 'C') then
            what='c'            !config
         else if (what == 'D') then
            what='d'            !display
         else if (what == 'E') then
            what='e'            !energy
         else if (what == 'F') then
            what='f'            !fermi
         else if (what == 'G') then
            what='g'            !green
         else if (what == 'H') then
            what='h'            !hopping integral
         else if (what == 'I') then
            what='i'            !iteration
         else if (what == 'K') then
            what = 'k'          !kougaku dendoudo
         else if (what == 'L') then
            what='l'            !loard
         else if (what == 'M') then
            what='m'            !meanfield
         else if (what == 'O') then
            what='o'            !orbital
         else if (what == 'P') then
            what='p'            !parameters
         else if (what == 'R') then
            what='r'            !resistivity
         else if (what == 'S') then
            what='s'            !suscept
         else if (what == 'V') then
            what='v'            !meanfield
         else if (what == 'W') then
            what='w'            !write
         else if (what == 'X') then
            what='x'            !RIXS L-edge
         end if

c## main menu]
         call menu(what)
         
         if (what == 'c') then
            call deallocation()
            call allocation()
         end if

      end do

 999  continue

       end


c####################################################################


c#########################################################
c#### menu:
c#########################################################

      subroutine menu(what)
      use common

      implicit none

      character, intent(inout) :: what*1

      character :: chara*15, charalong*40, filename*80
      integer :: ierr, iflag, mode
      integer :: idummy, i, iom, ii, ikmd
      real*8 :: dummy, theta, phi



c===========================================[
c== parameter ==============================
c==
      if (what == 'p') then
         write(*,*) '  U = ',REAL(U)
         write(*,*) 'Enter : U'
         read(5,*,err=110) U
         write(*,*) '===>   U =', U
         DJ = U * 0.2d0         !aaaaa
 110     write(*,*) '  DJ = ',REAL(DJ)
         write(*,*) 'Enter : DJ'
         read(5,*,err=120) DJ
         write(*,*) '===>   DJ =', DJ
 120     write(*,*) ' Dtemp=',REAL(Dtemp)
         write(*,*) 'Enter : Dtemp'
         read(5,*,err=130) Dtemp
 130     write(*,*) '===>   Dtemp =', Dtemp
         
c===========================================]

c===========================================[
c== config =================================
c==
      else if (what == 'c') then
         write(*,*) 'Enter #kx ', Nkx
         read(5,*,err=999) Nkx
         if (MOD(Nkx,4) /= 0) then
            Nkx = Nkx - MOD(Nkx,4)
            write(*,*) 'Nkx is set to', Nkx
         end if
         write(*,*) '===>   Nkx =', Nkx

         Nky = Nkx
c         write(*,*) 'Enter #ky ', Nky
c         read(5,*,err=999) Nky
c         if (MOD(Nky,4) /= 0) then
c            Nky = Nky - MOD(Nky,4)
c            write(*,*) 'Nky is set to', Nky
c         end if

         write(*,*) 'Enter Nqx ', Nqx
         read(5,*,err=999) Nqx
         write(*,*) '===>   Nqx =', Nqx
c         write(*,*) 'Enter Nqy ', Nqy
c         read(5,*,err=999) Nqy
         if (Nqx <= 3) Nqy = 1
         if (Nqx == 4) Nqy = 4
         
         Nredx = 4 / Nqx
         Nredy = 4 / Nqy
         
         write(*,*) 'Enter #electrons per site (<=', 2 * Ns, ')'
         write(*,*) 'Current setting: Dne=', Dne
         read(5,*,err=999) Dne
         do while ((Dne < 0.0d0).or.(Dne > 10.0d0))
            write(*,*) 'Enter 0 < n < 10'
            read(5,*,err=999) Dne
         end do
         
         write(*,*) '===>   #electrons =', Dne
      
c===========================================]

c===========================================[
c== display ================================
c==
      else if (what == 'd') then
         call display
c===========================================]

c===========================================[
c== energy-range ===========================
c==
      else if (what == 'e') then

         call menu_Eparam(Erng, Maxom, Deta, Ominit, Minom)
c         write(*,*) ' Enter energy range? ', REAL(Erng)
c         read (5,*,err=400) Erng
c         write(*,*) ' Enter #slices? ', Maxom
c         read (5,*,err=400) Maxom
c         write(*,*) ' Deta? ', REAL(Deta)
c         read (5,*,err=400) Deta
c         write(*,*) ' initom? ', REAL(Ominit)
c         read (5,*,err=400) Ominit
c         write(*,*) ' initial step of omega? ', Minom
c         read (5,*,err=400) Minom
         
 400     write(*,*) 'max energy',REAL(Erng),'  #slices',Maxom
         write(*,*) 'Deta=',REAL(Deta), 'initial energy',REAL(ominit)
c===========================================]


c===========================================[
c== fermi ================================
c==

      else if (what == 'f') then
         write(*,*)'drawing Fermi surface...'
c         call drawfermi()
         call get_param_int('size?',idummy)
         chara = 'v'
         call get_param_ch('kF only or vFw/kF (k or v)? ',chara)
         call fermi_KfandVf(0.0d0, i, idummy, idummy, chara(1:1))
c===========================================]

c===========================================[
c== hopping ===============================
c==
      else if (what == 'h') then
         write(*,*) 'hopping data read from:',fname(4)
         write(*,*) 'Do you want to change?'
         read(5,*,err=999) what
         if (what == 'y') then
            call get_param_ch('filename? ',fname(4))
            call set_hopping()
         end if
         what = 'h'
c===========================================]

c===========================================[
c== iteration ==============================
c==
      else if (what == 'i') then
         write(*,*) 'Enter #iterations (>1)'
         read (5,*,err=300) Maxit
         Maxit = MAX(Maxit,1)
         write(*,*) 'Enter convergence level (>1.0e-15)'
         read (5,*,err=300) Conv
         Conv = MAX(Conv,1.0e-15)
         
 300     write(*,*) 'max. #iterations=',Maxit,'  conv. level=',Conv
c===========================================]

c===========================================[
c== Optical conductivity ===================
c==

      else if (what == 'k') then
         mode = 0
         do while (mode == 0)
            write(*,*)'optical conductivity: mode?'
            write(*,*)'1: interband  2: Drude'
            write(*,*)'For Drude, vf is required'
            read(5,*,err=999) mode
            write(*,*) mode
         end do
         chara = '?'
         do while ( chara(1:1) == '?')
            write(*,*)'polarization?  (x, y, a, b)'
            read(5,*,err=999) chara
            write(*,*) chara(1:1)
         end do
c         idummy = 0
c         do while ( idummy == 0)
c            write(*,*)'number of energy step? '
c            read(5,*,err=999) idummy
c            write(*,*) idummy
c            idummy = ABS(idummy)
c            if (idummy > 10000) idummy = 10000
c         end do
         if (mode == 1) then
            i = 0
            do while ( i == 0)
               write(*,*)'size? '
               read(5,*,err=999) i
               write(*,*) i
               i = ABS(i)
               if (i > 100000) i = 100000
            end do
            
            write(*,*)'calculating optical conductivity...'
            call optcond4(chara(1:1),i)
         else
            call Drude4(chara(1:),dummy,'#')
         end if
c===========================================]

c===========================================[
c== load ===================================
c==
      else if (what == 'l') then
         write(*,*)'Enter the filename to load (<40 characters)'
         read(5,*,err=999) charalong
         call loaddata(charalong,ierr)
c===========================================]

c===========================================[
c== meanfield ==============================
c==
      else if (what == 'm') then
         call work(what,ierr)
c===========================================]

c===========================================[
c== initial conf ===========================
c==
      else if (what == '1') then
         call work(what,ierr)
c===========================================]

c===========================================[
c== write ===================================
c==
      else if (what == 'w') then
         write(*,*)'Enter filename to save'
         read(5,*,err=999) charalong
         call savedata(charalong)
c===========================================]

c===========================================[
c== interband excitation ===================
c==
      else if (what == 'b') then
         call menu_interband(what)
c===========================================]

c===========================================[
c== orbital ================================
c==
      else if (what == 'o') then
         call orbital()
c===========================================]


c===========================================[
c== resistivity ============================
c==

      else if (what == 'r') then
         chara = '?'
         do while ( chara(1:1) == '?')
            write(*,*)'polarization?  (x, y, a, b)'
            read(5,*,err=999) chara
            write(*,*) chara(1:1)
         end do
         write(*,*)'calculating rho...'
         call resistivity(chara(1:1),dummy,'#')
         write(*,'(" rho_",a1,"= ",g20.10)')chara(1:1),dummy 
c===========================================]

c===========================================[
c== suscept ================================
c==
      else if (what == 's') then
         call menu_suscept(what)
c===========================================]

c===========================================[
c== fermi velocity =========================
c==

      else if (what == 'v') then
         write(*,*)'calculating Fermi velocity...'
         call fermi_velocity2()
c===========================================]

c===========================================[
c== RIXS ===================================
c==

      else if (what == 'x') then
         call menu_RIXS(what)
      end if

 999  return


 2000 format(32A1)
      end

c####################################################################


c#########################################################
c#### work:
c#########################################################

      subroutine work(what,ierr)
      use common

      implicit none

      character, intent(in) :: what*1
      integer, intent(out) :: ierr
      character :: chara*40
      integer :: is, nk0, np
      real*8 :: dummy, ene

      
      call inicon(ierr)
      
      if(ierr /= 1) then        !ierr=1 : input error
         
         if(ierr == 0) then 
            write(*,*)'set the new config.' !ierr=-1: new configuration
         else if(ierr == -1) then
            write(*,*)'using the old config.'
         end if
            
c## save the initial state[
         chara=fname(3)
c     if(ierr == -1) chara='itemp.dat'
         call savedata(chara)
c## save the initial state]
         
         if (what /= '1') then
c            write(*,*) 'Enter nk0', nk0
c            read(5,*,err=999) nk0
            call mfcalc()
            call display
            chara=fname(2)
            call savedata(chara)
            write(*,*) '---------------------------------'
            write(*,*) 'saving the eigenenergies...'
            call eigenout()


            write(*,*) '---------------------------------'
            write(*,*) 'calculating fermi surface...'
            chara = 'out.fermi.dat'
            ene = -0.0d0
            call fermi2(chara,ene,np)
            chara = 'out.fs.dat'
            call fermi_seek_grid(chara,ene,np)
            if (what == 's') return
            nk0 = 100
c            call fermi3(chara,nk0,np)
            write(*,*) '---------------------------------'
            write(*,*) 'calculating DOS...'
            call calcDOS3(nk0*2)! 状態密度の計算 > optcond3.f
      write(*,*)'aaaaaaaaa'
         end if
      end if
         
c## check[
      dummy = ABS(Dne1 - Dne)
      if(dummy >= 1.0d-8)then
         write(*,*) '---------------------------------'
         write(*,*)'### calculation may be wrong ##'
         write(*,*)'Dne=',Dne
         write(*,*)'Dne1=',Dne1
         write(*,*)'|Dne1-Dne|=',dummy
      end if
c## check]

 999  return
      end
c####################################################################



c###########################################################
c##   deallocation:
c###########################################################
      subroutine deallocation()

      use common
      implicit none

      logical :: torf
      
      torf = allocated(Eall)
      if (torf .EQV. .true.) deallocate(Eall)
      torf = allocated(Zpsiall)
      if (torf .EQV. .true.) deallocate(Zpsiall)

      return
      end
c###########################################################

c###########################################################
c##   allocation:
c###########################################################
      subroutine allocation()

      use common
      implicit none

      logical :: torf
      
c      allocate( Eall(0:Nkx,0:Nky,Ns*Nqx*2) )
      allocate( Eall(0:Nkx,0:Nky,Ns*Nqx,2) )
      allocate( Zpsiall(0:Nkx,0:Nky,Ns*Nqx,Ns*Nqx,2) )



      return
      end
c###########################################################

c###########################################################
c##   display:
c###########################################################
      subroutine display()

      use common, only : Nkx, Nky, Ns, Dne, U, DJ, Dtemp, 
     &     Nqx, Fen, Erng, Maxom, Dnuu, Dens, Dnuuxyz, Densxyz
      implicit none

      
      write(*,*) '       ========= Current Setting ========'
      write(*,*) '  #kx=', Nkx, '  #ky=', Nky, '  #electrons=',REAL(Dne)
      write(*,*) '  U=',REAL(U), '  DJ=',REAL(DJ)
      write(*,*) '  Dtemp = ', REAL(Dtemp), '  Nqx =', Nqx 
      write(*,'(a,x,5f10.3)') 'op n(XYZ)=',Dnuu(1:Ns)
      write(*,'(a,x,5f10.3)') 'op n(xyz)=',Dnuuxyz(1:Ns)
      write(*,'(a,x,5f10.3)') 'dens1(XYZ)=',Dens(1:Ns,1)
      write(*,'(a,x,5f10.3)') 'dens2(XYZ)=',Dens(1:Ns,2)
      write(*,'(a,x,5f10.3)') 'dens1(xyz)=',Densxyz(1:Ns,1)
      write(*,'(a,x,5f10.3)') 'dens2(xyz)=',Densxyz(1:Ns,2)
      write(*,'(a,x,f12.8)') 'Fen=', Fen
      write(*,'(a,x,f12.6,a,x,i5)')'  Erng=', Erng, '  Maxom=', Maxom


      return
      end
c###########################################################


c#########################################################
c#### pathget:
c#########################################################

      subroutine pathget(chara)
      use common

      implicit none

      character, intent(inout) :: chara*5
      character :: ch*1
      integer :: i
      
      do i = 1, 5
         if (i == 3) cycle
         ch = chara(i:i)
         call check_int(ch)
         if (ch == 'F') then
            chara = 'sssss'
            exit
         end if
      end do

      return
      end
c###########################################################

c#########################################################
c#### pathget:
c#########################################################
      subroutine check_int(chara)
      use common

      implicit none

      character, intent(inout) :: chara*1
      character :: ch*1
      integer :: min, max
      integer :: i
      

      do i = 0, 9
         write (ch,'(i1)') i
         if (chara == ch) return
      end do

      chara = 'F'

      return
      end
c###########################################################

c#########################################################
c####  drawfermi:
c#########################################################

      subroutine drawfermi()
      use common

      implicit none

      character :: chara*15
      integer :: is, np
      real*8 :: ene, de, dp, rangek,dkx,dky

      integer :: nx1, nx2, ny1, ny2, mkx, mky, ispin

      mkx = 10000
      mky = 1000
      ispin = 1
      dp =  0.1319313d0            !dirac point
      rangek = 0.05d0

c      dp = 0.0d0
c      rangek = 0.09d0

      dkx = dp - INT(dp*mkx) / DBLE(mkx)
      dky = 0.0d0


      nx1 = INT((dp - rangek) * mkx)
      nx2 = INT((dp + rangek) * mkx)
      ny1 = -INT(rangek * mky)
      ny2 =  INT(rangek * mky)

c      nx1 = 0.0d0 * mkx
c      ny1 = 0.39d0 * mky
c      nx2 =  0.125d0 * mkx
c      ny2 = 0.61d0 * mky

      call Seigenk2(nx1, nx2, ny1, ny2,mkx, mky, ispin,dkx,dky)

      return

      de = 0.01d0
      
      chara = 'out.fermi.dat'
      ene = -0.0d0
      call fermi2(chara,ene,np)
      do is = -3, 3
         ene = DBLE(is) * de
         if(is == -3) chara(14:15) = '-3'
         if(is == -2) chara(14:15) = '-2'
         if(is == -1) chara(14:15) = '-1'
         if(is ==  0) chara(14:15) = '+0'
         if(is ==  1) chara(14:15) = '+1'
         if(is ==  2) chara(14:15) = '+2'
         if(is ==  3) chara(14:15) = '+3'
         call fermi2(chara,ene,np)
      end do
      
      

      return
      end
c###########################################################

c#########################################################
c####  Sraman:
c#########################################################

      subroutine Sraman()
      use common

      implicit none

      character :: chara*15
      integer :: is, ne0, np
      real*8 :: dummy, ene

      ne0 = 1000
c      call raman_v(ne0,'aa')
c         write(*,*) 'aaaaaaaaaaaa'
c      call raman_v(ne0,'ab')
c         write(*,*) 'aaaaaaaaaaaa'
c      call raman_v(ne0,'xx')
c      call raman_v(ne0,'yy')
c      call raman_v(ne0,'xy')
c         write(*,*) 'aaaaaaaaaaaa'
      

      return
      end
c###########################################################

c#########################################################
c####  wordcomb:
c#########################################################

      character(LEN=40) function wordcomb(ch1,ch2)
      use common

      implicit none

      character(LEN=*),intent(in)  :: ch1, ch2

      wordcomb = TRIM(ADJUSTL(ch1))//TRIM(ADJUSTL(ch2))

      return
      end
c###########################################################

c#########################################################
c#### menu_interband:
c#########################################################

      subroutine menu_interband(what)
      use common, only : Nkx, Maxom, Erng, Deta

      implicit none

      character, intent(inout) :: what*1
      integer :: iflag, mk
      integer :: idummy
      real*8 :: dummy

c===========================================[
c== interband excitation ===================
c==
  
      do while ( what == 'b')
         write(*,*)'1: optical conductivity  [x: menu]'
         write(*,*)'2: component (data required)'
         write(*,*)'3: both'
         write(*,*)'4: DOS'
         read(5,*) what
         if ((what == '1') .or. (what == 'o')) then
            iflag = 1
         else if ((what == '2') .or. (what == 'c')) then
            iflag = 2
         else if (what == '3') then
            iflag = 3
         else if (what == '4') then
            iflag = 4
         else if (what == 'x') then
            what = 'd'
            return
         else
            what = 'b'
         end if
      end do
      what = 'b'
      
      if ((iflag == 1) .or. (iflag == 3)) then
         write(*,*) "Do you wanna change the enrgy range? (y/n)"
         read(5,*) what
         if (what == 'y') then
            call menu_Eparam(Erng, Maxom, Deta, dummy, idummy)
            write(*,*) 'energy range',-REAL(Erng),'--',REAL(Erng)
            write(*,*) '  #slices',Maxom
            write(*,*) 'Deta=',REAL(Deta)
         end if
         mk = Nkx
         call get_param_int('system size?',mk)
         call get_param_ch('direction(x,y,a,b)?',what)
         call optcond3(what,mk)
         what = 'b'
      end if
      
      if ((iflag == 2) .or. (iflag == 3)) then
         call set_hopping()
         call interband()
      end if

      if (iflag == 4) then
         write(*,*) "Do you wanna change the energy range? (y/n)"
         read(5,*) what
         if (what == 'y') then
            call menu_Eparam(Erng, Maxom, Deta, dummy, idummy)
            write(*,*) 'energy range',-REAL(Erng),'--',REAL(Erng)
            write(*,*) '  #slices',Maxom
            write(*,*) 'Deta=',REAL(Deta)
         end if
         what = 'b'
         
         mk = Nkx
         call get_param_int('system size?',mk)
c         write(*,*) "system size?"
c         read(5,*) mk
         
         call calcDOS5(Maxom,Erng,Deta,mk)! 状態密度の計算 > optcond3.f
      end if
c===========================================]
      return


      end

c####################################################################


c#########################################################
c#### menu_suscept:  感受率の計算
c#########################################################

      subroutine menu_suscept(what)
      use common

      implicit none

      character, intent(inout) :: what*1

      character :: chara*15
      integer :: ierr, iflag
      integer :: idummy, i
      real*8 :: dummy


c===========================================[
c== suscept ================================
c==
      do while ( what == 's')
         write(*,*)'(l)ong, (t)rans or (b)oth?  [x: menu]'
         write(*,*)'chi(0)                      '
         read(5,*) what
         
         if ((what == 'l') .or. (what == 'L')) then
            iflag = 1
         else if ((what == 't') .or. (what == 'T')) then
            iflag = 2
         else if ((what == 'b') .or. (what == 'B')) then
            iflag = 3
         else if (what == '0') then
            iflag = 0
         else if (what == 'x') then
            what = 'd'
            return
         else
            what = 's'
         end if
      end do
      what = 's'
      
      chara = 's'
      do while ( chara(1:1) == 's')
         write(*,*)'Path?     (Pi = 1)', chara
         write(*,*)'a: 00-11, b: 00-10, c: 10-01  [x: menu]'
         write(*,*)'k: kmap, i: input by hand'
         read(5,*) chara
         
         if ((chara(1:1) == 'a') .or. (chara(1:1) == 'A')) then
            chara = '00-11'
         else if ((chara(1:1) == 'b') .or. (chara(1:1) == 'B')) then
            chara = '00-10'
         else if ((chara(1:1) == 'c') .or. (chara(1:1) == 'C')) then
            chara = '10-01'
         else if ((chara(1:1) == 'd') .or. (chara(1:1) == 'D')) then
            chara = '00-01'
         else if ((chara(1:1) == 'e') .or. (chara(1:1) == 'E')) then
            chara = '10-11'
         else if ((chara(1:1) == 'k') .or. (chara(1:1) == 'K')) then
            chara = 'k-map'
         else if ((chara(1:1) == 'i') .or. (chara(1:1) == 'I')) then
            chara = 'input'
         else if (chara(1:1) == 'x') then
            what = 'd'
            return
         else
            call pathget(chara(1:5))
         end if

         if (chara(1:1) == 'x') then
            what = 'd'
            return
         end if
      end do

      write(*,*) '---------------------------------'
      call work(what,ierr)
      if(ierr == 1) then
         write(*,*) 'calculation failed'
         return
      endif
      write(*,*) '---------------------------------'
!(MF)-> chi

      call suscept(iflag,chara)
c===========================================]

 999  return


 2000 format(32A1)
      end

c####################################################################

c#########################################################
c#### menu_RIXS:
c#########################################################

      subroutine menu_RIXS(what)
      use common

      implicit none

      character, intent(inout) :: what*1

      character :: chara*15, charalong*40, filename*80
      integer :: ierr, iflag, mode
      integer :: idummy, i, iom, ii, ikmd
      real*8 :: dummy, theta, phi


c===========================================[
c== RIXS ===================================
c==

      mode = 0
      do while (mode == 0)
         write(*,*)'mode?'
         write(*,*)'1: chi + RIXS  2: chi  3: RIXS'
         read(5,*,err=999) mode
         write(*,*) mode
      end do

      if (mode == 1) then
         write(*,*)'theta?'
         read(5,*) theta
         theta = theta / 180.0 * Pi
         write(*,*)'phi?'
         read(5,*) phi
         phi = phi / 180.0 * Pi
         
         chara = 's'
         do while ( chara(1:1) == 's')
            write(*,*)'Path?     (Pi = 1)', chara
            write(*,*)'a: 00-11, b: 00-10, c: 10-01  [x: menu]'
            write(*,*)'k: kmap, i: input by hand'
            read(5,*) chara
            
            if ((chara(1:1) == 'a')
     &           .or. (chara(1:1) == 'A')) then
               chara = '00-11'
            else if ((chara(1:1) == 'b') 
     &              .or. (chara(1:1) == 'B')) then
               chara = '00-10'
            else if ((chara(1:1) == 'c') 
     &              .or. (chara(1:1) == 'C')) then
               chara = '10-01'
            else if ((chara(1:1) == 'd') 
     &              .or. (chara(1:1) == 'D')) then
               chara = '00-01'
            else if ((chara(1:1) == 'e') 
     &              .or. (chara(1:1) == 'E')) then
               chara = '10-11'
            else if ((chara(1:1) == 'k') 
     &              .or. (chara(1:1) == 'K')) then
               chara = 'k-map'
            else if ((chara(1:1) == 'i') 
     &              .or. (chara(1:1) == 'I')) then
               chara = 'input'
            else if (chara(1:1) == 'x') then
c     chara(1:1) = 's'
               return
            else
               call pathget(chara(1:5))
            end if

            if (chara(1:1) == 'x') then
               chara(1:1) = 's'
               return
            end if
         end do
         chara(6:6) = '?'
         do while ( chara(6:6) == '?')
            write(*,*)'basis?', chara
            write(*,*)'1: Fe-As(XYZ), 2: Fe-Fe(xyz)  [x: menu]'
            read(5,*) chara(6:6)
            
            if (chara(6:6) == '1') then
               chara(6:10) = 'Fe-Fe'
            else if (chara(6:6) == '2') then
               chara(6:10) = 'Fe-As'
            else if (chara(6:6) == 'x') then
               what = 'd'
               return
            else
               chara(6:6) = '?'
            end if
         end do

         Muselect(:) = 0
         write(*,*)'select orbital'
         read(5,*) Muselect(1), Muselect(2),
     &        Muselect(3), Muselect(4)

         Ledgemode = 'both'
         write(*,*)'ledgemode?'
         read(5,*) Ledgemode

         Ledgeqdirection = '+'
         write(*,*)'ledgeqdirection?'
         read(5,*) Ledgeqdirection
         if (Ledgeqdirection /= '-')  Ledgeqdirection = '+'

         ModeJorbang = '---'
         do while (ModeJorbang(2:3) /= '/2')
            write(*,*)'orbaital angular momentum?'
            write(*,*)'1: 1/2   3: 3/2    [x: menu]'
            read(5,*,err=999) ModeJorbang(1:1)
            if (ModeJorbang(1:1) == 'x') then
               what = 'd'
               return
            else if (ModeJorbang(1:1) == '1') then
               ModeJorbang(2:3) = '/2'
            else if (ModeJorbang(1:1) == '3') then
               ModeJorbang(2:3) = '/2'
            end if
         end do
         write(*,*) 'ModeJorbang=',ModeJorbang

         write(*,*) '---------------------------------'
         call work(what,ierr)
         if(ierr == 1) then
            write(*,*) 'calculation failed'
            return
         endif
         write(*,*) '---------------------------------'
         
         write(*,*)'calculating RIXS spectra...'
         call Ledge(chara(1:5),theta,phi,chara(6:10))

      else if (mode == 2) then
         chara = 's'
         do while ( chara(1:1) == 's')
            write(*,*)'Path?     (Pi = 1)', chara
            write(*,*)'a: 00-11, b: 00-10, c: 10-01  [x: menu]'
            write(*,*)'k: kmap, i: input by hand'
            read(5,*) chara
            
            if ((chara(1:1) == 'a')
     &           .or. (chara(1:1) == 'A')) then
               chara = '00-11'
            else if ((chara(1:1) == 'b') 
     &              .or. (chara(1:1) == 'B')) then
               chara = '00-10'
            else if ((chara(1:1) == 'c') 
     &              .or. (chara(1:1) == 'C')) then
               chara = '10-01'
            else if ((chara(1:1) == 'd') 
     &              .or. (chara(1:1) == 'D')) then
               chara = '00-01'
            else if ((chara(1:1) == 'e') 
     &              .or. (chara(1:1) == 'E')) then
               chara = '10-11'
            else if ((chara(1:1) == 'k') 
     &              .or. (chara(1:1) == 'K')) then
               chara = 'k-map'
            else if ((chara(1:1) == 'i') 
     &              .or. (chara(1:1) == 'I')) then
               chara = 'input'
            else if (chara(1:1) == 'x') then
               what = 'd'
               return
            else
               call pathget(chara(1:5))
            end if

            if (chara(1:1) == 'x') then
               what = 'd'
               return
            end if
         end do

         chara(6:6) = '?'
         do while ( chara(6:6) == '?')
            write(*,*)'basis?', chara
            write(*,*)'1: Fe-As(XYZ), 2: Fe-Fe(xyz)  [x: menu]'
            read(5,*) chara(6:6)
            
            if (chara(6:6) == '1') then
               chara(6:10) = 'Fe-Fe'
            else if (chara(6:6) == '2') then
               chara(6:10) = 'Fe-As'
            else if (chara(6:6) == 'x') then
               what = 'd'
               return
            else
               chara(6:6) = '?'
            end if
         end do

         write(*,*) '---------------------------------'
         call work(what,ierr)
         if(ierr == 1) then
            write(*,*) 'calculation failed'
            return
         endif
         write(*,*) '---------------------------------'
         
         write(*,*)'calculating chi for RIXS spectra...'
         call Ledgechi(chara(1:5),chara(6:10))
      else if (mode == 3) then
         
         write(*,*)'filename?'
         read(5,*) filename

         write(*,*)'mode?'
         read(5,*) chara

         write(*,*)'theta?'
         read(5,*) theta
         theta = theta / 180.0 * Pi

         write(*,*)'phi?'
         read(5,*) phi
         phi = phi / 180.0 * Pi

         chara(6:6) = '?'
         do while ( chara(6:6) == '?')
            write(*,*)'basis?', chara
            write(*,*)'1: Fe-As(XYZ), 2: Fe-Fe(xyz)  [x: menu]'
            read(5,*) chara(6:6)
            
            if (chara(6:6) == '1') then
               chara(6:10) = 'Fe-Fe'
            else if (chara(6:6) == '2') then
               chara(6:10) = 'Fe-As'
            else if (chara(6:6) == 'x') then
               what = 'd'
               return
            else
               chara(6:6) = '?'
            end if
         end do

         Muselect(:) = 0
         write(*,*)'select orbital'
         read(5,*) Muselect(1), Muselect(2),
     &        Muselect(3), Muselect(4)

         Ledgemode = 'both'
         write(*,*)'ledgemode?'
         read(5,*) Ledgemode

         Ledgeqdirection = '+'
         write(*,*)'ledgeqdirection?'
         read(5,*) Ledgeqdirection
         if (Ledgeqdirection /= '-')  Ledgeqdirection = '+'

         ModeJorbang = '---'
         do while (ModeJorbang(2:3) /= '/2')
            write(*,*)'orbaital angular momentum?'
            write(*,*)'1: 1/2   3: 3/2  0: chitest  [x: menu]'
            read(5,*,err=999) ModeJorbang(1:1)
            if (ModeJorbang(1:1) == 'x') then
               what = 'd'
               return
            else if (ModeJorbang(1:1) == '1') then
               ModeJorbang(2:3) = '/2'
            else if (ModeJorbang(1:1) == '3') then
               ModeJorbang(2:3) = '/2'
            else if (ModeJorbang(1:1) == '0') then
               ModeJorbang(2:3) = '/2'
            else if (ModeJorbang(1:1) == 'v') then
               ModeJorbang(2:3) = '/2'
            end if
         end do
         write(*,*) 'ModeJorbang=',ModeJorbang

         write(*,*)'calculating RIXS spectra...'
         call Ledgespec(theta,phi,filename,chara(1:4),chara(6:10))
         call vertex
      end if
c===========================================]

 999  return


 2000 format(32A1)
      end

c####################################################################

c#########################################################
c#### menu_erange:
c#########################################################

      subroutine menu_Eparam(erng, maxom, deta, ominit, minom)
c      use common, only : Erng, Maxom, Deta, Ominit, Minom

      implicit none

c      integer, intent(inout) :: mode
      integer, intent(inout) :: maxom, minom
      real*8, intent(inout) :: erng, deta, ominit

      write(*,*) ' Enter energy range? ', REAL(erng)
      read (5,*,err=400) erng
      write(*,*) ' Enter #slices? ', maxom
      read (5,*,err=400) maxom

      deta = erng / DBLE(maxom) * 2.0d0
      write(*,*) ' Broadening? ', REAL(deta)
      read (5,*,err=400) deta
      write(*,*) ' initial energy? ', REAL(ominit)
      read (5,*,err=400) ominit
      write(*,*) ' initial step of omega? ', minom
      read (5,*,err=400) minom
         

 400  return
      end 
c####################################################################
      

c#########################################################
c#### get_param_int: set parameter
c#########################################################

      subroutine get_param_int(name,ival)

      implicit none

      character(LEN=*) :: name
      integer, intent(inout) :: ival
      character(LEN=1), parameter :: chend = '?'
      integer :: len_name

      len_name = INDEX(name,chend)

      write(*,*) 'Enter the value:'
      write(*,*) name(1:len_name),ival
      read (5,*,err=100) ival
      write(*,*) '==> Set to ', ival
      write(*,*)

 100  return
      end 
c####################################################################
      
c#########################################################
c#### menu_erange:
c#########################################################

      subroutine get_param_real(name,val)

      implicit none

      character(LEN=*) :: name
      real*8, intent(inout) :: val
      character(LEN=1), parameter :: chend = '?'
      integer :: len_name

      len_name = INDEX(name,chend)

      write(*,*) 'Enter the value:'
      write(*,*) name(1:len_name),val
      read (5,*,err=100) val
      write(*,*) '==> Set to ',val
      write(*,*)

 100  return
      end 
c####################################################################
      
c#########################################################
c#### menu_erange:
c#########################################################

      subroutine get_param_ch(name,val)

      implicit none

      character(LEN=*) :: name
      character(LEN=*), intent(inout) :: val
      character(LEN=1), parameter :: chend = '?'
      integer :: len_name

      len_name = INDEX(name,chend)

      write(*,*) 'Enter the value:'
      write(*,*) name(1:len_name),val
      read (5,*,err=100) val
      write(*,*) '==> Set to ',val
      write(*,*)

 100  return
      end 
c####################################################################
      

