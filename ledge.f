c################ 2010 Kaneshita ################

c#########################################################
c#### ledge:
c#########################################################

      Subroutine ledge(chara,theta,phi,basis)

      use common
      implicit none

c      integer, intent(inout):: iflag

      character, intent(in) :: chara*5, basis*5
      real*8, intent(in) :: theta, phi
      
      integer :: ikx, iky, kx1, kx2, ky1, ky2, ik, idummy
      integer :: kmax, kpi
      integer :: is, ir, ir1, ir2, k0x, k0y, iom, isl, ist
      integer :: mu1, mu2, mu3, mu4
      integer :: isE, isH, nsize, nQ
      integer :: kflag(0:Nkx,0:Nky)
      integer, external :: Ispinpair, Ivrtxinv, Ivrtx, Iorbcomb

      real*8 :: dummy1, dummy2, om
      real*8 :: ex1, ex2, ey1, ey2, ez1, ez2

      complex*16 :: vtex1(5,2,5,2),vtex2(5,2,5,2)

      real*8 :: spec1(0:Nkx,0:Nky), spec2(0:Nkx,0:Nky)

      complex*16 :: zdummy, zdummy1
      complex*16 :: zmtrx1(1:Ns**2*Nqx, 1:Ns**2*Nqx, 2)
      complex*16 :: zmtrx2(1:Ns**2*Nqx, 1:Ns**2*Nqx, 2)

      complex*16 :: zchil(1:Ns**2*Nqx, 1:Ns**2*Nqx)
      complex*16 :: zchit(1:Ns**2*Nqx, 1:Ns**2*Nqx)
      complex*16 :: zchiud(1:Ns**2*Nqx, 1:Ns**2*Nqx)

      real*8 :: ddtemp, ddmu, dummy

c      basis = 'Fe-As'

      isl = 1                   !spin configuration for longitudinal
      ist = 3                   !spin configuration for transverse

      nsize = Ns * Ns * Nqx
      Erange = MAXVAL(Eall) - MINVAL(Eall) + 1.0d0
      if (Erange < Erng) Erng = Erange

      if (basis == 'Fe-As') then
         call swfprod()
      else
         call Swfprod_orb()
      end if

      if (chara(1:1) == 'k') then
         write(*,'(a)')'path = k-map'
      else if (chara(1:1) == 'i') then
         write(*,'(a)')'input path'
         write(*,'(a)')'initial k?'
         read(5,*) kx1, ky1
         write(*,'(a)')'final k?'
         read(5,*) kx2, ky2
      else
         read(chara(1:1),*) kx1
         read(chara(2:2),*) ky1
         read(chara(4:4),*) kx2
         read(chara(5:5),*) ky2
         write(*,'(a,2i1,a,2i1)')'path = ', kx1, ky1, '-', kx2, ky2
      end if
      kx1 = kx1 * Nkx / 2     ! / 2
      ky1 = ky1 * Nky / 2     ! / 2

      kx2 = kx2 * Nkx / 2     !/ 2
      ky2 = ky2 * Nky / 2     !/ 2

      call ledgecoeff(theta,phi,vtex1,vtex2,kx2-kx1,ky2-ky1,basis)

c================================================================= 
      kmax = Nkx / Nqx / 2
      kmax = Nkx / Nqx
      kmax = MAX(ABS(kx2-kx1),ABS(ky2-ky1))
      
      open(30,file='out.ledge.dat',form='formatted')
      write(30,'(a)') '#Nkx__Nky'
      write(30,'(a,2i5)') '#', Nkx, Nky
      write(30,'(a)') '#U__J'
      write(30,'(a,2g15.8)') '#', U, DJ
      write(30,'(a)') '#Deta_Dtemp'
      write(30,'(a,2g15.8)') '#', Deta, Dtemp
      write(30,'(a)') '#Erange_Maxom'
      write(30,'(a,g20.10,i6)') '#', Erng, Maxom
      write(30,'(a,g20.10)') '#', Ominit
      
      if (chara /= 'k-map') then
         
         write(30,'(a,a)') '#path=', chara
         write(30,'(a)') 
     &        '#k om spec'
         
         do iom = Minom, Maxom
            kflag(:,:) = 0
            
            spec1(:,:) =0.0d0
            spec2(:,:) =0.0d0
            
            om = DBLE(iom) * Erng / DBLE(Maxom) + Ominit
            write(*,*)'chi0: om= ',iom,'/',maxom,' L-edge'
            
            do ik = 0, kmax
               write(*,'(10x,a,i5,a,i5,a)')'k=',ik,'/',kmax,' L-edge'
               if (ik == 0) then
                  k0x = kx1
                  k0y = ky1
               else
                  if (kx2 - kx1 == 0) then
                     k0x = kx1
                  else
                     k0x = k0x + (kx2 - kx1) / ABS(kx2 - kx1)
                  end if
                  if (ky2 - ky1 == 0) then
                     k0y = ky1
                  else
                     k0y = k0y + (ky2 - ky1) / ABS(ky2 - ky1)
                  end if
               end if
               
               dummy1 = 0.0d0
               ikx = k0x
               iky = k0y
               k0x = MOD(ikx,Nkx)
               k0y = MOD(iky,Nky)

      
c================== chizz & chinn ==========[
!     bare susceptibility
               
               do is = 1, 2
                  isE = Ispinpair(is,'E')

!$omp parallel do private(ir2,zdummy)
                  do ir1 = 1, nsize ; do ir2 = 1, nsize
                     call schi0_2(ir1,ir2,k0x,k0y,iom,is,zdummy)
                     zmtrx2(ir1,ir2,isE) = zdummy
                  end do ; end do
               end do
               
! susceptibility
               zmtrx1(:,:,:) = 0.0d0
               is = isl
               call schi_2(is,zmtrx2,zmtrx1,nsize)
               zchil(:,:) = zmtrx1(:,:,1) +  zmtrx1(:,:,2)
c================== chizz & chinn ==========]

c================== chix ===================[
               is = ist
               isE = Ispinpair(is,'E')
               isH = Ispinpair(is,'H')

!$omp parallel do private(ir2,zdummy)
               do ir1 = 1, nsize ; do ir2 = 1, nsize
                  call schi0_2(ir1,ir2,k0x,k0y,iom,is,zdummy)
                  zmtrx2(ir1,ir2,isE) = zdummy
               end do ; end do

               zmtrx1(:,:,:) = 0.0d0
               call schix2(is,zmtrx2,zmtrx1,nsize)
               zchit(:,:) = zmtrx1(:,:,isE)
c================== chix ===================]

               call calcledgespec(k0x,k0y,kmax,spec1,spec2,vtex1,vtex2,
     &              zchil,zchit,isl,ist,kflag)

            end do              ! DO LOOP for k

            do ik = 0, Nkx
               if (kx1 == kx2) then
                  k0x = kx1
               else if (kx1 < kx2) then
                  k0x = MOD(ik+kx1,Nkx)
               else if (kx1 > kx2) then
                  k0x = MOD(-ik+kx1+Nkx,Nkx)
               end if

               if (ky1 == ky2) then
                  k0y = ky1
               else if (ky1 < ky2) then
                  k0y = MOD(ik+ky1,Nky)
               else if (ky1 > ky2) then
                  k0y = MOD(-ik+ky1+Nky,Nky)
               end if
                                 
               write(30,'(i5,10g18.8)') ik, om,
     &              spec1(k0x,k0y), spec2(k0x,k0y)
            end do
            write(30,*)
         end do                 ! DO LOOP for Omega

      end if
c================== chizz & chinn ==========]

      close(30)
c============ chi(nn), chi(SzSz) ===========]

      close(40)

      return

c===========================================]

 999  return
      end

c########################################### ledge ####

c#########################################################
c#### ledgechi:
c#########################################################

      Subroutine ledgechi(chara,basis)

      use common
      implicit none

c      integer, intent(inout):: iflag

      character, intent(in) :: chara*5, basis*5
c      real*8, intent(in) :: theta

      logical :: torf
      character :: fileout*80, junk*80
      integer :: ikx, iky, kx1, kx2, ky1, ky2, ik, idummy
      integer :: kmax, kpi
      integer :: is, ir, ir1, ir2, k0x, k0y, iom
      integer :: mu1, mu2, mu3, mu4
      integer :: isE, isH, nsize, nQ
      integer :: kflag(0:Nkx,0:Nky)
      integer :: isl, ist
      integer, external :: Ispinpair, Ivrtxinv, Ivrtx, Iorbcomb

      real*8 :: dummy1, dummy2, om
      complex*16 :: zdummy, zdummy1
      complex*16, allocatable :: zmtrx1(:,:,:)
      complex*16, allocatable :: zmtrx2(:,:,:)

      complex*16, allocatable :: zchil(:,:)
      complex*16, allocatable :: zchit(:,:)
      complex*16, allocatable :: zchiud(:,:)

      real*8 :: ddtemp, ddmu, dummy

      isl = 1                   !spin configuration for longitudinal
      ist = 3                   !spin configuration for transverse

c==== allocation[
      torf = allocated(zmtrx1)
      if (torf .EQV. .true.) deallocate(zmtrx1)
      torf = allocated(zmtrx2)
      if (torf .EQV. .true.) deallocate(zmtrx2)
      torf = allocated(zchil)
      if (torf .EQV. .true.) deallocate(zchil)
      torf = allocated(zchit)
      if (torf .EQV. .true.) deallocate(zchit)
      torf = allocated(zchiud)
      if (torf .EQV. .true.) deallocate(zchiud)

      allocate(zmtrx1(Ns**2*Nqx, Ns**2*Nqx,2))
      allocate(zmtrx2(Ns**2*Nqx, Ns**2*Nqx,2))
      allocate(zchil(Ns**2*Nqx, Ns**2*Nqx))
      allocate(zchit(Ns**2*Nqx, Ns**2*Nqx))
      allocate(zchiud(Ns**2*Nqx, Ns**2*Nqx))
c==== allocation]

      nsize = Ns * Ns * Nqx
      Erange = MAXVAL(Eall) - MINVAL(Eall) + 1.0d0
      if (Erange < Erng) Erng = Erange

c      basis = 'Fe-As'
      if (basis == 'Fe-As') then
         call swfprod()         ! 2-Fe-per-unit-cell notation
      else
         call Swfprod_orb()     ! 1-Fe-per-unit-cell notation
      end if


      if (chara(1:1) == 'k') then
         write(*,'(a)')'path = k-map'
      else if (chara(1:1) == 'i') then
         write(*,'(a)')'input path'
         write(*,'(a)')'initial k?'
         read(5,*) kx1, ky1
         write(*,'(a)')'final k?'
         read(5,*) kx2, ky2
      else
         read(chara(1:1),*) kx1
         read(chara(2:2),*) ky1
         read(chara(4:4),*) kx2
         read(chara(5:5),*) ky2
         write(*,'(a,2i1,a,2i1)')'path = ', kx1, ky1, '-', kx2, ky2
      end if

      kx1 = kx1 * Nkx / 2     ! / 2
      ky1 = ky1 * Nky / 2     ! / 2

      kx2 = kx2 * Nkx / 2     !/ 2
      ky2 = ky2 * Nky / 2     !/ 2

c================================================================= 
      kmax = Nkx / Nqx / 2
      kmax = Nkx / Nqx
      kmax = MAX(ABS(kx2-kx1),ABS(ky2-ky1))
      
      open (unit=30, form='unformatted', status='unknown',
     &     file = 'out.chi4rixs.dat')
c      write(30) fileout

c     open(42,file='gomix.dat',form='formatted')
c     open(41,file='gomiy.dat',form='formatted')
c      open(30,file='out.ledge.dat',form='formatted')
      
      junk= '#Nkx_Nky_Ns_Nqx_U_J_Deta_Dtemp_'//
     &      'Erange_Minom_Maxom_Omini_kmax_isl_ist'
      write(30) junk
      write(30) Nkx, Nky, Ns, Nqx, U, DJ, Deta, Dtemp,
     &     Erng, Minom,Maxom, Ominit, kmax, isl, ist
      junk= '#path'
      write(30) junk
      write(30) chara(1:5)
      junk= '#basis'
      write(30) junk
      write(30) basis(1:5)
      
      do iom = Minom, Maxom
         kflag(:,:) = 0
         
         
         om = DBLE(iom) * Erng / DBLE(Maxom) + Ominit
         write(*,*)'chi0: om= ', iom, '/', maxom, ' L-edge'
         
         do ik = 0, kmax
            write(*,'(10x,a,i5,a,i5,a)')'k=',ik,'/',kmax,' L-edge'
            if (ik == 0) then
               k0x = kx1
               k0y = ky1
            else
               if (kx2 - kx1 == 0) then
                  k0x = kx1
               else
                  k0x = k0x + (kx2 - kx1) / ABS(kx2 - kx1)
               end if
               if (ky2 - ky1 == 0) then
                  k0y = ky1
               else
                  k0y = k0y + (ky2 - ky1) / ABS(ky2 - ky1)
               end if
            end if
            
            dummy1 = 0.0d0
            ikx = k0x
            iky = k0y
            k0x = MOD(ikx,Nkx)
            k0y = MOD(iky,Nky)
            
            
c==================chizz & chinn ==========[
!     bare susceptibility
            
            do is = 1, 2
               isE = Ispinpair(is,'E')
               
!     $omp parallel do private(ir2,zdummy)
               do ir1 = 1, nsize ; do ir2 = 1, nsize
                  call schi0_2(ir1,ir2,k0x,k0y,iom,is,zdummy)
                  zmtrx2(ir1,ir2,isE) = zdummy
               end do ; end do
            end do
            
!     susceptibility
            zmtrx1(:,:,:) = 0.0d0
            is = isl
            call schi_2(is,zmtrx2,zmtrx1,nsize)
            zchil(:,:) = zmtrx1(:,:,1)
            zchiud(:,:) = zmtrx1(:,:,2)
c            zchil(:,:) = zchil(:,:) + zchiud(:,:)
c==================chizz & chinn ==========]
            
c==================chix ===================[
            is = ist
            isE = Ispinpair(is,'E')
            isH = Ispinpair(is,'H')
            
!$omp parallel do private(ir2,zdummy)
            do ir1 = 1, nsize ; do ir2 = 1, nsize
               call schi0_2(ir1,ir2,k0x,k0y,iom,is,zdummy)
               zmtrx2(ir1,ir2,isE) = zdummy
            end do ; end do
            
            zmtrx1(:,:,:) = 0.0d0
            call schix2(is,zmtrx2,zmtrx1,nsize)
            zchit(:,:) = zmtrx1(:,:,isE)
c==================chix ===================]
            
c            junk= 'k0x__k0y__om'
c            write(30) junk
c            write(30) k0x,k0y,om
            
            junk= 'zchil'
            write(30) junk
            write(30) ((zchil(ir1,ir2), ir1=1,nsize), ir2=1,nsize)
            junk= 'zchit'
            write(30) junk
            write(30)  ((zchit(ir1,ir2), ir1=1,nsize), ir2=1,nsize)
            junk= 'zchiud'
            write(30) junk
            write(30)  ((zchiud(ir1,ir2), ir1=1,nsize), ir2=1,nsize)
            
         end do                 ! DO LOOP for k
      end do                    ! DO LOOP for Omega
            
c================== chizz & chinn ==========]

      close(30)
c============ chi(nn), chi(SzSz) ===========]


      return

c===========================================]

      return
      end

c########################################### ledgechi ####


c#########################################################
c#### ledgespec:
c#########################################################

      Subroutine ledgespec(theta0,phi0,filein,mode,basisin)

      use common
      implicit none

c      integer, intent(inout):: iflag

      real*8, intent(in) :: theta0, phi0
      character, intent(in) :: filein*80, mode*4, basisin*5

      logical :: torf
      character :: chara*5, basis*5
      character :: junk*80
      integer :: ikx, iky, kx1, kx2, ky1, ky2, ik, idummy
      integer :: kmax, kpi, kxf, kyf, i
      integer :: is, ir, ir1, ir2, k0x, k0y, iom, isl, ist
      integer :: mu1, mu2, mu3, mu4
      integer :: isE, isH, nsize, nQ
      integer, allocatable :: kflag(:,:)
      integer, external :: Ispinpair, Ivrtxinv, Ivrtx, Iorbcomb
      integer :: idirection, id

      real*8 :: theta, phi
      real*8 :: dummy1, dummy2, om
      real*8 :: ex1, ex2, ey1, ey2, ez1, ez2
      real*8 :: akx, aky, ak, akledge

      complex*16 :: ciqua(5,2,4),cidou(5,2,2)
      complex*16 :: coqua(5,2,4),codou(5,2,2)
      complex*16 :: vtex1(5,2,5,2),vtex2(5,2,5,2)

      real*8, allocatable :: spec1(:,:), spec2(:,:)
      real*8, allocatable :: sp1(:,:), sp2(:,:)

      complex*16 :: zdummy, zdummy1
      complex*16, allocatable :: zchil(:,:)
      complex*16, allocatable :: zchit(:,:)
      complex*16, allocatable :: zchiud(:,:)

      real*8 :: ddtemp, ddmu, dummy

c      call azimuth2 !aaaaaaa

      isl = 1                   !spin configuration for longitudinal
      ist = 3                   !spin configuration for transverse
c      phi = 90.0d0
c      theta=dfloat(90-i)/180.d0*pi
c      theta = 0.0d0

c=====.   phi: polarization
c=====. theta: incident angle 
c=====.    
c=====.   \theta| 
c=====.    \    |
c=====.     \   |
c=====.      \  |
c=====.       \ |
c=====. ----------------------

      theta = theta0
      phi = phi0
         
      open (unit=40, form='unformatted', status='OLD', 
     &     file=filein, err=998)

      read(40) junk
      if (junk(1:4) /= '#Nkx') then
         write(*,*) '#Nkx'
         write(*,*) junk
         goto 999
      end if
      read(40) Nkx, Nky, Ns, Nqx, U, DJ, Deta, Dtemp,
     &     Erng, Minom, Maxom, Ominit, kmax, isl, ist

      read(40) junk
      if (junk(1:5) /= '#path') then
         write(*,*) '#path'
         write(*,*) junk
         goto 999
      end if
      read(40) chara(1:5)

      read(40) junk
      if (junk(1:6) == '#basis') then
         read(40) basis
         if (basis /= basisin) then
            write(*,*) ' Error in ledgespec'
            write(*,*) ' inconsistent basis'
         end if
      else
         stop
         basis = basisin
         backspace(40)
      end if

      torf = allocated(zchil)
      if (torf .EQV. .true.) deallocate(zchil)
      torf = allocated(zchit)
      if (torf .EQV. .true.) deallocate(zchit)
      torf = allocated(zchiud)
      if (torf .EQV. .true.) deallocate(zchiud)
      torf = allocated(spec1)
      if (torf .EQV. .true.) deallocate(spec1)
      torf = allocated(spec2)
      if (torf .EQV. .true.) deallocate(spec2)
      torf = allocated(kflag)
      if (torf .EQV. .true.) deallocate(kflag)
      allocate(zchil(Ns**2*Nqx, Ns**2*Nqx))
      allocate(zchit(Ns**2*Nqx, Ns**2*Nqx))
      allocate(zchiud(Ns**2*Nqx, Ns**2*Nqx))
      allocate( sp1(-Nkx:Nkx,-Nky:Nky), sp2(-Nkx:Nkx,-Nky:Nky) )
      allocate( spec1(0:Nkx,0:Nky), spec2(0:Nkx,0:Nky) )
      allocate( kflag(0:Nkx,0:Nky) )
      nsize = Ns * Ns * Nqx

      read(chara(1:1),*) kx1
      read(chara(2:2),*) ky1
      read(chara(4:4),*) kx2
      read(chara(5:5),*) ky2
      write(*,'(a,2i1,a,2i1)')'path = ', kx1, ky1, '-', kx2, ky2
      
      kx1 = kx1 * Nkx / 2     ! / 2
      ky1 = ky1 * Nky / 2     ! / 2

      kx2 = kx2 * Nkx / 2     !/ 2
      ky2 = ky2 * Nky / 2     !/ 2

      if (mode(1:4) == 'kfix') then
         kxf = INT(theta * 180.0d0 * Nkx / Pi /2) 
         kyf = INT(phi * 180.0d0 * Nky / Pi / 2)
         idummy = 0
         if ( ( (kx2-kx1)*kxf == 0 )
     &        .AND.( kxf /= (kx2-kx1) ) ) idummy = 1
         if ( ( (ky2-ky1)*kyf == 0 )
     &        .AND.( kyf /= (ky2-ky1) ) ) idummy = 1
         if (idummy /= 0) then
            write(*,*) '---Error in ledgespec---'
            write(*,*) 'momentum transfer is illegal:'
            write(*,*) 'qx__qy = ', kxf, kyf
            if ((kxf == 0).AND.(kyf == 0)) then
               write(*,*) 'q needs to be a non-zero vector'
            else
               write(*,*) 'Input file has no chi data for this q'
            end if
            open(30,file='out.ledge.dat',form='formatted')
            write(30,*)'Error -- invalid q was input'
            close(30)
            stop
         end if
      end if
c================================================================= 
c      kmax = Nkx / Nqx / 2
c      kmax = Nkx / Nqx
c      kmax = MAX(ABS(kx2-kx1),ABS(ky2-ky1))

      open(30,file='out.ledge.dat',form='formatted')
      write(30,'(a)') '#Nkx__Nky'
      write(30,'(a,2i5)') '#', Nkx, Nky
      write(30,'(a)') '#U__J'
      write(30,'(a,2g15.8)') '#', U, DJ
      write(30,'(a)') '#Deta_Dtemp'
      write(30,'(a,2g15.8)') '#', Deta, Dtemp
      write(30,'(a)') '#Erange_Maxom'
      write(30,'(a,g20.10,i6)') '#', Erng, Maxom
      write(30,'(a,g20.10)') '#Ominit=', Ominit
      write(30,'(a,a)') '#path=', chara
      if (mode(1:4) /= 'kfix') then
         write(30,'(a,2f6.2)') '#theta_phi=', theta, phi
      else
         write(30,'(a,2i5)') '#kx_ky=', kxf, kyf
      end if
      write(30,'(a,x,a)')  '#excitation_mode:', Ledgemode
      write(30,'(a,4i5)')  '#orbital_selection_1234=', Muselect(1:4)
      write(30,'(a,x,a)')  '#basis:', basis
      write(30,'(a,a)') '#Jorbang=', ModeJorbang
      write(30,'(a)')  '#k om spec(plr1) spec(plr2)'
      
      write(*,'(a,2i5)') '#Nkx__Nky=', Nkx, Nky
      write(*,'(a,2g15.8)') '#U__J=', U, DJ
      write(*,'(a,2g15.8)') '#Deta_Dtemp=', Deta, Dtemp
      write(*,'(a,g20.10,i6)') '#Erange_Maxom=', Erng, Maxom
      write(*,'(a,g20.10)') '#Ominit=', Ominit
      write(*,'(a,a)') '#path=', chara(1:5)
      if (mode(1:4) /= 'kfix') then
         write(*,'(a,2f6.2)') '#theta_phi=', theta, phi
      else
         write(*,'(a,2i5,x,a,2(f10.5,a),a)') '#kx_ky=', kxf, kyf,
     &        '(', DBLE(2*kxf)/DBLE(Nkx),'pi',
     &        DBLE(2*kyf)/DBLE(Nky),'pi',')'
      end if
      write(*,'(a,x,a)')  '#excitation_mode:', Ledgemode
      write(*,'(a,4i5)')  '#orbital_selection_1234=', Muselect(1:4)
      write(*,'(a)')  '#k om spec(plr1) spec(plr2)'
      write(*,'(a,a)') '#Jorbang=', ModeJorbang

      do iom = Minom, Maxom
         kflag(:,:) = 0
         
         spec1(:,:) =0.0d0
         spec2(:,:) =0.0d0
         sp1(:,:) =0.0d0
         sp2(:,:) =0.0d0
         
         om = DBLE(iom) * Erng / DBLE(Maxom) + Ominit
         write(*,*)'chi0: om= ',iom,'/',maxom,' L-edge '
         
         do ik = 0, kmax
            write(*,'(10x,a,i5,a,i5,a)')'k=',ik,'/',kmax,' L-edge'
            if (ik == 0) then
               k0x = kx1
               k0y = ky1
            else
               if (kx2 - kx1 == 0) then
                  k0x = kx1
               else
                  k0x = k0x + (kx2 - kx1) / ABS(kx2 - kx1)
               end if
               if (ky2 - ky1 == 0) then
                  k0y = ky1
               else
                  k0y = k0y + (ky2 - ky1) / ABS(ky2 - ky1)
               end if
            end if
            akx = DBLE(k0x * 2) / DBLE(Nkx)
            aky = DBLE(k0y * 2) / DBLE(Nky)
            
            dummy1 = 0.0d0
            ikx = k0x
            iky = k0y
            k0x = MOD(ikx,Nkx)
            k0y = MOD(iky,Nky)
            
c================== chi ==========[
            read(40) junk

            if (junk(1:5) /= 'zchil') then
               write(*,*) 'zchil'
               write(*,*) junk
               do while (junk(1:5) /= 'zchil')
                  read(40) junk(1:5)
                  write(*,*) junk(1:5)
               end do
c               goto 999
            end if

c            read(40) junk
c               write(*,*) junk!aaaaaaaaaaaaaaaaaaaaaa
            read(40) ((zchil(ir1,ir2), ir1=1,nsize), ir2=1,nsize)

            read(40) junk
            if (junk(1:5) /= 'zchit') then
               write(*,*) 'zchit'
               write(*,*) junk
               goto 999
            end if
            read(40)  ((zchit(ir1,ir2), ir1=1,nsize), ir2=1,nsize)

            read(40) junk
            if (junk(1:5) /= 'zchiu') then
               write(*,*) 'zchiu'
               write(*,*) junk
c               do while (junk(1:5) /= 'zchiu')
c                  read(40) junk(1:5)
c                  write(*,*) junk(1:5)
c               end do
               goto 999
            end if
            read(40)  ((zchiud(ir1,ir2), ir1=1,nsize), ir2=1,nsize)
c================== chi ==========]

            if (mode(1:4) /= 'kfix') then !k-dependence
               akledge = 0.4526 ! kl=0.4526pi
               ak = SQRT(akx**2 + aky**2)
               if ((ABS(ak/akledge) > SQRT(2.0d0))
     &              .and.(ModeJorbang(1:1) /= '0')) cycle
               
               do idirection = 1, 2
                  theta = 0.25d0 * Pi
                  id = SIGN(1,(-1)**idirection)
c==== relation between momentum transfer & theta
c==== ak = -sqrt(2) * |k_in| * sin(theta - pi/4)
c     if (Ledgeqdirection == '-') then
c     theta = theta - ASIN(ak / akledge / SQRT(2.0d0))
c     else
c     theta = theta + ASIN(ak / akledge / SQRT(2.0d0))
c     end if
                  theta = theta + id * ASIN(ak / akledge / SQRT(2.0d0))
                  
                  call ledgecoeff(theta,phi,vtex1,vtex2,
     &                 kx2 - kx1,ky2 - ky1,basis)
                  call calcledgespec(k0x,k0y,kmax,spec1,spec2
     &                 ,vtex1,vtex2,zchil,zchit,zchiud,isl,ist,kflag)
                  sp1(id*k0x,id*k0y) = spec1(k0x,k0y)
                  sp2(id*k0x,id*k0y) = spec2(k0x,k0y)
                  if (id == -1) kflag(k0x,k0y) = 0
                  spec1(k0x,k0y) =0.0d0
                  spec2(k0x,k0y) =0.0d0
               end do
                              
            else                !momentum transfer fixed to (k0x,k0y)
               if ((k0x /= kxf).or.(k0y /= kyf)) cycle
               idummy = 0
               phi = 0.d0
               do i = 0, 90
                  kflag(k0x,k0y) = 0
                  spec1(k0x,k0y) =0.0d0
                  spec2(k0x,k0y) =0.0d0
                  theta = DBLE(90-i) / 180.d0 * pi
                  call ledgecoeff(theta,phi,vtex1,vtex2,k0x,k0y,basis)
                  call calcledgespec(k0x,k0y,kmax,spec1,spec2,
     &                 vtex1,vtex2,zchil,zchit,zchiud,isl,ist,kflag)
                  write(30,'(10g18.8)') DBLE(i + idummy)/90.0d0, om,
     &                 spec1(k0x,k0y),spec2(k0x,k0y)
               end do

               idummy = idummy + 90
               theta = 0.d0
               do i = 1, 90
                  kflag(k0x,k0y) = 0
                  spec1(k0x,k0y) =0.0d0
                  spec2(k0x,k0y) =0.0d0
                  phi = DBLE(i) / 180.d0 * pi
                  call ledgecoeff(theta,phi,vtex1,vtex2,k0x,k0y,basis)
                  call calcledgespec(k0x,k0y,kmax,spec1,spec2,
     &                 vtex1,vtex2,zchil,zchit,zchiud,isl,ist,kflag)
                  write(30,'(10g18.8)') DBLE(i + idummy)/90.0d0, om,
     &                 spec1(k0x,k0y),spec2(k0x,k0y)
               end do

               idummy = idummy + 90
               phi = 90.d0 / 180.d0 * pi
               do i = 1, 90
                  kflag(k0x,k0y) = 0
                  spec1(k0x,k0y) =0.0d0
                  spec2(k0x,k0y) =0.0d0
                  theta = DBLE(i) / 180.d0 * pi
                  call ledgecoeff(theta,phi,vtex1,vtex2,k0x,k0y,basis)
                  call calcledgespec(k0x,k0y,kmax,spec1,spec2,
     &                 vtex1,vtex2,zchil,zchit,zchiud,isl,ist,kflag)
                  write(30,'(10g18.8)') DBLE(i + idummy)/90.0d0, om,
     &                 spec1(k0x,k0y),spec2(k0x,k0y)
               end do
               
               idummy = idummy + 90
               theta = 90.d0 / 180.d0 * pi
               do i=1,90
                  kflag(k0x,k0y) = 0
                  spec1(k0x,k0y) =0.0d0
                  spec2(k0x,k0y) =0.0d0
                  phi = DBLE(90-i) / 180.d0 * pi
                  call ledgecoeff(theta,phi,vtex1,vtex2,k0x,k0y,basis)
                  call calcledgespec(k0x,k0y,kmax,spec1,spec2,
     &                 vtex1,vtex2,zchil,zchit,zchiud,isl,ist,kflag)
                  write(30,'(10g18.8)') DBLE(i + idummy)/90.0d0, om,
     &                 spec1(k0x,k0y),spec2(k0x,k0y)
               end do
               
            end if

         end do                 ! DO LOOP for k
         
         if (mode(1:4) /= 'kfix') then
            do ik = -Nkx, Nkx
               if (kx1 == kx2) then
                  k0x = kx1
               else if (kx1 < kx2) then
                  k0x = MOD(ik+kx1,Nkx)
               else if (kx1 > kx2) then
                  k0x = MOD(-ik+kx1+Nkx,Nkx)
               end if
               
               if (ky1 == ky2) then
                  k0y = ky1
               else if (ky1 < ky2) then
                  k0y = MOD(ik+ky1,Nky)
               else if (ky1 > ky2) then
                  k0y = MOD(-ik+ky1+Nky,Nky)
               end if
               
               write(30,'(i5,10g18.8)') ik, om,
     &              sp1(k0x,k0y), sp2(k0x,k0y)
            end do
         end if
         write(30,*)

      end do                    ! DO LOOP for Omega
      
      close(30)

      close(40)


      return

c===========================================]

      return

 998  write (*,*) 'Cannot open file'
      return

 999  write (*,*) 'ERROR while loading file'
      return

      end

c########################################### ledgespec ####

c#########################################################
c#### calcledgespec:
c#########################################################

      Subroutine calcledgespec(k0x,k0y,kmax,spec1,spec2,vtex1,vtex2,
     &     zchil,zchit,zchiud,isl,ist,kflag)

      use common, only : Nkx, Nky, Nqx, Ns, ZI,
     &     Muselect, Ledgemode, ModeJorbang
      implicit none

c      integer, intent(inout):: iflag

      real*8, intent(inout) :: spec1(0:Nkx,0:Nky), spec2(0:Nkx,0:Nky)
      complex*16, intent(in) :: zchil(1:Ns**2*Nqx, 1:Ns**2*Nqx)
      complex*16, intent(in) :: zchit(1:Ns**2*Nqx, 1:Ns**2*Nqx)
      complex*16, intent(in) :: zchiud(1:Ns**2*Nqx, 1:Ns**2*Nqx)
      complex*16, intent(in) :: vtex1(5,2,5,2),vtex2(5,2,5,2)
      integer, intent(in) :: k0x, k0y, isl, ist, kmax
      integer, intent(inout) :: kflag(0:Nkx,0:Nky)

      character :: mode*1, flagorb
      integer :: ikx, iky, ik, idummy
      integer :: is, ir, ir1, ir2
      integer :: mu1, mu2, mu3, mu4, nu1, nu2, nu3, nu4
      integer :: isE, isH, nQ, isE2, isH2

      character :: chitest*1
      character, external :: orbselect
      integer, external :: Ispinpair, Ivrtxinv, Ivrtx
      integer, external :: Ispinflip

      complex*16 :: zdummy, zdummy1, zdummy2, zzz=0.0d0
      real*8 :: test=0.0d0


      chitest = 'n'
      if (ModeJorbang(1:1) == '0') chitest = 'y'!aaaaaaaaaaa
      if (ModeJorbang(1:1) == 'v') chitest = 'v'!aaaaaaaaaaa
c      chitest = 'v'!aaaaaaaaaaa

      nu1 = Muselect(1)
      nu2 = Muselect(2)
      nu3 = Muselect(3)
      nu4 = Muselect(4)

      mode = Ledgemode(1:1)
      do nQ = 0, Nqx - 1
         ikx = k0x
         iky = k0y
         if (kmax <= Nkx / 2) ikx = k0x + nQ * Nkx / Nqx
         ikx = MOD(ikx,Nkx)
         
         if (kflag(ikx,iky) == 1) then
            cycle
         else if (kflag(ikx,iky) == 0) then
            kflag(ikx,iky) = 1
         else 
            write(*,*) 'error'
            stop
         end if
         do mu1 = 1, Ns ; do mu2 = 1, Ns
            ir1 = Ivrtx(mu1,mu2,nQ,'L')
            do mu3 = 1, Ns ; do mu4 = 1, Ns
               ir2 = Ivrtx(mu3,mu4,nQ,'R')
c               idummy = iorbcomb(mu1,mu2,mu3,mu4)
               flagorb = orbselect(mu1,mu2,mu3,mu4,nu1,nu2,nu3,nu4)
               if (flagorb == 'n') cycle
c                        if (mu1 /= 0) cycle
c                        if (idummy == 4) cycle 
c                        if (idummy == 3) cycle 

               if (mode /= 't') then
c                 is = isl
c                 isE = Ispinpair(is,'E')
c                 isH = Ispinpair(is,'H')
c                 isE2 = Ispinflip(isE)
c                 isH2 = Ispinflip(isH)
                  isE = 1
                  isH = 1
                  isE2 = 2
                  isH2 = 2

c==== :       axis:       x                 z
c==== :   L-vertex:  spin1 spin2 <==>  spin1 spin2 
c==== :   L-vertex:    up    up  <==>  (u+d) (u+d)
c==== :   L-vertex:    dn    dn  <==>  (u-d) (u-d)
c==== :   L-vertex:    up    dn  <==>  (u+d) (u-d)
c==== :   L-vertex:    dn    up  <==>  (u-d) (u+d)
c==== :   
c==== :   
c==== :    v(ch)12:  (uu + dd)/2 <==>  (uu + dd)/2
c==== :    v(sz)12:  (ud + du)/2 <==>  (uu - dd)/2
c==== :    v(sx)12:  (uu - dd)/2 <==>  (ud + du)/2


c==== polarization e1'
!     charge: [(12 + 12)^* (34 + 34)]/4 (chi_1234)
!     charge: [(uu + dd)^* (uu + dd)]/4 (uuuu + dddd + uudd + dduu)

                  zdummy1 =  vtex1(mu1,isH ,mu2,isE )
     &                 +     vtex1(mu1,isH2,mu2,isE2)
                  zdummy2 =  vtex1(mu4,isH ,mu3,isE )
     &                 +     vtex1(mu4,isH2,mu3,isE2)
                  zdummy = DCONJG(zdummy1) * zdummy2 * 0.25d0

                  if (chitest == 'y') zdummy = 1.0d0 !aaaaaaaaa

                  zdummy1 = (zchil(ir1,ir2) + zchiud(ir1,ir2)) * 2.0d0

                  if (chitest == 'v') zdummy1 = -ZI !aaaaaaaaa

                  zdummy = zdummy * zdummy1

                  if (mode /= 'l') then
                     spec1(ikx,iky) = spec1(ikx,iky) - DIMAG(zdummy)
                  end if

!     sz: [(12 + 12)^* (34 + 34)]/4 (chi_1234)
!     sz: [(ud + du)^* (du + ud)]/4 (uuuu + dddd - uudd - dduu)

                  zdummy1 =  vtex1(mu1,isH ,mu2,isE2)
     &                 +     vtex1(mu1,isH2,mu2,isE )
                  zdummy2 =  vtex1(mu4,isH ,mu3,isE2)
     &                 +     vtex1(mu4,isH2,mu3,isE )
                  zdummy = DCONJG(zdummy1) * zdummy2 * 0.25d0

                  if (chitest == 'y') zdummy = 1.0d0 !aaaaaaaaa

                  zdummy1 = (zchil(ir1,ir2) - zchiud(ir1,ir2)) * 2.0d0

                  if (chitest == 'v') zdummy1 = -ZI !aaaaaaaaa

                  zdummy = zdummy * zdummy1

                  if (mode /= 'c') then
                     spec1(ikx,iky) = spec1(ikx,iky) - DIMAG(zdummy)
                  end if

c==== polarization e2'
!     charge: [(12 + 12)^* (34 + 34)]/4 (chi_1234)
!     charge: [(uu + dd)^* (uu + dd)]/4 (uuuu + dddd + uudd + dduu)
                  zdummy1 =  vtex2(mu1,isH ,mu2,isE )
     &                 +     vtex2(mu1,isH2,mu2,isE2)
                  zdummy2 =  vtex2(mu4,isH ,mu3,isE )
     &                 +     vtex2(mu4,isH2,mu3,isE2)
                  zdummy = DCONJG(zdummy1) * zdummy2 * 0.25d0

                  if (chitest == 'y') zdummy = 1.0d0 !aaaaaaaaa

                  zdummy1 = (zchil(ir1,ir2) + zchiud(ir1,ir2)) * 2.0d0

                  if (chitest == 'v') zdummy1 = -ZI !aaaaaaaaa

                  zdummy = zdummy * zdummy1

                  if (mode /= 'l') then
                     spec2(ikx,iky) = spec2(ikx,iky) - DIMAG(zdummy)
                  end if

!     sz: [(12 + 12)^* (34 + 34)]/4 (chi_1234)
!     sz: [(ud + du)^* (du + ud)]/4 (uuuu + dddd - uudd - dduu)
                  zdummy1 =  vtex2(mu1,isH ,mu2,isE2)
     &                 +     vtex2(mu1,isH2,mu2,isE )
                  zdummy2 =  vtex2(mu4,isH ,mu3,isE2)
     &                 +     vtex2(mu4,isH2,mu3,isE )
                  zdummy = DCONJG(zdummy1) * zdummy2 * 0.25d0

                  if (chitest == 'y') zdummy = 1.0d0 !aaaaaaaaa

                  zdummy1 = (zchil(ir1,ir2) - zchiud(ir1,ir2)) * 2.0d0

                  if (chitest == 'v') zdummy1 = -ZI !aaaaaaaaa

                  zdummy = zdummy * zdummy1

                  if (mode /= 'c') then
                     spec2(ikx,iky) = spec2(ikx,iky) - DIMAG(zdummy)
                  end if
               end if

               if ((mode /= 'l').and.(mode /= 'c')) then
c                 is = ist
c                 isE = Ispinpair(is,'E')
c                 isH = Ispinpair(is,'H')
c                 isE2 = Ispinflip(isE)
c                 isH2 = Ispinflip(isH)
                  isE = 1
                  isH = 1
                  isE2 = 2
                  isH2 = 2
c==== polarization e1'
!     up-dn: (u+d)(u-d)
                  zdummy1 =  vtex1(mu1,isH ,mu2,isE )
     &                 -     vtex1(mu1,isH ,mu2,isE2)
     &                 +     vtex1(mu1,isH2,mu2,isE )
     &                 -     vtex1(mu1,isH2,mu2,isE2)
                  zdummy2 =  vtex1(mu4,isH ,mu3,isE )
     &                 -     vtex1(mu4,isH ,mu3,isE2)
     &                 +     vtex1(mu4,isH2,mu3,isE )
     &                 -     vtex1(mu4,isH2,mu3,isE2)
                  zdummy = DCONJG(zdummy1) * zdummy2 * 0.25d0

                  if (chitest == 'y') zdummy = 1.0d0 !aaaaaaaaa

                  zdummy1 = zchit(ir1,ir2)
                  if (chitest == 'v') zdummy1 = -ZI !aaaaaaaaa
                  zdummy = zdummy * zdummy1

                  spec1(ikx,iky) = spec1(ikx,iky) - DIMAG(zdummy)

!     dn-up: (u-d)(u+d)
                  zdummy1 =  vtex1(mu1,isH ,mu2,isE )
     &                 +     vtex1(mu1,isH ,mu2,isE2)
     &                 -     vtex1(mu1,isH2,mu2,isE )
     &                 -     vtex1(mu1,isH2,mu2,isE2)
                  zdummy2 =  vtex1(mu4,isH ,mu3,isE )
     &                 +     vtex1(mu4,isH ,mu3,isE2)
     &                 -     vtex1(mu4,isH2,mu3,isE )
     &                 -     vtex1(mu4,isH2,mu3,isE2)
                  zdummy = DCONJG(zdummy1) * zdummy2 * 0.25d0

                  if (chitest == 'y') zdummy = 1.0d0 !aaaaaaaaa

                  zdummy1 = zchit(ir1,ir2)
                  if (chitest == 'v') zdummy1 = -ZI !aaaaaaaaa
                  zdummy = zdummy * zdummy1

                  spec1(ikx,iky) = spec1(ikx,iky) - DIMAG(zdummy)

c==== polarization e2'
!     up-dn: (u+d)(u-d)
                  zdummy1 =  vtex2(mu1,isH ,mu2,isE )
     &                 -     vtex2(mu1,isH ,mu2,isE2)
     &                 +     vtex2(mu1,isH2,mu2,isE )
     &                 -     vtex2(mu1,isH2,mu2,isE2)
                  zdummy2 =  vtex2(mu4,isH ,mu3,isE )
     &                 -     vtex2(mu4,isH ,mu3,isE2)
     &                 +     vtex2(mu4,isH2,mu3,isE )
     &                 -     vtex2(mu4,isH2,mu3,isE2)
                  zdummy = DCONJG(zdummy1) * zdummy2 * 0.25d0

                  if (chitest == 'y') zdummy = 1.0d0 !aaaaaaaaa

                  zdummy1 = zchit(ir1,ir2)
                  if (chitest == 'v') zdummy1 = -ZI !aaaaaaaaa
                  zdummy = zdummy * zdummy1

                  spec2(ikx,iky) = spec2(ikx,iky) - DIMAG(zdummy)

!     dn-up: (u-d)(u+d)
                  zdummy1 =  vtex2(mu1,isH ,mu2,isE )
     &                 +     vtex2(mu1,isH ,mu2,isE2)
     &                 -     vtex2(mu1,isH2,mu2,isE )
     &                 -     vtex2(mu1,isH2,mu2,isE2)
                  zdummy2 =  vtex2(mu4,isH ,mu3,isE )
     &                 +     vtex2(mu4,isH ,mu3,isE2)
     &                 -     vtex2(mu4,isH2,mu3,isE )
     &                 -     vtex2(mu4,isH2,mu3,isE2)
                  zdummy = DCONJG(zdummy1) * zdummy2 * 0.25d0

                  if (chitest == 'y') zdummy = 1.0d0 !aaaaaaaaa

                  zdummy1 = zchit(ir1,ir2)
                  if (chitest == 'v') zdummy1 = -ZI !aaaaaaaaa
                  zdummy = zdummy * zdummy1

                  spec2(ikx,iky) = spec2(ikx,iky) - DIMAG(zdummy)

               end if
            end do ; end do
            
         end do ; end do
      end do                    ! nQ loop
            
c                  write(*,*)zzz,'test'
      return

      end

c########################################### calcledgespec ####

c#########################################################
c#### ledgecoeff:
c#########################################################

      Subroutine ledgecoeff(theta,phi,vtex1,vtex2,kqx,kqy,basis)

      use common, only : Pi, ModeJorbang
      implicit none

      character, intent(in) :: basis*5
      integer, intent(in) :: kqx, kqy
      real*8 :: theta, phi
      complex*16, intent(out) :: vtex1(5,2,5,2),vtex2(5,2,5,2)

      integer :: is1, is2, iorb, nJorbang
      real*8 :: ex1, ex2, ey1, ey2, ez1, ez2, angle
      real*8 :: dzero

      complex*16 :: ciqua(5,2,4),cidou(5,2,2)
      complex*16 :: coqua(5,2,4),codou(5,2,2)

      real*8 :: dummy, x, y

c==== . vtex1(mu,s;nu,s'): coefficient C^1_(mu,s;nu,s')
c==== .          (in) polarized to e
c==== .         (out) e1' polarization (parallel to incidnet k)
c==== . vtex2(mu,s;nu,s'): coefficient C^2_(mu,s;nu,s')
c==== .          (in) polarized to e
c==== .         (out) e2' polarization (perpendicular to incident k)
c==== .
c==== .          vtex(1,s';2,s) = vtex(out;in)
c==== .             out: beam outgoing + distract an electron in 3d
c==== .              in: beam incoming + create an electron in 3d
c==== .
c==== .
c==== .        L-vertex (1,2,Q)       |   R-vertex (1,2,Q)
c==== .                               |
c==== .               in              |              in
c==== .              2,s (k)+Q+(q) => |  => (k)+(q)+Q 3,s
c==== .   (k)+Q =>  <                 |                  > => (k)+Q
c==== .              1,s'  (q) <=     |      <= (q)   4,s'
c==== .               out             |              out
c==== . 
c==== .          vtex1(1,s';2,s)              [vtex1(4,s';3,s)]^*
c==== .          vtex2(1,s';2,s)              [vtex2(4,s';3,s)]^*
c==== . 
c==== .                   2  -->--  3
c==== .  chi(1,2,3,4) =    <///////>
c==== .                   1  --<--  4
c==== .
c==== .  I1 = Im{ vtex1(2,s;1,s') chi(1,2;3,4) [vtex1(3,s;4,s')]^* }
c==== .  I2 = Im{ vtex2(2,s;1,s') chi(1,2;3,4) [vtex2(3,s;4,s')]^* }
c==== .  I = I1 + I2
c==== .  
c==== .
c==== .  
c==== . / e_x' \    / cos(90-th)  0   -sin(90-th) \ / e_x \
c==== . | e_y' | =  |    0        1      0        | | e_y |
c==== . \ e_z' /    \ sin(90-th)  0    cos(90-th) / \ e_z /
c==== .  
c==== .             / sin(th)  0   -cos(th) \ / e_x \
c==== .          =  |    0     1      0     | | e_y |
c==== .             \ cos(th)  0    sin(th) / \ e_z /
c==== .  
c==== .  polarization
c==== .   incoming:  e_p = cos(phi) e_z' + sin(phi) e_y'
c==== .                  = -cos(phi) cos(th) e_x
c==== .                     + sin(phi) e_y
c==== .                     + cos(phi) sin(th) e_z
c==== .                  = ex1 * e_x + ey1 * e_y' + ez1 * e_z
c==== . 
c==== .   outgoing:  e_1' = e_x' (vtex1)
c==== .              e_2' = e_y' (vtex2)
c==== .  
c==== .   outgoing (45-degree rotated):
c==== .              e_1' = cos(45) * e_x' - sin(45) * e_y' 
c==== .              e_2' = sin(45) * e_x' + cos(45) * e_y'
c==== .  
c==== .  "angle" = the angle of q vector measured from the x (X) axis
c==== .             for the case of Fe-Fe (Fe-As)
c==== .  
c==== .  The polarization vector needs to be rotated by "angle",
c==== .    since it's always defined as perpendicular to q vector.

c      basis = "Fe-Fe"          !representation of orbital
c      basis = "Fe-As"

c==== direction of q vector (kqx, kqy)
      dzero = 1.0d-13
      if (basis(1:5) == 'Fe-Fe') then
         x = DBLE(kqx)          !qx
         y = DBLE(kqy)          !qy
      else                      !change the basis of q
         x = DBLE(kqx) -  DBLE(kqy) !qX 
         y = DBLE(kqx) + DBLE(kqy) !qY
      end if
      if (x > dzero) then
         angle = ATAN(y/x)
      else if (x < -dzero) then
         angle = Pi + ATAN(y/x)
      else if (y > dzero) then
         angle = Pi * 0.5d0    
      else if (y < -dzero) then
         angle = -Pi * 0.5d0
      else
         write(*,*) '---Error in ledgecoeff---'
         write(*,*) ' momentum transfer q is irregal:'
         write(*,*) 'qx_qy = ',x, y
         write(*,*) 'q needs to be non-zero'
         stop
      end if

c==== incoming polarization
c      dummy = Pi / 8.0d0
c      dummy = Pi / 4.0d0
      dummy = 0.0d0

c      phi = pi*0.5d0
      ex1 = -COS(phi) * COS(theta)
      ey1 =  SIN(phi)
      ez1 =  COS(phi) * SIN(theta)
c      write(*,'(a,3f8.5)') 'aaaaaaaaaangle=', ex1,ey1,ez1

      call rotate(ex1,ey1,ez1,'z',angle)


cc      call mkqua(ex1,ey1,ez1,ciqua) !kane
      if (ModeJorbang(1:3) == '1/2') then
         call mkdou(ex1,ey1,ez1,ciqua)
      else 
         call mkqua(ex1,ey1,ez1,ciqua)
      end if
c      write(*,'(a,3f8.5)') 'aaaaaaaaaangle=', ex1,ey1,ez1
c      write(*,*) 'aaaaaaaaaa'
c      stop

c==== outgoing polarization1
c      ex2 = SIN(theta)
c      ey2 = 0.d0
c      ez2 = COS(theta)
!aaaaaaaaaa
      ex2 =  SIN(theta) * COS(dummy)
      ey2 = -SIN(dummy)
      ez2 =  COS(theta) * COS(dummy)
      call rotate(ex2,ey2,ez2,'z',angle)
cc      call mkqua(ex2,ey2,ez2,coqua) !kane

      if (ModeJorbang(1:3) == '1/2') then
         call mkdou(ex2,ey2,ez2,coqua)
      else 
         call mkqua(ex2,ey2,ez2,coqua)
      end if

      nJorbang = 4
      if (ModeJorbang(1:3) == '1/2') nJorbang = 2

cc      call mkvtx(4,ciqua,coqua,vtex1) !kane
      call mkvtx(nJorbang,ciqua,coqua,vtex1)
!aaaaaaa
c      vtex1=0.0d0

c==== outgoing polarization2
c      ex2=0.d0
c      ey2=1.d0
c      ez2=0.d0
!aaaaaaaaaa
      ex2 = SIN(theta) * SIN(dummy)
      ey2 = COS(dummy)
      ez2 = COS(theta) * SIN(dummy)
      call rotate(ex2,ey2,ez2,'z',angle)
cc      call mkqua(ex2,ey2,ez2,coqua)
      if (ModeJorbang == '1/2') then
         call mkdou(ex2,ey2,ez2,coqua)
      else
         call mkqua(ex2,ey2,ez2,coqua)
      end if
      call mkvtx(nJorbang,ciqua,coqua,vtex2)

c==== change the orbital label
      do is1 = 1, 2 ; do is2 = 1, 2 ; do iorb = 1, 5
         call labelorbchange(vtex1(:,is1,iorb,is2))
         call labelorbchange(vtex2(:,is1,iorb,is2))
      end do ; end do ; end do
      do is1 = 1, 2 ; do is2 = 1, 2 ; do iorb = 1, 5
         call labelorbchange(vtex1(iorb,is1,:,is2))
         call labelorbchange(vtex2(iorb,is1,:,is2))
      end do ; end do ; end do


c        if(ABS(theta-pi*0.25)<1.0d-10)
c     &     write(*,*) 
c     &     vtex2,'vtex1aaaaaaaabbb'

      return
      end 

c#########################################################

c#########################################################
c#### labelorbchange:
c#########################################################

      Subroutine labelorbchange(vec)

      use common, only : Pi
      implicit none

      complex*16, intent(inout) :: vec(5)

      complex*16 :: vtemp(5)


      vtemp(:) = vec(:)

c==== 1: d3z^2-r^2
      vec(1) = vtemp(4)
c==== 2: zx
      vec(2) = vtemp(2)
c==== 3: yz
      vec(3) = vtemp(1)
c==== 4: x^2-y^2
      vec(4) = vtemp(5)
c==== 5: xy
      vec(5) = vtemp(3)
      return
      end 

c#########################################################

c
      subroutine mkvtx(nz,ci,co,vtex)
      implicit double precision(a-h,o-z)
      complex(kind=8) ci(5,2,*),co(5,2,*),vtex(5,2,5,2)
      vtex(1:5,1:2,1:5,1:2)=cmplx(0.d0,0.d0,kind=8)
      do n=1,nz
c      n=1
        do j2=1,2
          do j1=1,5
            do i2=1,2
              do i1=1,5
                vtex(i1,i2,j1,j2)=vtex(i1,i2,j1,j2)
     &             +co(i1,i2,n)*dconjg(ci(j1,j2,n))
              enddo
            enddo
          enddo
        enddo
      enddo
      return
      end
c
      subroutine mkqua(ex,ey,ez,cqua)
      implicit double precision(a-h,o-z)
      complex(kind=8) cqua(5,2,4),ai
c  cqua(d-orbital,spin,jz)
c
c d-orbital:
c 1: yz
c 2: zx
c 3: xy
c 4: 3z2-r2
c 5: x2-y2
c
c spin
c 1:up
c 2:down
c
c jz
c   cqua cdou
c 1: 3/2 1/2
c 2: 1/2 -1/2
c 3: -1/2
c 4: -3/2
c
c angle of the polarization: theta,phi
c      pi=datan(1.d0)*4.d0
      ai=cmplx(0.d0,1.d0,kind=8)
c      theta=30.d0/180.d0*pi
c      phi=0.d0/180.d0*pi
c      ex=dsin(theta)*dcos(phi)
c      ey=dsin(theta)*dsin(phi)
c      ez=dcos(theta)
c
c j=3/2
      rten=1.d0/dsqrt(10.d0)
      rthi=1.d0/dsqrt(30.d0)
      cqua(1:5,1:2,1:4)=cmplx(0.d0,0.d0,kind=8)
c jz=3/2
      cqua(1,1,1)=-ai*ez*rten
      cqua(2,1,1)=-ez*rten
      cqua(3,1,1)=(-ai*ex-ey)*rten
      cqua(4,1,1)=(ex+ai*ey)*rthi
      cqua(5,1,1)=(-ex+ai*ey)*rten
c jz=-3/2
      cqua(1,2,4)=-ai*ez*rten
      cqua(2,2,4)=ez*rten
      cqua(3,2,4)=(-ai*ex+ey)*rten
      cqua(4,2,4)=(-ex+ai*ey)*rthi
      cqua(5,2,4)=(ex+ai*ey)*rten
c jz=1/2
      cqua(1,1,2)=2*ey*rthi
      cqua(2,1,2)=2*ex*rthi
      cqua(4,1,2)=4/dsqrt(3.d0)*ez*rthi
      cqua(1,2,2)=-ai*ez*rthi
      cqua(2,2,2)=-ez*rthi
      cqua(3,2,2)=(-ai*ex-ey)*rthi
      cqua(4,2,2)=(ex+ai*ey)*rten/3.d0
      cqua(5,2,2)=(-ex+ai*ey)*rthi
c jz=-1/2
      cqua(1,1,3)=-ai*ez*rthi
      cqua(2,1,3)=ez*rthi
      cqua(3,1,3)=(-ai*ex+ey)*rthi
      cqua(4,1,3)=(-ex+ai*ey)*rten/3.d0
      cqua(5,1,3)=(ex+ai*ey)*rthi
      cqua(1,2,3)=2.d0*ey*rthi
      cqua(2,2,3)=2.d0*ex*rthi
      cqua(4,2,3)=4/dsqrt(3.d0)*ez*rthi
c
      return
      end
c
      subroutine mkdou(ex,ey,ez,cdou)
      implicit double precision(a-h,o-z)
      complex(kind=8) cdou(5,2,2),ai
c  cdou(d-orbital,spin,jz)
c
c
c spin
c 1:up
c 2:down
c
c jz
c   cqua cdou
c 1: 3/2 1/2
c 2: 1/2 -1/2
c 3: -1/2
c 4: -3/2
c
c angle of the polarization: theta,phi
c      pi=datan(1.d0)*4.d0
      ai=cmplx(0.d0,1.d0,kind=8)
c      theta=30.d0/180.d0*pi
c      phi=0.d0/180.d0*pi
c      ex=dsin(theta)*dcos(phi)
c      ey=dsin(theta)*dsin(phi)
c      ez=dcos(theta)
      cdou(1:5,1:2,1:2)=cmplx(0.d0,0.d0,kind=8)
c
c j=1/2
      rfif=1.d0/dsqrt(15.d0)
c jz=1/2
      cdou(1,1,1)=-ey*rfif
      cdou(2,1,1)=-ex*rfif
      cdou(4,1,1)=-2.d0*ez*rfif/dsqrt(3.d0)
      cdou(1,2,1)=-ai*ez*rfif
      cdou(2,2,1)=-ez*rfif
      cdou(3,2,1)=(-ai*ex-ey)*rfif
      cdou(4,2,1)=(ex+ai*ey)*rfif/dsqrt(3.d0)
      cdou(5,2,1)=(-ex+ai*ey)*rfif
c jz=-1/2
      cdou(1,1,2)=ai*ez*rfif
      cdou(2,1,2)=-ez*rfif
      cdou(3,1,2)=(ai*ex-ey)*rfif
      cdou(4,1,2)=(ex-ai*ey)*rfif/dsqrt(3.d0)
      cdou(5,1,2)=(-ex-ai*ey)*rfif
      cdou(1,2,2)=ey*rfif
      cdou(2,2,2)=ex*rfif
      cdou(4,2,2)=2.d0*ez*rfif/dsqrt(3.d0)
c
      return
      end
c###########################################################
      subroutine vertex
      implicit double precision(a-h,o-z)
      character :: basis*5
      complex(kind=8) ciqua(5,2,4),cidou(5,2,2)
      complex(kind=8) coqua(5,2,4),codou(5,2,2)
      complex(kind=8) vtex1(5,2,5,2),vtex2(5,2,5,2)
      dimension vin(5,2)

      integer :: kqx,kqy

c====. vtex1(l,k,lgs,lspn): vertex
c====. At the vertex,
c====. add an electron at (lgs, lspn)
c====. then create a hole at (l, k) 
c====. lgs, l: orbital, lspn, k: spin
c
c  cqua(d-orbital,spin,jz)
c  cdou(d-orbital,spin,jz)
c
c d-orbital:
c 1: 3z2-r2
c 2: zx
c 3: yz
c 4: x2-y2
c 5: xy
c
c spin
c 1:up
c 2:down
c
c jz
c   cqua cdou
c 1: 3/2 1/2
c 2: 1/2 -1/2
c 3: -1/2
c 4: -3/2
c

      basis = 'Fe-Fe'

      pi = datan(1.d0) * 4.d0

      kqx = 1
      kqy = 0
c
c
c
c      goto 100
      open(20,file = 'out.vertex.dat')

      idummy = 0

      lgs=3
      lspn=1
      write(20,1100) '# phi theta  up of yz zx xy 3z2 x2-y2'//
     &   ', down of yz zx xy 3z2 x2-y2'
 1100  format(a37,a30)
      phi=0.d0
      do i=0,90
        theta=dfloat(90-i)/180.d0*pi
c generate cqua
        call ledgecoeff(theta,phi,vtex1,vtex2,kqx,kqy,basis)
        do k=1,2
          do l=1,5
            vin(l,k)=vtex1(l,k,lgs,lspn)*dconjg(vtex1(l,k,lgs,lspn))
     &         +vtex2(l,k,lgs,lspn)*dconjg(vtex2(l,k,lgs,lspn))
c            vin(l,k)=(vtex1(l,k,5,1)+vtex2(l,k,5,1))
c     &         *dconjg(vtex1(l,k,5,1)+vtex2(l,k,5,1))
          enddo
        enddo
      write(20,1000) i + idummy
     &     ,(vin(l,1),l=1,5),(vin(l,2),l=1,5)
 1000   format(1i5,10e15.4)
      enddo
c
      idummy = idummy + 90
 100  continue
      theta=0.d0
      do i=1,90
        phi=dfloat(i)/180.d0*pi
        call ledgecoeff(theta,phi,vtex1,vtex2,kqx,kqy,basis)
        do k=1,2
          do l=1,5
            vin(l,k)=vtex1(l,k,lgs,lspn)*dconjg(vtex1(l,k,lgs,lspn))
     &         +vtex2(l,k,lgs,lspn)*dconjg(vtex2(l,k,lgs,lspn))
c            vin(l,k)=(vtex1(l,k,5,1)+vtex2(l,k,5,1))
c     &         *dconjg(vtex1(l,k,5,1)+vtex2(l,k,5,1))
          enddo
        enddo
      write(20,1000) i+idummy
     &     ,(vin(l,1),l=1,5),(vin(l,2),l=1,5)
      enddo
c
      idummy = idummy + 90
      phi=90.d0/180.d0*pi
      do i=1,90
        theta=dfloat(i)/180.d0*pi
        call ledgecoeff(theta,phi,vtex1,vtex2,kqx,kqy,basis)
        do k=1,2
          do l=1,5
            vin(l,k)=vtex1(l,k,lgs,lspn)*dconjg(vtex1(l,k,lgs,lspn))
     &         +vtex2(l,k,lgs,lspn)*dconjg(vtex2(l,k,lgs,lspn))
c            vin(l,k)=(vtex1(l,k,5,1)+vtex2(l,k,5,1))
c     &         *dconjg(vtex1(l,k,5,1)+vtex2(l,k,5,1))
          enddo
        enddo
      write(20,1000) i + idummy
     &     ,(vin(l,1),l=1,5),(vin(l,2),l=1,5)
      enddo
c
      idummy = idummy + 90
      theta=90.d0/180.d0*pi
      do i=1,90
        phi=dfloat(90-i)/180.d0*pi
        call ledgecoeff(theta,phi,vtex1,vtex2,kqx,kqy,basis)
        do k=1,2
          do l=1,5
            vin(l,k)=vtex1(l,k,lgs,lspn)*dconjg(vtex1(l,k,lgs,lspn))
     &         +vtex2(l,k,lgs,lspn)*dconjg(vtex2(l,k,lgs,lspn))
c            vin(l,k)=(vtex1(l,k,5,1)+vtex2(l,k,5,1))
c     &         *dconjg(vtex1(l,k,5,1)+vtex2(l,k,5,1))
          enddo
        enddo
      write(20,1000) i+idummy
     &     ,(vin(l,1),l=1,5),(vin(l,2),l=1,5)
      enddo
c
      return
      end


c#########################################################
c#### Iorbcomb : 
c#########################################################
      character*1 function orbselect(mu1,mu2,mu3,mu4,nu1,nu2,nu3,nu4)

      use common, only : Ns
      implicit none
      
      integer, intent(in) :: mu1, mu2, mu3, mu4, nu1, nu2, nu3, nu4

      orbselect = 'y'

         
      if ((mu1 /= nu1).and.(nu1 > 0)) then
         orbselect = 'n'
         return
      end if
      if ((mu1 == ABS(nu1)).and.(nu1 < 0)) then
         orbselect = 'n'
         return
      end if

      if ((mu2 /= nu2).and.(nu2 > 0)) then
         orbselect = 'n'
         return
      end if
      if ((mu2 == ABS(nu2)).and.(nu2 < 0)) then
         orbselect = 'n'
         return
      end if

      if ((mu3 /= nu3).and.(nu3 > 0)) then
         orbselect = 'n'
         return
      end if
      if ((mu3 == ABS(nu3)).and.(nu3 < 0)) then
         orbselect = 'n'
         return
      end if

      if ((mu4 /= nu4).and.(nu4 > 0)) then
         orbselect = 'n'
         return
      end if
      if ((mu4 == nu4).and.(nu4 < 0)) then
         orbselect = 'n'
         return
      end if


      return
      end
c#########################################################
      

c#########################################################
c#### rotate : 
c#########################################################
      subroutine rotate(vx,vy,vz,xyz,theta)

      use common, only : Ns
      implicit none
      
      character, intent(in) :: xyz*1 
      real*8, intent(in) :: theta
      real*8, intent(inout) :: vx, vy, vz
      real*8 :: ax, ay, az, bx, by, bz
c==== .
c==== .   rotate (ax,ay,az) to (bx,by,bz)
c==== .          by theta with respect to x, y, or z axis
c==== .       = (-theta) rotation of the basis (ax,ay) to (bx,by)
c==== .   ey
c==== .   |  /ey'           
c==== .   | /               
c==== .   |/                     b    theta from a to b  
c==== .   +-----ex             /|\ _                     
c==== .    \ )  -theta          |  /\ a                  
c==== .      \                  | /                      
c==== .        ex'              |/                         
c==== . 
c==== .  rotation of basis = rotation of vector
c==== .         by (-theta)            by theta

      if (xyz == 'x') then
         ax = vy
         ay = vz
         az = vx
      else if (xyz == 'y') then
         ax = vz
         ay = vx
         az = vy
      else if (xyz == 'z') then
         ax = vx
         ay = vy
         az = vz
      else
         write(*,*)'Error in rotate'
         stop
      end if
      
      bx =  ax * COS(theta) - ay * SIN(theta) 
      by =  ax * SIN(theta) + ay * COS(theta)
      bz =  az

      if (xyz == 'x') then
         vy = bx
         vz = by
         vx = bz
      else if (xyz == 'y') then
         vz = bx
         vx = by
         vy = bz
      else if (xyz == 'z') then
         vx = bx
         vy = by
         vz = bz
      else
         write(*,*)'Error in rotate'
         stop
      end if
      


      return
      end
c#########################################################
      

      integer function Ispinflip(is)

      implicit none
      
      integer, intent(in) :: is

      if (is == 1) then
         Ispinflip = 2
      else if (is == 2) then
         Ispinflip = 1
      else
         write(*,*) 'Error in Ispinflip'
         stop
      end if
      return
      end
c#########################################################
c
      subroutine mkctx(vtex,ctex)
      implicit double precision(a-h,o-z)
      complex(kind=8) vtex(5,2,5,2),ctex(5,5,4),ai
      ai=cmplx(0.d0,1.d0,kind=8)
c
c  sum_s,s' v(a,s,b,s') a^+_s b_s'
c = ( 1/2 v(a,u,b,u) + 1/2 v(a,d,b,d) )(a^+_u b_u + a^+_d b_d)
c + ( v(a,u,b,u) - v(a,d,b,d) ) 1/2(a^+_u b_u - a^+_d b_d)
c + ( v(a,u,b,d) + v(a,d,b,u) ) 1/2(a^+_u b_d + a^+_d b_u)
c + i( v(a,u,b,d) - v(a,d,b,u) ) 1/(2i)(a^+_u b_d - a^+_d b_u)
c
c = c(a,b,n) (a^+_u b_u + a^+_d b_d)
c + c(a,b,z) 1/2(a^+_u b_u - a^+_d b_d)
c + c(a,b,x) 1/2(a^+_u b_d + a^+_d b_u)
c + c(a,b,y) 1/(2i)(a^+_u b_d - a^+_d b_u)
c
      ctex(1:5,1:5,1:4)=cmplx(0.d0,0.d0,kind=8)
c
      do j=1,5
        do i=1,5
          ctex(i,j,1)=ctex(i,j,1)+0.5d0*(vtex(i,1,j,1)+vtex(i,2,j,2))
          ctex(i,j,2)=ctex(i,j,2)+vtex(i,1,j,1)-vtex(i,2,j,2)
          ctex(i,j,3)=ctex(i,j,3)+vtex(i,1,j,2)+vtex(i,2,j,1)
          ctex(i,j,4)=ctex(i,j,4)+ai*(vtex(i,1,j,2)-vtex(i,2,j,1))
        enddo
      enddo
      return
      end

      subroutine azimuth
c try to make coef. 
      implicit double precision(a-h,o-z)
      complex(kind=8) ciqua(5,2,4),cidou(5,2,2)
      complex(kind=8) coqua(5,2,4),codou(5,2,2)
      complex(kind=8) vtex1(5,2,5,2),vtex2(5,2,5,2)
      complex(kind=8) ctex1(5,5,4),ctex2(5,5,4)
      character filename*6
c
c ctex(alpha,beta,1):charge
c ctex(alpha,beta,2):sz
c ctex(alpha,beta,3):sx
c ctex(alpha,beta,4):sy
c
c  sum_s,s' v(a,s,b,s') a^+_s b_s'
c = ( 1/2 v(a,u,b,u) + 1/2 v(a,d,b,d) )(a^+_u b_u + a^+_d b_d)
c + ( v(a,u,b,u) - v(a,d,b,d) ) 1/2(a^+_u b_u - a^+_d b_d)
c + ( v(a,u,b,d) + v(a,d,b,u) ) 1/2(a^+_u b_d + a^+_d b_u)
c + i( v(a,u,b,d) - v(a,d,b,u) ) 1/(2i)(a^+_u b_d - a^+_d b_u)
c
c = c(a,b,n) (a^+_u b_u + a^+_d b_d)
c + c(a,b,z) 1/2(a^+_u b_u - a^+_d b_d)
c + c(a,b,x) 1/2(a^+_u b_d + a^+_d b_u)
c + c(a,b,y) 1/(2i)(a^+_u b_d - a^+_d b_u)
c
c
      dimension vin(5,2)
      dimension cin(5,5,4)
c
c  cqua(d-orbital,spin,jz)
c  cdou(d-orbital,spin,jz)
c
c d-orbital:
c 1: yz
c 2: zx
c 3: xy
c 4: 3z2-r2
c 5: x2-y2
c
c spin
c 1:up
c 2:down
c
c jz
c   cqua cdou
c 1: 3/2 1/2
c 2: 1/2 -1/2
c 3: -1/2
c 4: -3/2
c
      pi=datan(1.d0)*4.d0
      do m=1,5 ; do n = 0,3
         write(filename,'(a,i2)')'ctex',10*m+n
         open (unit=10*m+n, file=filename)
      end do ; end do
c
c
c
c      goto 100
      do m=1,5
      write(10*m,1100) '# theta psi charge of yz zx xy 3z2 x2-y2'
      write(10*m+1,1100) '# theta psi sz of yz zx xy 3z2 x2-y2'
      write(10*m+2,1100) '# theta psi sx of yz zx xy 3z2 x2-y2'
      write(10*m+3,1100) '# theta psi sy of yz zx xy 3z2 x2-y2'
      enddo
 1100  format(a37,a30)
c
c phi= 0:    pi-polarization
c phi=90: sigma-polarization
c
c      phi=0.d0/180.d0*pi
      phi=90.d0/180.d0*pi
      do m=1,5
      do k=0,3
        write(10*m+k,1200) '#  phi=',phi/pi*180.d0
      enddo
      enddo
 1200   format(a9,f20.13)
c
      do i=0,90
c        i=40
        theta=dfloat(i)/180.d0*pi
c
c  psi: azimuth angle from x-axis
c j=0: (pi,0)direction
c j=45: (pi,pi)direction
c      do j=0,360
        j=0
c        j=45
        psi=dfloat(j)/180.d0*pi
c generate cqua
        ex1=-dcos(phi)*dcos(theta)*dcos(psi)-dsin(theta)*dsin(psi)
        ey1=-dcos(phi)*dcos(theta)*dsin(psi)+dsin(phi)*dcos(psi)
        ez1=dcos(phi)*dsin(theta)
        call mkqua(ex1,ey1,ez1,ciqua)
        ex2=dsin(theta)*dcos(psi)
        ey2=dsin(theta)*dsin(psi)
        ez2=dcos(theta)
        call mkqua(ex2,ey2,ez2,coqua)
        call mkvtx(4,ciqua,coqua,vtex1)
        call mkctx(vtex1,ctex1)
        ex2=-dsin(psi)
        ey2=dcos(psi)
        ez2=0.d0
        call mkqua(ex2,ey2,ez2,coqua)
        call mkvtx(4,ciqua,coqua,vtex2)
        call mkctx(vtex2,ctex2)
c        do k=1,4
c          do l=1,5
c            vin(l,k)=ctex1(l,l,k)*dconjg(ctex1(l,l,k))
c     &         +ctex2(l,l,k)*dconjg(ctex2(l,l,k))
c          enddo
c        enddo
        do k=1,4
          do m=1,5
          do l=1,5
            cin(l,m,k)=ctex1(l,m,k)*dconjg(ctex1(l,m,k))
     &         +ctex2(l,m,k)*dconjg(ctex2(l,m,k))
          enddo
          enddo
        enddo
c
      do m=1,5
      write(10*m,1000) theta/pi*180.d0,psi/pi*180.d0
     &     ,(cin(l,m,1),l=1,5)
      write(10*m+1,1000) theta/pi*180.d0,psi/pi*180.d0
     &     ,(cin(l,m,2),l=1,5)
      write(10*m+2,1000) theta/pi*180.d0,psi/pi*180.d0
     &     ,(cin(l,m,3),l=1,5)
      write(10*m+3,1000) theta/pi*180.d0,psi/pi*180.d0
     &     ,(cin(l,m,4),l=1,5)
      enddo
 1000   format(2f6.1,20e15.4)
      enddo
c
      do m=1,5 ; do n = 1,5
         close (unit=10*m+n)
      end do ; end do

      return
      stop
      end


      subroutine azimuth2
c try to make coef. 
      implicit double precision(a-h,o-z)
      complex(kind=8) ciqua(5,2,4),cidou(5,2,2)
      complex(kind=8) coqua(5,2,4),codou(5,2,2)
      complex(kind=8) vtex1(5,2,5,2),vtex2(5,2,5,2)
      complex(kind=8) ctex1(5,5,4),ctex2(5,5,4)
      character filename*6
c
      dimension vin(5,2)
      dimension cin(5,5,4)

      pi=datan(1.d0)*4.d0
      do m=1,5 ; do n = 0,3
         write(filename,'(a,i2)')'ctex',10*m+n
         open (unit=10*m+n, file=filename)
      end do ; end do
c
c
c
c      goto 100
      do m=1,5
      write(10*m,1100) '# theta psi charge of yz zx xy 3z2 x2-y2'
      write(10*m+1,1100) '# theta psi sz of yz zx xy 3z2 x2-y2'
      write(10*m+2,1100) '# theta psi sx of yz zx xy 3z2 x2-y2'
      write(10*m+3,1100) '# theta psi sy of yz zx xy 3z2 x2-y2'
      enddo
 1100  format(a37,a30)
c
c phi= 0:    pi-polarization
c phi=90: sigma-polarization
c
      phi=0.d0/180.d0*pi
c      phi=90.d0/180.d0*pi
      do m=1,5
      do k=0,3
        write(10*m+k,1200) '#  phi=',phi/pi*180.d0
      enddo
      enddo
 1200   format(a9,f20.13)
c
      do i=0,90
c        i=40
        theta=dfloat(i)/180.d0*pi
c
c  psi: azimuth angle from x-axis
c j=0: (pi,0)direction
c j=45: (pi,pi)direction
c      do j=0,360
        j=0
c        j=45
        psi=dfloat(j)/180.d0*pi
c generate cqua

c==== incoming polarization
      ex1 = -COS(phi) * COS(theta)
      ey1 =  SIN(phi)
      ez1 =  COS(phi) * SIN(theta)
      call rotate(ex1,ey1,ez1,'z',psi)

c        ex1=-dcos(phi)*dcos(theta)*dcos(psi)-dsin(theta)*dsin(psi)
c        ey1=-dcos(phi)*dcos(theta)*dsin(psi)+dsin(phi)*dcos(psi)
c        ez1=dcos(phi)*dsin(theta)
        call mkqua(ex1,ey1,ez1,ciqua)

        dummy = 0.0d0
        ex2=dsin(theta)*dcos(psi)
        ey2=dsin(theta)*dsin(psi)
        ez2=dcos(theta)
        ex2 =  SIN(theta) * COS(dummy)
        ey2 = -SIN(dummy)
        ez2 =  COS(theta) * COS(dummy)
        call rotate(ex2,ey2,ez2,'z',psi)

        call mkqua(ex2,ey2,ez2,coqua)
        call mkvtx(4,ciqua,coqua,vtex1)


        ex2=-dsin(psi)
        ey2=dcos(psi)
        ez2=0.d0
      ex2 = SIN(theta) * SIN(dummy)
      ey2 = COS(dummy)
      ez2 = COS(theta) * SIN(dummy)
      call rotate(ex2,ey2,ez2,'z',psi)
        call mkqua(ex2,ey2,ez2,coqua)
        call mkvtx(4,ciqua,coqua,vtex2)


c==== change the orbital label
      do is1 = 1, 2 ; do is2 = 1, 2 ; do iorb = 1, 5
         call labelorbchange(vtex1(:,is1,iorb,is2))
         call labelorbchange(vtex2(:,is1,iorb,is2))
      end do ; end do ; end do
      do is1 = 1, 2 ; do is2 = 1, 2 ; do iorb = 1, 5
         call labelorbchange(vtex1(iorb,is1,:,is2))
         call labelorbchange(vtex2(iorb,is1,:,is2))
      end do ; end do ; end do

c        if(ABS(theta-pi*0.25)<1.0d-10)
c     &     write(*,*) 
c     &     vtex2,'vtex1aaaaaaaabbb'

        call mkctx(vtex1,ctex1)
        call mkctx(vtex2,ctex2)


        do k=1,4
          do m=1,5
          do l=1,5
            cin(l,m,k)=ctex1(l,m,k)*dconjg(ctex1(l,m,k))
     &         +ctex2(l,m,k)*dconjg(ctex2(l,m,k))
          enddo
          enddo
        enddo
        do k=1,4
          do m=1,5
          do l=1,5
             if ( ABS(cin(l,m,k)) < 1.0d-16) cin(l,m,k)=0.0d0
          enddo
          enddo
        enddo
c
      do m=1,5
      write(10*m,1000) theta/pi*180.d0,psi/pi*180.d0
     &     ,(cin(l,m,1),l=1,5)
      write(10*m+1,1000) theta/pi*180.d0,psi/pi*180.d0
     &     ,(cin(l,m,2),l=1,5)
      write(10*m+2,1000) theta/pi*180.d0,psi/pi*180.d0
     &     ,(cin(l,m,3),l=1,5)
      write(10*m+3,1000) theta/pi*180.d0,psi/pi*180.d0
     &     ,(cin(l,m,4),l=1,5)
      enddo
 1000   format(2f6.1,20e15.4)
      enddo
c
      do m=1,5 ; do n = 1,5
         close (unit=10*m+n)
      end do ; end do

      return
      stop
      end
