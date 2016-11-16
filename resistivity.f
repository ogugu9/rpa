c#########################################################
c#### resistivity                EK 2015
c#########################################################

      subroutine resistivity(xory,rho,finputvf)

      use common, only : Ns, Ns2, Pi,
     &     Eall, Zpsiall, Dmu, Dne, Dtemp, Nqx, Nqy, Dkshift,
     &     isite, isitex, isitey, ZI, t, lln, Erng, Deta, Maxom
      implicit none

      character(len=*), intent(in) :: finputvf
      character, intent(in) :: xory*1
      real*8, intent(out) :: rho

      character :: chara*80, filein*40
      integer :: nkf, ispin, nksize(2), i
      integer, allocatable :: ibandf(:)
      real*8 :: dummy(10), avec(2)
      real*8 :: ax, ay, dImMemfn, sigma0, drude
      real*8, allocatable :: dkf(:,:), dlenkf(:)
      real*8, allocatable :: vf(:,:), vfermi(:)
      real*8, parameter :: hbarinv = 1.0d0 !6.58244d16
      real*8, parameter :: ech = 1.0d0 !1.60217733d-19
      integer, external :: IcheckZone

c==== .    Drude と Memory fn から rho を計算．
c==== .    Drude の計算には vF が必要．             
c==== .    M の計算には kF の情報が必要．
c==== .    vF と kF の情報は fermi_KfandVf で計算したものを使う． 
c==== . 

c## polarization[
      call polarization(xory,avec)
c## polarization]
      
c## read k_fermi[
      filein = finputvf
      if(filein == '#') then
         write(*,*) 'Enter input filename (fermi velocity data)'
         read(5,*) filein
         if (filein(1:1) == '!') filein = 'v-out.fs.dat'
      end if

      open(50,file = filein)    !open data file
      write(*,*) 'vf data read from ',filein

      chara = '#####'
      do while (chara(1:5) /= '#size')
         read(50,*) chara
      end do
      backspace(50)
      read(50,*) chara, nksize(1:2)

      do while (chara(1:10) /= '#Beginning')
         read(50,*) chara
      end do
      nkf = 0
      do while (chara(1:4) /= '#EOF')
         read(50,*) dummy(1:6), ispin
         if (ispin == 1) nkf = nkf + 1
         read(50,*) chara
         backspace(50)
      end do
      rewind(50)
      write(*,*) 'nkf=',nkf

      allocate(dkf(nkf,2), vf(nkf,2), vfermi(nkf), dlenkf(nkf))
      allocate(ibandf(nkf))
      do while (chara(1:10) /= '#Beginning')
         read(50,*) chara
      end do

      do i = 1, nkf
         read(50,*) dkf(i,1:2), vf(i,1:2), 
     &        dlenkf(i),ibandf(i),dummy(1)
      end do
      vfermi(:) = SQRT(vf(:,1)**2 + vf(:,2)**2)

      close(50)
c## read k_fermi]

      write(*,*) 'calculating Drude...'
      call Drude4(xory,sigma0,filein)

      write(*,*) 'calculating memory function...'
      call memoryfn(avec,nksize,dkf,dlenkf,vfermi,ibandf,nkf,
     &     dImMemfn,drude)
      
c      rho = dImMemfn / drude**2 / 2.0d0 / hbarinv !error
      rho = dImMemfn / drude / 2.0d0 / hbarinv !corrected

      write(*,'(" sigma0_",a1,"= ",g20.10)') xory, sigma0
      write(*,'(" sigma0_",a1,"= ",g20.10)') xory, drude
      write(*,'(" ImM_",a1,"= ",g20.10)') xory, dImMemfn
      write(*,'(" rho_",a1,"= ",g20.10)') xory,  rho

      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### memoryfn: SUM[ Im{M} ]                EK 2015
c#########################################################

      subroutine memoryfn(avec,nksize,dkf,dl,vf,ibandf,nkf,
     &     dImMemfn, drude)

      use common, only : Ns, Nqx, Pi
      implicit none

      integer, intent(in) :: nksize(2)
      integer, intent(in) :: nkf
      integer, intent(in) :: ibandf(nkf)
      real*8, intent(in) :: dkf(nkf,2), vf(nkf), dl(nkf), avec(2)
      real*8, intent(out) :: dImMemfn, drude

      integer :: ie1, ie2, ik1, ik2, ispin, nnn
      real*8 :: adelta, akx, aky, pi2, cimp, dcir
      real*8 :: eig(Ns*Nqx)
      real*8 :: dImMemk(nkf)
      complex*16 :: zpsi1(Ns*Nqx,2), zpsi2(Ns*Nqx,2), zpsi3(Ns*Nqx,2)
      complex*16 :: zpsikf(Ns*Nqx,Ns*Nqx,nkf)
      complex*16 :: zAmtrx
      complex*16 :: zImtrx
      complex*16 :: zJmtrx(2)

      pi2 = Pi * 2.0d0
      nnn = 1 ! nksize(1) * nksize(2)
      dcir = SUM(dl(1:nkf))**2
      ispin = 1 
      cimp = 1.0d0              !0.01d0

      if (nkf < 1) return

!$omp parallel do private(ik1,akx,aky,eig) shared(zpsikf,dkf)
      do ik1 = 1, nkf
         akx = dkf(ik1,1) * pi2
         aky = dkf(ik1,2) * pi2
         call  diag(eig(:),zpsikf(:,:,ik1),akx,aky,ispin)
      end do
!$omp end parallel do

      dImMemk(:) = 0.0d0
      dImMemfn = 0.0d0
      drude = 0.0d0
c!$omp parallel do private(ik1,ik2,zpsi1,zpsi2,
c     &    zJmtrx,zImtrx,zAmtrx,ie1,ie2)
c     &     shared(dImMemfn,zpsikf,dkf,drude,ibandf,vf,dl,avec)
      do ik1 = 1, nkf
         ie1 = ibandf(ik1)
         zpsi1(:,1) = zpsikf(:,ie1,ik1)
         zpsi1(:,2) = zpsi1(:,1)
         call currentmtrx(avec,dkf(ik1,:),zpsi1(:,:),zJmtrx(1))

         drude = drude
     &        + ABS(zJmtrx(1))**2 / vf(ik1) * dl(ik1) / dcir * 2.0d0

         do ik2 = 1, nkf
            ie2 = ibandf(ik2)
            zpsi2(:,1) = zpsikf(:,ie2,ik2)
            zpsi2(:,2) = zpsi2(:,1)
            call currentmtrx(avec,dkf(ik2,:),zpsi2(:,:),zJmtrx(2))
            zpsi3(:,1) = zpsi1(:,1)
            zpsi3(:,2) = zpsi2(:,2)
            call impuritymtrx(zpsi3(:,:),zImtrx)
c            call calc_Amtrx(zJmtrx,zImtrx,zAmtrx,ie1,ie2)

            zAmtrx =  zJmtrx(1) * zImtrx - zImtrx * zJmtrx(2)

            dImMemk(ik1) = dImMemk(ik1)
     &           + ABS(zAmtrx)**2 / vf(ik2) * dl(ik2) / dcir
         end do 
      end do
c!$omp end parallel do
      dImMemfn = SUM(dImMemk(:) / vf(:) * dl(:) )
      dImMemfn = dImMemfn * Pi * cimp / dcir

      open(10,file='out.memfn.dat')
      do ik1 = 1, nkf
         write(10,'(3g20.10)')dkf(ik1,1),dkf(ik1,2),
     &        dImMemk(ik1)! / vf(ik1) * dl(ik1) / dcir
      end do
      close(10)

      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### calc_Amtrx: matrix elements of A (A = [j,H])
c####            by   EK 2015
c#########################################################

      subroutine calc_Amtrx(zJmtrx,zImtrx,zAmtrx,ibandf1,ibandf2)

      use common, only : Ns, Nqx
      implicit none

      integer, intent(in) :: ibandf1, ibandf2
      complex*16, intent(in) :: zJmtrx(2)
      complex*16, intent(in) :: zImtrx(Ns*Nqx,Ns*Nqx)
      complex*16, intent(out) :: zAmtrx(Ns*Nqx,Ns*Nqx)
      integer :: iband1, iband2

      zAmtrx(:,:) = 0.0d0
c      do iband1 = 1, Ns * Nqx ; do iband2 = 1, Ns * Nqx
c         if (iband1 /= ibandf1) cycle
c         if (iband2 /= ibandf2) cycle
         zAmtrx(ibandf1,ibandf2) = 
     &        zJmtrx(1) * zImtrx(ibandf2,ibandf2)
     &        - zImtrx(ibandf1,ibandf2) * zJmtrx(2)
c      end do ; end do

      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### impuritymtrx: matrix elements of Iimp (l=0)
c####           in band-quasiparticle representation
c####            by   EK 2015
c#########################################################

      subroutine impuritymtrx(zpsi,zImtrx)

      use common, only : Ns, Nqx
      implicit none

      complex*16, intent(in) :: zpsi(Ns*Nqx,2)
c      complex*16, intent(in) :: zpsi2(Ns*Nqx)
      complex*16, intent(out) :: zImtrx
      integer :: mu, nu, iQ1, iQ2
      integer :: iorb1, iorb2
      real*8 :: dIimp

c==== .    H' = I_imp SUM_{l,ispin,alpha} c^+ c             
c==== .       = (1/N) SUM Imtrx(k0,ie;k0',ie') gamma^+ gamma
c==== .    See Eq.(12) in PRB vol.90, 125157(2014)
c==== .
c==== .    Imtrx = I_imp SUM_{mu,iQ,iQ'} psi(k0)^* psi(k0') 
c==== .         (Here, we consider the l = 0 case)      
c==== .    See Eq.(13) in PRB vol.90, 125157(2014)        
c==== . 

      dIimp = 1.0d0
      zImtrx = 0.0d0
      do mu = 1, Ns ; do iQ1 = 0, Nqx-1; do iQ2 = 0, Nqx-1
         iorb1 = mu + iQ1 * Ns
         iorb2 = mu + iQ2 * Ns
         zImtrx = zImtrx + zpsi(iorb1,1) * CONJG(zpsi(iorb2,2))
      end do ; end do ; end do

      zImtrx = zImtrx * dIimp

      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### currentmtrx: matrix elements of J (ie=ie')
c####           in band-quasiparticle representation
c####            by   EK 2015
c#########################################################

      subroutine currentmtrx(avec,dkvec0,zpsi,zJmtrx)

      use common, only : Ns, Nqx, Nqy, t, lln, isitex, isitey, Pi, ZI
      implicit none

      real*8, intent(in) :: dkvec0(2), avec(2)
      complex*16, intent(in) :: zpsi(Ns*Nqx,2)
      complex*16, intent(out) :: zJmtrx

      integer :: iband, mu, nu, iQ, ix, iy, is
      integer :: iorb1, iorb2, ie
      real*8 :: akx, aky, pi2, adelta
      real*8, parameter :: ddd = 1.0d-15 ! 6.58244d16
      real*8, parameter :: hbarinv = 1.0d0 ! 6.58244d16
      real*8, parameter :: ech = 1.0d0 ! 1.60217733d-19
      complex*16 :: zcoeff

c==== .    j = SUM Jmtrx(ie;ie') gamma^+ gamma
c==== .    See Eq.(5) in PRB vol.90, 125157(2014)
c==== .
c==== .    Jmtrx = (i/N)(e/hbar) 
c==== .            * SUM[ Delta * t * EXP{-ZI*k*Delta} psi^* psi ] 
c==== .         (Here, we consider the ie = ie'case)      
c==== .    See Eq.(6) in PRB vol.90, 125157(2014)        
c==== . 


      pi2 = Pi * 2.0d0
      zJmtrx = 0.0d0

      do is = 1, lln - 1
         ix = isitex(is)
         iy = isitey(is)
         
         adelta = DBLE(ix) * avec(1) + DBLE(iy) * avec(2)
         if (ABS(adelta) < ddd) cycle  
         do mu = 1, Ns ; do nu = 1, Ns
            do iQ = 0, Nqx-1
               iorb1 = mu + iQ * Ns
               iorb2 = nu + iQ * Ns
               akx = (dkvec0(1) + DBLE(iQ)  / DBLE(Nqx))
               aky = (dkvec0(2) + DBLE(iQ)  / DBLE(Nqy))
               akx = akx * DBLE(ix) * pi2
               aky = aky * DBLE(iy) * pi2

               zJmtrx = zJmtrx
     &              + zpsi(iorb1,1) * CONJG(zpsi(iorb2,2))
     &              * EXP(-ZI * (akx + aky))
     &              * adelta * t(ix,iy,mu,nu)
            end do
         end do ; end do
      end do
      zcoeff = ZI
      zcoeff = zcoeff * ech * hbarinv
      zJmtrx = zJmtrx * zcoeff


      return
 1    format(i0)
      end
c######################################################################

c%%%%
c%%%% END_OF_FILE:resistivity.f
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
