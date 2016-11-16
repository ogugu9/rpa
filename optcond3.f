c#########################################################
c#### optcond4                EK 2015
c#########################################################

      subroutine optcond5(xory,mk)

      use common, only : Ns, Pi, Dmu, Dne, Dtemp, Nqx, Nqy,
     &     isite, isitex, isitey, ZI, t, lln, Erng, Deta, Maxom
      implicit none

      character, intent(in) :: xory*1
      integer, intent(in) :: mk
      logical :: torf
      character :: chara*80, what*1
      integer :: ne0
      integer :: i, mu, nu, ikx, iky, kx, ky, iom, mkx, mky
      integer :: nxbl, nxbu, nybl, nybu, nnn
      integer :: iep1, iep2, iorb1, iorb2
      integer :: npoles, npl0
      integer :: is, ix, iy, im, in, iQ, isign, iQloop
      integer :: kdelx, kdely
      integer :: ispin
c      integer :: npp, mpp, kx1, ky1, kx2, ky2
c      integer :: idummy, isn, ism
      integer, external :: inverse
      real*8, parameter :: ddd = 1.0d-10 
      real*8 :: coeff, pi2
      real*8 :: e0, w, etest
      real*8 :: avec(2)
      real*8 :: dummy, broadw, spectra, akx, aky
      real*8 :: sc, ss, cond, df1, df2, E1, E2
      real*8 :: adelta, ax, ay
      real*8 :: dkvec(2)
      real*8, allocatable :: sint(:)
      real*8, allocatable :: poles(:)
      real*8, allocatable :: eig(:)
      complex*16, allocatable :: zpsi(:,:)
      complex*16 :: zpsi12(Ns*Nqx,2)
      real*8, external :: Dffn
      integer, external :: IcheckZone
      complex*16 :: zzz1, zzz2, zzz3, zpsi1, zpsi2, zdummy


      pi2 = 2.0d0 * Pi

      mkx = mk
      mky = mk

      torf = allocated(sint)
      if (torf .EQV. .true.) deallocate(sint)
      torf = allocated(poles)
      if (torf .EQV. .true.) deallocate(poles)
      torf = allocated(eig)
      if (torf .EQV. .true.) deallocate(eig)
      torf = allocated(zpsi)
      if (torf .EQV. .true.) deallocate(zpsi)

      allocate( eig(Ns*Nqx) )
      allocate( zpsi(Ns*Nqx,Ns*Nqx) )
      allocate( sint(mkx*mky*Ns*Nqx*(Ns*Nqx+1)/2) )
      allocate( poles(mkx*mky*Ns*Nqx*(Ns*Nqx+1)/2) )


      e0 = Erng
      ne0 = Maxom
      broadw = Deta
      write(*,*)'e0=', e0
      write(*,*)'ne0=', ne0
      write(*,*)'broadw=', broadw
c      etest = 0.075
      etest = 0.0d0
 

      open(10,file = "out.poles.dat")
      open(20,file = "out.condall.dat")
c      open(30,file = "out.bandall.dat")
      write(10,'(a,3i6)') '#Nk= ', mkx
      write(10,'(a)') '#E  intensity   E_initial   E_final'
      write(10,'(a)') '#Beginning'
      write(20,'(a,3i6)') 'mkx_mky_Nqx= ', mkx, mky, Nqx
      write(20,'(a,3i6)') '#Nkx_Nky_Nqx= ', mkx, mky, Nqx

c      chara = '    '

c## polarization[
      call polarization(xory,avec)
      ax = avec(1)
      ay = avec(2)
c## polarization]

      nxbl = 0 
      nxbu = mkx / Nqx - 1
c      nxbu = mkx  - 1
c      nxbl = -nxbu
      nybl = 0
      nybu = mky   - 1
c      nybl = -nybu
      nnn = (nxbu - nxbl + 1) * (nybu - nybl + 1)

      write(20,'(a,4i6)') 'kx0_kxn_ky0_kyn=', nxbl, nxbu, nybl, nybu

      npoles = 0
      npl0 = 0
      i = 0
      do ispin = 1, 2
         write(20,'(a)') ' '
         write(20,'(a,i6)') '#spin=', ispin

         do ikx = 0, mkx-1 ; do iky = 0, mky-1
            kx = MOD(ikx + mkx, mkx)
            ky = MOD(iky + mky, mky)

            if (IcheckZone(kx,ky,nxbl,nybl,nxbu,nybu) == 0) cycle
            
            akx = DBLE(kx) / DBLE(mkx)
            aky = DBLE(ky) / DBLE(mky)
            dkvec = (/ akx, aky/)
            akx = akx * pi2
            aky = aky * pi2
            call  diag(eig(:),zpsi(:,:),akx,aky,ispin)
            
            do in = 1, Ns * Nqx ; do im = in, Ns * Nqx
               
               zzz2 = 0.0d0
               zzz3 = 0.0d0

               iep1 = in
               iep2 = im
               E1 = eig(iep1)
               E2 = eig(iep2)

               if (ABS(E1 - E2) > e0 + broadw * 4.0d0) cycle
               if (ABS(E1 - E2) <= ddd) cycle

               if (E1 < E2) then ! Label the staes 1 & 2 so that E1 > E2
                  iep1 = im
                  iep2 = in
                  E1 = eig(iep1)
                  E2 = eig(iep2)
               end if
               
               df1 = Dffn(E1 - Dmu, Dtemp)
               df2 = Dffn(E2 - Dmu, Dtemp)
               
               if (ABS(df1 - df2) <= ddd) cycle

               zpsi12(:,1) = zpsi(:,iep1)
               zpsi12(:,2) = zpsi(:,iep2)
               call currentmtrx(avec,dkvec,zpsi12(:,:),zzz2)

               coeff = (df1 - df2) / (E1 - E2)         
               cond = - ABS(zzz2)**2
               cond = cond * coeff
               cond = cond * Pi / DBLE(nnn * Nqx)
 
               
               if (ABS(cond) <= ddd) cycle

               write(10,'(4g25.15)') E1 - E2, cond, 
     &              E2 - Dmu, E1 - Dmu
               write(20,'(2i6,2g20.10,5i6)') ikx, iky, E1 - E2, cond,
     &              iep1, iep2, ispin
               npoles = npoles + 1
               
               
            end do ; end do
         end do ; end do
      end do
      
      write(10,'(a,3i20)') '#Npoles= ', npoles

      write(10,'(a)')'#EOF' 
      close(10)

      write(20,'(a)')'#EOF' 
      close(20)


      open(10,file = "out.poles.dat")
      chara(1:10)='###'
      what = '?'
      i = 0
      do while(chara(1:10) /= '#Beginning')
         read(10,*) chara
      end do

      read(10,*) what
      do while (what /= '#')
         backspace(10)
         if (i > mkx*mky*Ns*Nqx*(Ns*Nqx+1)) then
            write(*,*) 'ERROR in optcond'
            stop
         end if
         i = i + 1
         read(10,*) poles(i), sint(i)
         read(10,*) what
      end do
      close(10)
      npoles = i

      write(*,*) MINVAL(poles(1:npoles)),MAXVAL(poles(1:npoles))

      open(10,file = "out.cond.dat")
         write(10,'(a,f10.6,i10)')'# broadw,Npoles=',broadw,npoles
      do iom = 0, ne0
         spectra = 0.0d0
         w = DBLE(iom) / DBLE(ne0) * e0
         
         do i = 1, npoles
            spectra = spectra + sint(i) / ((poles(i)-w)**2 + broadw**2)
         end do

         write(10,*) w, spectra * broadw
      end do
      close(10)
 100  continue

      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### optcond4                EK 2015
c#########################################################

      subroutine optcond4(xory,mk)

      use common, only : Ns, Ns2, Pi,
     &     Eall, Zpsiall, Dmu, Dne, Dtemp, Nqx, Nqy, Dkshift,
     &     isite, isitex, isitey, ZI, t, lln, Erng, Deta, Maxom
      implicit none

      character, intent(in) :: xory*1
      integer, intent(in) :: mk
      logical :: torf
      character :: chara*80, what*1
      integer :: ne0
      integer :: i, mu, nu, ikx, iky, kx, ky, iom, mkx, mky
      integer :: nxbl, nxbu, nybl, nybu, nnn
      integer :: iep1, iep2, iorb1, iorb2
      integer :: npoles, npl0
      integer :: is, ix, iy, im, in, iQ, isign, iQloop
      integer :: kdelx, kdely
      integer :: ispin
c      integer :: npp, mpp, kx1, ky1, kx2, ky2
c      integer :: idummy, isn, ism
      integer, external :: inverse
      real*8, parameter :: ddd = 1.0d-10 
      real*8 :: coeff, pi2
      real*8 :: e0, w, etest
      real*8 :: avec(2)
      real*8 :: dummy, broadw, spectra, akx, aky
      real*8 :: sc, ss, cond, df1, df2, E1, E2
      real*8 :: adelta, ax, ay
      real*8, allocatable :: sint(:)
      real*8, allocatable :: poles(:)
      real*8, allocatable :: eig(:,:)
      complex*16, allocatable :: zpsi(:,:,:)
      real*8, external :: Dffn
      integer, external :: IcheckZone
      complex*16 :: zzz1, zzz2, zzz3, zpsi1, zpsi2, zdummy


      pi2 = 2.0d0 * Pi

      mkx = mk
      mky = mk

      torf = allocated(sint)
      if (torf .EQV. .true.) deallocate(sint)
      torf = allocated(poles)
      if (torf .EQV. .true.) deallocate(poles)
      torf = allocated(eig)
      if (torf .EQV. .true.) deallocate(eig)
      torf = allocated(zpsi)
      if (torf .EQV. .true.) deallocate(zpsi)

      allocate( eig(Ns*Nqx,2) )
      allocate( zpsi(Ns*Nqx,Ns*Nqx,2) )
      allocate( sint(mkx*mky*Ns*Nqx*(Ns*Nqx+1)/2) )
      allocate( poles(mkx*mky*Ns*Nqx*(Ns*Nqx+1)/2) )


      e0 = Erng
      ne0 = Maxom
      broadw = Deta
      write(*,*)'e0=', e0
      write(*,*)'ne0=', ne0
      write(*,*)'broadw=', broadw
c      etest = 0.075
      etest = 0.0d0
 

      open(10,file = "out.poles.dat")
      open(20,file = "out.condall.dat")
c      open(30,file = "out.bandall.dat")
      write(10,'(a,3i6)') '#Nk= ', mkx
      write(10,'(a)') '#E  intensity   E_initial   E_final'
      write(10,'(a)') '#Beginning'
      write(20,'(a,3i6)') 'mkx_mky_Nqx= ', mkx, mky, Nqx
      write(20,'(a,3i6)') '#Nkx_Nky_Nqx= ', mkx, mky, Nqx

c      chara = '    '

c## polarization[
      call polarization(xory,avec)
      ax = avec(1)
      ay = avec(2)
c## polarization]

      nxbl = 0 
      nxbu = mkx / Nqx - 1
c      nxbu = mkx  - 1
c      nxbl = -nxbu
      nybl = 0
      nybu = mky   - 1
c      nybl = -nybu
      nnn = (nxbu - nxbl + 1) * (nybu - nybl + 1)

      write(20,'(a,4i6)') 'kx0_kxn_ky0_kyn=', nxbl, nxbu, nybl, nybu

      npoles = 0
      npl0 = 0
      i = 0
      do ispin = 1, 2
         write(20,'(a)') ' '
         write(20,'(a,i6)') '#spin=', ispin

         do ikx = 0, mkx-1 ; do iky = 0, mky-1
            kx = MOD(ikx + mkx, mkx)
            ky = MOD(iky + mky, mky)

            if (IcheckZone(kx,ky,nxbl,nybl,nxbu,nybu) == 0) cycle
            
            akx = DBLE(kx) / DBLE(mkx) * pi2
            aky = DBLE(ky) / DBLE(mky) * pi2
            call  diag(eig(:,ispin),zpsi(:,:,ispin),akx,aky,ispin)
            
            do in = 1, Ns * Nqx ; do im = in, Ns * Nqx
               
               zzz2 = 0.0d0
               zzz3 = 0.0d0

               iep1 = in
               iep2 = im
               E1 = eig(iep1,ispin)
               E2 = eig(iep2,ispin)

               if (ABS(E1 - E2) > e0 + broadw * 4.0d0) cycle
               if (ABS(E1 - E2) <= ddd) cycle

               if (E1 < E2) then ! Label the staes 1 & 2 so that E1 > E2
                  iep1 = im
                  iep2 = in
                  E1 = eig(iep1,ispin)
                  E2 = eig(iep2,ispin)
               end if
               
               df1 = Dffn(E1 - Dmu, Dtemp)
               df2 = Dffn(E2 - Dmu, Dtemp)
               
               if (ABS(df1 - df2) <= ddd) cycle

               cond = 0.0d0
               do is = 1, lln - 1
                  ix = isitex(is)
                  iy = isitey(is)
                  
                  adelta = DBLE(ix) * ax + DBLE(iy) * ay
                  
                  if (ABS(adelta) < ddd) cycle  
                  do mu = 1, Ns ;  do nu = 1, Ns
                     
                     zzz1 = 0.0d0
                     do iQ = 0, Nqx-1
                        kdelx = (kx + iQ * mkx / Nqx) * ix
                        kdely = (ky + iQ * mky / Nqy) * iy
                        akx = Pi2 * DBLE(kdelx) / DBLE(mkx)
                        aky = Pi2 * DBLE(kdely) / DBLE(mky)
                        
                        do iQloop = 0, Nqx-1
                           iorb1 = mu + MOD(iQloop + iQ,Nqx) * Ns
                           iorb2 = nu + MOD(iQloop + iQ,Nqx) * Ns
                           zpsi1 = zpsi(iorb1,iep1,ispin)
                           zpsi2 = zpsi(iorb2,iep2,ispin)
                        
                           zdummy = zpsi1 * CONJG(zpsi2)
                           zdummy = zdummy * EXP(-ZI*(akx + aky))
                           zzz1 = zzz1 + zdummy
                        end do
                     end do
                     zzz2 = zzz2
     &                    + zzz1 * adelta * t(ix,iy,mu,nu)
                  end do ; end do

                  
                  
               end do
                  
               coeff = (df1 - df2) / (E1 - E2)         
               cond = - DBLE(zzz2 * DCONJG(zzz2))
               cond = cond * coeff
               cond = cond * Pi / DBLE(nnn * Nqx)
               
               if (ABS(cond) <= ddd) cycle

               write(10,'(4g25.15)') E1 - E2, cond, 
     &              E2 - Dmu, E1 - Dmu
               write(20,'(2i6,2g20.10,5i6)') ikx, iky, E1 - E2, cond,
     &              iep1, iep2, ispin
               npoles = npoles + 1
               
               
            end do ; end do
         end do ; end do
      end do
      
      write(10,'(a,3i20)') '#Npoles= ', npoles

      write(10,'(a)')'#EOF' 
      close(10)

      write(20,'(a)')'#EOF' 
      close(20)


      open(10,file = "out.poles.dat")
      chara(1:10)='###'
      what = '?'
      i = 0
      do while(chara(1:10) /= '#Beginning')
         read(10,*) chara
      end do

      read(10,*) what
      do while (what /= '#')
         backspace(10)
         if (i > mkx*mky*Ns*Nqx*(Ns*Nqx+1)) then
            write(*,*) 'ERROR in optcond'
            stop
         end if
         i = i + 1
         read(10,*) poles(i), sint(i)
         read(10,*) what
      end do
      close(10)
      npoles = i

      write(*,*) MINVAL(poles(1:npoles)),MAXVAL(poles(1:npoles))

      open(10,file = "out.cond.dat")
         write(10,'(a,f10.6,i10)')'# broadw,Npoles=',broadw,npoles
      do iom = 0, ne0
         spectra = 0.0d0
         w = DBLE(iom) / DBLE(ne0) * e0
         
         do i = 1, npoles
            spectra = spectra + sint(i) / ((poles(i)-w)**2 + broadw**2)
         end do

         write(10,*) w, spectra * broadw
      end do
      close(10)
 100  continue

      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### optcond2
c#########################################################

      subroutine optcond2(xory,mk)

      use common, only : Ns, Ns2, Pi,
     &     Eall, Zpsiall, Dmu, Dne, Dtemp, Nqx, Nqy, Dkshift,
     &     isite, isitex, isitey, ZI, t, lln, Erng, Deta, Maxom
      implicit none

      character, intent(in) :: xory*1
      integer, intent(in) :: mk
      logical :: torf
      character :: chara*80, what*1
      integer :: ne0
      integer :: i, mu, nu, ikx, iky, kx, ky, iom, mkx, mky
      integer :: nxbl, nxbu, nybl, nybu, nnn, nrgn
      integer :: iep1, iep2, iorb1, iorb2
      integer :: npoles, npl0
      integer :: is, ix, iy, im, in, iQ, isign
      integer :: kdelx, kdely
      integer :: ispin
c      integer :: npp, mpp, kx1, ky1, kx2, ky2
c      integer :: idummy, isn, ism
      integer, external :: inverse
      real*8, parameter :: ddd = 1.0d-15 
      real*8 :: avec(2)
      real*8 :: coeff, pi2
      real*8 :: e0, w, etest
      real*8 :: dummy, broadw, spectra, akx, aky
      real*8 :: sc, ss, cond, df1, df2, E1, E2
      real*8 :: adelta, ax, ay
      real*8, allocatable :: sint(:)
      real*8, allocatable :: poles(:)
      real*8, allocatable :: eig(:,:)
      complex*16, allocatable :: zpsi(:,:,:)
      real*8, external :: Dffn
      complex*16 :: zzz1, zzz2, zpsi1, zpsi2, zdummy


      pi2 = 2.0d0 * Pi

      mkx = mk
      mky = mk

      torf = allocated(sint)
      if (torf .EQV. .true.) deallocate(sint)
      torf = allocated(poles)
      if (torf .EQV. .true.) deallocate(poles)
      torf = allocated(eig)
      if (torf .EQV. .true.) deallocate(eig)
      torf = allocated(zpsi)
      if (torf .EQV. .true.) deallocate(zpsi)
      allocate( sint(mkx*mky*Ns*Nqx*(Ns*Nqx+1)/2) )
      allocate( poles(mkx*mky*Ns*Nqx*(Ns*Nqx+1)/2) )
      allocate( eig(Ns*Nqx,2) )
      allocate( zpsi(Ns*Nqx,Ns*Nqx,2) )


      e0 = Erng
      ne0 = Maxom
      broadw = Deta
      write(*,*)'e0=', e0
      write(*,*)'ne0=', ne0
      write(*,*)'broadw=', broadw
c      etest = 0.075
      etest = 0.0d0
 

      open(10,file = "out.poles.dat")
      open(20,file = "out.condall.dat")
c      open(30,file = "out.bandall.dat")
      write(10,'(a,3i6)') '#Nk= ', mkx
      write(10,'(a)') '#E  intensity   E_initial   E_final'
      write(10,'(a)') '#Beginning'
      write(20,'(a,3i6)') 'mkx_mky_Nqx= ', mkx, mky, Nqx
      write(20,'(a,3i6)') '#Nkx_Nky_Nqx= ', mkx, mky, Nqx

c      chara = '    '

c## polarization[
      call polarization(xory,avec)
      ax = avec(1)
      ay = avec(2)
c## polarization]

      nrgn = Nqx
      nxbl = 0 
      nxbu = mkx / Nqx * nrgn / 2   - 1
c      nxbl = -nxbu
      nybl = 0
      nybu = mky   - 1
c      nybl = -nybu
      nnn = (nxbu - nxbl + 1) * (nybu - nybl + 1)

      write(20,'(a,4i6)') 'kx0_kxn_ky0_kyn=', nxbl, nxbu, nybl, nybu

      npoles = 0
      npl0 = 0
      i = 0
      do ispin = 1, 2
         write(20,'(a)') ' '
         write(20,'(a,i6)') '#spin=', ispin

         do ikx = nxbl, nxbu ; do iky = nybl, nybu
            kx = MOD(ikx + mkx, mkx)
            ky = MOD(iky + mky, mky)
            
            akx = DBLE(kx) / DBLE(mkx) * pi2
            aky = DBLE(ky) / DBLE(mky) * pi2
            call  diag(eig(:,ispin),zpsi(:,:,ispin),akx,aky,ispin)
            
            do in = 1, Ns * Nqx ; do im = in, Ns * Nqx
               
               zzz2 = 0.0d0

               iep1 = in
               iep2 = im
               E1 = eig(iep1,ispin)
               E2 = eig(iep2,ispin)

               if (E1 < E2) then ! Label the staes 1 & 2 so that E1 > E2
                  iep1 = im
                  iep2 = in
                  E1 = eig(iep1,ispin)
                  E2 = eig(iep2,ispin)
               end if
               
               if (ABS(E1 - E2) > e0 + broadw * 4.0d0) cycle
               if (ABS(E1 - E2) <= ddd) cycle
               df1 = Dffn(E1 - Dmu, Dtemp)
               df2 = Dffn(E2 - Dmu, Dtemp)
               
               ! zero temperature
c               if (E1 <= Dmu) cycle
c               if (E2 >= Dmu) cycle

               cond = 0.0d0
               if (ABS(df1 - df2) <= ddd) cycle
               do is = 1, lln - 1
                  ix = isitex(is)
                  iy = isitey(is)
                  
                     !aaaaaaaaaaaaaaaaaa
c     if((ABS(ix) == 2).and.(ABS(iy) == 2)) cycle
                  
                  adelta = DBLE(ix) * ax + DBLE(iy) * ay
                     
                  if (ABS(adelta) <= ddd) cycle      
                  do mu = 1, Ns ;  do nu = 1, Ns
                     isign = inverse(mu,nu)
                     sc = DBLE(1 - isign)
                     ss = DBLE(1 + isign)
                     zzz1 = 0.0d0
                     do iQ = 0, Nqx-1
                        kdelx = (kx + iQ * mkx / Nqx) * ix
                        kdely = (ky + iQ * mky / Nqy) * iy
                        akx = Pi2 * DBLE(kdelx) / DBLE(mkx)
                        aky = Pi2 * DBLE(kdely) / DBLE(mky)
                        
                        iorb1 = mu ! + iQ * Ns
                        iorb2 = nu ! + iQ * Ns
                        zpsi1 = zpsi(iorb1,iep1,ispin)
                        zpsi2 = zpsi(iorb2,iep2,ispin)
                        
                        zdummy = zpsi1 * CONJG(zpsi2)
                        zdummy = zdummy
     &                       * (    sc * COS(akx + aky)
     &                       - ZI * ss * SIN(akx + aky) )
                        zzz1 = zzz1 + zdummy
                     end do
                     zzz2 = zzz2
     &                    + zzz1 * adelta * t(ix,iy,mu,nu)
                  end do ; end do
            
                     
               end do
                  
               coeff = (df1 - df2) / (E1 - E2)
               cond = -DBLE(zzz2 * DCONJG(zzz2))
               cond = cond * coeff
               cond = cond * Pi / DBLE(nnn * Nqx)
               if (ABS(cond) <= ddd) cycle
               
               write(10,'(4g25.15)') E1 - E2, cond, 
     &              E2 - Dmu, E1 - Dmu
               write(20,'(2i6,2g20.10,5i6)') ikx, iky, E1 - E2, cond,
     &              iep1, iep2, ispin
               npoles = npoles + 1

                              
            end do ; end do
         end do ; end do
      end do
      
      write(10,'(a,3i20)') '#Npoles= ', npoles

      write(10,'(a)')'#EOF' 
      close(10)

      write(20,'(a)')'#EOF' 
      close(20)


c      write(30,'(a)')'#EOF' 
c      close(30)

c      call optcondout(ne0,xory)
c      return


c      goto 100
      open(10,file = "out.poles.dat")
      chara(1:10)='###'
      what = '?'
      i = 0
      do while(chara(1:10) /= '#Beginning')
         read(10,*) chara
      end do

      read(10,*) what
      do while (what /= '#')
         backspace(10)
         if (i > mkx*mky*Ns*(Ns*Nqx+1) * nrgn) then
            write(*,*) 'ERROR in optcond'
            stop
         end if
         i = i + 1
         read(10,*) poles(i), sint(i)
         read(10,*) what
      end do
      close(10)
      npoles = i

      open(10,file = "out.cond.dat")
      do iom = 0, ne0
         spectra = 0.0d0
         w = DBLE(iom) / DBLE(ne0) * e0
         
         do i = 1, npoles
            spectra = spectra + sint(i) / ((poles(i)-w)**2 + broadw**2)
         end do

         write(10,*) w, spectra * broadw
      end do
      close(10)
 100  continue

      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### optcond3
c#########################################################

      subroutine optcond3(xory,mk)

      use common, only : Nkx, Nky, Ns, Ns2, Pi,
     &     Eall, Zpsiall, Dmu, Dne, Dtemp, Nqx, Nqy,
     &     isite, isitex, isitey, ZI, t, lln, Erng, Deta, Maxom
      implicit none

      character, intent(in) :: xory*1
      integer, intent(in) :: mk
      logical :: torf
      character :: chara*80, what*1
      integer :: ne0
      integer :: i, mu, nu, ikx, iky, kx, ky, iom, mkx, mky
      integer :: nxbl, nxbu, nybl, nybu, nnn, nrgn
      integer :: iep1, iep2, iorb1, iorb2
      integer :: npoles, npl0
      integer :: is, ix, iy, im, in, iQ, isign
      integer :: kdelx, kdely
      integer :: ispin
c      integer :: npp, mpp, kx1, ky1, kx2, ky2
c      integer :: idummy, isn, ism
      integer, external :: inverse
      real*8, parameter :: ddd = 1.0d-15 
      real*8 :: avec(2)
      real*8 :: coeff, pi2
      real*8 :: e0, w, etest
      real*8 :: dummy, broadw, spectra, akx, aky
      real*8 :: sc, ss, cond, df1, df2, E1, E2
      real*8 :: adelta, ax, ay
      complex*16 :: zcond(0:lln-1)
      real*8, allocatable :: sint(:)
      real*8, allocatable :: poles(:)
      real*8, allocatable :: eig(:,:)
      complex*16, allocatable :: zpsi(:,:,:)
      real*8, external :: Dffn
      complex*16 :: zzz1, zzz2, zpsi1, zpsi2, zdummy


      pi2 = 2.0d0 * Pi

      mkx = mk
      mky = mk

      torf = allocated(sint)
      if (torf .EQV. .true.) deallocate(sint)
      torf = allocated(poles)
      if (torf .EQV. .true.) deallocate(poles)
      torf = allocated(eig)
      if (torf .EQV. .true.) deallocate(eig)
      torf = allocated(zpsi)
      if (torf .EQV. .true.) deallocate(zpsi)
      allocate( sint(mkx*mky*Ns2*(Ns2+1)/2) )
      allocate( poles(mkx*mky*Ns2*(Ns2+1)/2) )
      allocate( eig(Ns*Nqx,2) )
      allocate( zpsi(Ns*Nqx,Ns*Nqx,2) )


      e0 = Erng
      ne0 = Maxom
      broadw = Deta
      write(*,*)'e0=', e0
      write(*,*)'ne0=', ne0
      write(*,*)'broadw=', broadw
c      etest = 0.075
      etest = 0.0d0
 

      open(10,file = "out.poles.dat")
      open(20,file = "out.condall.dat")
c      open(30,file = "out.bandall.dat")
      write(10,'(a,3i6)') '#Nk= ', mkx
      write(10,'(a)') '#E  intensity   E_initial   E_final'
      write(20,'(a,3i6)') 'mkx_mky_Nqx= ', mkx, mky, Nqx
      write(20,'(a,3i6)') '#Nkx_Nky_Nqx= ', mkx, mky, Nqx

c      chara = '    '

c## polarization[
      call polarization(xory,avec)
      ax = avec(1)
      ay = avec(2)
c## polarization]

      nrgn = 2
      nxbl = 0 
      nxbu = mkx / Nqx * nrgn / 2   - 1
c      nxbl = -nxbu
      nybl = 0
      nybu = mky   - 1
c      nybl = -nybu
      nnn = (nxbu - nxbl + 1) * (nybu - nybl + 1)

      write(20,'(a,4i6)') 'kx0_kxn_ky0_kyn=', nxbl, nxbu, nybl, nybu

      npoles = 0
      npl0 = 0
      i = 0
      do ispin = 1, 2
         write(20,'(a)') ' '
         write(20,'(a,i6)') '#spin=', ispin

         do ikx = nxbl, nxbu ; do iky = nybl, nybu
            kx = MOD(ikx + mkx, mkx)
            ky = MOD(iky + mky, mky)
            
            akx = DBLE(kx) / DBLE(mkx) * 2.0d0 * Pi
            aky = DBLE(ky) / DBLE(mky) * 2.0d0 * Pi
            call  diag(eig(:,ispin),zpsi(:,:,ispin),akx,aky,ispin)
            
            do in = 1, Ns * Nqx ; do im = in, Ns * Nqx
               
               zzz2 = 0.0d0

               iep1 = in
               iep2 = im
               E1 = eig(iep1,ispin)
               E2 = eig(iep2,ispin)

               if (E1 < E2) then ! Label the staes 1 & 2 so that E1 > E2
                  iep1 = im
                  iep2 = in
                  E1 = eig(iep1,ispin)
                  E2 = eig(iep2,ispin)
               end if
               
               df1 = Dffn(E1 - Dmu, Dtemp)
               df2 = Dffn(E2 - Dmu, Dtemp)
               
               ! zero temperature
c               if (E1 <= Dmu) cycle
c               if (E2 >= Dmu) cycle

               cond = 0.0d0
               if ( (ABS(E1 - E2) > ddd)
     &              .and.(ABS(df1 - df2) > ddd)) then

!$omp parallel do private(is) shared(zpsi,zcond)
                  do is = 1, lln - 1
                     call calcoptcond3(is,mkx,mky,kx,ky,ax,ay,
     &                    iep1,iep2,ispin,zpsi,zcond(is),ddd)
                  end do
!$omp end parallel do
                  zzz2 = SUM(zcond(0:lln-1))
                  
                  coeff = 0.0d0
                  if ( (ABS(E1 - E2) > ddd)
     &                 .and.(ABS(df1 - df2) > ddd) ) then
                     coeff = (df1 - df2) / (E1 - E2)
                  else if (ABS(E1 - E2) < ddd) then !Drude
                  end if
                  
                  cond = -DBLE(zzz2 * DCONJG(zzz2))
                  cond = cond * coeff
                  cond = cond * Pi / DBLE(nnn * Nqx)
                  if (ABS(cond) < ddd) cond = 0.0d0

                  if (ABS(cond) >= ddd) then
                     write(10,'(4g25.15)') E1 - E2, cond, 
     &                    E2 - Dmu, E1 - Dmu
                  end if

               end if
               
               if (cond > ddd) then
                  write(20,'(2i6,2g20.10,5i6)') ikx, iky, E1 - E2, cond,
     &                 iep1, iep2, ispin
                  npoles = npoles + 1
               end if

            end do ; end do
         end do ; end do
      end do
      
      write(10,'(a,3i20)') '#Npoles= ', npoles

      write(10,'(a)')'#EOF' 
      close(10)

      write(20,'(a)')'#EOF' 
      close(20)


c      write(30,'(a)')'#EOF' 
c      close(30)

c      call optcondout(ne0,xory)
c      return


c      goto 100
      open(10,file = "out.poles.dat")
      what = '?'
      i = 0
      do while (what /= '#')
         read(10,*) what
         if(what /= '#') then
            backspace(10)
            if (i > mkx*mky*Ns*(Ns2+1) * nrgn) then
               write(*,*) 'ERROR in optcond'
               stop
            end if
            i = i + 1
            read(10,*) poles(i), sint(i)
         end if
      end do
      close(10)
      npoles = i

      open(10,file = "out.cond.dat")
      write(10,'(a,f10.6,i10)')'# broadw,Npoles=',broadw,npoles
      do iom = 1, ne0
         spectra = 0.0d0
         w = DBLE(iom) / DBLE(ne0) * e0
         
         do i = 1, npoles
            spectra = spectra + sint(i) / ((poles(i)-w)**2 + broadw**2)
         end do

         write(10,*) w, spectra * broadw
      end do
      close(10)
 100  continue

      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### calcoptcond3
c#########################################################

      subroutine calcoptcond3(is,mkx,mky,kx,ky,ax,ay,iep1,iep2,ispin,
     &     zpsi,zcond,ddd)

      use common, only : Nkx, Nky, Ns, Pi, Nqx, Nqy, 
     &     isite, isitex, isitey, ZI, t
      implicit none

      integer, intent(in) :: mkx, mky, kx, ky
      integer, intent(in) :: iep1, iep2
      integer, intent(in) :: ispin, is
      real*8, intent(in) :: ddd
      real*8, intent(in) :: ax, ay
      complex*16, intent(in) :: zpsi(Ns*Nqx,Ns*Nqx,2)
      complex*16, intent(out) :: zcond

      integer :: mu, nu
      integer :: ix, iy, im, in, iQ, isign
      integer :: iorb1, iorb2
      integer :: kdelx, kdely
      real*8 :: akx, aky, sc, ss
      real*8 :: adelta
      complex*16 :: zzz1, zpsi1, zpsi2, zdummy
      integer, external :: inverse


      zcond = 0.0d0


      ix = isitex(is)
      iy = isitey(is)
      
!aaaaaaaaaaaaaaaaaa
c      if (ABS(ix) /= 2) return
c      if (ABS(iy) /= 1) return

      adelta = DBLE(ix) * ax + DBLE(iy) * ay
      
      if (ABS(adelta) > ddd) then                  
         do mu = 1, Ns ;  do nu = 1, Ns
            isign = inverse(mu,nu)
            sc = DBLE(1 - isign)
            ss = DBLE(1 + isign)
            zzz1 = 0.0d0
            do iQ = 0, Nqx-1
               kdelx = (kx + iQ * mkx / Nqx) * ix
               kdely = (ky + iQ * mky / Nqy) * iy
               akx = 2.0d0 * Pi * DBLE(kdelx) / DBLE(mkx)
               aky = 2.0d0 * Pi * DBLE(kdely) / DBLE(mky)
               
               iorb1 = mu + iQ * Ns
               iorb2 = nu + iQ * Ns
               zpsi1 = zpsi(iorb1,iep1,ispin)
               zpsi2 = zpsi(iorb2,iep2,ispin)
               
               zdummy = zpsi1 * CONJG(zpsi2)
               zdummy = zdummy
     &              * (    sc * COS(akx + aky)
     &              - ZI * ss * SIN(akx + aky) )
               zzz1 = zzz1 + zdummy
            end do
            zcond = zcond
     &           + zzz1 * adelta * t(ix,iy,mu,nu)
            
         end do ; end do
      end if
         
      return
 1    format(i0)
      end
c######################################################################


c#########################################################
c#### Drude:         modified by EK 2015
c#########################################################

      subroutine Drude4(xory,sigma0,finputvf)

      use common, only : Nkx, Nky, Ns, Ns2, Pi,
     &     Eall, Zpsiall, Dmu, Dne, Dtemp, Nqx, Nqy,
     &     isite, isitex, isitey, ZI, t, lln, fname,
     &     Erng, Maxom, Deta, Dkshift
      implicit none

      character, intent(in) :: xory*1
      character(len=*), intent(in) :: finputvf
      real*8, intent(out) :: sigma0

      integer :: ne0

      logical :: torf
      character :: chara*80
      integer :: mu, nu, iom, iband
      integer :: mkx, mky, ngrid
      integer :: iep1, iorb1, iorb2
      integer :: is, ix, iy, iQ, isign
      integer :: ispin, isp

      character(len=40) :: filein
      integer :: i
      integer, external :: inverse
      real*8 :: avec(2)
      real*8 :: coeff, dl, dcir
      real*8 :: e0, w, pi2
      real*8 :: dummy, broadw, spectra, akx, aky
      real*8 :: sc, ss, cond
      real*8 :: adelta, ax, ay
      real*8, allocatable :: eig(:)
      real*8 :: vfx, vfy, vf, drudew, gamma, tau
      real*8 :: dkx, dky, dkxp, dkyp
      real*8 :: veltot
      real*8, external :: Dffn
      complex*16 :: zzz1, zzz2, zpsi1, zpsi2, zdummy
      complex*16, allocatable :: zpsi(:,:)

c      real*8 :: akx, aky, dkx, dky, vx, vy, vk, qphoton, px, py
c      real*8, allocatable :: poles(:), sint(:)

c==== .
c==== .         b  y  a                           
c==== .          \ | /                  
c==== .           \|/__x                                   
c==== .                                 
c==== .                                        
c==== . 
      
      pi2 = 2.0d0 * Pi
      open(20,file = "out.drudeall.dat")

      torf = allocated(eig)
      if (torf .EQV. .true.) deallocate(eig)
      torf = allocated(zpsi)
      if (torf .EQV. .true.) deallocate(zpsi)
      allocate(eig(Ns*Nqx))
      allocate(zpsi(Ns*Nqx,Ns*Nqx))

      e0 = Erng
      ne0 = Maxom
      broadw = Deta
      write(*,*)'e0=', e0
      write(*,*)'ne0=', ne0
      write(*,*)'broadw=', broadw

      gamma = broadw
      tau = 0.5d0 / gamma 
      
c## polarization[
      call polarization(xory,avec)
      ax = avec(1)
      ay = avec(2)
c## polarization]

      filein = finputvf
      if(filein(1:1) == '#') then
         write(*,*) 'Enter input filename (fermi velocity data)'
         read(5,*) filein
         if (filein(1:1) == '!') filein = 'v-out.fs.dat'
      end if
      open(50,file = filein)    !open data file
      write(*,*) 'vf data read from ',filein

      chara = '###'
      do while (chara(1:10) /= '#Beginning')
         read(50,*) chara
      end do
      dcir = 0.0d0
      do while (chara(1:4) /= '#EOF')
         read(50,*) dkx, dky, vfx, vfy, dl, iband, ispin
         if (ispin == 1) dcir = dcir + dl
         read(50,*) chara
         backspace(50)
      end do
      rewind(50)
      write(*,*)'dcir=',dcir

      chara = '#'
      do while (chara(1:10) /= '#Beginning')
         read(50,*) chara
      end do

      veltot = 0.0d0 
      drudew = 0.0d0
      do while (chara(1:4) /= '#EOF')

         read(50,*) dkx, dky, vfx, vfy, dl, iband, ispin
         read(50,*) chara
         backspace(50)

         iep1 = iband
         isp = ispin
         vf = SQRT(vfx**2 + vfy**2)

         dkxp = dkx * pi2
         dkyp = dky * pi2

         call  diag(eig(:),zpsi(:,:),dkxp,dkyp,ispin)

         zzz2 = 0.0d0
         do is = 1, lln - 1
            
            ix = isitex(is)
            iy = isitey(is) 
            adelta = DBLE(ix) * ax + DBLE(iy) * ay
            
            if (ABS(adelta) > 0.0d0) then
               do mu = 1, Ns ;  do nu = 1, Ns

                  zzz1 = 0.0d0
                  do iQ = 0, Nqx-1
                     akx = (dkx + iQ / DBLE(Nqx)) * pi2 * DBLE(ix)
                     aky = (dky + iQ / DBLE(Nqy)) * pi2 * DBLE(iy)
c                     akx = dkx * pi2 * ix
c                     aky = dky * pi2 * iy

                     iorb1 = mu + iQ * Ns
                     iorb2 = nu + iQ * Ns
                     zpsi1 = zpsi(iorb1,iep1)
                     zpsi2 = zpsi(iorb2,iep1) !iep2 = iep1 for drude
                     
                     zdummy = zpsi1 * CONJG(zpsi2)
                     zdummy = zdummy * EXP(- ZI * (akx + aky))

                     zzz1 = zzz1 + zdummy
                  end do
                  zzz2 = zzz2
     &                 + zzz1 * adelta * t(ix,iy,mu,nu)
               end do ; end do
               
            end if
         end do
               
         veltot = veltot + zzz2 
         cond = ABS(zzz2)**2

         coeff = 1.0d0
         coeff = coeff * (dl / dcir)
         coeff = coeff / vf

         cond = cond * coeff

         if (ABS(cond) > 1.0d-15) drudew = drudew + cond
         write(20,'(3g20.10)') dkx, dky, cond

      end do 
      close(50)
      close(20)

      sigma0 = drudew

      write(*,*) 'veltot=', veltot
      open(10,file = "out.drude.dat")
      write(10,'(2(a,x),2(f8.3,x))') '#polarization=', xory, ax, ay
      write(10,'(a,x,f15.8,x,2(a,x,f8.3,x))') '#sigma0=', sigma0,
     &     'tau=', tau, 'gamma=', gamma
      do iom = 0, ne0
         spectra = 0.0d0
         w = DBLE(iom) / DBLE(ne0) * e0
         
         spectra = sigma0 * tau / ((tau * w)**2 + 1.0d0)/ Pi
         
         write(10,'(4(f13.8,x))') w, spectra, tau, drudew
      end do

      close(10)

      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### interband:
c#########################################################

      subroutine interband()

      use common, only : Nkx, Nky, Ns, Ns2, Pi,
     &     Eall, Zpsiall, Dmu, Dne, Dtemp, Nqx, Nqy,
     &     isite, isitex, isitey, ZI, t, lln
      implicit none

      character :: chara*15, what*1
      integer :: ikx, iky, kx, ky, nnx, nny, nnq
      integer :: iep1, iep2
      integer :: kx0, kxn, ky0, kyn
      real*8 :: w, emin, emax, cond, ddd
      real*8 :: sint(0:200,0:200)

      ddd = 1.0d-3

      write(*,*) 'Enter energy range:'
      write(*,*) 'Lower bound?'
      read(5,*) emin
      write(*,*) 'Upper bound?'
      read(5,*) emax
      write(*,*) 'Minimal intensity?'
      read(5,*) ddd
      
      open(20,file = "out.condall.dat")
      read(20,*) chara, nnx, nny, nnq
      write(*,*) chara, nnx, nny, nnq
      read(20,*) chara, kx0, kxn, ky0, kyn
      write(*,*) chara, kx0, kxn, ky0, kyn

      open(30,file = "out.bandall.dat")
      write(30,'(a,3i6)') '## Nkx_Nky_Nqx= ', nnx, nny, nnq
      write(30,'(a)') '## E1-E2, kx, ky,'//
     &     ' iep1, iep2, cond'

      sint(:,:)=0.0d0
      do while(chara(1:4) /='#EOF')
         read(20,*) chara
         if (chara(1:4) /='#EOF') then
            backspace(20)
            read(20,*) kx, ky, w, cond, iep1, iep2
            if ((w >= emin).and.(w <= emax))then
               sint(kx,ky) = sint(kx,ky) + cond
               if (ABS(cond) > ddd) then 
                  if (ABS(cond) > 1.0d-1) then
                     what = '*'
                  else
                     what = ' '
                  end if
                  write(30,'(g20.10,2i6,2x,2i6,4g20.10,a)') 
     &                 w, kx, ky,
     &                 iep1, iep2, cond, what
               end if        
            end if
         end if
      end do
      close(20)
      write(30,'(a)') '#EOF'
      close(30)

      open(10,file = "out.inter.dat")
      do iky = ky0 , kyn 
         do ikx = kx0, kxn
            kx = ikx
            ky = iky

            kx = MOD(kx, nnx/nnq)
            ky = MOD(ky, nny)

            write(10,*) DBLE(2*ikx)/DBLE(nnx), DBLE(2*iky)/DBLE(nny)
     &           , sint(kx,ky)
         end do
         write(10,*)
      end do
      close(10)


      return
      end
c#########################################################


c#########################################################
c#### calcDOS3:
c#########################################################

      subroutine calcDOS3(ne0)

      use common, only : Nkx, Nky, Ns, Ns2, Pi,
     &     Eall, Zpsiall, Dmu, Dne, Dtemp, Nqx, Nqy
      implicit none

      integer, intent(in) :: ne0
      character :: unitcell*5
      integer :: i, mu, nu, ikx, iky, kx, ky, iom
      integer :: nxbl, nxbu, nybl, nybu
      integer :: npoles
      integer :: ispin
      real*8 :: e0
      real*8 :: dummy, broadw, spectra
      real*8 :: sign
      real*8 :: tot, w

      real*8, allocatable :: dos(:)
      real*8, allocatable :: poles(:), sint(:,:)

c      real*8 :: dos(Ns2+1)
c      real*8 :: poles(Ns*Nqx*Nkx*Nky), sint(Ns*Nqx*Nkx*Nky,Ns*Nqx+1)
      complex*16 :: zdummy1, zdummy2
      complex*16, external :: zwvfn


c==== : sint(n,mu):  intensity of the n-th pole for mu orbital 
c==== : mu = [1,Ns] for the majority spin
c==== : mu = [Ns+1,Ns2] for the minority spin
c==== : mu = N2+1 for total

      unitcell = 'Fe-Fe'
c      unitcell = 'Fe-As'

      broadw = 0.1d-1
      e0 = 2.5d0


      allocate( dos(Ns2+1))
      allocate(poles(Ns*Nqx*Nkx*Nky), sint(Ns*Nqx*Nkx*Nky,Ns2+1))


      open(10,file = "out.dospl.dat")

      nxbl = 0
      nxbu = Nkx / 2 / Nqx - 1
      nybl = 0
      nybu = Nky /2  - 1

      npoles = 0
      do ikx = nxbl, nxbu ; do iky = nybl, nybu
         kx = MOD(ikx + Nkx, Nkx)
         kx = MOD(ikx + 2 * Nkx - Nkx / 2, Nkx)
         ky = MOD(iky + Nky, Nky)
         do ispin = 1, 2 ; do i = 1, Ns * Nqx
            sign = DBLE(Ns + 0.5d0 - i)
            sign = sign / ABS(sign)
            dummy = Eall(kx,ky,i,ispin) - Dmu

            npoles = npoles + 1
            poles(npoles) = dummy
            sint(npoles,Ns2+1) = 0.0d0
            do mu = 1, Ns
               nu = MOD(mu + Ns, Ns2)
               if (nu == 0) nu = Ns2
               if (Nqx == 1) nu = mu
               
               zdummy1 = DCONJG( Zpsiall(kx,ky,mu,i,ispin))
               zdummy2 = Zpsiall(kx,ky,nu,i,ispin)

               spectra = ABS(zdummy1)**2 + ABS(zdummy2)**2
               sint(npoles,Ns2+1) = sint(npoles,Ns2+1) + spectra
               
               sint(npoles,mu) = spectra
     &              + 2.0d0 * DBLE(zdummy1 * zdummy2) * sign
               sint(npoles,mu+Ns) = spectra
     &              - 2.0d0 * DBLE(zdummy1 * zdummy2) * sign
            end do
            write(10,*) dummy, sint(npoles,Ns2+1)

         end do ; end do
      end do ; end do

      tot = SUM(sint(1:npoles,1:Ns2)) * 0.5d0 
      write(*,*) 'sum of DOS =', tot, npoles
      tot = DBLE(npoles)

      write(10,'(a)') '#'
      write(10,'(a,g20.10)') '#total=',tot
      close(10)


      open(10,file = "out.spec.dat")
         write(10,'(a,a)') '### unitcell=' , unitcell
         if (unitcell == 'Fe-Fe') then
            write(10,'(a)') '###'//
     &           '  orbital 1:  3z^2-r^2'//
     &           '  2:  zx'//
     &           '  3:  yz'//
     &           '  4:  x^2-y^2'//
     &           '  5:  xy'
         else
            write(10,'(a)') '###'//
     &           '  orbital 1:  3Z^2-R^2'//
     &           '  2:  ZX'//
     &           '  3:  YZ'//
     &           '  4:  X^2-Y^2'//
     &           '  5:  XY'
         end if

         write(10,'(a)') '### omega,'//
     &     ' majority-spin (mu=1-5), minority-spin (mu=1-5),'//
     &     ' total' 
      do iom = -ne0, ne0
         w = DBLE(iom) / DBLE(ne0) * e0
         
         do mu = 1, Ns2
            dos(mu) = 0.0d0
            do i = 1, npoles
               dos(mu) = dos(mu) 
     &              + sint(i,mu) / ((poles(i)-w)**2 + broadw**2)
            end do
         end do
         dos(Ns2+1) = 0.0d0
         do i = 1, npoles
            dos(Ns2+1) = dos(Ns2+1) 
     &           + sint(i,Ns2+1) / ((poles(i)-w)**2 + broadw**2)
         end do

         write(10,'(12g20.10)') w, 
     &        ((dos(mu) * broadw / tot), mu=1,Ns2+1)
      end do


      close(10)

      deallocate( dos,poles,sint)

      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### transform of wave function
c#########################################################

      complex*16 function zwvfn(kx,ky,iorb,in,ispn)

      use common, only : Zpsiall, Ns, Ns2, Nqx
      implicit none

c#### The calculation is perfomed in reduced unit cell (one Fe per cell),
c#### but the orbital notation used in the calculation is
c#### corresponding to the original one (two Fe per cell).
c#### (X, Y, Z) used in calucation (notation in original unit cell).
c#### (x, y, z) (notation in reduced unit cell).
c#### (x,y,z) is 45-degree rotated from (X,Y,z)
c#### with respect to z-axis countercclockwise.
c#### (x, y, z) = (X+Y,-X+Y, z).
c#### 
c####      y Y
c####      |/_x
c####       \
c####        X


      integer, intent(in) :: kx, ky, iorb, in
      integer, intent(in), optional :: ispn
      integer :: i, is
      real*8 :: fac

      if (present(ispn) .eqv. .True.) then
         is = ispn
      else
         is = 1
      end if

      fac = 1.0d0 / SQRT(2.0d0)

      i = MOD(iorb,Ns)
      
      if (i == 1) then       ! 3z^2-r^2 = 3Z^2-R^2
         zwvfn = Zpsiall(kx,ky,iorb,in,is)
      else if (i == 2) then  ! zx = (ZX + YZ) / SQRT(2)
         zwvfn = Zpsiall(kx,ky,iorb,in,is)
     &        + Zpsiall(kx,ky,iorb+1,in,is)
         zwvfn = zwvfn * fac
      else if (i == 3) then  ! yz = (-ZX + YZ) / SQRT(2)
         zwvfn = -Zpsiall(kx,ky,iorb-1,in,is)
     &        + Zpsiall(kx,ky,iorb,in,is)
         zwvfn = zwvfn * fac
      else if (i == 4) then  ! x^2-y^2 = XY
         zwvfn = Zpsiall(kx,ky,iorb+1,in,is)
      else if (i == 0) then  ! xy = -X^2 + Y^2
         zwvfn = -Zpsiall(kx,ky,iorb-1,in,is)
      end if

      if (iorb > Ns * Nqx) then
         write(*,*) 'Error in fnction zwvfn'
         write(*,*) kx,ky,iorb,in
         stop
      end if

      return
      end
c######################################################################

c#########################################################
c#### calcDOS4:
c#########################################################

      subroutine calcDOS5(ne0,e0,broadw,mk)

      use common, only : Nkx, Nky, Ns, Ns2, Pi,
     &     Dmu, Dne, Dtemp, Nqx, Nqy, Deta
      implicit none

      logical :: torf
      integer, intent(in) :: ne0, mk
      real*8, intent(in) :: broadw, e0
      character :: unitcell*5
      integer :: i, mu, nu, ikx, iky, kx, ky, iom, mkx, mky
      integer :: nxbl, nxbu, nybl, nybu
      integer :: npoles, npoles1, npoles2
      integer :: is, nfile
      real*8 :: dummy, spectra
      real*8 :: sign
      real*8 :: tot, w
      real*8 :: dos(Ns2+1)
      real*8 :: akx, aky

      real*8 :: sint(Ns*mk**2,Ns*Nqx+1,2)
      real*8 :: poles(Ns*mk**2,2)
      real*8 :: eig(Ns*Nqx,2)
      complex*16 :: zpsi(Ns*Nqx,Ns*Nqx,2)

      complex*16 :: zdummy1, zdummy2
      complex*16, external :: zwvfn

      mkx = mk
      mky = mk

c==== : sint(n,mu):  intensity of the n-th pole for mu orbital 
c==== : mu = [1,Ns] for the majority spin
c==== : mu = [Ns+1,Ns2] for the minority spin
c==== : mu = N2+1 for total

      unitcell = 'Fe-Fe'
c      unitcell = 'Fe-As'

c      broadw = 0.1d-1
c      e0 = 2.5d0

      open(10,file = "out.dospl1.dat")
      open(20,file = "out.dospl2.dat")


      nxbl = 0
      nxbu = mkx / 2 / Nqx - 1
      nybl = 0
      nybu = mky /2  - 1

      do is = 1, 2
         npoles = 0
         do ikx = nxbl, nxbu ; do iky = nybl, nybu
            kx = MOD(ikx + mkx, mkx)
            ky = MOD(iky + mky, mky)
            
            akx = DBLE(kx) / DBLE(mkx) * 2.0d0 * Pi
            aky = DBLE(ky) / DBLE(mky) * 2.0d0 * Pi
            
            nfile = is * 10 
            call diag(eig(:,is),zpsi(:,:,is),akx,aky,is)

            do i = 1, Ns * Nqx
               npoles = npoles + 1
               poles(npoles,is) = eig(i,is) - Dmu
               sint(npoles,Ns*Nqx+1,is) = 0.0d0
               call lDOS(zpsi(:,i,is),sint(npoles,:,is))
               write(nfile,'(100g20.10)')
     &              poles(npoles,is), sint(npoles,:,is)
            end do
         end do ; end do
         if (is ==1) npoles1 = npoles 
         if (is ==2) npoles2 = npoles 
      end do

      if (npoles1 /= npoles2) then
         write(*,*) "# of poles are different for up and down"
      end if

      tot = SUM(sint(1:npoles1,1:Ns*Nqx,1)) / DBLE(Nqx)
      write(*,*) 'sum of DOS =', tot, npoles1
c      tot = DBLE(npoles)
      write(10,'(a)') '#'
      write(10,'(a,g20.10)') '#total=',tot
      write(10,'(a,g20.10)') '#total=',npoles1

      tot = SUM(sint(1:npoles2,1:Ns*Nqx,2)) / DBLE(Nqx)
      write(*,*) 'sum of DOS =', tot, npoles2
      write(20,'(a)') '#'
      write(20,'(a,g20.10)') '#total=',tot
      write(20,'(a,g20.10)') '#total=',npoles2

      close(10)
      close(20)

      open(10,file = "out.spec1.dat")
      open(20,file = "out.spec2.dat")
      do is = 1, 2
         nfile = is * 10
         write(nfile,'(a,a)') '### unitcell=' , unitcell
         if (unitcell == 'Fe-Fe') then
            write(10,'(a)') '###'//
     &           '  orbital 1:  3z^2-r^2'//
     &           '  2:  zx'//
     &           '  3:  yz'//
     &           '  4:  x^2-y^2'//
     &           '  5:  xy'
         else
            write(nfile,'(a)') '###'//
     &           '  orbital 1:  3Z^2-R^2'//
     &           '  2:  ZX'//
     &           '  3:  YZ'//
     &           '  4:  X^2-Y^2'//
     &           '  5:  XY'
         end if

         write(10,'(a)') '### omega,'//
     &     ' majority-spin (mu=1-5), minority-spin (mu=1-5),'//
     &     ' total' 

         npoles = npoles1
         if (is == 2) npoles = npoles2
         call spec(nfile,ne0,npoles,
     &           poles(1:npoles,is),sint(1:npoles,:,is),e0,Deta)
         end do

      close(10)


      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### lDOS:
c#########################################################

      subroutine lDOS(zpsi,sint)

      use common, only : Nkx, Nky, Ns, Zi, Pi, Nqx, Nqy
      implicit none

      real*8, intent(out) :: sint(Ns*Nqx)
      complex*16, intent(in) :: zpsi(Ns*Nqx)

      integer :: iorb, ix, iy, iQ, mu
      real*8 :: qx, qy, qr
      complex*16 :: zzz(Ns*Nqx)
      
      iy = 1
      do ix = 1, Nqx
         zzz(:) = 0.0d0
         do iQ = 0, Nqx -1
            qx = 2.0d0 * Pi * iQ / DBLE(Nqx) 
            qy = 2.0d0 * Pi * iQ / DBLE(Nqy)
            qr = qx * ix + qy * iy
            do mu = 1, Ns
               iorb = mu + Ns * iQ
               iorb = MOD(iorb, Ns * Nqx)
               if (iorb == 0) iorb = Ns * Nqx
               if (Nqx == 1) iorb = mu
               zzz(mu) = zzz(mu) +
     &              zpsi(iorb) * EXP(-Zi * qr)               
            end do
         end do
         do mu = 1, Ns
            sint(mu + (ix-1) * Ns) = ABS(zzz(mu))**2
         end do
      end do

      return
      end
c######################################################################

c#########################################################
c#### spec:
c#########################################################

      subroutine spec(nfile,ne0,npoles,poles,sint,e0,broadw)


      use common, only : Ns, Nqx
      implicit none

      logical :: torf
      integer, intent(in) :: ne0, nfile, npoles
      real*8, intent(in) :: e0, broadw
      real*8, intent(in) :: sint(npoles,Ns*Nqx)
      real*8, intent(in) :: poles(npoles)

      integer :: i, mu, iom
      real*8 :: w, tot
      real*8 :: dos(Ns * Nqx)

      tot = DBLE(npoles)
      do iom = -ne0, ne0
         w = DBLE(iom) / DBLE(ne0) * e0
         
         do mu = 1, Ns * Nqx
            dos(mu) = 0.0d0
            do i = 1, npoles
               dos(mu) = dos(mu) 
     &              + sint(i,mu) / ((poles(i)-w)**2 + broadw**2)
            end do
         end do
         
         write(nfile,'(25g20.10)') w, 
     &        ((dos(mu) * broadw / tot), mu = 1, Ns * Nqx)
      end do



      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### 
c#########################################################

      integer function IcheckZone(kx,ky,kx1,ky1,kx2,ky2)

      implicit none

      integer, intent(in) :: kx, ky, kx1, ky1, kx2, ky2

      IcheckZone = 0

      if(kx < kx1) return
      if(ky < ky1) return
      if(kx > kx2) return
      if(ky > ky2) return

      IcheckZone = 1

      return
      end
c######################################################################

c#########################################################
c#### 
c#########################################################

      subroutine polarization(xory,avec)

      implicit none

      character, intent(in) :: xory*1
      real*8, intent(out):: avec(2)
      real*8 :: dummy, ax, ay

      if (xory == 'x') then
         ax = 1.0d0
         ay = 0.0d0
      else if (xory == 'y') then
         ax = 0.d0
         ay = 1.0d0
      else if (xory == 'a') then
         ax = 1.0d0
         ay = -1.0d0
      else if (xory == 'b') then
         ax = 1.0d0
         ay = 1.0d0
      else
         stop
      end if

      dummy = SQRT(ax**2 + ay**2)
      avec = (/ ax, ay /) / dummy

      return
      end
c######################################################################


c%%%%
c%%%% END_OF_FILE:optcond2.f
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
