c#########################################################
c#### savedata: save MF data
c#########################################################
      subroutine savedata(chara)
      
      use common
      implicit none

      character , intent(in) :: chara*40
      character :: chara1*40
      integer :: iorb, mu, nu, is
      real*8 :: sz, ch, dummy
      real*8 :: drl, dim

      open(10,file=chara)

c## parameters[
      write(10,*) 'Nkx__Nky'
      write(10,*) Nkx, Nky
      write(10,*) 'Nqx__Nqy'
      write(10,*) Nqx, Nqy
      write(10,*) 'U_DJ'
      write(10,*) U,DJ
      write(10,*) 'Dtemp_Dne'
      write(10,*) Dtemp, Dne
      write(10,*) 'Erng__Maxom__Deta'
      write(10,*) Erng, Maxom, Deta
c## parameters]

      write(10,*) 'Erange'
      write(10,*) Erange
      write(10,*) 'Fen_Dmu'
      write(10,*) Fen,Dmu

      drl = 0.0d0
      dim = 0.0d0
      write(10,*)'uu'
      do mu = 1, Ns
         write(10,*) mu
         drl = DBLE(Zop(mu,mu,1))
         dim = DIMAG(Zop(mu,mu,1))
         if (ABS(drl) < 1.0d-16) drl = 0.0d0
         if (ABS(dim) < 1.0d-16) dim = 0.0d0
         write(10,'(2g24.15)') drl, dim
      end do

      write(10,*)'order_parameters'
      do is = 1,2 ; do mu = 1, Ns ; do nu = 1, Ns
         drl = DBLE(Zop(mu,nu,is))
         dim = DIMAG(Zop(mu,nu,is))
         if (ABS(drl) < 1.0d-16) drl = 0.0d0
         if (ABS(dim) < 1.0d-16) dim = 0.0d0
         write(10,'(3i5)') mu, nu, is 
         write(10,'(4g24.15)') drl, dim
c         Zop(mu,nu,1) = cmplx(drl,dim)
      end do ; end do ; end do

      write(10,*)'order_parameters2'
      do is = 1,2 ; do mu = 1, Ns ; do nu = 1, Ns
         drl = DBLE(Zop2(mu,nu,is))
         dim = DIMAG(Zop2(mu,nu,is))
         if (ABS(drl) < 1.0d-16) drl = 0.0d0
         if (ABS(dim) < 1.0d-16) dim = 0.0d0
         write(10,'(3i5)') mu, nu, is
         write(10,'(4g24.15)') drl, dim
c         Zop(mu,nu,1) = cmplx(drl,dim)
      end do ; end do ; end do

      write(10,*)'order_parameters0'
      do is = 1,2 ; do mu = 1, Ns ; do nu = 1, Ns
         drl = DBLE(Zop0(mu,nu,is))
         dim = DIMAG(Zop0(mu,nu,is))
         if (ABS(drl) < 1.0d-16) drl = 0.0d0
         if (ABS(dim) < 1.0d-16) dim = 0.0d0
         write(10,'(3i5)') mu, nu, is
         write(10,'(4g24.15)') drl, dim
c         Zop(mu,nu,1) = cmplx(drl,dim)
      end do ; end do ; end do



      write(10,*)'density'
      write(10,*)'orbital__sz__n'
      do iorb = 1, Ns
         sz = 2.0d0 * DBLE(Zop(iorb,iorb,1))
         ch = Dens(iorb,1) + Dens(iorb,2)
         if (iorb == 1) chara1 = "3z^2-r^3:"
         if (iorb == 2) chara1 = "xz______:"
         if (iorb == 3) chara1 = "yz______:"
         if (iorb == 4) chara1 = "x^2-y^2_:"
         if (iorb == 5) chara1 = "xy______:"
         write(10,'(a9,2x,2g24.15)') chara1(1:9), sz, ch
      end do

      write(10,*)
      write(10,*)'conrs=',conrs

      dummy = ABS(Dne1 - Dne)
      write(10,*) 'Dne1=', Dne1
      write(10,*) '|Dne1-Dne|=', dummy
      write(10,*)' temperature=', Dtemp

      close(10)
      write(*,'(2a16)')' data saved ==> ',chara
      return

      end
c###################################################################

c#########################################################
c#### loaddata: load MF data (saved by savedata)
c#########################################################

      subroutine loaddata(chara,ierr)
      
      use common
      implicit none

      character, intent(inout) :: chara*40
      integer, intent(out) :: ierr
      character :: chara1*40
      integer :: i, j, is, js, mu, nu
      real*8 :: dummy

      ierr=0

      open(10,file=chara,status='old', err=998)

      write(*,*)'MF data: file=',chara
c## parameters[
      write(*,*)'load parameters'
c      read(10,*,err=997) chara1
c      chara1=ADJUSTL(chara1)
c      lador2d = chara1
      read(10,*,err=997) chara1
      read(10,*,err=997) Nkx, Nky
      read(10,*,err=997) chara1
      read(10,*,err=997) Nqx, Nqy
      read(10,*,err=997) chara1
      read(10,*,err=997) U, DJ
      read(10,*,err=997) chara1
      read(10,*,err=997) Dtemp
      read(10,*,err=997) chara1
      read(10,*,err=997) Erng, Maxom, Deta

      call deallocation()
      call allocation()

      Nredx = 2 / Nqx
      Nredy = 2 / Nqy
c## parameters]

      read(10,*,err=997) chara1
      read(10,*,err=997) Erange
      read(10,*,err=997) chara1
      read(10,*,err=997) Fen, Dmu

      write(*,*)'load Dnuu'
      read(10,*,err=997)chara1
      do i = 1, Ns
         read(10,*) mu
         read(10,'(2g24.15)',err=997) Zop(mu,mu,1)
      end do
      
      write(*,*)'load order parameters'
      read(10,*,err=997)chara1
      do js = 1, 2 ; do i = 1, Ns ; do j = 1, Ns
         read(10,'(3i5)',err=997,end=997) mu, nu, is
         read(10,'(4g24.15)',err=997) Zop(mu,nu,is)
      end do ; end do ; end do

      read(10,*,err=997)chara1
      do js = 1, 2 ; do i = 1, Ns ; do j = 1, Ns
         read(10,'(3i5)',err=997,end=997) mu, nu, is
         read(10,'(4g24.15)',err=997) Zop2(mu,nu,is)
      end do ; end do ; end do

      read(10,*,err=997)chara1
      do js = 1, 2 ; do i = 1, Ns ; do j = 1, Ns
         read(10,'(3i5)',err=997,end=997) mu, nu, is
         read(10,'(4g24.15)',err=997) Zop0(mu,nu,is)
      end do ; end do ; end do

      Dne1 = 0.0d0
      do i = 1, Ns
         Dnuu(i) = DBLE(Zop(i,i,1))
         Dndd(i) = DBLE(Zop(i,i,2))
         Dne1 = Dne1 + DBLE(Zop0(i,i,1) + Zop0(i,i,2))
         Dens(i,:) = DBLE(Zop0(i,i,:))
      end do

      close (10)
      
      write(*,*) 
      write(*,*) '################'
      write(*,*) 'Loading complete'
      write(*,*) '################'
      return
      
      
 997  write(*,*) 'Cannot read data in file ', chara
      if (Xmode(1:1) /= 'm') stop
      write(*,*) 'return to Main Menu? (enter "y")'
      read(*,*)chara1
      ierr=-1
      return      

 998  write (6,*) 'Cannot open file ', chara
      if (Xmode(1:1) /= 'm') stop
      write(*,*) 'return to Main Menu? (enter "y")'
      read(*,*)chara1
      ierr=-1
      return      

      end
c###################################################################
      
c#########################################################
c#### eigenout: output eigen energy & eigen function
c#########################################################

      subroutine eigenout()

      use common, only : Nkx, Nky, Ns, Ns2, Nqx, Nqy,
     &     Eall, Zpsiall, Dmu, Dne, Dtemp
      implicit none

      integer :: i, j, k, ikx, iky, kx, ky, ispin
      integer :: lmtl, lmtu, mx, my
      real*8 :: p1, p2, p3, p4, p5, p6, emin, emax, e0, drl, dim
      real*8 :: dummy1, dummy2


      goto 100
      open(10,file = "eigen.out")
      write(10,*) Nkx, Nky
      write(10,*) 'Dne=', Dne
      write(10,*) 'Dmu=', Dmu
      write(10,*) 'Dtemp=', Dtemp

c===============
      write(10,*)'WAVE_FUNCTION'
      do kx = 0, Nkx-1 ; do ky = 0, Nky - 1
         do i = 1, Ns*Nqx ; do j = 1, Ns*Nqx
            write(10,'(4i5)') i, j, kx, ky
c            write(10,'(8g20.10,a)') Eall(kx,ky,j)
c     &           ,Zpsiall(kx,ky,i,j,1),Zpsiall(kx,ky,i,j,2)
c     &           ,Zpsiall(ky,kx,i,j,1)
            do ispin = 1,2
               write(10,'(8g20.10,a)') Eall(kx,ky,j,ispin)
     &              ,Zpsiall(kx,ky,i,j,ispin), Zpsiall(ky,kx,i,j,ispin)
            end do

c            if (Nqx == 1) then 
c               ikx = kx + Nkx / 2
c               write(10,'(4i5)') i, j, ikx, ky
c               write(10,'(4g20.10)')Eall(ikx,ky,j), 
c     &              Zpsiall(ikx,ky,i,j,1)
c            end if
         end do ; end do
      end do ; end do
c===============
      write(10,*) 'ENERGY_EIGEN_VALUE'
      emax = -100.
      emin = 100.
      do ispin = 1,2 ; do i = 1, Ns * Nqx
         do kx = 0, Nkx - 1 ; do ky = 0, Nky - 1
               write(10,'(4i5,2g20.10)') kx, ky, i, ispin,
     &           Eall(kx,ky,i,ispin)
               emax = MAX(Eall(kx,ky,i,ispin),emax)
               emin = MIN(Eall(kx,ky,i,ispin),emin)
         end do ; end do
      end do ; end do
      close(10)


 100   continue

      open(10,file="out.eplot3d.dat")
      dummy1 = 0.0d0
      do kx = 0, Nkx
         do ky = 0, Nky
            write(10,'(44g20.10)') kx/DBLE(Nkx), ky/DBLE(Nky), 
     &           ((Eall(kx,ky,i,ispin)-Dmu,i=1,Ns*Nqx),ispin=1,2),dummy1
         end do
         write(10,*)
      end do
      close(10)


      open(10,file="out.ene.dat")
      write(10,'(a)') '#ENERGY_EIGEN_VALUE'
      do ispin = 1,2 ; do i = 1, Ns * Nqx
         do kx = 0, Nkx - 1 ; do ky = 0, Nky - 1
               write(10,'(g20.10,i5)') Eall(kx,ky,i,ispin)-Dmu, 1
         end do ; end do
      end do ; end do
      close(10)

      if(Ns*Nqx >= 15) then
      open(10,file="out.ene.dat")
      write(10,'(a)') '#ENERGY_EIGEN_VALUE'
c      do ispin = 1,2 
      do kx = 0, Nkx-1
         do ky = 0, Nky - 1
               write(10,'(2g20.10,i5)') 
     &           ky/DBLE(Nky),Eall(kx,ky,15,1)-Dmu
         end do 
         write(10,*)
      end do
      close(10)
      end if
               emax = MAXVAL(Eall(:,:,1:Ns*Nqx,:))
               emin = MINVAL(Eall(:,:,1:Ns*Nqx,:))

c============================[
c== output data for gnuplot
c==
      open(20,file="out.eplot.dat")

      mx = Nkx / 2
      my = Nky / 2
      if (mx == 0) mx = 1
      if (my == 0) my = 1

c      e0 = 0.0d0
      e0 = Dmu
      do ispin = 1,2 ; do k = 1, Ns * Nqx
         do ikx = 0, mx - 1
            kx = ikx
            ky = 0
            p1 = ikx / REAL(mx)
            dummy1 =  Eall(kx,ky,k,ispin) - e0
            kx = kx + Nkx / 2
            ky = ky + Nky / 2
            kx = MOD(kx,Nkx)
            ky = MOD(ky,Nky)
            dummy2 =  Eall(kx,ky,k,ispin) - e0
            write(20,*) p1, dummy1, dummy2
         end do
         p1 = 1.0d0
         do iky = 0, my - 1
            kx = mx
            if(kx == Nkx) kx = 0
            ky = iky
            p2 = p1 + iky / REAL(my)
            dummy1 =  Eall(kx,ky,k,ispin) - e0
            kx = kx + Nkx / 2
            ky = ky + Nky / 2
            kx = MOD(kx,Nkx)
            ky = MOD(ky,Nky)
            dummy2 =  Eall(kx,ky,k,ispin) - e0
            write(20,*) p2, dummy1, dummy2
         end do
         p2 = p1 + 1.0d0
         do ikx = 0, mx - 1
            kx = mx - ikx
            if(kx == Nkx) kx = 0
            ky = kx
            if(Nkx /= Nky) ky = 0
            p3 = p2 + ikx / REAL(mx) * SQRT(2.0d0)
            dummy1 =  Eall(kx,ky,k,ispin) - e0
            kx = kx + Nkx / 2
            ky = ky + Nky / 2
            kx = MOD(kx,Nkx)
            ky = MOD(ky,Nky)
            dummy2 =  Eall(kx,ky,k,ispin) - e0
            write(20,*) p3, dummy1, dummy2
         end do
         p3 = p2 + SQRT(2.0d0)
         do iky = 0, my-1
            kx = 0
            ky = iky
            ky = MOD(ky + Nky,Nky)
            p4 = p3 + iky / REAL(my)
            dummy1 =  Eall(kx,ky,k,ispin) - e0
            kx = kx + Nkx / 2
            ky = ky + Nky / 2
            kx = MOD(kx,Nkx)
            ky = MOD(ky,Nky)
            dummy2 =  Eall(kx,ky,k,ispin) - e0
            write(20,*) p4, dummy1, dummy2
         end do
         p4 = p3 + 1.0d0
         do ikx = 0, mx
            kx = ikx
            ky = my - ikx
            ky = MOD(ky + Nky,Nky)
            if(Nkx /= Nky) ky = Nky
            p5 = p4 + ikx / REAL(mx) * SQRT(2.0)
            dummy1 =  Eall(kx,ky,k,ispin) - e0
            kx = kx + Nkx / 2
            ky = ky + Nky / 2
            kx = MOD(kx,Nkx)
            ky = MOD(ky,Nky)
            dummy2 =  Eall(kx,ky,k,ispin) - e0
            write(20,*) p5, dummy1, dummy2
         end do
         p5 = p4 + SQRT(2.0d0)
         write(20,1)
      end do ; end do
      
      lmtu = INT(emax - e0) + 1
      lmtl = INT(emin - e0) - 1
      
      write(20,*) "? ? ?", 0.0d0, lmtu
      write(20,*) "? ? ?", 0.0d0, lmtl
      write(20,1)
      write(20,*) "? ? ?", p1, lmtu
      write(20,*) "? ? ?", p1, lmtl
      write(20,1)
      write(20,*) "? ? ?", p2, lmtu
      write(20,*) "? ? ?", p2, lmtl
      write(20,1)
      write(20,*) "? ? ?", p3, lmtu
      write(20,*) "? ? ?", p3, lmtl
      write(20,1)
      write(20,*) "? ? ?", p4, lmtu
      write(20,*) "? ? ?", p4, lmtl
      write(20,1)
      write(20,*) "? ? ?", p5, lmtu
      write(20,*) "? ? ?", p5, lmtl
      write(20,1)
      close(20)

c================================]


      return
 1    format(i0)
      end
c######################################################################
c#########################################################
c#### Fermi surface
c#########################################################

      subroutine fermi2(filename, ene, np)

      use common, only : Nkx, Nky, Ns, Ns2, Pi, Nredx, Nredy, Nqx,
     &     Eall, Zpsiall, Dmu, Dne, Dtemp
      implicit none


      character, intent(in) :: filename*15
      integer, intent(out) :: np
      real*8, intent(in) :: ene

      character :: chara*80
      integer :: i, j, ikx, iky, kx, ky, kx1, ky1, kx2, ky2, ispin
      real*8 :: dummy1, dummy2
      real*8 :: ddd, dkx, dky, akx, aky

      dkx = 1.0d0 / DBLE(Nkx)
      dky = 1.0d0 / DBLE(Nky)

      open(10,file = filename)
      write(10,'(a)') '#fermi_surface'

      np = 0
      chara = '    '
      do ispin = 1, 2 ; do i = 1, Ns*Nqx
         j = i * 4
         do ikx = -Nkx/2, Nkx/2 ; do iky = -Nky/2, Nky/2
            akx = DBLE(ikx) / DBLE(Nkx)
            aky = DBLE(iky) / DBLE(Nky)

            kx1 = MOD(ikx+Nkx,Nkx)
            ky1 = MOD(iky+Nky,Nky)

            kx2 = MOD(kx1+1,Nkx)
            ky2 = ky1
            dummy1 =  Eall(kx1,ky1,i,ispin) - Dmu - ene
            dummy2 =  Eall(kx2,ky2,i,ispin) - Dmu - ene
            if (dummy1*dummy2 <= 0.0d0) then
               ddd = ABS(dummy1 / (dummy1 - dummy2)) * dkx
               write(10,'(a,2g20.10)') chara(1:j), akx + ddd, aky
               np = np + 1
            end if
            
            kx2 = kx1
            ky2 = MOD(ky1+1,Nky)
            dummy1 =  Eall(kx1,ky1,i,ispin) - Dmu - ene
            dummy2 =  Eall(kx2,ky2,i,ispin) - Dmu - ene
            if (dummy1*dummy2 <= 0.0d0) then
               ddd = ABS(dummy1 / (dummy1 - dummy2)) * dky
               write(10,'(a,2g20.10)') chara(1:i*4), akx, aky + ddd
               np = np + 1
            end if
         end do ; end do
         chara = chara(1:j)//'? ? '
      end do ; end do

      write(10,'(2a)') '10 10 10 10 10 10 10 10 10 10',
     &     ' 10 10 10 10 10 10 10 10 10 10'

      close(10)


      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### Fermi surface
c#########################################################

      subroutine fermi2k(filename, ene, eig,nx1,nx2,ny1,ny2,mkx,mky)

      use common, only :  Ns, Ns2, Pi, Nredx, Nredy, Nqx,
     &     Dmu, Dne, Dtemp
      implicit none


      character, intent(in) :: filename*15
      integer, intent(in) :: nx1, nx2, ny1, ny2, mkx, mky
      real*8, intent(in) :: ene
      real*8, intent(in) :: eig(nx1:nx2,ny1:ny2,Ns2*2)

      character :: chara*80
      integer :: i, j, ikx, iky, kx, ky, kx1, ky1, kx2, ky2
      real*8 :: dummy1, dummy2
      real*8 :: ddd, dkx, dky, akx, aky

      dkx = 1.0d0 / DBLE(mkx)
      dky = 1.0d0 / DBLE(mky)

      open(10,file = filename)
      write(10,'(a)') '#fermi_surface'
      write(10,'(a,g20.10)') '#Dmu=',Dmu



      chara = '    '
      do i = 1, Ns*Nqx
         j = i * 4
         do ikx = nx1, nx2 ; do iky = ny1, ny2
            akx = DBLE(ikx) / DBLE(mkx)
            aky = DBLE(iky) / DBLE(mky)

            kx1 = ikx
            ky1 = iky

            kx2 = kx1 + 1
            ky2 = ky1

            if ( (kx2 <= nx2).and.(ky2 <= ny2)) then
               
               dummy1 =  eig(kx1,ky1,i) - Dmu - ene
               dummy2 =  eig(kx2,ky2,i) - Dmu - ene
               if (dummy1*dummy2 <= 0.0d0) then
                  ddd = ABS(dummy1 / (dummy1 - dummy2)) * dkx
                  write(10,'(a,2g20.10)') chara(1:j), akx + ddd, aky
               end if
            end if
            
            kx2 = kx1
            ky2 = ky1 + 1
            if ( (kx2 <= nx2).and.(ky2 <= ny2)) then
               dummy1 =   eig(kx1,ky1,i) - Dmu - ene
               dummy2 =  eig(kx2,ky2,i) - Dmu - ene
               if (dummy1*dummy2 <= 0.0d0) then
                  ddd = ABS(dummy1 / (dummy1 - dummy2)) * dky
                  write(10,'(a,2g20.10)') chara(1:i*4), akx, aky + ddd
               end if
            end if

         end do ; end do
         chara = chara(1:j)//'? ? '
      end do

      write(10,'(a)') '#EOF'
      write(10,'(2a)') '10 10 10 10 10 10 10 10 10 10',
     &     ' 10 10 10 10 10 10 10 10 10 10'

      close(10)


      return
 1    format(i0)
      end
c######################################################################

c#########################################################
c#### Fermi surface
c#########################################################

      subroutine fermi3(filename,nk0,np)

      use common, only : Nkx, Nky, Ns, Ns2, Pi, Nredx, Nredy,
     &     Eall, Zpsiall, Dmu, Dne, Dtemp
      implicit none

      character, intent(in) :: filename*15
      integer, intent(in) :: nk0, np
      character :: chara*5
      integer :: i, kx0, ky0
      real*8 :: broadk, spectra, aax, aay
      real*8 :: akx(np), aky(np)
      real*8 :: dummyk


      broadk = 0.1d-2


      open(20,file = filename)

      i = 1
      do while(i /= 0)
         read(20,*)chara
         if (chara(1:1) /= '#') i = 0
      end do
      backspace(20)
      i = 0
      do i = 1, np
         read(20,*) akx(i), aky(i)
      end do
      close(20)

      open(10,file = "out.fermi2.dat")
      do ky0 = -nk0/2, nk0/2
         do kx0 = -nk0/2, nk0/2
            spectra = 0.0d0
            aax = DBLE(kx0) / DBLE(nk0)
            aay = DBLE(ky0) / DBLE(nk0)

            do i = 1, np
               dummyk = (akx(i)-aax)**2 + (aky(i)-aay)**2 + broadk**2
               spectra = spectra + 1.0d0 / dummyk
            end do
               
            write(10,*) aax, aay, spectra
         end do 
         write(10,1)
      end do

      close(10)


      return
 1    format(i0)
      end
c######################################################################
