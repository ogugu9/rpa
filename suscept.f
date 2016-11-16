c################ 2009 Kaneshita ################

c#########################################################
c#### suscept:
c#########################################################

      Subroutine suscept(iflag,chara)

      use common
      implicit none

      integer, intent(inout):: iflag

      character :: chara*5

      integer :: ikx, iky, kx1, kx2, ky1, ky2, ik, idummy
      integer :: kmax, kpi
      integer :: is, ir, ir1, ir2, k0x, k0y, iom
      integer :: mu1, mu2, mu3, mu4
      integer :: isE, nsize, nQ
      integer :: kflag(0:Nkx,0:Nky)
      integer, external :: Ispinpair, Ivrtxinv, Ivrtx, Iorbcomb
      integer, external :: Iorbselect

      real*8 :: dummy1, dummy2, om
      real*8 :: chixinter(0:Nkx,0:Nky)
      real*8 :: chixintra(0:Nkx,0:Nky)
      real*8 :: chix(0:Nkx,0:Nky,0:Ns,0:Ns)
      real*8 :: chixinter0(0:Nkx,0:Nky)
      real*8 :: chixintra0(0:Nkx,0:Nky)
      real*8 :: chix0(0:Nkx,0:Nky,0:Ns,0:Ns)

      real*8 :: chiz(0:Nkx,0:Nky)
      real*8 :: chizinter (0:Nkx,0:Nky)
      real*8 :: chizintra(0:Nkx,0:Nky)
      real*8 :: chin(0:Nkx,0:Nky)
      real*8 :: chininter(0:Nkx,0:Nky)
      real*8 :: chinintra(0:Nkx,0:Nky)
      complex*16 :: chi0(0:Nkx,0:Nky)

      complex*16 :: zdummy, zdummy1
c      complex*16 :: zchi(Nkx,Nky,2)
c      complex*16 :: zchi2(0:Nkx,0:Ns,0:Ns)
      complex*16 :: zmtrx1(1:Ns**2*Nqx, 1:Ns**2*Nqx, 2)
      complex*16 :: zmtrx2(1:Ns**2*Nqx, 1:Ns**2*Nqx, 2)

      real*8 :: ddtemp, ddmu, dummy


c#
c# output| chi data in zmtrx1
c#       | chi0 data in zmtrx2
c#
  
c==================== chemical potential ===[
c      call smuval(Dmu)
c      write(*,*)'Dmu(T=0)=',Dmu
c
c      ddtemp=0.0001d0
c      ddmu=1.0d-8
c
c      write(*,*)'calculating the chemical potential'
c      call Smuval2(ddmu,ddtemp)
c==================== chemical potential ===]


      nsize = Ns * Ns * Nqx
      Erange = MAXVAL(Eall) - MINVAL(Eall) + 1.0d0
      if (Erange < Erng) Erng = Erange


      call swfprod()
c===========================================[
c== chi(nn), chi(SzSz) =====================
c==
c# 
c# chi(nn) = chi(uu,uu) + chi(dd,dd) + chi(uu,dd) + chi(dd,uu) -> zmtrx2(1)
c# chi(zz) = chi(uu,uu) + chi(dd,dd) - chi(uu,dd) - chi(dd,uu) -> zmtrx2(2)
c#

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
      kx1 = kx1 * Nkx / Nqx     ! / 2
      ky1 = ky1 * Nky / Nqy     ! / 2

      kx2 = kx2 * Nkx / Nqx     !/ 2
      ky2 = ky2 * Nky / Nqy     !/ 2

      


      if (iflag == 0) then
         write(*,*) 'calculating susceptibility...'

         open(30,file='out.chi.dat',form='formatted')
         write(30,'(a)') '#Nkx__Nky'
         write(30,'(a,2i5)') '#', Nkx, Nky
         write(30,'(a)') '#Deta_Dtemp'
         write(30,'(a,2g15.8)') '#', Deta, Dtemp
         write(30,'(a)') '#Erange_Maxom'
         write(30,'(a,g20.10,i6)') '#', Erng, Maxom
         write(30,'(a,g20.10)') '#', Ominit
         write(30,'(a)') '#---chi0---'
         write(30,'(a,a)') '#path=', chara
         
         if (chara(1:1) == 'k') then
            iom = 1
            is = 1
            om = DBLE(iom) * Erng / DBLE(Maxom) + Ominit
            write(30,'(a,g20.10)') '#omega=', om
            do k0y = 0, Nky
               write(*,*) 'ky = ', k0y, '/', Nky/2
               do k0x = 0, Nkx
                  dummy1 = 0.0d0
                  ikx = k0x + Nkx / 2
                  iky = k0y + Nky / 2
                  ikx = MOD(ikx,Nkx)
                  iky = MOD(iky,Nky)
                  do mu1 = 1, Ns ; do nQ = 0, Nqx-1
                     ir1 = Ivrtx(mu1,mu1,nQ,'L')
                     call schi0_2(ir1,ir1,ikx,iky,iom,is,zdummy)
                     dummy1 = dummy1 - DIMAG(zdummy)
                  end do ; end do
                  write(30,*) k0x, k0y, dummy1
               end do
               write(30,*)
            end do
            close(30)
            write(*,*)'chi0 k-map calculated at E= ', om
         else
            is = 1
            do iom = 0, Maxom
               write(*,*)'chi0: om= ', iom,'/', Maxom,' chi0'
               om = DBLE(iom) * Erng / DBLE(Maxom) + Ominit
               do ik = 0, MAX(ABS(kx2-kx1),ABS(ky2-ky1))
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
                  ikx = MOD(ikx,Nkx)
                  iky = MOD(iky,Nky)
                  do mu1 = 1, Ns ; do nQ = 0, Nqx-1
                     ir1 = Ivrtx(mu1,mu1,nQ,'L')
                     call schi0_2(ir1,ir1,ikx,iky,iom,is,zdummy)
                     dummy1 = dummy1 - DIMAG(zdummy)
                  end do ; end do
                  
c                  if (ABS(dummy1) < 1.0d0-14) dummy1 = 0.0d0
                  write(30,*) ik, REAL(om), dummy1
               end do
               write(30,*) 
            end do
            close(30)
            write(*,*)
            write(*,*)'chi0 calculated along ', chara
         end if
         return
      end if
c=================================================================      
      if (iflag /= 2) then      ! longitudinal
         kmax = Nkx / Nqx / 2
         kmax = Nkx / Nqx
         kmax = MAX(ABS(kx2-kx1),ABS(ky2-ky1))

c         open(42,file='gomix.dat',form='formatted')
c         open(41,file='gomiy.dat',form='formatted')
         open(30,file='out.chil.dat',form='formatted')
         write(30,'(a)') '#Nkx__Nky'
         write(30,'(a,2i5)') '#', Nkx, Nky
         write(40,'(a)') '#U__J'
         write(40,'(a)') '#', U, DJ
         write(30,'(a)') '#Deta_Dtemp'
         write(30,'(a,2g15.8)') '#', Deta, Dtemp
         write(30,'(a)') '#Erange_Maxom'
         write(30,'(a,g20.10,i6)') '#', Erng, Maxom
         write(30,'(a,g20.10)') '#', Ominit

         if (chara /= 'k-map') then

         write(30,'(a,a)') '#path=', chara
         write(30,'(a)') 
     &        '#k om chiz chizintra chizinter chin chinintra chininter'
         
         do iom = Minom, Maxom
            kflag(:,:) = 0
            chiz(:,:) = 0.0d0
            chizintra(:,:) = 0.0d0
            chizinter(:,:) = 0.0d0
            chin(:,:) = 0.0d0
            chinintra(:,:) = 0.0d0
            chininter(:,:) = 0.0d0
            chi0(:,:) = 0.0d0

            om = DBLE(iom) * Erng / DBLE(Maxom) + Ominit
            write(*,*)'chi0: om= ',iom,'/',maxom,' chiz & chin'
            
            do ik = 0, kmax
               write(*,'(10x,a,i5,a,i5,a)')'k=',ik,'/',kmax,' chiz'
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
!$omp end parallel do
               end do
               
! susceptibility
               zmtrx1(:,:,:) = 0.0d0
               do is = 1, 2
                  call schi_2(is,zmtrx2,zmtrx1,nsize)
               end do
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
                        idummy = iorbcomb(mu1,mu2,mu3,mu4)
                        idummy = iorbselect(idummy)
                        if (idummy == 0) cycle

                        dummy1 = -DIMAG(zmtrx1(ir1,ir2,1)) 
     &                       +    DIMAG(zmtrx1(ir1,ir2,2))  
                        if (ABS(dummy1) < 1.0d-10) dummy1 = 0.0d0
                        chiz(ikx,iky) = chiz(ikx,iky) + dummy1
                        if (ir1 == ir2) then
                           chizintra(ikx,iky) = chizintra(ikx,iky) 
     &                          + dummy1
                        else
                           chizinter(ikx,iky) = chizinter(ikx,iky)
     &                          + dummy1
                        end if
                           
                        dummy1 = -DIMAG(zmtrx1(ir1,ir2,1)) 
     &                       -    DIMAG(zmtrx1(ir1,ir2,2))
                        if (ABS(dummy1) < 1.0d-10) dummy1 = 0.0d0
                           
                        chin(ikx,iky) = chin(ikx,iky) + dummy1
                        if (ir1 == ir2) then
                           chinintra(ikx,iky) = chinintra(ikx,iky) 
     &                          + dummy1
                        else
                           chininter(ikx,iky) = chininter(ikx,iky)
     &                          + dummy1
                        end if

                        dummy1 = -DBLE(zmtrx2(ir1,ir2,1)) 
     &                       -    DBLE(zmtrx2(ir1,ir2,2))
                        if (ABS(dummy1) < 1.0d-10) dummy1 = 0.0d0

                        dummy2 = -DIMAG(zmtrx2(ir1,ir2,1)) 
     &                       -    DIMAG(zmtrx2(ir1,ir2,2))
                        if (ABS(dummy2) < 1.0d-10) dummy2 = 0.0d0
                           
                        chi0(ikx,iky) = chi0(ikx,iky)
     &                       + DCMPLX(dummy1,dummy2)
                           
                     end do ; end do
                     
                  end do ; end do
               end do           ! nQ loop


               write(*,'(10x,a,6g18.8)')' chi =',
     &              chiz(k0x,k0y),
     &              chiz(MOD(k0x+Nkx/Nqx,Nkx),k0y),
     &              chin(k0x,k0y),
     &              chin(MOD(k0x+Nkx/Nqx,Nkx),k0y),
     &              DBLE(chi0(k0x,k0y)), DIMAG(chi0(k0x,k0y))

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
                                 
               write(30,'(i5,10g18.8)') ik, om, chiz(k0x,k0y),
     &              chizintra(k0x,k0y), chizinter(k0x,k0y), 
     &              chin(k0x,k0y),
     &              chinintra(k0x,k0y), chininter(k0x,k0y),
     &              DBLE(chi0(k0x,k0y)), DIMAG(chi0(k0x,k0y))
            end do
            write(30,*)
         end do                 ! DO LOOP for Omega

c--------------------------------------------------------------
         else if (chara == 'k-map') then
            kflag(:,:) = 0
            chiz(:,:) = 0.0d0
            chizintra(:,:) = 0.0d0
            chizinter(:,:) = 0.0d0
            chin(:,:) = 0.0d0
            chinintra(:,:) = 0.0d0
            chininter(:,:) = 0.0d0
            chi0(:,:) = 0.0d0
            iom = Minom
            om = DBLE(iom) * Erng / DBLE(Maxom) + Ominit

            write(30,'(a,x,f16.6)') '#omega=', om
            write(30,'(a)') 
     &       '#k om chiz chizintra chizinter chin chinintra chininter'
         
            do k0x = 0, Nkx/2 ; do k0y = 0, Nky/2
               write(*,'(10x,a,2i5,a,2i5)')'k=',k0x,k0y,'/',Nkx,Nky

               do is = 1, 2
                  isE = Ispinpair(is,'E')

!$omp parallel do private(ir2,zdummy)
                  do ir1 = 1, nsize ; do ir2 = 1, nsize
                     call schi0_2(ir1,ir2,k0x,k0y,iom,is,zdummy)
                     zmtrx2(ir1,ir2,isE) = zdummy
                  end do ; end do
!$omp end parallel do
               end do
               
! susceptibility
               zmtrx1(:,:,:) = 0.0d0
               do is = 1, 2
                  call schi_2(is,zmtrx2,zmtrx1,nsize)
               end do
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
                        idummy = iorbcomb(mu1,mu2,mu3,mu4)
                        idummy = iorbselect(idummy)
                        if (idummy == 0) cycle

                        dummy1 = -DIMAG(zmtrx1(ir1,ir2,1)) 
     &                       +    DIMAG(zmtrx1(ir1,ir2,2))  
                        if (ABS(dummy1) < 1.0d-10) dummy1 = 0.0d0
                        chiz(ikx,iky) = chiz(ikx,iky) + dummy1
                        if (ir1 == ir2) then
                           chizintra(ikx,iky) = chizintra(ikx,iky) 
     &                          + dummy1
                        else
                           chizinter(ikx,iky) = chizinter(ikx,iky)
     &                          + dummy1
                        end if
                           
                        dummy1 = -DIMAG(zmtrx1(ir1,ir2,1)) 
     &                       -    DIMAG(zmtrx1(ir1,ir2,2))
                        if (ABS(dummy1) < 1.0d-10) dummy1 = 0.0d0
                        chin(ikx,iky) = chin(ikx,iky) + dummy1
                        if (ir1 == ir2) then
                           chinintra(ikx,iky) = chinintra(ikx,iky) 
     &                          + dummy1
                        else
                           chininter(ikx,iky) = chininter(ikx,iky)
     &                          + dummy1
                        end if

                        dummy1 = -DBLE(zmtrx2(ir1,ir2,1)) 
     &                       -    DBLE(zmtrx2(ir1,ir2,2))
                        if (ABS(dummy1) < 1.0d-10) dummy1 = 0.0d0

                        dummy2 = -DIMAG(zmtrx2(ir1,ir2,1)) 
     &                       -    DIMAG(zmtrx2(ir1,ir2,2))
                        if (ABS(dummy2) < 1.0d-10) dummy2 = 0.0d0
                           
                        chi0(ikx,iky) = chi0(ikx,iky)
     &                       + DCMPLX(dummy1,dummy2)
                           
                     end do ; end do !mu3, mu4 loop                    
                  end do ; end do !mu1, mu2 loop
               end do           ! nQ loop

               if (ABS(chiz(k0x,k0y)) < 1.0d-10) then
                  chiz(k0x,k0y)=0.0d0
               end if

               write(30,'(i5,10g18.8,i5)') k0x, k0y, chiz(k0x,k0y),
     &              chizintra(k0x,k0y), chizinter(k0x,k0y), 
     &              chin(k0x,k0y),
     &              chinintra(k0x,k0y), chininter(k0x,k0y),
     &              DBLE(chi0(k0x,k0y)), DIMAG(chi0(k0x,k0y))
            end do ; write(30,*) ; end do
         end if
c================== chizz & chinn ==========]

         close(30)
c============ chi(nn), chi(SzSz) ===========]

      end if
c====================================================================
      if (iflag /= 1) then      ! transverse

c===========================================[
c== chi(xx) ================================
c==
c# 
c# chi(xx) = chi(ud,du) + chi(du,ud) -> zmtrx2(1)
c#

         kmax = Nkx / Nqx / 2
         kmax = Nkx / Nqx
         kmax = MAX(ABS(kx2-kx1),ABS(ky2-ky1))

         open(42,file='gomix.dat',form='formatted')
         open(41,file='gomiy.dat',form='formatted')
         open(40,file='out.chix.dat',form='formatted')
         write(40,'(a)') '#Nkx__Nky'
         write(40,'(a,2i5)') '#', Nkx, Nky
         write(40,'(a)') '#U__J'
         write(40,'(a,2g18.8)') '#', U, DJ
         write(40,'(a)') '#Deta_Dtemp'
         write(30,'(a,2g15.8)') '#', Deta, Dtemp
         write(40,'(a)') '#Erange_Maxom'
         write(40,'(a,g20.10,i6)') '#', Erng, Maxom
         write(40,'(a,g20.10)') '#', Ominit
c         write(40,'(a)') '#---chixx--'
c         write(40,'(a,a)') '#path=', chara
c         write(40,'(a)') 
c     &        '#k om chix chixintra chixinter chix0 intra0 inter0'
         

         if (chara /= 'k-map') then

         write(40,'(a,a)') '#path=', chara
         write(40,'(a,2i4,a,2i4,a)') '#(',kx1,ky1,')-(',kx2,ky2,')'
         write(40,'(a)') 
     &        '#k om chix chixintra chixinter chix0 intra0 inter0'
         
         do iom = Minom, Maxom
            kflag(:,:) = 0
            chix(:,:,:,:) = 0.0d0
            chixintra(:,:) = 0.0d0
            chixinter(:,:) = 0.0d0
            chix0(:,:,:,:) = 0.0d0
            chixintra0(:,:) = 0.0d0
            chixinter0(:,:) = 0.0d0
            chi0(:,:) = 0.0d0
c            zmtrx1(:,:,:) = 0.0d0
c            zmtrx2(:,:,:) = 0.0d0
            om = DBLE(iom) * Erng / DBLE(Maxom) + Ominit
            write(*,*)'chi0: om= ',iom,'/',maxom,' chixx'
            

            do ik = 0, kmax
               write(*,'(10x,a,i5,a,i5,a)')'k=',ik,'/',kmax,' chixx'
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

c            k0y = 0
c            do k0x = 0, kmax-1
c            write(*,'(10x,a,i5,a,i5,a)')'k0x=',k0x,'/',kmax,' chixx'
               
c==================chi for chixx ==========[
!     bare susceptibility
               
               do is = 3, 3
                  isE = Ispinpair(is,'E')

!$omp parallel do private(ir2,zdummy)
                  do ir1 = 1, nsize ; do ir2 = 1, nsize
                     call schi0_2(ir1,ir2,k0x,k0y,iom,is,zdummy)
                     zmtrx2(ir1,ir2,isE) = zdummy
                  end do ; end do
!$omp end parallel do

                  call schix2(is,zmtrx2,zmtrx1,nsize)

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
                           idummy = iorbcomb(mu1,mu2,mu3,mu4)
                           idummy = iorbselect(idummy)
                           if (idummy == 0) cycle

                           dummy1 = -DIMAG(zmtrx1(ir1,ir2,isE))  
                           if (ABS(dummy1) < 1.0d-10) dummy1 = 0.0d0
                           
                           chix(ikx,iky,0,mu2) = chix(ikx,iky,0,mu2)
     &                          + dummy1
                           chix(ikx,iky,mu1,0) = chix(ikx,iky,mu1,0)
     &                          + dummy1
                           chix(ikx,iky,0,0) = chix(ikx,iky,0,0)
     &                          + dummy1
                           if (ir1 == ir2) then
                              chixintra(ikx,iky) = chixintra(ikx,iky) 
     &                             + dummy1
                           else
                              chixinter(ikx,iky) = chixinter(ikx,iky)
     &                             + dummy1
                           end if
                           
                           dummy1 = -DIMAG(zmtrx2(ir1,ir2,isE))
                           if (ABS(dummy1) < 1.0d-10) dummy1 = 0.0d0
                           
                           chix0(ikx,iky,0,mu2) = chix0(ikx,iky,0,mu2)
     &                          + dummy1
                           chix0(ikx,iky,mu1,0) = chix0(ikx,iky,mu1,0)
     &                          + dummy1
                           chix0(ikx,iky,0,0) = chix0(ikx,iky,0,0)
     &                          + dummy1
                           if (ir1 == ir2) then
                              chixintra0(ikx,iky) = chixintra0(ikx,iky) 
     &                             + dummy1
                           else
                              chixinter0(ikx,iky) = chixinter0(ikx,iky)
     &                             + dummy1
                           end if

                        dummy1 = -DBLE(zmtrx2(ir1,ir2,isE)) 
                        if (ABS(dummy1) < 1.0d-14) dummy1 = 0.0d0

                        dummy2 = -DIMAG(zmtrx2(ir1,ir2,isE))
                        if (ABS(dummy2) < 1.0d-14) dummy2 = 0.0d0
                           
                        chi0(ikx,iky) = chi0(ikx,iky)
     &                       + DCMPLX(dummy1,dummy2)

                        end do ; end do
                        
                     end do ; end do
                  end do        ! nQ loop
               end do           ! is loop


               write(*,'(10x,a,4g18.8)')' chi =',
     &              chix(k0x,k0y,0,0), chix0(k0x,k0y,0,0),
     &              DBLE(chi0(k0x,k0y)), DIMAG(chi0(k0x,k0y))

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
                                 
               write(40,'(i5,10g18.8)') ik, om, chix(k0x,k0y,0,0),
     &              chixintra(k0x,k0y), chixinter(k0x,k0y), 
     &              chix0(k0x,k0y,0,0),
     &              chixintra0(k0x,k0y), chixinter0(k0x,k0y),
     &              DBLE(chi0(k0x,k0y)), DIMAG(chi0(k0x,k0y))
            end do
            write(40,*)
         end do
c--------------------------------------------------------------
         else if (chara == 'k-map') then
            kflag(:,:) = 0
            chix(:,:,:,:) = 0.0d0
            chixintra(:,:) = 0.0d0
            chixinter(:,:) = 0.0d0
            chix0(:,:,:,:) = 0.0d0
            chixintra0(:,:) = 0.0d0
            chixinter0(:,:) = 0.0d0
            chi0(:,:) = 0.0d0
            iom = Minom
            om = DBLE(iom) * Erng / DBLE(Maxom) + Ominit

            write(40,'(a,x,f16.6)') '#omega=', om
            write(40,'(a)') 
     &           '#kx ky chix chixintra chixinter chix0 intra0 inter0'

            
            do k0x = 0, Nkx/2 ; do k0y = 0, Nky/2
               write(*,'(10x,a,2i5,a,2i5)')'k=',k0x,k0y,'/',Nkx,Nky
               do is = 3, 3
                  isE = Ispinpair(is,'E')
                     
!$omp parallel do private(ir2,zdummy)
                  do ir1 = 1, nsize ; do ir2 = 1, nsize
                     call schi0_2(ir1,ir2,k0x,k0y,iom,is,zdummy)
                     zmtrx2(ir1,ir2,isE) = zdummy
                  end do ; end do
!$omp end parallel do
                  
                  call schix2(is,zmtrx2,zmtrx1,nsize)
                  
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
                           idummy = iorbcomb(mu1,mu2,mu3,mu4)
                           idummy = iorbselect(idummy)
                           if (idummy == 0) cycle

                           dummy1 = -DIMAG(zmtrx1(ir1,ir2,isE))  
                           if (ABS(dummy1) < 1.0d-14) dummy1 = 0.0d0
                           
                           chix(ikx,iky,0,mu2) = chix(ikx,iky,0,mu2)
     &                          + dummy1
                           chix(ikx,iky,mu1,0) = chix(ikx,iky,mu1,0)
     &                          + dummy1
                           chix(ikx,iky,0,0) = chix(ikx,iky,0,0)
     &                          + dummy1
                           if (ir1 == ir2) then
                              chixintra(ikx,iky) = chixintra(ikx,iky) 
     &                             + dummy1
                           else
                              chixinter(ikx,iky) = chixinter(ikx,iky)
     &                             + dummy1
                           end if
                           
                           dummy1 = -DIMAG(zmtrx2(ir1,ir2,isE))
                           if (ABS(dummy1) < 1.0d-14) dummy1 = 0.0d0
                           
                           chix0(ikx,iky,0,mu2) = chix0(ikx,iky,0,mu2)
     &                          + dummy1
                           chix0(ikx,iky,mu1,0) = chix0(ikx,iky,mu1,0)
     &                          + dummy1
                           chix0(ikx,iky,0,0) = chix0(ikx,iky,0,0)
     &                          + dummy1
                           if (ir1 == ir2) then
                              chixintra0(ikx,iky) = chixintra0(ikx,iky) 
     &                             + dummy1
                           else
                              chixinter0(ikx,iky) = chixinter0(ikx,iky)
     &                             + dummy1
                           end if

                        dummy1 = -DBLE(zmtrx2(ir1,ir2,isE)) 
                        if (ABS(dummy1) < 1.0d-14) dummy1 = 0.0d0

                        dummy2 = -DIMAG(zmtrx2(ir1,ir2,isE))
                        if (ABS(dummy2) < 1.0d-14) dummy2 = 0.0d0
                           
                        chi0(ikx,iky) = chi0(ikx,iky)
     &                       + DCMPLX(dummy1,dummy2)



                        end do ; end do
                        
                     end do ; end do
                  end do        ! nQ loop
               end do           ! is loop

               if (ABS(chix(k0x,k0y,0,0)) < 1.0d-14) then
                  chix(k0x,k0y,0,0)=0.0d0
               end if

               write(40,'(2i5,10g18.8)') k0x, k0y, chix(k0x,k0y,0,0),
     &              chixintra(k0x,k0y), chixinter(k0x,k0y), 
     &              chix0(k0x,k0y,0,0),
     &              chixintra0(k0x,k0y), chixinter0(k0x,k0y),
     &              DBLE(chi0(k0x,k0y)), DIMAG(chi0(k0x,k0y))
            end do ; write(40,*) ; end do
         end if
c================== chi for chixx ==========]
c               end do ; end do


         
         close(40)
         close(41)
         close(42)
      end if

      return

c===========================================]

 999  return
      end

c########################################### suscept ####

c#########################################################
c#### old
c#### schi0: chi0
c####        output bare susceptibility -> Zmtrx2
c####
c####        Zpsiall, Eall are required
c#########################################################

      Subroutine schi0(ir1,ir2,k0x,k0y,iom,is,zmtrx2)

      use common, only : Ns, Ns2, Nkx, Nky, Zpsiall, Eall, ZI,
     &     Maxom, Ominit, Erng, Deta, Dmu, Dtemp, Pi, Nqx
      implicit none

      integer, intent(in) :: ir1, ir2, k0x, k0y, iom, is
      complex*16, intent(out):: zmtrx2(1:ns2,1:ns2,2)

      integer :: isE, isH
      integer :: iorb1, iorb2, iorb3, iorb4, iq
      integer :: jom, kk0x, kk0y, ik0x, ik0y, n, m
      integer, external :: Ispinpair

      real*8 :: akx, aky
      real*8 :: om
      real*8, external :: Dffn

      complex*16 :: zfacx, zfacy

                                ! Zmtrx2: chi0 (output)
c#
c#  is = 1:up-up, 2:dn-dn, 3:up-dn, 4:dn-up
c#  chi0(uu):1, chi0(dd):2 for <<n;n>>, <<Sz;Sz>>
c#  chi0(ud):1, chi0(du):2 for <<Sx;Sx>>, <<Sy;Sy>>
c#

      isE = Ispinpair(is,'E')
      isH = Ispinpair(is,'H')

      zmtrx2(ir1,ir2,isE) = DCMPLX(0.0d0,0.0d0)

      akx = 2.0d0 * Pi / DBLE(Nkx)
      aky = 2.0d0 * Pi / DBLE(Nky)

                                !omega
      om = DBLE(iom) * Erng / DBLE(Maxom) + Ominit

      do kk0x = 0,Nkx-1 ; do kk0y = 0,Nky-1
         ik0x = k0x + kk0x
         ik0y = k0y + kk0y


         zfacx = 1.0d0
         zfacy = 1.0d0
c         if (ik0x >= Nkx) then
c            zfacx= EXP( Zi * akx * (icx(ir1) - icx(ir2)) )
c         end if
c         if (ik0y >= Nky) then
c            zfacy= EXP( Zi * aky * (icy(ir1) - icy(ir2)) )
c         end if

         ik0x = MOD(ik0x,Nkx)
         ik0y = MOD(ik0y,Nky)

         do n = 1, Ns * Nqx ; do m = 1, Ns * Nqx
            do iq = 1, Nqx
               iorb1 = ir1
               iorb2 = MOD(ir2 + (iq-1) * Ns, Ns * Nqx)
               if (iorb2 == 0) iorb2 = Ns * Nqx
               iorb3 = MOD(MOD(ir2,Ns) +  (iq-1) * Ns, Ns * Nqx)
               if (iorb3 == 0) iorb3 = Ns * Nqx
               iorb4 = MOD(ir1,Ns)
               if (iorb4 == 0) iorb4 = Ns * Nqx
               zmtrx2(ir1,ir2,isE)
     &              = zmtrx2(ir1,ir2,isE)
     &              -        Zpsiall(ik0x,ik0y,iorb2,n,isE)
     &              * DCONJG(Zpsiall(ik0x,ik0y,iorb1,n,isE))
     &              *        Zpsiall(kk0x,kk0y,iorb4,m,isH)
     &              * DCONJG(Zpsiall(kk0x,kk0y,iorb3,m,isH))
c     &           * zfacx * zfacy
     &              * ( Dffn(Eall(ik0x,ik0y,n,isE)-Dmu,Dtemp)
     &              -   Dffn(Eall(kk0x,kk0y,m,isH)-Dmu,Dtemp) )
     &              / ( om - (Eall(ik0x,ik0y,n,isE)
     &              -         Eall(kk0x,kk0y,m,isH))
     &              + Zi * Deta )
     &              /DBLE(Nkx * Nky)
            end do
         end do ; end do
      end do ; end do

      return
      end
c############################################## schi0 ####


c#########################################################
c#### schi02: chi0
c####        output bare susceptibility -> Zmtrx2
c####
c####        Zpsiall, Eall are required
c#########################################################

      Subroutine schi0_2(ir12,ir34,k0x,k0y,iom,is,zmtrx2)

      use common, only : Ns, Ns2, Nkx, Nky, Zpsiall, Eall, ZI,
     &     Maxom, Ominit, Erng, Deta, Dmu, Dtemp, Pi, Nqx,
     &     Zwfprod
      implicit none

      integer, intent(in) :: ir12, ir34
      integer, intent(in) :: k0x, k0y, iom, is
      complex*16, intent(out) :: zmtrx2

      integer :: isE, isH
      integer :: ip0x, ip0y, ik0x, ik0y, n, m
c      integer :: nx, ny, nperi, nss
      integer :: mu1, mu2, mu3, mu4, nQ1, nQ2, mQ
      integer :: iorb1, iorb2, iorb3, iorb4, iq, npQ
      integer, external :: Ispinpair, Ivrtxinv

      real*8 :: om, omin
c      real*8 :: ene(0:Nkx,0:Nky,Ns2), ffn(0:Nkx,0:Nky,Ns2)
c      real*8 :: temp, eta, omax
      real*8, external :: Dffn
      
      complex*16 :: zdummy, zzz
c      complex*16 :: zimg


                                ! Zmtrx2: chi0 (output)
c#
c#  is = 1:up-up, 2:dn-dn, 3:up-dn, 4:dn-up
c#  chi0(uu):1, chi0(dd):2 for <<n;n>>, <<Sz;Sz>>
c#  chi0(ud):1, chi0(du):2 for <<Sx;Sx>>, <<Sy;Sy>>
c#
c==== .   chi0(k+Q1, k+Q2) 
c==== .                
c==== .                  k+p+Q1  k+p+Q2+Q'         n, s
c==== .                   2  -->--  3
c==== . chi0(1,2,3,4) =    <       >        (SUM: p, Q', n, m)
c==== .                   1  --<--  4
c==== .                    p     p+Q'             m, s'
c==== . 
c==== . 1, 2, ...:  orbital, i.e., mu1, mu2,...
c==== .      n, m:  band index 
c==== .     s, s':  spin index
c==== .         k:  chi0(k+Q1, k+Q2)
c==== .
c==== .     p = p0 + Qp (Qp = npQ * )
c==== . 
c==== .             k+p0+Q1+Qp   k+p0+Q2+Q'+Qp         n, s
c==== .                   2  -->--  3
c==== . chi0(1,2,3,4) =    <      >        (SUM: p, Q', n, m)
c==== .                   1  --<--  4
c==== .                p0+Qp      p0+Q'+Qp             m, s'
c==== .

      mu1 = Ivrtxinv(ir12,1,'L')
      mu2 = Ivrtxinv(ir12,2,'L')
      nQ1 = Ivrtxinv(ir12,3,'L')

      mu3 = Ivrtxinv(ir34,1,'R')
      mu4 = Ivrtxinv(ir34,2,'R')
      nQ2 = Ivrtxinv(ir34,3,'R')


      isE = Ispinpair(is,'E')     !s
      isH = Ispinpair(is,'H')     !s' 

      zzz = DCMPLX(0.0d0,0.0d0)
           
                                !omega
      om = DBLE(iom) * Erng / DBLE(Maxom) + Ominit


c!$omp parallel do reduction(+:zzz) private(ik0x,ik0y,zdummy,
c!$omp1  iorb1, iorb2, iorb3, iorb4, mQ, n, m, iq, npQ, ip0y)
      do ip0x = 0, Nkx / Nqx - 1 ; do ip0y = 0, Nky - 1 ! p0
         ik0x = k0x + ip0x      ! k+p0
         ik0y = k0y + ip0y      ! k+p0
         
         ik0x = MOD(ik0x, Nkx)
         ik0y = MOD(ik0y, Nky)
         
         do n = 1, Ns * Nqx ; do m = 1, Ns * Nqx
            zdummy = ( Dffn(Eall(ik0x,ik0y,n,isE)-Dmu,Dtemp)
     &           -   Dffn(Eall(ip0x,ip0y,m,isH)-Dmu,Dtemp) )
     &           / ( om + Zi * Deta
     &           - (Eall(ik0x,ik0y,n,isE) - Eall(ip0x,ip0y,m,isH)))

c            if ((ir12 == 1).and.(ir34 == 2))
c     &                 write(*,'(2i5,2g18.5,a)')n,m,zdummy,
c     &                 'aaaaaaaaaa'
            if (ABS(zdummy) < 1.0d-13) cycle

            do npQ = 0, Nqx - 1     !Qp

               mQ = npQ         !Qp
               iorb1 = mu1 + mQ * Ns
               
               mQ = nQ1 + npQ   !Q1+Qp
               iorb2 = mu2 + MOD(mQ, Nqx) * Ns
               
               do iq = 0, Nqx - 1 !Q'
                  mQ = nQ2 + iq + npQ !Q2+Q'+Qp
                  iorb3 = mu3 + MOD(mQ, Nqx) * Ns

                  mQ = iq + npQ !Q'+Qp
                  iorb4 = mu4 + MOD(mQ, Nqx) * Ns

                  zzz = zzz
     &                 - Zwfprod(ik0x,ik0y,iorb2,iorb3,n,isE)
     &                 * Zwfprod(ip0x,ip0y,iorb4,iorb1,m,isH)
     &                 * zdummy

c            if ((ir12 == 1).and.(ir34 == 2))
c     &                 write(*,*) Zwfprod(ik0x,ik0y,iorb2,iorb3,n,isE)
c     &                 * Zwfprod(ip0x,ip0y,iorb4,iorb1,m,isH),
c     &                 'aaaaaaaaaa'
               end do
            end do
         end do ; end do
      end do ; end do
      
c            if ((ir12 == 1).and.(ir34 == 2))write(*,*)zzz,
c     &     'aaa++++++++'
      zmtrx2 = zzz / DBLE(Nkx * Nky)
      if(ABS(zmtrx2) < 1.0d-13) zmtrx2 = 0.0d0

      return
      end
c############################################## schi0 ####


c#########################################################
c#### schi: chi -> zmtrx1
c#########################################################
      subroutine schi_2(is,zmtrx2,zmtrx1,nsize)

      use common, only : Ns, Ns2, U
      implicit none
      
      integer, intent(in) :: is, nsize
      complex*16, intent(inout):: zmtrx1(1:nsize,1:nsize,2)
      complex*16, intent(in):: zmtrx2(1:nsize,1:nsize,2)

      integer :: ir12, ir34, ir56, is1, is2
      integer :: info, ipiv(nsize*3)

      real*8, external :: Wint
      complex*16 :: za(nsize*2,nsize*2), zb(nsize*2,nsize)
           
c==== .   
c==== . chiuu(1,2;3,4)= chi0u(1,2;3,4)
c==== .                + chi0u(1,2;5,6) W(5,6;7,8) chiuu(7,8;3,4)
c==== .                + chi0u(1,2;5,6) W(5,6;7,8) chidu(7,8;3,4)
c==== . 
c==== . chidu(1,2;3,4)=  chi0d(1,2;5,6) W(5,6;7,8) chiuu(7,8;3,4)
c==== .                + chi0d(1,2;5,6) W(5,6;7,8) chidu(7,8;3,4)
c==== .
c==== . chi0u(1,2;3,4)= chi0u(1,2;3,4)       
c==== . 
c==== . /       \   /                          \  /       \   /       \
c==== . | chiuu |   | chi0u*W(uu)  chi0u*W(ud) |  | chiuu |   | chi0u |
c==== . |       | = |                          |  |       | + |       |
c==== . | chidu |   | chi0d*W(du)  chi0d*W(dd) |  | chidu |   |   0   |
c==== . \       /   \                          /  \       /   \       /
c==== . 
c==== .        x  = M x  + b 
c==== .    
c==== .        (1-M)x = Ax = b 
c==== . 
c==== .                                                         
c==== .       /               \                                 
c==== .       |    up   up    |                       up   up    
c==== .       |  2  -->--  3  |                     2  -->--  3  
c==== .       |   <///////>   |   chiuu(1,2;3,4) =   <///////>   
c==== .       |  1  --<--  4  |                     1  --<--  4  
c==== .       |    up   up    |                       up   up    
c==== .  x =  |               |                                    
c==== .       |    dn   up    |                       dn   up     
c==== .       |  2  -->--  3  |                     2  -->--  3   
c==== .       |   <///////>   |   chidu(1,2;3,4) =   <///////>    
c==== .       |  1  --<--  4  |                     1  --<--  4   
c==== .       |    dn   up    |                       dn   up     
c==== .       \               / 
c==== .                                       
c==== .       /               \ 
c==== .       |      up       |                           up       
c==== .       |  2  -->--  3  |                       2  -->--  3  
c==== .       |   <       >   |     chi0u(1,2;3,4) =   <       >   
c==== .       |  1  --<--  4  |                       1  --<--  4  
c==== .       |      up       |                            up      
c==== .  b =  |               | 
c==== .       |               |  
c==== .       |               |  
c==== .       |       0       |  
c==== .       |               |  
c==== .       |               |  
c==== .       \               / 
c==== .                                    up
c==== .                                 2 -->-- 3      6 
c==== .  chi0u(1,2,3,4) W(3,4;5,6)  =  <         >~~~~~  
c==== .                                 1 --<-- 4      5 
c==== .                                    up
c==== .   
c==== .               up                    dn 
c==== .             1-->--   W(1,2;3,4)   -->--4
c==== .                   \              /
c==== .  W(ud)  =          >~~~~~~~~~~~~<  
c==== .                   /              \
c==== .             2--<--                --<--3
c==== .               up                    dn  
c==== .

      if(is == 1)then
         is1 = 1; is2 = 2       ! (up,up)---(dn,dn)
      else if(is == 2)then
         is1 = 2; is2 = 1       ! (dn,dn)---(up,up)
      else
         write(*,*)"error: is(=spin) should be 1 or 2:",is
         stop
      end if

      za(:,:) = DCMPLX(0.0d0,0.0d0)

      do ir12 = 1, nsize 
         za(ir12      ,ir12      ) =  1.0d0
         za(ir12+nsize,ir12+nsize) =  1.0d0
         do ir56 = 1, nsize
            
            do ir34 = 1, nsize
               za(ir12      ,ir56      ) = za(ir12      ,ir56      )
     &              - zmtrx2(ir12,ir34,is1)
     &              * Wint(ir34,ir56,is1,is1,is1,is1)

               za(ir12      ,ir56+nsize) = za(ir12      ,ir56+nsize)
     &              - zmtrx2(ir12,ir34,is1)
     &              * Wint(ir34,ir56,is1,is1,is2,is2)

               za(ir12+nsize,ir56      ) = za(ir12+nsize,ir56      )
     &              - zmtrx2(ir12,ir34,is2)
     &              * Wint(ir34,ir56,is2,is2,is1,is1)

               za(ir12+nsize,ir56+nsize) = za(ir12+nsize,ir56+nsize)
     &              - zmtrx2(ir12,ir34,is2)
     &              * Wint(ir34,ir56,is2,is2,is2,is2)
            end do
                        
         end do 
      end do

      zb(      1:nsize  ,1:nsize) =  zmtrx2(1:nsize  ,1:nsize,is1)
      zb(nsize+1:nsize*2,1:nsize) =  0.0d0
      call ZGESV(nsize*2, nsize, za, nsize*2, ipiv, zb, nsize*2, INFO)


      if(INFO.ne.0) then
         write(*,*) 'Lapack ZGESV: Info=',Info
         stop
      end if

      zmtrx1(1:nsize,1:nsize,1) = zmtrx1(1:nsize,1:nsize,1)
     &     + zb(      1:nsize  ,1:nsize)
      zmtrx1(1:nsize,1:nsize,2) =zmtrx1(1:nsize,1:nsize,2)
     &     +  zb(nsize+1:nsize*2,1:nsize)

      return
      end
c################################################## schi #


c#########################################################
c#### Swfprod: productions of wave functions
c#########################################################

      Subroutine Swfprod()

      use common, only : Ns, Nkx, Nky, Zpsiall, Zwfprod, Nqx
      implicit none

      logical :: torf
      integer :: ikx, iky, nb
      integer :: iorb1, iorb2, is, n

      complex*16 :: zdummy

c==== .   psi(1)^dagger * psi(2) 
c==== .            
c==== .        1-->--2,  2--<--1


      nb = Ns * Nqx
      torf = allocated(Zwfprod)
      if (torf .EQV. .true.) deallocate(Zwfprod)
      allocate( Zwfprod(0:Nkx,0:Nky,nb,nb,nb,2) )

c!$omp parallel do private(iky, iorb1, iorb2)
      do ikx = 0, Nkx  ; do iky = 0, Nky
            do iorb1 = 1, nb ; do iorb2 = 1, nb

            do n = 1, nb ; do is = 1, 2
                  zdummy
     &                 =        Zpsiall(ikx,iky,iorb2,n,is)
     &                 * DCONJG(Zpsiall(ikx,iky,iorb1,n,is))
                  if (ABS(DBLE(zdummy)) < 1.0d-14) 
     &                 zdummy = DCMPLX(0.0d0,DIMAG(zdummy))
                  if (ABS(DIMAG(zdummy)) < 1.0d-14) 
     &                 zdummy = DCMPLX(DBLE(zdummy),0.0d0)

                  Zwfprod(ikx,iky,iorb1,iorb2,n,is) = zdummy
               end do ; end do

c                  Zwfprod(ikx,iky,iorb1,iorb2,1:nb,1:2)
c     &                 =        Zpsiall(ikx,iky,iorb2,1:nb,1:2)
c     &                 * DCONJG(Zpsiall(ikx,iky,iorb1,1:nb,1:2))
         end do ; end do
      end do ; end do

      return
      end
c############################################## schi0 ####

c#########################################################
c#### Swfprod: productions of wave functions
c#########################################################

      Subroutine Swfprod_orb()

      use common, only : Ns, Nkx, Nky, Zpsiall, Zwfprod, Nqx
      implicit none

      logical :: torf
      integer :: ikx, iky, nb
      integer :: iorb1, iorb2, is, n

      complex*16 :: zdummy
      complex*16, external :: zwvfn
      

c==== .   psi(1)^dagger * psi(2) 
c==== .            
c==== .        1-->--2,  2--<--1


      nb = Ns * Nqx
      torf = allocated(Zwfprod)
      if (torf .EQV. .true.) deallocate(Zwfprod)
      allocate( Zwfprod(0:Nkx,0:Nky,nb,nb,nb,2) )

c!$omp parallel do private(iky, iorb1, iorb2)
      do ikx = 0, Nkx  ; do iky = 0, Nky
            do iorb1 = 1, nb ; do iorb2 = 1, nb

            do n = 1, nb ; do is = 1, 2
                  zdummy
     &                 =        zwvfn(ikx,iky,iorb2,n,is)
     &                 * DCONJG(zwvfn(ikx,iky,iorb1,n,is))
                  if (ABS(DBLE(zdummy)) < 1.0d-14) 
     &                 zdummy = DCMPLX(0.0d0,DIMAG(zdummy))
                  if (ABS(DIMAG(zdummy)) < 1.0d-14) 
     &                 zdummy = DCMPLX(DBLE(zdummy),0.0d0)

                  Zwfprod(ikx,iky,iorb1,iorb2,n,is) = zdummy
               end do ; end do

c                  Zwfprod(ikx,iky,iorb1,iorb2,1:nb,1:2)
c     &                 =        Zpsiall(ikx,iky,iorb2,1:nb,1:2)
c     &                 * DCONJG(Zpsiall(ikx,iky,iorb1,1:nb,1:2))
         end do ; end do
      end do ; end do

      return
      end
c############################################## schi0 ####

c#########################################################
c#### schix: chix2
c#########################################################
      subroutine schix2(is,zmtrx2,zmtrx1,nsize)

      use common, only : Ns, U, Nqx
      implicit none
      
      integer, intent(in) :: is, nsize
      complex*16, intent(out):: zmtrx1(nsize,nsize,2)
      complex*16, intent(in):: zmtrx2(nsize,nsize,2)

      integer :: ir12, ir34, ir56
      integer :: isE, isH, is1, is2, is3, is4, is5, is6  
      integer :: mu1, mu2, mu5, mu6, nQ1, nQ2
      integer :: info, ipiv(nsize)
      integer, external :: Ispinpair, Ivrtxinv

      real*8 :: det
      real*8, external :: Wint

      complex*16 :: za(nsize,nsize), zb(nsize,nsize)
      complex*16 :: zdummy
      
c==== .                   2  -->--  3
c==== . chi0(1,2,3,4) =    <       >
c==== .                   1  --<--  4
c==== .   
c==== .   
c==== . chi0(1,2,3,4) W(3,4;5,6) chi0(5,6,7,8)
c==== .
c==== .        2 -->-- 3       6 -->-- 7
c==== .    =  <         >~~~~~<         >
c==== .        1 --<-- 4       5 --<-- 8
c==== . 
c==== .   
c==== .  [1 - chi0(1,2,3,4) W(3,4;5,6)] chi(5,6,7,8) = chi0(1,2,7,8)
c==== .
c==== .   /       2 -->-- 3      6 \    6 -->-- 7       2 -->-- 7   
c==== .  |  1 -  <         >~~~~~   |  </////////>  =  <         >
c==== .   \       1 --<-- 4      5 /    5 --<-- 8       1 --<-- 8   
c==== . 
c==== . 

      isE = Ispinpair(is,'E')
      isH = Ispinpair(is,'H')

      is1 = isH
      is2 = isE

      is3 = isE
      is4 = isH

      is5 = isH
      is6 = isE

      za(:,:) = DCMPLX(0.0d0,0.0d0)
      do ir12 = 1, nsize 
         za(ir12,ir12) = 1.0d0
         do ir56 = 1, nsize
            
c            zdummy = 0.0d0
c            if (ir12 == ir56) then
c               zdummy = 1.0d0
c            end if
            
            do ir34 = 1, nsize
               za(ir12,ir56) = za(ir12,ir56) - zmtrx2(ir12,ir34,isE)
     &              * Wint(ir34,ir56,is3,is4,is5,is6)
            end do
                        
         end do 
      end do
      
      zb(:,:) = zmtrx2(:,:,isE)

caaaaaaaaaaaaaaaaaaaaaaaa
c      call zgetrf(nsize,nsize,za,nsize,ipiv,info)
c      zdummy = 1.0d0
c      do ir12 = 0, nsize-2
c         if (ipiv(ir12) /= ir12 + 1) zdummy = -zdummy
c      end do
c      do ir12 = 1, nsize
c         zdummy = zdummy * za(ir12,ir12)
c      end do
c      zmtrx1(1,1.isE) = zdummy
c      return
caaaaaaaaaaaaaaaaaaaa

      call ZGESV(nsize, nsize, za, nsize, ipiv, zb, nsize, INFO)

      if(INFO /= 0) then
         write(*,*) 'Lapack ZGESV: Info=',Info
         stop
      end if

      zmtrx1(:,:,isE) = zb(:,:)
      

      return
      end
c################################################# schix #

c#########################################################
c#### Wint : interaction
c#########################################################
      real*8 function Wint(ir1,ir2,is1,is2,is3,is4)

      use common, only : U, DJ
      implicit none
      
      integer, intent(in) :: ir1, ir2
      integer, intent(in) :: is1, is2, is3, is4
      integer :: mu1, mu2, mu3, mu4, nQ1, nQ2
      integer :: iorbs, ispns
      integer, external :: Ivrtxinv

      real *8 :: wU, wJ

      
c==== .                
c==== .        1-->--   W(1,2;3,4)   -->--4
c==== .              \              /
c==== .               >~~~~~~~~~~~~<  
c==== .              /              \
c==== .        2--<--                --<--3
c====   
c==== .                
c==== .        1-->--         -->--4         1-->----->--4   
c==== .              \       /                     :          
c==== .   =           >.....<           +          :            
c==== .              /       \                     :          
c==== .        2--<--         --<--3         2--<-----<--3   
c==== .                

c==== .                     |  ispns1    | ispns2     |  ispns3
c==== .              \ spin |  1=2 = 3=4 | 1=2 /= 3=4 | 1=4 /= 2=3
c==== .       orbital \     |            |            |
c==== .  ---------------------------------------------------- 
c==== .  iorbs1: 1=2  = 3=4 |            |    U       |    -U
c==== .  iprbs2: 1=3 /= 2=4 |            |    J       |    -J
c==== .  iorbs3: 1=2 /= 3=4 |   U - 3J   |  U - 2J    |    -J
c==== .  iorbs4: 1=4 /= 2=3 |  -U + 3J   |    J       |  -U + 2J
c==== .                

      Wint = 0.0d0

      nQ1 = Ivrtxinv(ir1,3,'R')
      nQ2 = Ivrtxinv(ir2,3,'L')
      if (nQ1 /= nQ2) return

      mu1 = Ivrtxinv(ir1,1,'R')
      mu2 = Ivrtxinv(ir1,2,'R')
      mu3 = Ivrtxinv(ir2,1,'L')
      mu4 = Ivrtxinv(ir2,2,'L')

      wU = 0.0d0
      wJ = 0.0d0

      wU = U
      wJ = DJ


c==== orbitals[
      if (mu1 == mu2) then
         if (mu3 /= mu4) return 
         if (mu2 == mu3) then 
            iorbs = 1           !1=2 = 3=4
         else
            iorbs = 3           !1=2 /= 3=4
         end if
      else if ((mu1 == mu3).and.(mu2 == mu4)) then
         iorbs = 2              !1=3 /= 2=4
      else if ((mu1 == mu4).and.(mu2 == mu3)) then
         iorbs = 4              !1=4 /= 2=3
      else
         return
      end if
c==== orbitals]

c==== spins[
      if (is1 == is2) then
         if (is3 /= is4) return
         if (is2 == is3) then
            ispns = 1           !1=2 = 3=4
         else
            ispns = 2           !1=2 /= 3=4
         end if
      else if ((is1 == is4).and.(is2 == is3)) then
         ispns = 3              !1=4 /= 2=3
      else
         return
      end if
c==== spins]


      if (iorbs == 1) then
         if (ispns == 1) then
            return
         else if (ispns == 2) then !iorbs1 ispns2
            Wint = wU
         else if (ispns == 3) then  !iorbs1 ispns3
            Wint = -wU
         end if
      else if (iorbs == 2) then
         if (ispns == 1) then 
            return
         else if (ispns == 2) then !iorbs2 ispns2
            Wint = wJ
         else if (ispns == 3) then !iorbs2 ispns3
            Wint = -wJ
         end if
      else if (iorbs == 3) then
         if (ispns == 1) then   !iorbs3 ispns1
            Wint = wU - 3.0d0 * wJ
         else if (ispns == 2) then !iorbs3 ispns2
            Wint = wU - 2.0d0 * wJ
         else if (ispns == 3) then !iorbs3 ispns3
            Wint = -wJ
         end if
      else if (iorbs == 4) then
         if (ispns == 1) then   !iorbs4 ispns1
            Wint = -wU + 3.0d0 * wJ
         else if (ispns == 2) then !iorbs4 ispns2
            Wint = wJ
         else if (ispns == 3) then !iorbs4 ispns3
            Wint = -wU + 2.0d0 * wJ
         end if
      end if


      return
      end
c#########################################################

c#########################################################
c#### Ivrtx : label vertex
c#########################################################
      integer function Ivrtx(mu1,mu2,nQ,LorR)

      use common, only : Ns, Nqx
      implicit none
      
      character, intent(in) :: LorR*1
      integer, intent(in) :: mu1, mu2, nQ

c==== .  LorR = 'L':  the left side vertex of the bubble
c==== .  LorR = 'R':  the right side vertex of the bubble
c==== .
c==== .           (line 1: in, line 2: out)
c==== .
c==== .        L-vertex (1,2,Q)      |   R-vertex (1,2,Q)
c==== .                              | 
c==== .              2 (k)+Q+(q) =>  |  => (k)+(q)+Q 1
c==== .   (k)+Q =>  <                |                > => (k)+Q
c==== .              1  (q) <=       |      <= (q)   2
c==== .                              |
c==== .
c==== . Label vertices so that Ivrtx(i,j,Q,'L') = Ivrtx(j,i,Q,'R')
c==== .                  j->-j
c==== .   i.e.,    Q => <     > => Q    is diagonal
c==== .                  i-<-i   
c==== .  
c==== .

      if (LorR == 'L') then
         Ivrtx = mu1 + (mu2-1) * Ns + nQ * Ns * Ns
      else if (LorR == 'R') then
         Ivrtx = mu2 + (mu1-1) * Ns + nQ * Ns * Ns
      else
         write(*,*) ' -- ERROR in Ivrtx -- '
         stop
      end if

      if ((mu1 > Ns).or.(mu2 > Ns).or.(nQ > Nqx-1)) then
         write(*,*) ' -- ERROR in Ivrtx -- '
         stop
      end if


      return
      end
c#########################################################

c#########################################################
c#### Ivrtxinv : vertex
c#########################################################
      integer function Ivrtxinv(ivrtx,n,LorR)

      use common, only : Ns, Nqx
      implicit none
      
      character, intent(in) :: LorR*1
      integer, intent(in) :: ivrtx, n
      integer :: mu1, mu2, nQ, iv
      integer :: n1, n2, n3, n4

c==== .
c==== .           (line 1: in, line 2: out)
c==== .
c==== .        L-vertex (1,2,Q)      |   R-vertex (1,2,Q)
c==== .                              | 
c==== .              2 (k)+Q+(q) =>  |  => (k)+(q)+Q 1
c==== .   (k)+Q =>  <                |                > => (k)+Q
c==== .              1  (q) <=       |      <= (q)   2
c==== .                              |
c==== .
c==== . Vertices labeled so that Ivrtx(i,j,Q,'L') = Ivrtx(j,i,Q,'R')
c==== .
c==== .  n = 1 :  line coming in = 1st index of Ivrtx(1,2,Q,LorR)
c==== .  n = 2 :  line going out = 2nd index of Ivrtx(1,2,Q,LorR)
c==== .  n = 3 :  For L, Q = k2-k1 (k1 = q,  k2 = k+Q+q)
c==== .           For R, Q = k1-k2 (k1 = k+Q+q, k2 = q)

      iv = ivrtx - 1

      if (LorR == 'L') then
         n1 = 1 ; n2 = 2 ; n3 = 3
      else if (LorR == 'R') then
         n1 = 2 ; n2 = 1 ; n3 = 3
      else
         write(*,*) ' -- ERROR in Ivrtx -- '
         stop
      end if
      
      Ivrtxinv = iv / Ns / Ns   !Q
      if ((Nqx == 1).and.(Ivrtxinv /= 0)) then
         write(*,*) ' -- Error in Ivrtxinv --'
         write(*,*) ' Q need to be 0 for Nqx=1 '
      end if
      if (n == n3) return
      
      nQ = Ivrtxinv
      Ivrtxinv = (iv - nQ * Ns**2) / Ns + 1 !line out
      if (n == n2) return
      
      mu2 = Ivrtxinv
      Ivrtxinv =  iv - (mu2-1) * Ns - nQ * Ns * Ns + 1 !line in
      if (n == n1) return


      write(*,*) ' -- ERROR in Ivrtx -- '
      stop

      return
      end
c#########################################################


c#########################################################
c#### Ispin : 
c#########################################################
      integer function Ispinpair(is,EorH)

      use common, only : Ns
      implicit none
      
      character, intent(in) :: EorH*1
      integer, intent(in) :: is

c==== .        isE -->-- isE   
c==== .       <             >
c==== .        isH --<-- isH  
c==== .
c==== . A pair of spin (isE, isH) for an electron-hole pair.
c==== . isE = Ispinpair(is,'E'),  isH = Ispinpair(is,'H')
c==== .      
c==== .         |  isE = 1         isE = 2
c==== . --------+-------------------------- 
c==== . isH = 1 |   is = 1          is = 4
c==== .         |                    
c==== . isH = 2 |   is = 3          is = 2
c==== .

      Ispinpair = 0
      if (EorH == 'E') then
         if(is == 1)then
            Ispinpair = 1
         else if(is == 2)then
            Ispinpair = 2
         else if(is == 3)then
            ispinpair = 1
         else if(is == 4)then
            ispinpair = 2
         else
            write(*,*)"error in Ispinpair: 'is' should be 1 2 3 or 4"
            write(*,*)"                  : is = ",is
            stop
         end if
      else if (EorH == 'H') then
         if(is == 1)then
            Ispinpair = 1 
         else if(is == 2)then
            Ispinpair = 2 
         else if(is == 3)then
            Ispinpair = 2 
         else if(is == 4)then
            Ispinpair = 1 
         else
            write(*,*)"error in Ispinpair: 'is' should be 1 2 3 or 4"
            write(*,*)"                  : is = ",is
            stop
         end if
      else
         write(*,*)"error: 'n' should be 1 or 2:"
         stop
      end if

      
      return
      end
c#########################################################
      

c#########################################################
c#### Iorbcomb : 
c#########################################################
      integer function Iorbcomb(mu1,mu2,mu3,mu4)

      use common, only : Ns
      implicit none
      
      integer, intent(in) :: mu1, mu2, mu3, mu4
      integer :: iorb(Ns), num(Ns), i, imax, iflag

c==== .  (a, b, c, d)  => 0
c==== .  (a,a,a,a) => 1
c==== .  (a,a,b,b) => 2
c==== .  (a,b,b,a) => 3
c==== .  (a,b,a,b) => 4


      iorbcomb = 0

      if ((mu1 == mu2).and.(mu3 == mu4)) then
         iorbcomb = 2           !1=2 /= 3=4
         if (mu1 == mu3) iorbcomb = 1 !1=2=3=4
         return
      else if((mu1 == mu4).and.(mu2 == mu3)) then
         iorbcomb = 3           !1=4 /= 2=3
         return
      else if((mu1 == mu3).and.(mu2 == mu4)) then
         iorbcomb = 4           ! 1=3 /= 2=4
         return
      end if


      return
      end
c#########################################################
      
c#########################################################
c#### Iorbcomb : 
c#########################################################
      integer function Iorbselect(iorbs)

      implicit none
      
      integer, intent(in) :: iorbs

c==== .  (a, b, c, d)  => 0
c==== .  (a,a,a,a) => 1
c==== .  (a,a,b,b) => 2
c==== .  (a,b,b,a) => 3
c==== .  (a,b,a,b) => 4
      
      iorbselect = 0

      if (iorbs == 0) return 
      if (iorbs == 4) return
      if (iorbs == 2) return 
         
      iorbselect = 1
      if (iorbs == 1) return 
      if (iorbs == 3) return
      
      write(*,*) 'Error in iorbselect'
      write(*,*) 'check input'
      stop


      return
      end
c#########################################################
      
