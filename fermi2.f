c#########################################################
c#### Fermi seek: find fermi surface    EK 2015 
c#########################################################

      subroutine fermi_seek_grid(filename, ene, np)

      use common, only : Nkx, Nky, Ns, Ns2, Pi, Nredx, Nredy, Nqx,
     &     Eall, Zpsiall, Dmu, Dne, Dtemp
      implicit none

      character, intent(in) :: filename*15
      integer, intent(out) :: np
      real*8, intent(in) :: ene

      character :: chara*80, ch
      integer :: i, j, ikx, iky, kx, ky, ispin, iflag
      integer :: kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4
      integer :: izero
      real*8 :: dummy1, dummy2, dummy3, dummy4
      real*8 :: ddd, dkx, dky, akx, aky
      real*8 :: akfx(1:4), akfy(1:4)

c==== .    ゾーンを格子に切ってゼロエネルギーの点を探す．             
c==== .    ひとつの格子について，下のように番号付け，
c==== .    バンドがどの辺を横切るか調べる．
c==== .
c==== .             3_________4                           
c==== .             |         |          
c==== .             |         |                               
c==== .             |         |          
c==== .            1|_________|2                         
c==== . 

      dkx = 1.0d0 / DBLE(Nkx)
      dky = 1.0d0 / DBLE(Nky)

      open(10,file = filename)
      write(10,'(a)') '#fermi_surface'
      write(10,'(a)') '#data: ikx,iky, iband,ispin, '//
     &     'kFx1,kFy1,kFx2,kFy2, izero, E1,E2,E3,E4'

      write(10,'(a)') '#Beginning'

      np = 0
      chara = '    '
      do ispin = 1, 2 ; do i = 1, Ns*Nqx
         iflag = 0
         do ikx = -Nkx/2, Nkx/2-1 ; do iky = -Nky/2, Nky/2-1
            akx = DBLE(ikx) / DBLE(Nkx)
            aky = DBLE(iky) / DBLE(Nky)

            kx1 = MOD(ikx+Nkx,Nkx)
            ky1 = MOD(iky+Nky,Nky)
            dummy1 =  Eall(kx1,ky1,i,ispin) - Dmu - ene

            kx2 = MOD(kx1+1,Nkx)
            ky2 = ky1
            dummy2 =  Eall(kx2,ky2,i,ispin) - Dmu - ene

            kx3 = kx1
            ky3 = MOD(ky1+1,Nky)
            dummy3 =  Eall(kx3,ky3,i,ispin) - Dmu - ene

            kx4 = MOD(kx1+1,Nkx)
            ky4 = MOD(ky1+1,Nky)
            dummy4 =  Eall(kx4,ky4,i,ispin) - Dmu - ene

            izero = 0
            if (dummy1 * dummy2 < 0.0d0) then
               izero = izero + 1
               ddd = ABS(dummy1 / (dummy1 - dummy2)) * dkx
               akfx(izero) = akx + ddd
               akfy(izero) = aky
            end if
            if (dummy1 * dummy3 < 0.0d0) then
               izero = izero + 1
               ddd = ABS(dummy1 / (dummy1 - dummy3)) * dky
               akfx(izero) = akx
               akfy(izero) = aky + ddd
            end if
            if (dummy2 * dummy4 < 0.0d0) then
               izero = izero + 1
               ddd = ABS(dummy2 / (dummy2 - dummy4)) * dky
               akfx(izero) = akx + dkx
               akfy(izero) = aky + ddd
            end if
            if (dummy3 * dummy4 < 0.0d0) then
               izero = izero + 1
               ddd = ABS(dummy3 / (dummy3 - dummy4)) * dkx
               akfx(izero) = akx + ddd
               akfy(izero) = aky + dky
            end if
            
            if (izero == 0) cycle
            if (izero == 1) write(*,*) 'WARNING: Grid size may be large'
            if (izero > 2) write(*,*) 'WARNING: Grid size may be large'
            if (izero /= 2) ch ='x'
c            write(10,'(4i5,4g20.10,i5)') 
c     &           ikx, iky, i, ispin,
c     &           akfx(1), akfy(1), akfx(2), akfy(2), izero
            write(10,'(4i5,4g20.10,i5,4g20.10)') 
     &           ikx, iky, i, ispin,
     &           akfx(1), akfy(1), akfx(2), akfy(2),izero,
     &           dummy1,dummy2,dummy3,dummy4
            np = np + 1
            iflag = 1
         end do ; end do
         if(iflag == 1) write(10,'(a)')' '
      end do ; end do

      write(10,'(a)') '#EOF'

      close(10)


      return
 1    format(i0)
      end
c######################################################################

c###################################################################
c##### Fermi_velocity:
c###################################################################

      subroutine fermi_velocity2()

      use common, only : Ns, Nqx, Pi, Zi
      implicit none

      character(len=32) :: filein, fileout
      real*8 :: akx, aky, akx1, aky1, akx2, aky2, dkx, dky, dl
      real*8 :: vx(Ns*Nqx), vy(Ns*Nqx)

      integer :: i, j, ispin
      character :: chara*74
      
      integer :: idummy(1:10), iband

      write(*,*) 'Enter input filename (fermi surface data)'
      read(5,*) filein

      if (filein(1:1) == '!') filein = 'out.fs.dat'
      fileout = 'v-'//filein
      open(10,file = filein)
      open(20,file = fileout)
      chara = '#'
      write(20,'(a)') '#fermi_velocity' 
      write(20,'(a)') '#data: akx, aky, Vfx, Vfy, iband, ispin' 
      write(20,'(a)') '#Beginning'
 
      do while(chara(1:10) /= '#Beginning')
         read(10,*) chara
      end do


      do while (chara(1:4) /= '#EOF')
         
         read(10,*) idummy(1:2), iband, ispin,
     &        akx1, aky1, akx2, aky2, chara

         akx = (akx1 + akx2) * 0.5d0
         aky = (aky1 + aky2) * 0.5d0
         dl = ABS( (akx1 + aky1 * Zi) - (akx2 + aky2 * Zi) )

         dkx = akx * 2.0d0 * Pi
         dky = aky * 2.0d0 * Pi

         call Svelo(dkx, dky, vx, vy, ispin)
         
         write(20,'(5f12.8, 2i5)') akx, aky, vx(iband), vy(iband), 
     &        dl, iband, ispin
         
         read(10,*,err=999) chara
         backspace(10)
      end do
      write(20,'(a)') chara


      close(10)
      close(20)

 999  return
      end

c###################################################################

c#########################################################
c#### fermi_KfandVf: find Kf and calculate Vf  EK 2015 
c#########################################################

      subroutine fermi_KfandVf(ene, np, mkx, mky, mode)

      use common, only : Ns, Ns2, Pi, Nredx, Nredy, Nqx,
     &     Eall, Zpsiall, Dmu, Dne, Dtemp, Zi
      implicit none

      character, intent(in) :: mode
      integer, intent(in) :: mkx, mky
      integer, intent(out) :: np
      real*8, intent(in) :: ene

      character :: chara*80, ch, modevel*1
      character(len=32) :: filekfout, filevfout
      integer :: i, j, ikx, iky, kx, ky, ispin, iflag
      integer :: kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4
      integer :: izero
      real*8 :: dummy1, dummy2, dummy3, dummy4
      real*8 :: ddd, dkx, dky, akx, aky, dkx0, dky0, dl
      real*8 :: dkxm, dkym, delkx, delky
      real*8 :: dkfx(1:4), dkfy(1:4)
      real*8 :: edata(4,Ns*Nqx), eig(Ns*Nqx)
      real*8 :: vx(Ns*Nqx), vy(Ns*Nqx)

      complex*16 :: zpsi(Ns*Nqx,Ns*Nqx)

c==== .    ゾーンを格子に切ってゼロエネルギーの点を探す．             
c==== .    ひとつの格子について，下のように番号付け，
c==== .    バンドがどの辺を横切るか調べる．
c==== .
c==== .             3_________4                           
c==== .             |         |          
c==== .             |         |                               
c==== .             |         |          
c==== .            1|_________|2                         
c==== . 
c==== .    modevel = 'v' のときは，フェルミ速度も計算する．             
c==== .    辺上の 2 交点の中点での速度を計算する．      
c==== . 

      modevel = mode(1:1)
      if (mode(1:1) == 'V') modevel = 'v'

      write(*,*) 'Enter output filename (fermi surface data)'
      read(5,*) filekfout

      if (filekfout(1:1) == '!') filekfout = 'out.fs.dat'
      write(*,*) 'output filename set to ',filekfout


      delkx = 1.0d0 / DBLE(mkx)
      delky = 1.0d0 / DBLE(mky)

      open(10,file = filekfout)
      write(10,'(a)') '#fermi_surface'
      write(10,'(a,g20.10)') '#Energy=', ene
      write(10,'(a,2i7)') '#size=', mkx, mky
      write(10,'(a)') '#data: ikx,iky, iband,ispin, '//
     &     'kFx1,kFy1,kFx2,kFy2, izero, E1,E2,E3,E4'

      write(10,'(a)') '#Beginning'

      if (modevel(1:1) == 'v') then
         filevfout = 'v-'//filekfout
         open(20,file = filevfout)
         chara = '#'
         write(20,'(a)') '#fermi_velocity' 
         write(20,'(a,g20.10)') '#Energy=', ene
         write(20,'(a,2i7)') '#size=', mkx, mky
         write(20,'(a)') '#data: dkx, dky, Vfx, Vfy, dl, iband, ispin' 
         write(20,'(a)') '#Beginning'
      end if

      np = 0
      chara = '    '
      do ispin = 1, 2
         if(ispin == 2) cycle
         do ikx = -mkx/2, mkx/2-1 ; do iky = -mky/2, mky/2-1
            dkx0 = DBLE(ikx) / DBLE(mkx)
            dky0 = DBLE(iky) / DBLE(mky)

            kx1 = MOD(ikx+mkx,mkx)
            ky1 = MOD(iky+mky,mky)
            akx = DBLE(kx1) / DBLE(mkx) * Pi * 2.0d0
            aky = DBLE(ky1) / DBLE(mky) * Pi * 2.0d0
            call  diag(edata(1,:),zpsi(:,:),akx,aky,ispin)

            kx2 = MOD(kx1+1,mkx)
            ky2 = ky1
            akx = DBLE(kx2) / DBLE(mkx) * Pi * 2.0d0
            aky = DBLE(ky2) / DBLE(mky) * Pi * 2.0d0
            call  diag(edata(2,:),zpsi(:,:),akx,aky,ispin)

            kx3 = kx1
            ky3 = MOD(ky1+1,mky)
            akx = DBLE(kx3) / DBLE(mkx) * Pi * 2.0d0
            aky = DBLE(ky3) / DBLE(mky) * Pi * 2.0d0
            call  diag(edata(3,:),zpsi(:,:),akx,aky,ispin)

            kx4 = MOD(kx1+1,mkx)
            ky4 = MOD(ky1+1,mky)
            akx = DBLE(kx4) / DBLE(mkx) * Pi * 2.0d0
            aky = DBLE(ky4) / DBLE(mky) * Pi * 2.0d0
            call  diag(edata(4,:),zpsi(:,:),akx,aky,ispin)

            do i = 1, Ns*Nqx

               dummy1 =  edata(1,i) - Dmu - ene
               dummy2 =  edata(2,i) - Dmu - ene
               dummy3 =  edata(3,i) - Dmu - ene
               dummy4 =  edata(4,i) - Dmu - ene
               
               izero = 0
               if (dummy1 * dummy2 < 0.0d0) then
                  izero = izero + 1
                  ddd = ABS(dummy1 / (dummy1 - dummy2)) * delkx
                  dkfx(izero) = dkx0 + ddd
                  dkfy(izero) = dky0
               end if
               if (dummy1 * dummy3 < 0.0d0) then
                  izero = izero + 1
                  ddd = ABS(dummy1 / (dummy1 - dummy3)) * delky
                  dkfx(izero) = dkx0
                  dkfy(izero) = dky0 + ddd
               end if
               if (dummy2 * dummy4 < 0.0d0) then
                  izero = izero + 1
                  ddd = ABS(dummy2 / (dummy2 - dummy4)) * delky
                  dkfx(izero) = dkx0 + delkx
                  dkfy(izero) = dky0 + ddd
               end if
               if (dummy3 * dummy4 < 0.0d0) then
                  izero = izero + 1
                  ddd = ABS(dummy3 / (dummy3 - dummy4)) * delkx
                  dkfx(izero) = dkx0 + ddd
                  dkfy(izero) = dky0 + delky
               end if
               
               if (izero == 0) cycle
               if (izero == 1) write(*,*)
     &              'WARNING: Grid size may be large'
               if (izero > 2) write(*,*) 
     &              'WARNING: Grid size may be large'
               if (izero /= 2) cycle

               write(10,'(4i5,4g20.10,i5,4g20.10)') 
     &              ikx, iky, i, ispin,
     &              dkfx(1), dkfy(1), dkfx(2), dkfy(2),izero,
     &              dummy1,dummy2,dummy3,dummy4
               np = np + 1

c#### velocity
               if (modevel(1:1) /= 'v')  cycle
               dkxm = (dkfx(1) + dkfx(2)) * 0.5d0
               dkym = (dkfy(1) + dkfy(2)) * 0.5d0
               dl = ABS( (dkfx(1) + dkfy(1) * Zi)
     &              -    (dkfx(2) + dkfy(2) * Zi) )
               
               akx = dkxm * 2.0d0 * Pi
               aky = dkym * 2.0d0 * Pi
               
               call diag(eig(:),zpsi(:,:),akx,aky,ispin)
               call velocity(vx,zpsi(:,:),akx,aky,'x')
               call velocity(vy,zpsi(:,:),akx,aky,'y')
      
               write(20,'(5f12.8, 2i5)') dkxm, dkym, 
     &              vx(i), vy(i), dl, i, ispin
            end do

         end do ; end do
      end do 

      write(10,'(a)') '#EOF'

      close(10)

      if (modevel(1:1) == 'v') then
         write(20,'(a)') '#EOF'
         close(20)
      end if

      return
 1    format(i0)
      end
c######################################################################



c%%%%
c%%%% END_OF_FILE:fermi2.f
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
