c###################################################################
c##### mfcalc: meanfield calculation
c###################################################################
      subroutine mfcalc()
      
      use common
      implicit none

      character :: chara*40

      character :: cdwmode*3, n0mode*3

      integer :: i, j, ix, iy, it, is, js, js0, js2
      integer :: ispin, ii
      integer :: kx, ky, kkx, kky, mu, nu
      integer :: idummy, nzone, nsQ
      
      real*8 :: volk, dummy
      real*8 :: ratio
      real*8 :: dnewnuu(Ns),dnewndd(Ns)
      real*8 :: ddtemp, ddmu
      real*8 :: dummy1, dummy2, dnop0
      complex*16 :: zdummy
      complex*16 :: zop2temp(Ns,Ns,2), zop0temp(Ns,Ns,2)

      real*8, external :: Dffn, Dntot, Dneq

      cdwmode = 'off'
      n0mode = 'on'

c==== cdw, n0 モード off の場合,hamiltonian 行列にこれらの項を含めない．
      if (cdwmode(1:2) /= 'on') then
          Zop2(:,:,:) = 0.0d0
      end if
      if (n0mode(1:2) /= 'on') then
          Zop0(:,:,:) = 0.0d0
      end if

c      call set_hopping()

      dnewnuu(:) = 0.0d0
      dnewndd(:) = 0.0d0

      write(*,*)'MF calculation'
      write(*,*)'Maxit = ',Maxit
      write(*,*)'Conv = ',Conv

      nsQ = Ns * Nqx
      nzone = Nkx / Nqx
c      nzone = Nkx

      write(*,'(g30.16)')Dmu
      zpsiall(:,:,:,:,:)=0.0d0
      chara = 'init'
c      call Hcheck(0.25d0,0.250d0,1,chara)!aaaaaaaaaaaaa
c      call Hcheckall(Nkx,Nky,1,chara)!aaaaaaaaaaaaa
c      call seigen2(nzone, Nky,1)!aaaaaaaaaaaaa
c      call eigenfncheck(1,1,1,chara)!aaaaaaaaaaaaa

c## ITERATION[
      it = 1 ; conrs = 1000.0d0
      write(*,*)
      do while ((it < Maxit).and.(conrs >= Conv))

         call seigen2(nzone, Nky,1)
         call seigen2(nzone, Nky,2)
c## CHEMICAL POTENTIAL[
         call smuval(Dmu)
         if (it == 1) then
            ddtemp = 0.00001d0
            ddmu = 1.0d-8
c            call Smuval2(ddmu,ddtemp)
         end if

c## CHEMICAL POTENTIAL[


c## NEW ORDER PARAMETER[
         Dens(:,:) = 0.0d0
         dummy = 0.0d0
         if (Nqx == 1) then
            Zopnew(:,:,:) = 0.0d0
            Zopnew0(:,:,:) = 0.0d0
            Zopnew2(:,:,:) = 0.0d0
         else
            do mu = 1, Ns ; do nu = 1, Ns
               zdummy = 0.0d0
               do kx = 0, Nkx - 1 ; do ky = 0, Nky - 1
                  is = nu
                  js0 = mu
                  js = mu + Ns
                  js2 = mu + Ns * 2
                  if (Nqx == 2) js2 = js0
                  do i = 1, nsQ

                     Zopnew(mu,nu,:) = Zopnew(mu,nu,:)
     &                    + DCONJG(Zpsiall(kx,ky, js, i,:))
     &                    *        Zpsiall(kx,ky, is, i,:)
     &                    * Dffn(Eall(kx,ky,i,:)-Dmu,Dtemp)
                     
                     Zopnew0(mu,nu,:) = Zopnew0(mu,nu,:)
     &                    + DCONJG(Zpsiall(kx,ky, js0, i,:))
     &                    *        Zpsiall(kx,ky,  is, i,:)
     &                    * Dffn(Eall(kx,ky,i,:)-Dmu,Dtemp)
                     
                     Zopnew2(mu,nu,:) = Zopnew2(mu,nu,:)
     &                    + DCONJG(Zpsiall(kx,ky, js2, i,:))
     &                    *        Zpsiall(kx,ky,  is, i,:)
     &                    * Dffn(Eall(kx,ky,i,:)-Dmu,Dtemp)
                     if (mu == nu) then
                        Dens(mu,:) = Dens(mu,:)
     &                       + DCONJG(Zpsiall(kx,ky,is,i,:))
     &                       *       Zpsiall(kx,ky,is,i,:)
     &                       * Dffn(Eall(kx,ky,i,:)-Dmu,Dtemp)
                     end if
                  end do
               end do ; end do
            end do ; end do
            volk = DBLE(Nkx * Nky)
            Zopnew(:,:,:) = Zopnew(:,:,:) / volk
            Zopnew0(:,:,:) = Zopnew0(:,:,:) / volk
            Zopnew2(:,:,:) = Zopnew2(:,:,:) / volk
         end if
         dnop0 = 0.0d0
         do mu = 1, Ns
            Dnuu(mu) = Zopnew(mu,mu,1)
c            Dndd(mu) = -Dnuu(mu)
            Dndd(mu) = Zopnew(mu,mu,2)
            dnop0 = dnop0 + Zopnew0(mu,mu,1)+Zopnew0(mu,mu,2)
         end do
         Dens(:,:) = Dens(:,:) / volk
         Dne1 = Dntot(Dmu) * DBLE(Ns)
c         Dne1 = SUM(Dens)
         write(*,*) 'density =', dnop0
c## NEW ORDER PARAMETER]


c## SELF-CONSISTENCY[
         conrs = -1.0d0

         conrs = MAXVAL(ABS(Zopnew(:,:,:)-Zop(:,:,:)))
c         write(*,*) MAXLOC(ABS(Zopnew(:,:,:)-Zop(:,:,:)))
c         write(*,*) Zopnew(1,5,1),'==',Zop(1,5,1)

         if (cdwmode(1:2) == 'on') then
            dummy = MAXVAL(ABS(Zopnew2(:,:,:)-Zop2(:,:,:)))
            conrs = MAX(conrs,dummy)
c            write(*,*) MAXLOC(ABS(Zopnew2(:,:,:)-Zop2(:,:,:))),dummy
c            write(*,*) Zopnew2(5,5,1),Zop2(5,5,1)
c            write(*,*) Zopnew2(5,5,2),Zop2(5,5,2)
         end if

         if (n0mode(1:2) == 'on') then
            dummy = MAXVAL(ABS(Zopnew0(:,:,:)-Zop0(:,:,:)))
            conrs = MAX(conrs,dummy)
c            write(*,*) MAXLOC(ABS(Zopnew0(:,:,:)-Zop0(:,:,:))),dummy
         end if

c         dummy = 0.0d0
c         do i = 1, Ns
c            dummy = dummy + Zopnew0(i,i,1)
c         end do
c         write(*,*) dummy

         write(*,'(a,i6,g15.8,a,g15.8,a,g15.8,a)')
     &        ' it,conrs=  ', it,conrs, '  Dmu= ', Dmu,
     &        ' M =  ', SUM(Dnuu(:)-Dndd(:)),'muB'
c     &        ' Sz =  ', ABS(Zop(2,2,1)-Zop(3,3,1))

         ratio = 0.3d0
         Zop(:,:,:) = Zopnew(:,:,:) * (1.0d0 - ratio)
     &        + ratio * Zop(:,:,:)

         ratio = 0.0d0
         Zop0(:,:,:) = Zopnew0(:,:,:) * (1.0d0 - ratio)
     &        + ratio * Zop0(:,:,:)

         Zop2(:,:,:) = Zopnew2(:,:,:) * (1.0d0 - ratio)
     &        + ratio * Zop2(:,:,:)

         Zopnew(:,:,:) = 0.0d0
         Zopnew0(:,:,:) = 0.0d0
         Zopnew2(:,:,:) = 0.0d0
         do i = 1, Ns
            Dnuu(i) = DBLE(Zopnew(i,i,1))+DIMAG(Zopnew(i,i,1))
            Dndd(i) = DBLE(Zopnew(i,i,2))+DIMAG(Zopnew(i,i,2))
         end do

c         Zop0(:,:,:) = 0.0d0
         zop2temp(:,:,:) = Zop2(:,:,:)
         if (cdwmode(1:2) /= 'on') then
             Zop2(:,:,:) = 0.0d0
         end if
         zop0temp(:,:,:) = Zop0(:,:,:)
         if (n0mode(1:2) /= 'on') then
             Zop0(:,:,:) = 0.0d0
         end if

c         chara = 'temp'
c         call Hcheckall(Nkx,Nky,1,chara) !aaaaaaaaaaaaa
c         stop!aaaaaaaaaaa

c         do i = 1, Ns; do j = 1,Ns
c            if (i /= j) then
c               zop2temp(i,j,:) = 0.0d0
c               zop0temp(i,j,:) = 0.0d0
c               Zop2(i,j,:) = 0.0d0
c               Zop0(i,j,:) = 0.0d0
c            end if
c         end do ; end do

c         do i = 1, Ns
c            Zop2(i,i,:) = 0.0d0
c         end do
c         do i = 1, Ns
c            Zop(i,i,:) = CMPLX(DBLE(Zop(i,i,:)),DBLE(Zop(i,i,:)))
c         end do

c            Zop(:,:,:) = CMPLX(DBLE(Zop(:,:,:)),DBLE(Zop(:,:,:)))

c         do i = 1, Ns
c            Zop(i,i,:) = CMPLX(DBLE(Zop(i,i,:)),0.0d0)
c         end do

c         do i = 1, Ns-1 ; do j = i+1, Ns
c            if ((i == 2).and.(j == 3))  cycle
c            if ((i == 3).and.(j == 2))  cycle
c            Zop(i,j,:) = 0.0d0
c            Zop(j,i,:) = 0.0d0
c         end do ; end do
c            Zop(2,3,:) = 0.0d0
c            Zop(3,2,:) = 0.0d0
            
c         Zop(3,3,:) = Zop(2,2,:)
         
c         Zop(1,2,:) = 0.0d0
c         Zop(1,3,:) = 0.0d0
c         Zop(1,4,:) = 0.0d0
c         Zop(2,1,:) = 0.0d0
c         Zop(2,4,:) = 0.0d0
c         Zop(2,5,:) = 0.0d0
c         Zop(3,1,:) = 0.0d0
c         Zop(3,4,:) = 0.0d0
c         Zop(3,5,:) = 0.0d0
c         Zop(4,1,:) = 0.0d0
c         Zop(4,2,:) = 0.0d0
c         Zop(4,3,:) = 0.0d0
c         Zop(4,5,:) = 0.0d0
c         Zop(5,2,:) = 0.0d0
c         Zop(5,3,:) = 0.0d0
c         Zop(5,4,:) = 0.0d0

c         Zop(3,2,:) = 0.0d0
c         Zop(2,3,:) = 0.0d0


c## SELF-CONSISTENCY]

c## output temporary result[
         if(MOD(it,10) == 0)then
            Zop2(:,:,:) = zop2temp(:,:,:)
            Zop0(:,:,:) = zop0temp(:,:,:)
            chara = 'temp.dat'
            call savedata(chara)
            if (cdwmode(1:2) /= 'on') then
               if (conrs >= Conv) Zop2(:,:,:) = 0.0d0
            end if
            if (n0mode(1:2) /= 'on') then
               if (conrs >= Conv) Zop0(:,:,:) = 0.0d0
            end if

            if (ABS(Dne1-Dne) >= 1.0d-4)then
               write(*,*) '### calculation may be wrong ##'
               write(*,*) 'Dne=', Dne
               write(*,*) 'Dne1=', Dne1
               write(*,*) '|Dne1-Dne|=',ABS(Dne1-Dne)
            end if
         end if
c## output temporary result]
         it = it +1
      end do
c## ITERATION]

c==== cdw, n0 モード off の場合,hamiltonian 行列にこれらの項を含めない．
      if (cdwmode(1:2) /= 'on') then
          Zop2(:,:,:) = 0.0d0
      end if
      if (n0mode(1:2) /= 'on') then
          Zop0(:,:,:) = 0.0d0
      end if

c      call energy_check()

      call density()
      write(*,*)  'energy =', Fen

      if (conrs >= Conv) then
         write(*,*)
         write(*,*)'#### NOT CONVERGED YET #####'
         write(*,*)
         return
      end if


c      chara = 'final'
c      call Hcheck(0.25d0,0.250d0,1,chara)
c      call eigenfncheck(1,1,1,chara)

      return

      end
c###################################################################

c###################################################################
c##### boundary:                                    
c###################################################################

      subroutine boundary2(mkx,mky,nn,eig,zpsi)

      implicit none
      
      integer, intent(in) :: mkx, mky, nn
      real*8 :: eig(0:mkx,0:mky,nn)
      complex*16 :: zpsi(0:mkx,0:mky,nn,nn)

      eig(  mkx, 0:mky, :) = eig(0    , 0:mky, :)
      eig(0:mkx,   mky, :) = eig(0:mkx, 0    , :)
      eig(mkx,mky,:) = eig(0,0,:)
      
      zpsi(  mkx, 0:mky,:,:) = zpsi(0    , 0:mky,:,:)
      zpsi(0:mkx,   mky,:,:) = zpsi(0:mkx, 0    ,:,:)
      zpsi(  mkx,   mky,:,:) = zpsi(0    , 0    ,:,:)

      return
      end 
      
c###################################################################


c###################################################################
c##### boundary:                                    
c###################################################################

      subroutine boundary()


      use common, only : Nkx, Nky, Eall, Zpsiall
      implicit none

      Eall(  Nkx, 0:Nky, :,:) = Eall(0    , 0:Nky, :,:)
      Eall(0:Nkx,   Nky, :,:) = Eall(0:Nkx, 0    , :,:)
      Eall(Nkx,Nky,:,:) = Eall(0,0,:,:)
      
      Zpsiall(  Nkx, 0:Nky,:,:,:) = Zpsiall(0    , 0:Nky,:,:,:)
      Zpsiall(0:Nkx,   Nky,:,:,:) = Zpsiall(0:Nkx, 0    ,:,:,:)
      Zpsiall(  Nkx,   Nky,:,:,:) = Zpsiall(0    , 0    ,:,:,:)

      return
      end 
      
c###################################################################

c###################################################################
c##### density:                                    
c###################################################################

      subroutine density()


      use common, only : Nkx, Nky, Ns, Ns2, Nqx, 
     &     Eall, Zpsiall, Zdens, Dmu, Dtemp, Dens, Zop
      implicit none

      integer :: i, j, kx, ky, mu, nu, ispin
      real*8 :: dummy, dn, volk
      complex*16 :: zdummy
      real*8, external :: Dffn



      volk = DBLE(Nkx*Nky) * 2 !2 for spin
      Zdens(:,:) = 0.0d0
      do mu = 1, Ns ; do nu = 1, Ns
         do kx = 0, Nkx - 1 ; do ky = 0, Nky - 1
            do i = 1, Ns*Nqx ; do ispin = 1, 2
               Zdens(mu,nu) = Zdens(mu,nu)
     &              + DCONJG(Zpsiall(kx,ky, mu, i,ispin))
     &              *        Zpsiall(kx,ky, nu, i,ispin)
     &              * Dffn(Eall(kx,ky,i,ispin)-Dmu,Dtemp)
            end do ; end do
         end do ; end do
      end do ; end do
      Zdens(1:Ns,1:Ns) = Zdens(1:Ns,1:Ns) / volk
      do mu = 1, Ns 
         Dens(mu,1) = DBLE(Zdens(mu,mu) + Zop(mu,mu,1))
         Dens(mu,2) = DBLE(Zdens(mu,mu) - Zop(mu,mu,1))
      end do

      dn = 0.0d0
      do mu = 1, Ns 
         dummy = DIMAG(Zdens(mu,mu))
         if (dummy > 1.0d-13) then
            write(*,*) 'ERROR in density'
            write(*,*) mu, dummy
            read(5,*)
         end if
         dummy = Dens(mu,1) + Dens(mu,2) ! n(mu,up) + n(mu,dn)
         dn = dn + dummy
      end do
      write(*,*) 'density =', dn


      return
      end 
      
c###################################################################

c###################################################################
c##### Senergy: calculate energy                                   
c###################################################################
      subroutine Senergy()

      use common
      implicit none

      integer :: j, mu, nu, ispin
      integer :: k0x, k0y
      real*8 :: vol, volk0
      real*8 :: df, eav, eav0, entt, elat
      real*8 :: dummy, cu, cj
      complex*16 :: zdummy

      real*8, external :: Dffn, Dentrp


      volk0 = DBLE(Nkx*Nky*2) !2 for spin

c==== H_0 
      eav = 0.d0
      entt = 0.d0

      do k0x = 0, Nkx/Nqx - 1
         do k0y = 0, Nky - 1
            do j = 1, Ns*Nqx ; do ispin = 1,2
               df = Dffn(Eall(k0x,k0y,j,ispin)-Dmu,Dtemp)
               eav = eav + df * Eall(k0x,k0y,j,ispin)
               entt = entt + Dentrp(Df)
            end do ; end do
         end do
      end do

      eav = eav * 2.0d0 ! 2.0 for spin
      eav = eav / volk0

      entt = entt * 2.0d0 ! 2.0 for spin
      entt = entt / volk0 * Dtemp ! entt = S*T

      write(*,*) ' <H_0>=', eav

c==== constant from H_U
      cu = SUM( ABS(Zop(:,:,1))**2 ) - SUM( ABS(Zdens(:,:))**2 )
      cu = cu + Dne1**2 * 0.5d0
      cu = cu * U
      write(*,*) '   C_U=', cu

c==== constant from H_J
      zdummy = 0.0d0
      do mu = 1, Ns
         zdummy = zdummy + Zop(mu,mu,1)
      end do
      cj = ABS(zdummy)**2

      zdummy = 0.0d0
      do mu = 1, Ns ; do nu = 1, Ns
         zdummy = zdummy + Zop(mu,nu,1) * DCONJG(Zop(nu,mu,1))
      end do ; end do
      cj = cj + ABS(zdummy)**2

      cj = cj - 2.0d0 * SUM( ABS(Zop(:,:,1))**2 )
      cj = cj + 4.0d0 * SUM( ABS(Zdens(:,:))**2 )
      cj = cj + SUM( Zdens(:,:)**2 )
      cj = cj - Dne1**2 * 5.0d0 * 0.25d0
      cj = cj * DJ
      write(*,*) '   C_J=', cj

c==== total energy
      eav = eav + cu + cj

      Fen = eav - entt          ! F = E - TS
      write(*,*) '   Fen=', Fen

      return
      end
c####################################################################

      double precision function Dentrp(df)

      use common
      implicit none

      real*8 :: daa1, daa2, dlimit, dg, df


      dlimit = 1.0d-80

      if(DABS(Df) < dlimit) then 
         daa1=0.d0
       else
         daa1= - df * DLOG(df)
      end if

      dg = 1.d0 - df
      if(DABS(dg) < dlimit) then
         Daa2 = 0.d0
       else
         Daa2 = - dg * DLOG(dg)
      endif

      Dentrp = daa1 + daa2
      return
      end
c####################################################################

c###################################################################
c##### orbital:
c###################################################################
      subroutine orbital()
      
      use common
      implicit none

      character :: chara*15

      integer :: i, j, ix, iy, it, is, js
      integer :: kx, ky, kkx, kky, mu, nu
      integer :: idummy, nzone
      
      real*8 :: volk, dummy
      real*8 :: ratio
      real*8 :: dnewnuu(Ns),dnewndd(Ns)
      real*8 :: ddtemp, ddmu
      real*8 :: akx, aky
      real*8 :: dummy1, dummy2
      real*8 :: w(5)
      complex*16 :: zdummy

      complex*16, external :: Zwvfn

      write(*,*) 'label site'
      call labelsite()
      write(*,*) 'hopping integral'
      call hopping_integral(fname(4))

      dnewnuu(:) = 0.0d0
      dnewndd(:) = 0.0d0

      write(*,*)'MF calculation'
      write(*,*)'Maxit = ',Maxit
      write(*,*)'Conv = ',Conv

      nzone = Nkx / Nqx
c      nzone = Nkx

      zpsiall(:,:,:,:,1)=0.0d0

c## ITERATION[
c      it = 1 ; conrs = 1000.0d0
c      write(*,*)
c      do while ((it < Maxit).and.(conrs >= Conv))

      kx = 0
      do while (kx >= 0)
 100     write(*,*) 'Enter  kx,ky (Nkx,Nky=', Nkx, Nky, ')'

c         read(5,*) kx, ky
         read(5,*,err=100) dummy1, dummy2
         kx = INT(dummy1 * Nkx)
         ky = INT(dummy2 * Nky)
         
         write(*,*) 'kx,ky = ', kx, ky
         if (kx*ky < 0) cycle
c         call hamiltonian(kx,ky,1)

         akx = 2.0d0 * Pi * DBLE(kx) / DBLE(Nkx)
         aky = 2.0d0 * Pi * DBLE(ky) / DBLE(Nky)
         call  diag(Eall(kx,ky,:,1),Zpsiall(kx,ky,:,:,1),akx,aky,1)
         call  diag(Eall(kx,ky,:,2),Zpsiall(kx,ky,:,:,2),akx,aky,2)


         write(*,'(a,a)')'  band    ',
     &        '  3z^2-r^2       zx         yz        x^2-y^2       xy'
         do i = 6, 7
            dummy1 = 0.0d0
            do mu = 1, Ns
               nu = MOD(mu + Ns, Ns2)
               if (nu == 0) nu = Ns2
               if (nqx == 1) nu = mu
               w(mu) = ABS(Zwvfn(kx,ky,mu,i))**2
     &              +ABS(Zwvfn(kx,ky,nu,i))**2
               dummy1 = dummy1 + w(mu)
            end do
            write(*,'(a,i5,5f12.6)') 'i= ',i, (w(is), is=1, Ns)
         end do
      end do

      return
      end
c###################################################################

c###################################################################
c##### Seigenk2:  calculate eigen energies at a specific k point
c###################################################################

      subroutine Seigenk2(nx1, nx2, ny1, ny2, mkx, mky, ispin, dkx,dky)

      use common, only : Nkx, Nky, Nqx, Nqy, Dne, Ns, Ns2, Dmu, 
     &     Dnuu, Dndd, Pi
      implicit none

      integer, intent(in) :: nx1, nx2, ny1, ny2, ispin, mkx, mky
      real*8, intent(in) :: dkx, dky

      integer :: kkx, kky, kx, ky
      real*8 :: eig(nx1:nx2,ny1:ny2,Ns*Nqx,2)
      complex*16 :: zpsi(Ns*Nqx,Ns*Nqx,2)

      character :: chara*15
      integer :: is, np, i
      real*8 :: ene, de, dummy, tilt, akx, aky, dp


!$omp parallel do private(kky, kx, ky) shared(eig,zpsi)
      do kkx = nx1, nx2 ; do kky = ny1, ny2
         kx = kkx
         ky = kky
         akx = 2.0d0 * Pi * DBLE(kx) / DBLE(mkx) + dkx
         aky = 2.0d0 * Pi * DBLE(ky) / DBLE(mky) + dky
         call  diag(eig(kx,ky,:,ispin),zpsi(:,:,ispin),akx,aky,ispin)
      end do ; end do
!$omp end parallel do

      open(10,file="out.eplot3d.dat")
      dummy = 0.0d0
      np = 0
      do kx = nx1, nx2
         do ky = ny1, ny2
            write(10,'(13g20.10)') kx/DBLE(mkx)+dkx, ky/DBLE(mky), 
     &           (eig(kx,ky,i,1)-Dmu,i=1,Ns2),dummy
            if ((eig(kx,ky,6,1) > Dmu).or.(eig(kx,ky,7,1) < Dmu))
     &           np = np + 1
         end do
         write(10,*)
      end do
      close(10)

      de = 0.01d0
      
      chara = 'out.fermid.dat'
      ene = -0.0d0
      call fermi2k(chara, ene, eig,nx1,nx2,ny1,ny2,mkx,mky)
      do is = -3, 3
         ene = DBLE(is) * de
         if(is == -3) chara(14:15) = '-3'
         if(is == -2) chara(14:15) = '-2'
         if(is == -1) chara(14:15) = '-1'
         if(is ==  0) chara(14:15) = '+0'
         if(is ==  1) chara(14:15) = '+1'
         if(is ==  2) chara(14:15) = '+2'
         if(is ==  3) chara(14:15) = '+3'
         call fermi2k(chara, ene, eig,nx1,nx2,ny1,ny2,mkx,mky)
      end do
      
      ky = 0
      if (ny1 > 0) ky = ny1
      if (ny2 < 0) ky = ny1

      open(20,file="out.dirac.dat")
      do kx = nx1,nx2
         dummy = kx / DBLE(mkx) + dkx
         write(20,'(13g20.10)') dummy,
     &        (eig(kx,ky,i,1) - Dmu, i=1,Ns*Nqx)
      end do
      close(20)
      

      call Diracposition(nx1,nx2,mkx,dkx,de,dp,eig(:,0,:,1))

      call density()
      write(*,'(a,f12.5)') 'Dirac point energy=',de-Dmu
      write(*,'(a,f12.7)') 'Dirac point position=',  dp
      write(*,'(a)') 'density, energy, position, area, M'
      write(*,'(2f12.5,f12.7,2f12.5)') Dne, de - Dmu, dp,
     &     100.0d0 * np / DBLE(mkx*mky),
     &     SUM(Dnuu(:)-Dndd(:))
      

      
      return
      end

c###################################################################

c###################################################################
c##### Fermi_velocity:
c###################################################################

      subroutine fermi_velocity()

      use common, only : Ns2, Pi
      implicit none

      character(len=32) :: filein, fileout
      real*8 :: akx, aky, dkx, dky, ddmu
      real*8 :: vx(Ns2), vy(Ns2)

      integer :: i, j, ispin, ngomi
      character :: chara*74, chgomi(Ns2)
      

      write(*,*) 'Enter input filename (fermi surface data)'
      read(5,*) filein
      if (filein(1:1) == '!') filein = 'out.fermid.dat'
      fileout = 'v'//filein
      open(10,file=filein)
      open(20,file=fileout)
      chara = '#'

      do while (chara(1:1) == '#')
         read(10,*) chara
         if( chara(1:4) == '#Dmu') then
            backspace(10)
            read(10,*) chara, ddmu
         end if
      end do
      write(20,'(a,f12.8)') '#Dmu=', ddmu
      
      ngomi = 0
      do while(chara(1:1) == '?')
         if (chara(1:1) == '?') then 
            ngomi = ngomi + 1
            backspace(10)
            read(10,*) (chgomi(i), i=1,ngomi),  chara
         end if
      end do

      do while (chara(1:4) /= '#EOF')

         backspace(10)
         if (ngomi == 0) then
            read(10,*) akx, aky
         else
            read(10,*) (chgomi(i), i=1,ngomi), akx, aky
         end if

         dkx = akx * 2.0d0 * Pi
         dky = aky * 2.0d0 * Pi
         call Svelo(dkx, dky, vx, vy, 1)
         
         if (ngomi == 0) then
            write(20,'(20f12.8)') akx, aky, vx(6), vy(6), vx(7), vy(7)
         else
            write(20,'(20f12.8)') akx, aky, 
     &           vx(ngomi/2+1), vy(ngomi/2+1)
         end if

         read(10,*,err=999) chara
      end do
      write(20,'(a)') chara


      close(10)
      close(20)

 999  return
      end

c###################################################################

c###################################################################
c##### Svelo:
c###################################################################

      subroutine Svelo(akx, aky, vx, vy, ispn)

      use common, only : Nqx, Nqy, Ns, Ns2
      implicit none

      integer, intent(in) :: ispn
      real*8, intent(in) :: akx, aky
      real*8, intent(out) :: vx(Ns*Nqx), vy(Ns*Nqx)

      integer :: kkx, kky, kx, ky
      real*8 :: eig(Ns*Nqx,2)
      complex*16 :: zpsi(Ns*Nqx,Ns*Nqx,2)

      character :: chara*15
      integer :: i, j, ispin

       ispin = ispn
      call diag(eig(:,ispin),zpsi(:,:,ispin),akx,aky,ispin)
      call velocity(vx,zpsi(:,:,ispin),akx,aky,'x')
      call velocity(vy,zpsi(:,:,ispin),akx,aky,'y')
      

      return
      end

c###################################################################

c###################################################################
c##### tiltd:  
c###################################################################

      subroutine tiltd()

      use common, only : Nkx, Nky, Nqx, Nqy, Dne, Ns, Ns2, Dmu
      implicit none



c      b1 = bx - ax
c      b2 = by - ay
c      d1 = dx - ax
c      d2 = dy - ay
c
c      a = d1**2 + d2**2
c      a = a / ((b1 - d1)**2 + (b2 - d2)**2)
c      

c      beta = b2 / b1

c      dummy = - a * (a - 1.0) * (1.0 + beta**2) * (b1**2 + b2**2)
c      dummy = dummy + a**2 * (b1 + b2 * beta)
c      if (dummy < 0.0d0) then
c         write(*,*) 'error'
c         stop
c      else
c         dummy = SQRT(dummy)
c      end if

c      x = (a * (b1 + b2 * beta) + dummy)/ (a-1.0) / (1.0 + beta**2)
c      if ((x <= 0).or.(x >= b1)) then
c         x = (a * (b1 + b2 * beta) - dummy)/ (a-1.0) / (1.0 + beta**2)
c      end if
c      y = beta * x
            
      
      return
      end

c###################################################################

c###################################################################
c##### Sarea:  
c###################################################################

      subroutine Sarea()

      use common, only : Nkx, Nky, Nqx, Nqy, Dne, Ns, Ns2, Dmu
      implicit none



c      do i = 1, Ns*Nqx
c         do ikx = nx1, nx2 ; do iky = ny1, ny2
c            dummy = eig(ikx,iky,i)
c            if (eorh == 'h') then
c               if (dummy > Dmu) np = np + 1
c            else
c               if (dummy < Dmu) np = np + 1
c            end if
c         end do ; end do
c      end do

      return
      end

c###################################################################

c###################################################################
c##### Diracposition:  
c###################################################################

      subroutine Diracposition(nx1,nx2,mkx,dkx,de,dp,eig)

      use common, only : Nkx, Nky, Nqx, Nqy, Dne, Ns, Ns2, Dmu, 
     &     Dnuu, Dndd, Pi
      implicit none

      integer, intent(in) :: nx1, nx2, mkx
      real*8, intent(in) :: dkx
      real*8, intent(out) :: dp, de

      integer :: kkx, kx
      real*8 :: eig(nx1:nx2,Ns2*2)

      character :: chara*15
      integer :: is, np, i
      real*8 :: ene, dummy, tilt, akx, aky


      de = 1.0d10
      do kx = nx1, nx2
         dummy = kx / DBLE(mkx) + dkx
         dummy = ABS(eig(kx,6)-eig(kx,7))
         if (dummy >= de) then
            kkx = kx
            if (dummy < 0.1d0) exit
         end if
         de = dummy
      end do

      kx = kkx - 1
      kkx = kkx + 1

      de = eig(kx,7)*eig(kkx,7)-eig(kx,6)*eig(kkx,6)
      de = de /
     *     (ABS(eig(kx,7)-eig(kx,6))
     &     +  ABS(eig(kkx,7)-eig(kkx,6)))

      dp = (de-eig(kx,6)) / (eig(kx,6) - eig(kkx,7))
      dp = dp * DBLE(kx-kkx) + DBLE(kx)
      dp = dp / DBLE(mkx) + dkx


      return
      end

c###################################################################

c###################################################################
c##### Seigen2:                                    
c###################################################################

      subroutine Seigen2(nx, ny, ispin)

      use common, only : Nkx, Nky, Nqx, Nqy, Eall, Zpsiall, Ns, Ns2,
     &     Pi
      implicit none

      integer, intent(in) :: nx, ny, ispin
      integer :: kkx, kky, kx, ky
      integer :: iqq, nsQ
      integer :: nn0, nn1, iQ


      real*8 :: akx, aky


!$omp parallel do private(akx,aky,kkx, kky, kx, ky) shared(Eall,Zpsiall)
      do kkx = 0, nx-1 ; do kky = 0, ny-1
         kx = kkx
         ky = kky
         akx = 2.0d0 * Pi * DBLE(kx) / DBLE(Nkx)
         aky = 2.0d0 * Pi * DBLE(ky) / DBLE(Nky)
c         call  h_test(Eall(kx,ky,:,:),Zpsiall(kx,ky,:,:,:),
c     &        akx,aky,ispin)
c     call  hamiltoniank_any3(Eall(kx,ky,:,:),Zpsiall(kx,ky,:,:,:),
c     &        akx,aky,ispin)
         call  diag(Eall(kx,ky,:,ispin),Zpsiall(kx,ky,:,:,ispin)
     &        ,akx,aky,ispin)
      end do ; end do
!$omp end parallel do
      
      if ((nx /= Nkx).or.(ny /= Nky).or.(Nqx /= 1)) then
         call Seigen_expand(Nkx,Nky,nx,ny,
     &        Eall(:,:,:,ispin),Zpsiall(:,:,:,:,ispin))
      end if 

      call boundary()


      return
      end

c###################################################################

c###################################################################
c##### Seigen_expand:                                    
c###################################################################

      subroutine Seigen_expand(mx,my,mmx,mmy,eig,zpsi)

      use common, only : Nqx, Nqy, Ns
      implicit none

      integer, intent(in) :: mx, my, mmx, mmy
      real*8, intent(inout) :: eig(0:mx,0:my,Ns*Nqx)
      complex*16, intent(inout) :: zpsi(0:mx,0:my,Ns*Nqx,Ns*Nqx)

      integer :: kkx, kky, kx, ky
      integer :: nn0, nn1, iQ, iqq, nsQ

      nsQ = Ns * Nqx

      do iQ = 1, Nqx-1
         do kx = 0, mmx - 1 ; do ky = 0, mmy - 1
            kkx = MOD(kx + iQ * mx / Nqx, mx)
            kky = MOD(ky + iQ * my / Nqy, my)
            
            do iqq = 0, Nqx-1
               nn0 = Ns * iqq
               nn1 = Ns * MOD(iqq+iQ,Nqx)
               eig(       kkx, kky, :) 
     &              = eig( kx,  ky, :)
               
               zpsi(    kkx, kky, 1+nn0:Ns+nn0, :)= 
     &              zpsi(kx,  ky, 1+nn1:Ns+nn1, :)
            end do
         end do ; end do
      end do

      return
      end

c###################################################################

c###################################################################
c##### set_hopping:
c###################################################################
      subroutine set_hopping()
      
      use common, only : fname
      implicit none

      write(*,*) 'label site'
      call labelsite()
      write(*,*) 'hopping integral'
      call hopping_integralxyz(fname(4)) !xyz


      return
      end
c###################################################################


c%%%%
c%%%% END_OF_FILE:mfcalc2.f
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
