c###################################################################
c##### calcMeanField: meanfield calculation
c###################################################################
      subroutine calcMeanField()

      use common
      implicit none

      character :: cdwmode*3, n0mode*3

      integer :: i, it, mu, ispin

      real*8 :: dnewnuu(Nband),dnewndd(Nband)
      real*8 :: ddmu, dnop0
      complex*16 :: zop2temp(Nband,Nband,2), zop0temp(Nband,Nband,2)

      real*8, external :: fermiDirac, Dntotal, diffDens

      cdwmode = 'off'
      n0mode = 'on'

c==== cdw, n0 モード off の場合,hamiltonian 行列にこれらの項を含めない．
      if (cdwmode(1:2) /= 'on') then
          Zop2(:,:,:) = 0.0d0
      end if
      if (n0mode(1:2) /= 'on') then
          Zop0(:,:,:) = 0.0d0
      end if

      dnewnuu(:) = 0.0d0
      dnewndd(:) = 0.0d0

      write(*,*)' calculating Meanfield ...'
      write(*,*)
      print("('  MaxIteration = ', I5)"), maxIter
      print("('  Convergence  = ', F5.3)"), conv
      print("('  Initial Dmu  = ', F5.3)"), Dmu

      write(*,*)
      write(*,*)'=========== Iterations ==========='
      zpsiall(:,:,:,:,:)=0.0d0
c## ITERATION[
      it = 1 ; conrs = 1000.0d0
      do while ((it < maxIter).and.(conrs >= Conv))

         write(*,*)
         call calcEigenvalue(Nkx/Nqx, Nky,1)
         call calcEigenvalue(Nkx/Nqx, Nky,2)
         write(*,*)' calculating Eigenvalue is done ...'

         call calcChemical(Dmu)
         write(*,*) ' estimating chemical potential is done ...'

         write(*,*)

         if (it == 1) then
            ddmu = 1.0d-8
         end if

c## NEW ORDER PARAMETER[
         Dens(:,:) = 0.0d0

         if (Nqx == 1) then
            !** Para Magnetic Mode
            Zopnew(:,:,:) = 0.0d0
            Zopnew0(:,:,:) = 0.0d0
            Zopnew2(:,:,:) = 0.0d0
         else
            !** Not Para Magnetic mode
            dnop0 = 0.0d0
            call calcModifiedDens(dnop0)
         end if

         Zopnew(:,:,:) = Zopnew(:,:,:) / DBLE(Nkx * Nky)
         Zopnew0(:,:,:) = Zopnew0(:,:,:) / DBLE(Nkx * Nky)
         Zopnew2(:,:,:) = Zopnew2(:,:,:) / DBLE(Nkx * Nky)

         do mu = 1, Nband
            Dnuu(mu) = Zopnew(mu,mu,1)
            Dndd(mu) = Zopnew(mu,mu,2)
            dnop0 = dnop0 + Zopnew0(mu,mu,1)+Zopnew0(mu,mu,2)
         end do

         Dens(:,:) = Dens(:,:) / DBLE(Nkx * Nky)
         Dne1 = Dntotal(Dmu) * DBLE(Nband)

         if (Nqx /= 1) then
            print("('  Density   = ', F5.3)"), dnop0
         endif
c## NEW ORDER PARAMETER]

         call calcSelfConsistent(zop2temp(:,:,:), zop0temp(:,:,:), it)
         it = it +1

      end do
c## ITERATION]

      write(*,*)
      write(*,*)'=========== Iterations ==========='
      write(*,*)

c==== cdw, n0 モード off の場合,hamiltonian 行列にこれらの項を含めない．
      if (cdwmode /= 'on') then
          Zop2(:,:,:) = 0.0d0
      end if
      if (n0mode /= 'on') then
          Zop0(:,:,:) = 0.0d0
      end if

      call calcDens()

      if (conrs >= Conv) then
         write(*,*)
         write(*,*)'#### NOT CONVERGED YET #####'
         write(*,*)
         return
      end if

      return
      end
c###################################################################

c###################################################################
c##### calcSelfConsistent:
c###################################################################

      subroutine calcSelfConsistent(zop0temp,zop2temp,it)

      use common
      implicit none

      character :: fnameTemp*40
      character :: cdwmode*3, n0mode*3
      real*8 :: diff, ratio
      integer :: i
      integer,intent(in) :: it
      complex*16,intent(out) :: zop0temp(Nband,Nband,2)
      complex*16,intent(out) :: zop2temp(Nband,Nband,2)

      !** Commonに入れたい
      cdwmode = 'off'
      n0mode = 'on'

      diff = 0.0d0
      conrs = -1.0d0
      conrs = MAXVAL(ABS(Zopnew(:,:,:)-Zop(:,:,:)))

      if (cdwmode == 'on') then
         diff = MAXVAL(ABS(Zopnew2(:,:,:)-Zop2(:,:,:)))
         conrs = MAX(conrs,diff)
      end if
      if (n0mode == 'on') then
         diff = MAXVAL(ABS(Zopnew0(:,:,:)-Zop0(:,:,:)))
         conrs = MAX(conrs,diff)
      end if

      print("('  Iteration = ', I5)"), it
      print("('  Conrs     = ', F5.3)"), conrs
      print("('  Dmu = ', F5.3)"), Dmu
      print("('  M   = ', F5.3, ' muB')"), SUM(Dnuu(:)-Dndd(:))

      ratio = 0.3d0
      Zop(:,:,:) = Zopnew(:,:,:) * (1.0d0 - ratio)
     &      + ratio * Zop(:,:,:)

      ratio = 0.0d0
      Zop0(:,:,:) = Zopnew0(:,:,:) * (1.0d0 - ratio)
     &      + ratio * Zop0(:,:,:)

      Zop2(:,:,:) = Zopnew2(:,:,:) * (1.0d0 - ratio)
     &      + ratio * Zop2(:,:,:)

      Zopnew(:,:,:) = 0.0d0
      Zopnew0(:,:,:) = 0.0d0
      Zopnew2(:,:,:) = 0.0d0
      do i = 1, Nband
         Dnuu(i) = DBLE(Zopnew(i,i,1))+DIMAG(Zopnew(i,i,1))
         Dndd(i) = DBLE(Zopnew(i,i,2))+DIMAG(Zopnew(i,i,2))
      end do

      zop2temp(:,:,:) = Zop2(:,:,:)
      if (cdwmode /= 'on') then
          Zop2(:,:,:) = 0.0d0
      end if
      zop0temp(:,:,:) = Zop0(:,:,:)
      if (n0mode /= 'on') then
          Zop0(:,:,:) = 0.0d0
      end if

c## output temporary result[
      if(MOD(it,10) == 0)then
         Zop2(:,:,:) = zop2temp(:,:,:)
         Zop0(:,:,:) = zop0temp(:,:,:)
         fnameTemp = 'temp.dat'
         !** call saveData(fnameTemp)
         if (cdwmode /= 'on') then
            if (conrs >= Conv) Zop2(:,:,:) = 0.0d0
         end if
         if (n0mode /= 'on') then
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

      return
      end

c###################################################################

c###################################################################
c##### calcModifiedDensities: calc Densities using
c###################################################################

      subroutine calcModifiedDens(dnop0)

      use common
      implicit none

      real*8, intent(inout) :: dnop0

      integer :: i, is, js, js0, js2
      integer :: kx, ky, mu, nu

      real*8, external :: fermiDirac, Dntotal

      do mu = 1, Nband ; do nu = 1, Nband
         do kx = 0, Nkx - 1 ; do ky = 0, Nky - 1
            is = nu
            js0 = mu
            js = mu + Nband
            js2 = mu + Nband * 2
            if (Nqx == 2) js2 = js0
            do i = 1, Nband * Nqx

               Zopnew(mu,nu,:) = Zopnew(mu,nu,:)
     &               + DCONJG(Zpsiall(kx,ky, js, i,:))
     &               *        Zpsiall(kx,ky, is, i,:)
     &               * fermiDirac(Eall(kx,ky,i,:)-Dmu,kT)

               Zopnew0(mu,nu,:) = Zopnew0(mu,nu,:)
     &               + DCONJG(Zpsiall(kx,ky, js0, i,:))
     &               *        Zpsiall(kx,ky,  is, i,:)
     &               * fermiDirac(Eall(kx,ky,i,:)-Dmu,kT)

               Zopnew2(mu,nu,:) = Zopnew2(mu,nu,:)
     &               + DCONJG(Zpsiall(kx,ky, js2, i,:))
     &               *        Zpsiall(kx,ky,  is, i,:)
     &               * fermiDirac(Eall(kx,ky,i,:)-Dmu,kT)
               if (mu == nu) then
                  Dens(mu,:) = Dens(mu,:)
     &                  + DCONJG(Zpsiall(kx,ky,is,i,:))
     &                  *       Zpsiall(kx,ky,is,i,:)
     &                  * fermiDirac(Eall(kx,ky,i,:)-Dmu,kT)
               end if
            end do
         end do ; end do
      end do ; end do

      return
      end
c###################################################################

c###################################################################
c##### density:
c###################################################################

      subroutine calcDens()

      use common, only : Nkx, Nky, Nband, Nqx,
     &     Eall, Zpsiall, Zdens, Dmu, kT, Dens, Zop
      implicit none

      integer :: i, j, kx, ky, mu, nu, ispin
      real*8 :: sum, dn
      complex*16 :: zdummy
      real*8, external :: fermiDirac

      Zdens(:,:) = 0.0d0
      do mu = 1, Nband ; do nu = 1, Nband
         do kx = 0, Nkx - 1 ; do ky = 0, Nky - 1
            do i = 1, Nband*Nqx ; do ispin = 1, 2
               Zdens(mu,nu) = Zdens(mu,nu)
     &              + DCONJG(Zpsiall(kx,ky, mu, i,ispin))
     &              *        Zpsiall(kx,ky, nu, i,ispin)
     &              * fermiDirac(Eall(kx,ky,i,ispin)-Dmu,kT)
            end do ; end do
         end do ; end do
      end do ; end do
      Zdens(1:Nband,1:Nband) = Zdens(1:Nband,1:Nband) / DBLE(Nkx*Nky*2)
      ! /2 for spin
      do mu = 1, Nband
         Dens(mu,1) = DBLE(Zdens(mu,mu) + Zop(mu,mu,1))
         Dens(mu,2) = DBLE(Zdens(mu,mu) - Zop(mu,mu,1))
      end do

      dn = 0.0d0
      do mu = 1, Nband
         sum = DIMAG(Zdens(mu,mu))
         if (sum > 1.0d-13) then
            write(*,*) 'ERROR in density'
            write(*,*) mu, sum
            read(5,*)
         end if
         sum = Dens(mu,1) + Dens(mu,2) ! n(mu,up) + n(mu,dn)
         dn = dn + sum
      end do
      write(*,*) 'density =', dn


      return
      end

c###################################################################

c%%%%
c%%%% END_OF_FILE:mfcalc2.f
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
