c######################################################################
c#### estimateChemicalPotential:
c####
c####    - 化学ポテンシャル mu を +- 1.0 ずらす
c####    - n_total(mu) と ne の差分を調べる
c####    - 差分が 0 になれば終わり、ならなければ更にずらして調べる
c####
c######################################################################

      subroutine calcChemical(ddmu)
      implicit none

      real*8, intent(inout) :: ddmu
      real*8 :: dx1, dx2, dxacc, d1, d2
      real*8, external :: Drtbis, diffDens

      if (diffDens(ddmu) == 0.0d0) then
         return
      else
         dx1 = ddmu - 1.0d0
         dx2 = ddmu + 1.0d0
         if (ddmu == 0.0d0) then
            dx1 = -10.d0
            dx2 = 10.d0
         end if
         dxacc = 1.0d-10

         d1 = diffDens(dx1)
         d2 = diffDens(dx2)

         do while (d1 * d2 >= 0.0d0)
            if (d1 == 0.d0) then
               ddmu = dx1
               return
            else if(d2 == 0.0d0) then
               ddmu = dx2
               return
            end if
            dx1 = dx1 - 1.0d0
            dx2 = dx2 + 1.0d0
            d1 = diffDens(dx1)
            d2 = diffDens(dx2)
         end do

         ddmu = Drtbis(diffDens,dx1,dx2,Dxacc)
      end if

      return
      end
c######################################################################
c     Difference Between Ntotal(mu) And Nelectlon
c **********************************
      double precision function diffDens(ddmu)

      use common, only: Dne, Nband
      implicit none

      real*8, intent(in) :: ddmu
      real*8, external :: Dntotal

      diffDens = Dntotal(ddmu) * DBLE(Nband) - Dne
      return
      end
c **********************************
      double precision function Dntotal(ddmu)

      use common, only: Nkx, Nky, Nkz, Nband, Eall, kT, Nqx
      implicit none

      real*8, intent(in) :: ddmu
      integer :: kx, ky, kz, i, ispin
      real*8 :: sum, Vk !Vk: Volume of k-space
      real*8, external :: fermiDirac

      sum = 0.0d0

      Vk = DBLE(Nkx * Nky * Nkz)
      do kx = 0, Nkx-1 ; do ky = 0, Nky-1 ; do kz = 0, Nkz-1
         do i = 1, Nband * Nqx ; do ispin = 1, 2
            sum = sum + fermiDirac(Eall(kx,ky,kz,i,ispin) - ddmu, kT)
         end do ; end do
      end do ; end do ; end do
      sum = sum

      Dntotal = sum  / DBLE(Nband * Nqx) / Vk
c      Dntot = Dntot * 2.0d0 ! 2.0 for spin

      return
      end
c*******************************************************************
      double precision function fermiDirac(Degn,DkT)
c******* Fermi function  (Degn:energy, DkT:Temperature)

      implicit none

      real*8, intent(in) :: Degn, DkT
      real*8 :: De, ddd

      ddd = 1.0d-15
        if (DkT >= 1.0d-15) then         ! finite temperature
          De = Degn / DkT
          if (De > 80.d0)  fermiDirac = 0.d0
          if (De < -80.d0) fermiDirac = 1.d0
          if (DABS(De) <= 80.d0) then
             fermiDirac = 1.d0 / (1.d0 + DEXP(De))
          endif
        else                          ! zero temperature
          if (Degn > 0.0d0) fermiDirac = 0.d0
          if (Degn < 0.0d0) fermiDirac = 1.d0
c          if (Degn == 0.0d0) Dffn = 1.d0
          if (ABS(Degn) <= ddd) fermiDirac = 0.5d0
        end if

        return
        end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c ** ここの処理が理解できてない
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      double precision FUNCTION Drtbis(diffDens,x1,x2,xacc)
      INTEGER JMAX
      REAL*8 x1,x2,xacc,diffDens
      EXTERNAL diffDens
      PARAMETER (JMAX=40)
      INTEGER j
      REAL*8 dx,f,fmid,xmid

      fmid=diffDens(x2)
      f=diffDens(x1)
      if(f*fmid.ge.0.d0) stop 'root must be bracketed in Drtbis'
      if(f.lt.0.d0)then
        Drtbis=x1
        dx=x2-x1
      else
        Drtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*0.5d0
        xmid=Drtbis+dx
        fmid=diffDens(xmid)
        if(fmid.le.0.d0)Drtbis=xmid
        if(dabs(dx).lt.xacc .or. fmid.eq.0.d0) return
 11   continue
      stop 'too many bisections in Drtbis'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .2#.
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
