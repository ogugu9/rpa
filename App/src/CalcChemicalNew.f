c######################################################################
c#### Smuval: estimate the chemical potential
c######################################################################

      subroutine calcChemical(ddmu)
      implicit none

      real*8, intent(inout) :: ddmu
      real*8 :: dx1, dx2, dxacc, d1, d2
      real*8, external :: Drtbis, Dneq

      if (Dneq(ddmu) == 0.0d0) then
         return
      else
         dx1 = ddmu - 1.0d0
         dx2 = ddmu + 1.0d0
         if (ddmu == 0.0d0) then
            dx1 = -10.d0
            dx2 = 10.d0
         end if
         dxacc = 1.0d-10

         d1 = Dneq(Dx1)
         d2 = Dneq(Dx2)

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
            d1 = Dneq(dx1)
            d2 = Dneq(dx2)
         end do

         ddmu = Drtbis(Dneq,dx1,dx2,Dxacc)
      end if

      return
      end
c######################################################################

c **********************************
      double precision function Dneq(ddmu)

      use common, only: Dne, Nband
      implicit none

      real*8, intent(in) :: ddmu
      real*8 :: vol
      real*8, external :: Dntotal

      vol = DBLE(Nband)
      Dneq = Dntotal(ddmu) * vol - Dne
      return
      end
c **********************************
      double precision function Dntotal(ddmu)

      use common, only: Nkx, Nky, Nband, Eall, kT, Nqx
      implicit none

      real*8, intent(in) :: ddmu
      integer :: kx, ky, i, ispin
      real*8 :: dummy, volk0
      real*8, external :: fermiDirac

      dummy = 0.0d0

      volk0 = DBLE(Nkx * Nky)
      do kx = 0, Nkx - 1 ; do ky = 0, Nky - 1
         do i = 1, Nband * Nqx ; do ispin = 1, 2
            dummy = dummy + Dffn(Eall(kx,ky,i,ispin) - ddmu, kT)
         end do ; end do
      end do ; end do
      dummy = dummy

      Dntotal = dummy  / DBLE(Nband * Nqx) / volk0
c      Dntot = Dntot * 2.0d0 ! 2.0 for spin

      return
      end
c*******************************************************************
      double precision function fermiDirac(Degn,Dtemp)
c*** *** Fermi function  (Degn:energy, Dtemp:Temperature)

c      use common, only : Deltae
      implicit none

      real*8, intent(in) :: Degn, Dtemp
      real*8 :: De, ddd

      ddd = 1.0d-15
        if (Dtemp >= 1.0d-15) then         ! finite temperature
          De = Degn / Dtemp
          if (De > 80.d0)  fermiDirac = 0.d0
          if (De < -80.d0) fermiDirac = 1.d0
          if (DABS(De) <= 80.d0) fermiDirac = 1.d0 / (1.d0 + DEXP(De))
         else                          ! zero temperature
          if (Degn > 0.0d0) fermiDirac = 0.d0
          if (Degn < 0.0d0) fermiDirac = 1.d0
c          if (Degn == 0.0d0) Dffn = 1.d0
          if (ABS(Degn) <= ddd) fermiDirac = 0.5d0
        end if

        return
        end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      double precision FUNCTION Drtbis(Dneq,x1,x2,xacc)
      INTEGER JMAX
      REAL*8 x1,x2,xacc,Dneq
      EXTERNAL Dneq
      PARAMETER (JMAX=40)
      INTEGER j
      REAL*8 dx,f,fmid,xmid

      fmid=Dneq(x2)
      f=Dneq(x1)
      if(f*fmid.ge.0.d0) pause 'root must be bracketed in Drtbis'
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
        fmid=Dneq(xmid)
        if(fmid.le.0.d0)Drtbis=xmid
        if(dabs(dx).lt.xacc .or. fmid.eq.0.d0) return
 11   continue
      pause 'too many bisections in Drtbis'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .2#.
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
