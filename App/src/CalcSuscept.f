c#########################################################
c#### calcSuscepts
c####        output bare susceptibility -> val
c#########################################################

      subroutine calcSuscepts(iflag)

      use common, only : Nksize, Nkpath, Nkmesh, Nkx, Nky, Nkz, Nqx,
     &   chi0, maxOmega, initOmega, Erange, Eta, kT, Nband, Zwfprod

      implicit none
      integer, intent(in) :: iflag

      integer :: iq, mu, nQ, is, iom, ir
      integer, external :: Ivrtx
      real*8 :: Omega
      complex*16 :: val
      logical :: is_allocated

      write(*,*) 'caliculating bare-susceptibility...'


      chi0(:) = 0.0d0
      is = 1      ! Non-spin-polarization
      iom = 0     ! Omega = 0

      Omega = DBLE(iom)*Erange/DBLE(maxOmega) + initOmega
      write(*,*) 'Omega=', Omega

      is_allocated = allocated(Zwfprod)
      if (is_allocated .EQV. .true.) deallocate(Zwfprod)
      allocate( Zwfprod(Nksize,0:Nkpath*Nkmesh,Nband,Nband,Nband,2) )

      ! Gamma点だけ除いて計算する。Gamma点で発散するのが自明なため。
      ! iq = 0 or Nkpath*Nkmesh : Gamma点
      do iq = 0, Nkpath*Nkmesh

         ! Caluculate Eigenvalue at vector k+q
         call calcEigenvalue(Nkx/Nqx, Nky, Nkz,iq,1)
         call calcEigenvalue(Nkx/Nqx, Nky, Nkz,iq,2)
         ! Calculate Productions of WaveFunction
         call productWaveFunction(iq)

         do mu = 1, Nband ; do nQ = 0, Nqx-1
            ir = Ivrtx(mu,mu,nQ,'L')
            call calcChi0(ir,ir,iq,iom,is,val)
            chi0(iq) = chi0(iq) - DIMAG(val)
         end do ; end do

         write(*,'(I8,ES10.3)'), iq, chi0(iq)

      end do

      call outSuscept()

      return
      end

c#########################################################

c#########################################################
c#### calcChi0
c####        output bare susceptibility -> val
c#########################################################

      subroutine calcChi0(ir12,ir34,iq,iom,is,chi)

      use common, only : Nband, Nkpt, Nqx, Zi, Pi,
     &     maxOmega, initOmega, Erange, Eta, kT,
     &     Zwfprod, Zpsiall, Eall
      implicit none

      integer, intent(in) :: ir12, ir34, iq, iom, is
      complex*16, intent(out) :: chi

      integer :: ik, mu, nu, nQ, npQ
      integer :: isE, isH
      integer :: mu1, mu2, mu3, mu4, nQ1, nQ2, mQ
      integer :: iorb1, iorb2, iorb3, iorb4
      complex*16 :: val
      integer, external :: Ispinpair, Ivrtxinv

      real*8 :: Omega
      real*8, external :: Dffn

      mu1 = Ivrtxinv(ir12,1,'L')
      mu2 = Ivrtxinv(ir12,2,'L')
      nQ1 = Ivrtxinv(ir12,3,'L')

      mu3 = Ivrtxinv(ir34,1,'R')
      mu4 = Ivrtxinv(ir34,2,'R')
      nQ2 = Ivrtxinv(ir34,3,'R')

      isE = Ispinpair(is,'E')     !s
      isH = Ispinpair(is,'H')     !s'

      val = DCMPLX(0.0d0,0.0d0)
      chi = DCMPLX(0.0d0,0.0d0)

      Omega = DBLE(iom)*Erange/DBLE(maxOmega) + initOmega

      do ik = 1, Nkpt
         do mu = 1, Nband*Nqx ; do nu = 1, Nband*nqx

            val =
     &      ( Dffn(Eall(ik,iq,mu,isE),kT)
     &         - Dffn(Eall(ik,0,nu,isH),kT) )
     &      / ( Omega + Zi * Eta
     &         - (Eall(ik,iq,mu,isE) - Eall(ik,0,nu,isH)))

c            if((Dffn(Eall(ik,0,mu,isE),kT)
c     &      -Dffn(Eall(ik,iq,nu,isH),kT)) .ne. 0.00d0) then
c               write(*,*) (Dffn(Eall(ik,0,mu,isE),kT)
c     &            - Dffn(Eall(ik,iq,nu,isH),kT) )
c               write(*,*) ( Omega + Zi * Eta
c     &            - (Eall(ik,0,mu,isE) - Eall(ik,iq,nu,isH)))
c               write(*,*) val
c            endif

            if (ABS(val) < 1.0d-13) cycle

c            write(*,*) 'val > 1.0d-14'
            do npQ = 0, Nqx - 1     !Qp
               mQ = npQ         !Qp
               iorb1 = mu1 + mQ * Nband
               mQ = nQ1 + npQ   !Q1+Qp
               iorb2 = mu2 + MOD(mQ, Nqx) * Nband
               do nQ = 0, Nqx - 1 !Q'
                  mQ = nQ2 + nQ + npQ !Q2+Q'+Qp
                  iorb3 = mu3 + MOD(mQ, Nqx) * Nband
                  mQ = nQ + npQ !Q'+Qp
                  iorb4 = mu4 + MOD(mQ, Nqx) * Nband

c                  if(iq .ne. 0) then
c                     write(*,*)
c     &               Zwfprod(ik, 0,iorb2,iorb3,mu,isE)
c                  endif

                  chi = chi
     &                 - Zwfprod(ik, 0,iorb2,iorb3,nu,isE)
     &                 * Zwfprod(ik,iq,iorb4,iorb1,mu,isH)
     &                 * val

               end do
            end do

         end do ; end do
      end do

      chi = chi / DBLE(Nkpt)
c      if(ABS(chi) < 1.0d-13) chi = 0.0d0

      return
      end
c#########################################################
