
c#########################################################
c#### eigenout: output eigen energy & eigen function
c#########################################################

      subroutine outEigenvalue()

      use common, only : Nkx, Nky, Nkz, Nband, Nqx, Nqy, Nqz,
     &      Nkpath, Nkmesh, Eband, Dmu, kpath, dkfrac,
     &      recipLat, EF
      implicit none

      integer :: ik, mu, ip, im, ispin
      real*8 :: dk, pathlen
      real*8 :: pathvec(1:3)

      open(40,file="out.band.dn.dat")
      do mu = 1, Nband * Nqx
         dk = 0.0d0
         ik = 0
         do ip = 1, Nkpath
            pathvec(:) = kpath(ip)%iniPosition
     &       - kpath(ip)%finPosition
            pathvec(:) = MATMUL(pathvec, recipLat)
            pathlen = DOT_PRODUCT(pathVec,pathvec)
            do im = 1, Nkmesh
               write(40,'(F8.3,F8.3)'),
     &         dk, Eband(ik,mu*Nqx,2)-EF
               dk = dk + pathlen/DBLE(Nkmesh)
               ik = ik + 1
            end do
         end do
         write(40,'(F8.3,F8.3)'),
     &         dk, Eband(0,mu*Nqx,2)-EF
         write(40,1)
      end do
      close(40)

      return
 1    format(i0)
      end

c#########################################################

c#########################################################
c#### eigenout: output eigen energy & eigen function
c#########################################################

      subroutine outSuscept()

      use common, only : Nkx, Nky, Nkz, Nqx, Nqy, Nqz,
     &      Nkpath, Nkmesh, chi0, kpath, dkfrac, recipLat
      implicit none

      integer :: iq, mu, ip, im, ispin
      real*8 :: dk, pathlen
      real*8 :: pathvec(1:3)

      open(40,file="out.chi0.dat")
      dk = 0.0d0
      iq = 0
      do ip = 1, Nkpath
         pathvec(:) = kpath(ip)%iniPosition
     &   - kpath(ip)%finPosition
         pathvec(:) = MATMUL(pathvec, recipLat)
         pathlen = DOT_PRODUCT(pathvec,pathvec)
         do im = 1, Nkmesh
            write(40,'(F8.3,ES12.3)') dk, chi0(iq)
            dk = dk + pathlen/DBLE(Nkmesh)
            iq = iq + 1
         end do
      end do
      write(40,'(F8.3,2ES12.3)') dk, chi0(0) ! Gamma点
      write(40,1)
      close(40)

      return
 1    format(i0)
      end

c#########################################################
