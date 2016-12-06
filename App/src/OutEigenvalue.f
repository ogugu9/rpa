
c#########################################################
c#### eigenout: output eigen energy & eigen function
c#########################################################

      subroutine outEigenvalue()

      use common, only : Nkx, Nky, Nkz, Nband, Nqx, Nqy, Nqz,
     &      Nkpath, Nkmesh, Eband, Dmu, kpathList, kpointList,
     &      recipLat, EF
      implicit none

      integer :: ik, mu, nkpt, nkms, ispin
      real*8 :: dk, pathLength
      real*8 :: pathVec(1:3)

c============================[
c== output data for gnuplot
c==

      open(30,file="out.band.up.dat")
      do mu = 1, Nband * Nqx
         dk = 0.0d0
         do nkpt = 1, Nkpath
            ik = 0
            pathVec(:) = kpathList(nkpt)%iniPosition
     &       - kpathList(nkpt)%finPosition
            pathLength = DOT_PRODUCT(pathVec,pathVec)
            do nkms = 1, Nkmesh
               write(30,'(F8.3,F8.3)'),
     &         dk, Eband(ik,mu*Nqx,1)-EF
               dk = dk + pathLength/DBLE(Nkmesh)
               ik = ik + 1
            end do
         end do
         write(30,1)
      end do
      close(30)

      open(40,file="out.band.dn.dat")
      do mu = 1, Nband * Nqx
         dk = 0.0d0
         ik = 0
         do nkpt = 1, Nkpath
            pathVec(:) = kpathList(nkpt)%iniPosition
     &       - kpathList(nkpt)%finPosition
            pathVec(:) = MATMUL(pathVec, recipLat)
            pathLength = DOT_PRODUCT(pathVec,pathVec)
            do nkms = 1, Nkmesh
               write(40,'(F8.3,F8.3)'),
     &         dk, Eband(ik,mu*Nqx,2)-EF
               dk = dk + pathLength/DBLE(Nkmesh)
               ik = ik + 1
            end do
         end do
         write(40,'(F8.3,F8.3)'),
     &         dk, Eband(0,mu*Nqx,2)-EF
         write(40,1)
      end do
      close(40)

c================================]

      return
 1    format(i0)
      end
