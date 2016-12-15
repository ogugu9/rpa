c###################################################################
c##### calcEigen: calc Eigen vector and value in all k-points
c###################################################################

      subroutine calcEigenvalue(nx, ny, nz, iq, ispin)

      use common, only : Nkx, Nky, Nkz, Nkpt, Nqx, Nqy, Nqz,
     &      Eall, Zpsiall, Pi, dk1, dk2, dk3, dkfrac
      implicit none

      integer, intent(in) :: nx, ny, nz, iq, ispin
      integer :: ik
      complex*16 :: ddk1, ddk2, ddk3

      do ik = 1, Nkpt
         ddk1 = dk1(ik) + dkfrac(iq,1)
         ddk2 = dk2(ik) + dkfrac(iq,2)
         ddk3 = dk3(ik) + dkfrac(iq,3)
         call  diagH(Eall(ik,iq,:,ispin),Zpsiall(ik,iq,:,:,ispin),
     &   ddk1,ddk2,ddk3,ispin)
c         write(*,'(5F8.3)'), Eall(ik,iq,:,ispin)
      end do

      return
      end

c###################################################################

c###################################################################
c##### calcEigenBand: calc Eigen vector and value for bandplot
c###################################################################

      subroutine calcBandplot(nx, ny, nz, ispin)

      use common, only : Nkx, Nky, Nkz, Nqx, Nqy, Nqz,
     &      Eband, Zpsiband, Pi, Nkpath, Nkmesh, dkfrac
      implicit none

      integer, intent(in) :: nx, ny, nz, ispin
      integer :: ik

      real*8 :: dk1, dk2, dk3

      do ik = 0, Nkpath*Nkmesh
         dk1 = dkfrac(ik,1)
         dk2 = dkfrac(ik,2)
         dk3 = dkfrac(ik,3)
         call  diagH(Eband(ik,:,ispin),Zpsiband(ik,:,:,ispin),
     &   dk1,dk2,dk3,ispin)
         !write(*,'(8F8.3)'), dkx, dky, dkz, Eband(nkpt,:,ispin)
      end do

      return
      end

c###################################################################
