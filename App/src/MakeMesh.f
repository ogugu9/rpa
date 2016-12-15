      subroutine makeMesh()

      use common, only: Nkx, Nky, Nkz, Nkpt, dk1, dk2, dk3, Pi
      implicit none

      integer :: ik, ikx, iky, ikz
      real*8 :: ds, dt, rate, la, lb, lc
      real*8 :: dkx, dky, dkz

      la = 3.745d0
      lb = 9.700d0
      lc = 9.815d0
      rate = lb**2/la**2

      write(*,*) 'making k-mesh ...'

      open(10,file="kmesh_cart.dat")
      open(20,file="kmesh_frac.dat")

      ik = 1
      do ikx = -Nkx,Nkx ; do iky = -Nky,Nky ; do ikz = -Nkz,Nkz
         dkx = DBLE(ikx)/DBLE(Nkx)
         dky = DBLE(iky)/DBLE(Nky)
         dkz = DBLE(ikz)/DBLE(Nkz)
         if ((ABS(dky)<=-(rate-0.01)*ABS(dkx)+(1.0+rate)/2.0)) then
            dk1(ik) = (dkx-dky)/2.0d0
            dk2(ik) = (dkx+dky)/2.0d0
            dk3(ik) = dkz
            write(10,'(I5, 3F8.3)')
     &      ik, dkx*2.0d0*Pi/la, dky*2.0d0*Pi/lb, dky*2.0d0*Pi/lc
            write(20,'(I5, 3F8.3)') ik, dk1(ik), dk2(ik), dk3(ik)
            write(*,'(I5, 3F8.3)') ik, dk1(ik), dk2(ik), dk3(ik)
            ik = ik + 1
         endif
      end do ; end do ; end do

      close(10)
      close(20)

      Nkpt = ik - 1
      write(*,*) '#k-points = ', Nkpt
      write(*,*) 'making k-mesh is done ...'

c        symmetric k-points
      open(20,file="sympts.dat")
      write(20,'(2F8.3)') 0.58d0*2.0d0*Pi/la, 0.00d0*2.0d0*Pi/lb !D
      write(20,'(2F8.3)') 0.51d0*2.0d0*Pi/la, -0.51d0*2.0d0*Pi/lb !S
      write(20,'(2F8.3)') 0.42d0*2.0d0*Pi/la, -1.00d0*2.0d0*Pi/lb !R
      write(20,'(2F8.3)') 0.00d0*2.0d0*Pi/la, -1.00d0*2.0d0*Pi/lb !T
      close(20)

      return
      end
