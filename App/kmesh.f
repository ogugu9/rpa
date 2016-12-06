      program main
         implicit none

         integer :: ikx, iky, Nkx
         real*8 :: ds, dt, rate, dkx, dky, la, lb, Pi

         Pi = 3.14d0
         Nkx = 100
         la = 3.745d0
         lb = 9.700d0
         rate = lb**2/la**2

         open(10,file="kmesh.dat")
         do ikx = -Nkx,Nkx ; do iky = -Nkx,Nkx
c            if (DBLE(ABS(-ikx+iky))/DBLE(Nkx).lt.1) then
c               dkx = DBLE(ikx+iky)/DBLE(Nkx)*2.0d0*Pi/la
c               dky = DBLE(-ikx+iky)/DBLE(Nkx)*2.0d0*Pi/lb
c               write(10,'(2F8.3)') dkx,
            dkx = DBLE(ikx)/DBLE(Nkx)
            dky = DBLE(iky)/DBLE(Nkx)
            if ((ABS(dky)<=-rate*ABS(dkx)+(1.0+rate)/2.0)) then
               write(10,'(6F8.3)')
     &         ikx, iky, dkx*2.0d0*Pi/la, dky*2.0d0*Pi/lb
            endif
c            endif
         end do ; end do
         close(10)

c        symmetric k-points
         open(20,file="sympt.dat")
         write(20,'(2F8.3)') 0.58d0*2.0d0*Pi/la, 0.00d0*2.0d0*Pi/lb !D
         write(20,'(2F8.3)') 0.51d0*2.0d0*Pi/la, -0.51d0*2.0d0*Pi/lb !S
         write(20,'(2F8.3)') 0.42d0*2.0d0*Pi/la, -1.00d0*2.0d0*Pi/lb !R
         write(20,'(2F8.3)') 0.00d0*2.0d0*Pi/la, -1.00d0*2.0d0*Pi/lb !T
         close(20)
      end
