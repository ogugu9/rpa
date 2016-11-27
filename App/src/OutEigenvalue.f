
c#########################################################
c#### eigenout: output eigen energy & eigen function
c#########################################################

      subroutine outEigenvalue()

      use common, only : Nkx, Nky, Nband, Nqx, Nqy,
      &     Eall, Zpsiall, Dmu, Dne, kT
      implicit none

      integer :: i, j, k, ikx, iky, kx, ky, ispin
      integer :: lmtl, lmtu, mx, my
      real*8 :: p1, p2, p3, p4, p5, p6, emin, emax, e0, drl, dim
      real*8 :: dummy1, dummy2

      goto 100
      open(10,file = "eigen.out")
      write(10,*) Nkx, Nky
      write(10,*) 'Dne=', Dne
      write(10,*) 'Dmu=', Dmu
      write(10,*) 'kT=', kT

c============================[
c== output data for gnuplot
c==
      open(20,file="out.eplot.dat")

      if (mx == 0) mx = 1
      if (my == 0) my = 1

      do ispin = 1,2 ; do mu = 1, Nband * Nqx
         do nkpt = 1, Nkpath*Nkmesh

            write(20,*) p1, dummy1, dummy2
         end do
         write(20,1)
      end do ; end do

      lmtu = INT(emax - e0) + 1
      lmtl = INT(emin - e0) - 1

      write(20,*) "? ? ?", 0.0d0, lmtu
      write(20,*) "? ? ?", 0.0d0, lmtl
      write(20,1)
      write(20,*) "? ? ?", p1, lmtu
      write(20,*) "? ? ?", p1, lmtl
      write(20,1)
      write(20,*) "? ? ?", p2, lmtu
      write(20,*) "? ? ?", p2, lmtl
      write(20,1)
      write(20,*) "? ? ?", p3, lmtu
      write(20,*) "? ? ?", p3, lmtl
      write(20,1)
      write(20,*) "? ? ?", p4, lmtu
      write(20,*) "? ? ?", p4, lmtl
      write(20,1)
      write(20,*) "? ? ?", p5, lmtu
      write(20,*) "? ? ?", p5, lmtl
      write(20,1)
      close(20)

c================================]

      return
 1    format(i0)
      end
