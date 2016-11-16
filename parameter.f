
c###################################################
c#### parameter
c###################################################

      module parameters
      implicit none

      integer, parameter :: llx = 2, lly = 2, llorb = 5
      integer, parameter :: lln = (llx * 2 + 1) * (lly * 2 + 1)
      real*8, parameter :: Deltae = 1.0d-8

      end

c#####:  llorb: number of orbitals (= 5 for d orbital)
c#####:  Neighbors in [-llx:llx] and [-lly:lly] are involved
c#####:   in the hopping from (0,0) site.
c#####:  The number of such neighbors is equal to (lln - 1).

c########################################################
