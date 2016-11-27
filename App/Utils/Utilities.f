
c###########################################################
c##   deallocation:
c###########################################################
      subroutine deallocation()

      use common
      implicit none

      logical :: isAllocated

      isAllocated = allocated(Eall)
      if (isAllocated .EQV. .true.) deallocate(Eall)
      isAllocated = allocated(Zpsiall)
      if (isAllocated .EQV. .true.) deallocate(Zpsiall)

      return
      end
c###########################################################

c###########################################################
c##   allocation:
c###########################################################
      subroutine allocation()

      use common
      implicit none

      allocate( Eall(0:Nkx,0:Nky,Nband*Nqx,2) )
      allocate( Zpsiall(0:Nkx,0:Nky,Nband*Nqx,Nband*Nqx,2) )

      return
      end
c###########################################################
