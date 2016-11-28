
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
      isAllocated = allocated(t)
      if (isAllocated .EQV. .true.) deallocate(t)
      isAllocated = allocated(Eorbit)
      if (isAllocated .EQV. .true.) deallocate(Eorbit)
      isAllocated = allocated(Isitex)
      if (isAllocated .EQV. .true.) deallocate(Isitex)
      isAllocated = allocated(Isitey)
      if (isAllocated .EQV. .true.) deallocate(Isitey)
      isAllocated = allocated(Dnuu)
      if (isAllocated .EQV. .true.) deallocate(Dnuu)
      isAllocated = allocated(Dndd)
      if (isAllocated .EQV. .true.) deallocate(Dndd)
      isAllocated = allocated(Dens)
      if (isAllocated .EQV. .true.) deallocate(Dens)
      isAllocated = allocated(Zwfprod)
      if (isAllocated .EQV. .true.) deallocate(Zwfprod)
      isAllocated = allocated(Zop)
      if (isAllocated .EQV. .true.) deallocate(Zop)
      isAllocated = allocated(Zop0)
      if (isAllocated .EQV. .true.) deallocate(Zop0)
      isAllocated = allocated(Zop2)
      if (isAllocated .EQV. .true.) deallocate(Zop2)
      isAllocated = allocated(Zopnew0)
      if (isAllocated .EQV. .true.) deallocate(Zopnew0)
      isAllocated = allocated(Zopnew2)
      if (isAllocated .EQV. .true.) deallocate(Zopnew2)
      isAllocated = allocated(Zdens)
      if (isAllocated .EQV. .true.) deallocate(Zdens)

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
      allocate( t(Nsite,Nband,Nband) )
      allocate( Eorbit(Nband) )
      allocate( Isitex(0:Nsite-1) )
      allocate( Isitey(0:Nsite-1) )
      allocate( Dnuu(Nband) )
      allocate( Dndd(Nband) )
      allocate( Dens(Nband,2) )
      allocate( Zwfprod(0:Nkx,0:Nky,Nband,Nband,Nband,2) )
      allocate( Zop(Nband,Nband,2) )
      allocate( Zop0(Nband,Nband,2) )
      allocate( Zop2(Nband,Nband,2) )
      allocate( Zopnew0(Nband,Nband,2) )
      allocate( Zopnew2(Nband,Nband,2) )
      allocate( Zdens(Nband,Nband) )

      return
      end
c###########################################################
