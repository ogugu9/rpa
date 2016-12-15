
c###################################################################
c####   deallocation:
c###################################################################
      subroutine deallocation()

      use common
      implicit none

      logical :: isAllocated

      isAllocated = allocated(t)
      if (isAllocated .EQV. .true.) deallocate(t)
      isAllocated = allocated(Eorbit)
      if (isAllocated .EQV. .true.) deallocate(Eorbit)
      isAllocated = allocated(Isitex)
      if (isAllocated .EQV. .true.) deallocate(Isitex)
      isAllocated = allocated(Isitey)
      if (isAllocated .EQV. .true.) deallocate(Isitey)
      isAllocated = allocated(Isitez)
      if (isAllocated .EQV. .true.) deallocate(Isitez)
      isAllocated = allocated(Dnuu)
      if (isAllocated .EQV. .true.) deallocate(Dnuu)
      isAllocated = allocated(Dndd)
      if (isAllocated .EQV. .true.) deallocate(Dndd)
      isAllocated = allocated(Dens)
      if (isAllocated .EQV. .true.) deallocate(Dens)
      isAllocated = allocated(Zop)
      if (isAllocated .EQV. .true.) deallocate(Zop)
      isAllocated = allocated(Zop0)
      if (isAllocated .EQV. .true.) deallocate(Zop0)
      isAllocated = allocated(Zop2)
      if (isAllocated .EQV. .true.) deallocate(Zop2)
      isAllocated = allocated(Zopnew)
      if (isAllocated .EQV. .true.) deallocate(Zopnew)
      isAllocated = allocated(Zopnew0)
      if (isAllocated .EQV. .true.) deallocate(Zopnew0)
      isAllocated = allocated(Zopnew2)
      if (isAllocated .EQV. .true.) deallocate(Zopnew2)
      isAllocated = allocated(Zdens)
      if (isAllocated .EQV. .true.) deallocate(Zdens)
      isAllocated = allocated(dk1)
      if (isAllocated .EQV. .true.) deallocate(dk1)
      isAllocated = allocated(dk2)
      if (isAllocated .EQV. .true.) deallocate(dk2)
      isAllocated = allocated(dk3)
      if (isAllocated .EQV. .true.) deallocate(dk3)

      return
      end
c###################################################################

c###################################################################
c####   allocation:
c###################################################################
      subroutine allocation()

      use common
      implicit none

      allocate( t(Nsite,Nband,Nband) )
      allocate( Eorbit(Nband) )
      allocate( Isitex(Nsite) )
      allocate( Isitey(Nsite) )
      allocate( Isitez(Nsite) )
      allocate( Dnuu(Nband) )
      allocate( Dndd(Nband) )
      allocate( Dens(Nband,2) )
      allocate( Zop(Nband,Nband,2) )
      allocate( Zop0(Nband,Nband,2) )
      allocate( Zop2(Nband,Nband,2) )
      allocate( Zopnew(Nband,Nband,2) )
      allocate( Zopnew0(Nband,Nband,2) )
      allocate( Zopnew2(Nband,Nband,2) )
      allocate( Zdens(Nband,Nband) )
      allocate( dk1(Nksize) )
      allocate( dk2(Nksize) )
      allocate( dk3(Nksize) )

      return
      end
c###################################################################

c###################################################################
c####   deallocBandPlot:
c###################################################################
      subroutine deallocPath()

      use common
      implicit none

      logical :: isAllocated

      isAllocated = allocated(Eall)
      if (isAllocated .EQV. .true.) deallocate(Eall)
      isAllocated = allocated(Zpsiall)
      if (isAllocated .EQV. .true.) deallocate(Zpsiall)
      isAllocated = allocated(kpath)
      if (isAllocated .EQV. .true.) deallocate(kpath)
      isAllocated = allocated(dkfrac)
      if (isAllocated .EQV. .true.) deallocate(dkfrac)
      isAllocated = allocated(Eband)
      if (isAllocated .EQV. .true.) deallocate(Eband)
      isAllocated = allocated(Zpsiband)
      if (isAllocated .EQV. .true.) deallocate(Zpsiband)
      isAllocated = allocated(Zwfprod)
      if (isAllocated .EQV. .true.) deallocate(Zwfprod)
      isAllocated = allocated(chi0)
      if (isAllocated .EQV. .true.) deallocate(chi0)

      return
      end
c###################################################################

c###################################################################
c####   allocBandPlot:
c###################################################################
      subroutine allocPath()

      use common
      implicit none

      allocate( Eall(Nksize,0:Nkpath*Nkmesh,Nband*Nqx,2) )
      allocate( Zpsiall(Nksize,0:Nkpath*Nkmesh,Nband*Nqx,Nband*Nqx,2) )
      allocate( Zwfprod(Nksize,0:Nkpath*Nkmesh,Nband,Nband,Nband,2) )
      allocate( kpath(Nkpath) )
      allocate( dkfrac(0:Nkpath*Nkmesh,1:3) )
      allocate( Eband(0:Nkpath*Nkmesh,Nband*Nqx,2) )
      allocate( Zpsiband(0:Nkpath*Nkmesh,Nband*Nqx,Nband*Nqx,2) )
      allocate( chi0(0:Nkpath*Nkmesh) )

      return
      end
c###################################################################
