c###################################################
c#### common statements
c###################################################

      module common

      implicit none

      ! Main
      character(len=1) :: req
      character(len=20), public :: fnameHop, fnameOut, fnameInit
      character(len=1), public :: Check
      character(len=6), public :: Xmode
      real*8, public :: Pi
      complex*16, public :: Zi

c     ## setParameter ##
      real*8, public :: U, J, kT, Dne
c     ## setConfig ##
      integer, public :: Nkx, Nky, Nkz, Nksize, Nband !Nband2, Nrs, Nb
      integer, public :: Nqx, Nqy, Nqz, Nredx, Nredy, Nredz !** Nredxとは?
      real*8, public :: EF
      real*8, public :: la, lb, lc
      real*8, public :: recipLat(1:3,1:3)
c     ## makeMesh ##
      integer, public :: Nkpt
      real*8, allocatable, public :: dk1(:), dk2(:), dk3(:)
c     ## setIteration calcSelfConsistent ##
      integer, public :: maxIter
      real*8, public  :: conv, conrs
c     ## setEnergyRange ##
      integer, public :: maxOmega, minOmega
      real*8, public  :: Erange, Erng, initOmega, eta !!Erange, Erngの違い
c     ## loadHopping ##
      real*8, allocatable, public :: t(:,:,:), Eorbit(:)
      integer, public :: Nsite
      integer, allocatable, public :: Isitex(:), Isitey(:), Isitez(:)

c     ## calcChemicalPotential ##
      real*8, public :: Dne1, Dmu
      real*8, allocatable, public :: Eall(:,:,:,:)
      real*8, allocatable, public :: Dnuu(:), Dndd(:)
      real*8, allocatable, public :: Dens(:,:)
      complex*16, allocatable, public :: Zpsiall(:,:,:,:,:)
c      complex*16, allocatable, public :: Zwfprod(:,:,:,:,:,:)
      complex*16, allocatable, public :: Zop(:,:,:)
      complex*16, allocatable, public :: Zdens(:,:)
      complex*16, allocatable, public :: Zopnew(:,:,:)
      complex*16, allocatable, public :: Zop0(:,:,:)
      complex*16, allocatable, public :: Zopnew0(:,:,:)
      complex*16, allocatable, public :: Zop2(:,:,:)
      complex*16, allocatable, public :: Zopnew2(:,:,:)

c     ## BandPlot ##
      integer, public :: Nkmesh, Nkpath
      real*8, allocatable, public :: dkfrac(:,:)
      real*8, allocatable, public :: Eband(:,:,:)
      complex*16, allocatable, public :: Zpsiband(:,:,:,:)
      type kpathType
         character*16 :: iniName, finName
         real*8 :: iniPosition(1:3)
         real*8 :: finPosition(1:3)
      end type kpathType
      type(kpathType),allocatable :: kpath(:)

c     ## CalcSuscept ##
      real*8, allocatable, public :: chi0(:)


c      integer, public :: Muselect(4)
c      character(len=4), public :: Ledgemode
c      character(len=1), public :: Ledgeqdirection
c      character(len=3), public :: ModeJorbang

c     ## Input/Output Files ##

      end

c## order parameters ##
c     Dnuu: density of electrons with up spin
c     Dndd: density of electrons with dn spin
c     Dnewnuu: new Dnuu
c     Dnewndd: new Dndd
c
c## hamiltonian ##
c     Dhm: total hamiltonian
c     Dhmu: hamiltonian for up spin
c     Dhmd: hamiltonian for dn spin
c     E: all eigen values
c     Eu: eigen values for up spin
c     Ed: eigen values for dn spin
c     Dpsi: all eigen functions
c
c##

c########################################################
