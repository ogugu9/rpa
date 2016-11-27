c###################################################
c#### common statements
c###################################################

      module common

      implicit none

      ! Main
      character(len=1) :: request
      character(len=20), dimension(10), public :: fnameHop, fnameOut, fnameInit
      character(len=1), public :: Check
      character(len=6), public :: Xmode
      real*8, public :: Pi
      complex*16, public :: Zi

c     ## setParameter ##
      real*8, public :: U, J, kT, Dne
c     ## setConfig ##
      integer, public :: Nkx, Nky, Nband !Nband2, Nrs, Nb
      integer, public :: Nqx, Nqy, Nredx, Nredy !** Nredxとはなにか
c      ## setIteration calcSelfConsistent ##
      integer, public :: maxIter
      real*8, public  :: conv, conrs
c     ## setEnergyRange ##
      integer, public :: maxOmega, minOmega
      real*8, public  :: Erange, Erng, initOmega, eta !!Erange, Erngの違い

c     ## loadHopping ##
      real*8, public :: t(-2:2,-2:2,Nband,Nband), Eorbit(Nband)
      integer, public :: Nsite
      integer, public :: Isitex(0:Nsite-1), Isitey(0:Nsite-1)

c     ## calcChemicalPotential ##
      real*8, public :: Dne1, Dmu, Dne
      real*8, allocatable, public :: Eall(:,:,:,:)
      real*8, allocatable, public :: Dnuu(:), Dndd(:)
      real*8, allocatable, public :: Dens(:,:)
      complex*16, allocatable, public :: Zpsiall(:,:,:,:,:)
      complex*16, allocatable, public :: Zwfprod(:,:,:,:,:,:)
      complex*16, public :: Zop(Nband,Nband,2), Zdens(Nband,Nband)
      complex*16, public :: Zopnew(Nband,Nband,2)
      complex*16, public :: Zop0(Nband,Nband,2), Zopnew0(Nband,Nband,2)
      complex*16, public :: Zop2(Nband,Nband,2), Zopnew2(Nband,Nband,2)

c      integer, public :: Muselect(4)
c      character(len=4), public :: Ledgemode
c      character(len=1), public :: Ledgeqdirection
c      character(len=3), public :: ModeJorbang

c     ## Input/Output Files ##
      fnameHop = 'hop.in'
      fnameInit = 'init.out'
      fnameOut = 'lin90.out'

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
