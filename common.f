c###################################################
c#### common statements
c###################################################

      module common

      use parameters

      implicit none

      character(len=20), dimension(10), public :: fname
      character(len=1), public :: Check
      character(len=6), public :: Xmode
      character(len=4), public :: Ledgemode
      character(len=1), public :: Ledgeqdirection
      character(len=3), public :: ModeJorbang

      integer, public :: Maxit
      integer, public :: Nkx, Nky, Ns, Ns2, Nrs, Nb
      integer, public :: Nqx, Nqy, Nredx, Nredy
      integer, public :: Maxom, Minom
      integer, public :: Isite(-llx:llx,-lly:lly)
      integer, public :: Isitex(0:lln-1), Isitey(0:lln-1)
      integer, public :: Muselect(4)

      real*8, public :: Conv, Dtemp, conrs
      real*8, public :: t(-2:2,-2:2,llorb,llorb)
      real*8, public :: U, DJ
      real*8, public :: Pi
      real*8, public :: Fen, Dne1, Dmu, Dne
      real*8, public :: Erange, Erng, Ominit, Deta
      real*8, allocatable, public :: Eall(:,:,:,:)
      real*8, allocatable, public :: Dnuu(:), Dndd(:), Dnuuxyz(:)
      real*8, allocatable, public :: Dens(:,:), Densxyz(:,:)
      complex*16, public :: Zi
      complex*16, allocatable, public :: Zpsiall(:,:,:,:,:)
      complex*16, allocatable, public :: Zwfprod(:,:,:,:,:,:)
      complex*16, public :: Zop(llorb,llorb,2), Zdens(llorb,llorb)
      complex*16, public :: Zopnew(llorb,llorb,2)
      complex*16, public :: Zop0(llorb,llorb,2), Zopnew0(llorb,llorb,2)
      complex*16, public :: Zop2(llorb,llorb,2), Zopnew2(llorb,llorb,2)
      real*8, public :: Dop(llorb,llorb)
      real*8, public :: Dopnew(llorb,llorb)
      real*8, public :: Eneorb(llorb)
      real*8, public :: Dkshift(2)


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
