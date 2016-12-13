c###########################################################
c##   showStartMenu:
c###########################################################

      subroutine showMenu()

      use common, only : Xmode

      write(*,*)
      if (Xmode == 'manual') then
         write(*,*) '             -- MAIN MENU -- '
         write(*,*) '           ENTER THE CHARACTER'
         write(*,*)
         write(*,*) 'SETUP: (P)aram  (C)onfig  (I)teration'
         write(*,*) '     : (1)st-config  (D)isplay'
         write(*,*) '     : (E)nergy-range (#)system size'
         write(*,*) '     : (2)phase factor'
         write(*,*)
         write(*,*) ' FILE: (L)oad   (W)rite'
         write(*,*)
         write(*,*) ' CALC: (M)F  (S)uscept  inter(B)and'
         write(*,*) '     : (F)ermi (V)elocity'
         write(*,*) '     : (O)rbital'
         write(*,*) '     : (Q)uit'
      else
         write(*,*) '             -- AUTO MODE -- '
      end if

      return
      end

c###########################################################

c###########################################################
c##   toLowerCase:
c###########################################################

      subroutine toLowerCase()

      use common, only: req
      implicit none

      if (req == 'Q') then
         req='q'            !quit
      else if (req == 'B') then
         req='b'            !bandplot
      else if (req == 'C') then
         req='c'            !config
      else if (req == 'D') then
         req='d'            !display
      else if (req == 'E') then
         req='e'            !energy
      else if (req == 'F') then
         req='f'            !fermi
      else if (req == 'G') then
         req='g'            !green
      else if (req == 'H') then
         req='h'            !hopping integral
      else if (req == 'I') then
         req='i'            !iteration
      else if (req == 'K') then
         req = 'k'          !kougaku dendoudo
      else if (req == 'L') then
         req='l'            !loard
      else if (req == 'M') then
         req='m'            !meanfield
      else if (req == 'O') then
         req='o'            !orbital
      else if (req == 'P') then
         req='p'            !parameters
      else if (req == 'R') then
         req='r'            !resistivity
      else if (req == 'S') then
         req='s'            !suscept
      else if (req == 'V') then
         req='v'            !meanfield
      else if (req == 'W') then
         req='w'            !write
      else if (req == 'X') then
         req='x'            !RIXS L-edge
      end if

      return
      end

c###########################################################

c#########################################################
c#### getFilename:
c#########################################################

      subroutine getFilename(name,val)

      implicit none

      character(LEN=*) :: name
      character(LEN=*), intent(inout) :: val
      character(LEN=1), parameter :: endChara = '?'
      integer :: nameLength

      nameLength = INDEX(name,endChara)

      write(*,*) 'Enter the value:'
      write(*,*) name(1:nameLength),val
      read (5,*,err=100) val
      write(*,*) '==> Set to ',val
      write(*,*)

 100  return
      end
c####################################################################
