c###########################################################
c##   showStartMenu:
c###########################################################

      subroutine showMenu()

      use common, only : Xmode

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

      use common, only: request
      implicit none

      if (request == 'Q') then
         request='q'            !quit
      else if (request == 'B') then
         request='b'            !interband excitation
      else if (request == 'C') then
         request='c'            !config
      else if (request == 'D') then
         request='d'            !display
      else if (request == 'E') then
         request='e'            !energy
      else if (request == 'F') then
         request='f'            !fermi
      else if (request == 'G') then
         request='g'            !green
      else if (request == 'H') then
         request='h'            !hopping integral
      else if (request == 'I') then
         request='i'            !iteration
      else if (request == 'K') then
         request = 'k'          !kougaku dendoudo
      else if (request == 'L') then
         request='l'            !loard
      else if (request == 'M') then
         request='m'            !meanfield
      else if (request == 'O') then
         request='o'            !orbital
      else if (request == 'P') then
         request='p'            !parameters
      else if (request == 'R') then
         request='r'            !resistivity
      else if (request == 'S') then
         request='s'            !suscept
      else if (request == 'V') then
         request='v'            !meanfield
      else if (request == 'W') then
         request='w'            !write
      else if (request == 'X') then
         request='x'            !RIXS L-edge
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
