c#########################################################
c#### Ivrtx : label vertex
c#########################################################
      integer function Ivrtx(mu1,mu2,nQ,LorR)

      use common, only : Nband, Nqx
      implicit none

      character, intent(in) :: LorR*1
      integer, intent(in) :: mu1, mu2, nQ

c==== .  LorR = 'L':  the left side vertex of the bubble
c==== .  LorR = 'R':  the right side vertex of the bubble
c==== .
c==== .           (line 1: in, line 2: out)
c==== .
c==== .        L-vertex (1,2,Q)      |   R-vertex (1,2,Q)
c==== .                              |
c==== .              2 (k)+Q+(q) =>  |  => (k)+(q)+Q 1
c==== .   (k)+Q =>  <                |                > => (k)+Q
c==== .              1  (q) <=       |      <= (q)   2
c==== .                              |
c==== .
c==== . Label vertices so that Ivrtx(i,j,Q,'L') = Ivrtx(j,i,Q,'R')
c==== .                  j->-j
c==== .   i.e.,    Q => <     > => Q    is diagonal
c==== .                  i-<-i
c==== .
c==== .

      if (LorR == 'L') then
         Ivrtx = mu1 + (mu2-1) * Nband + nQ * Nband * Nband
      else if (LorR == 'R') then
         Ivrtx = mu2 + (mu1-1) * Nband + nQ * Nband * Nband
      else
         write(*,*) ' -- ERROR in Ivrtx -- '
         stop
      end if

      if ((mu1 > Nband).or.(mu2 > Nband).or.(nQ > Nqx-1)) then
         write(*,*) ' -- ERROR in Ivrtx -- '
         stop
      end if


      return
      end
c#########################################################

c#########################################################
c#### Ivrtxinv : vertex
c#########################################################
      integer function Ivrtxinv(ivrtx,n,LorR)

      use common, only : Nband, Nqx
      implicit none

      character, intent(in) :: LorR*1
      integer, intent(in) :: ivrtx, n
      integer :: mu1, mu2, nQ, iv
      integer :: n1, n2, n3, n4

c==== .
c==== .           (line 1: in, line 2: out)
c==== .
c==== .        L-vertex (1,2,Q)      |   R-vertex (1,2,Q)
c==== .                              |
c==== .              2 (k)+Q+(q) =>  |  => (k)+(q)+Q 1
c==== .   (k)+Q =>  <                |                > => (k)+Q
c==== .              1  (q) <=       |      <= (q)   2
c==== .                              |
c==== .
c==== . Vertices labeled so that Ivrtx(i,j,Q,'L') = Ivrtx(j,i,Q,'R')
c==== .
c==== .  n = 1 :  line coming in = 1st index of Ivrtx(1,2,Q,LorR)
c==== .  n = 2 :  line going out = 2nd index of Ivrtx(1,2,Q,LorR)
c==== .  n = 3 :  For L, Q = k2-k1 (k1 = q,  k2 = k+Q+q)
c==== .           For R, Q = k1-k2 (k1 = k+Q+q, k2 = q)

      iv = ivrtx - 1

      if (LorR == 'L') then
         n1 = 1 ; n2 = 2 ; n3 = 3
      else if (LorR == 'R') then
         n1 = 2 ; n2 = 1 ; n3 = 3
      else
         write(*,*) ' -- ERROR in Ivrtx -- '
         stop
      end if

      Ivrtxinv = iv / Nband / Nband   !Q
      if ((Nqx == 1).and.(Ivrtxinv /= 0)) then
         write(*,*) ' -- Error in Ivrtxinv --'
         write(*,*) ' Q need to be 0 for Nqx=1 '
      end if
      if (n == n3) return

      nQ = Ivrtxinv
      Ivrtxinv = (iv - nQ * Nband**2) / Nband + 1 !line out
      if (n == n2) return

      mu2 = Ivrtxinv
      Ivrtxinv =  iv - (mu2-1) * Nband
     &         - nQ * Nband * Nband + 1 !line in
      if (n == n1) return


      write(*,*) ' -- ERROR in Ivrtx -- '
      stop

      return
      end
c#########################################################


c#########################################################
c#### Ispin :
c#########################################################
      integer function Ispinpair(is,EorH)

      use common, only : Nband
      implicit none

      character, intent(in) :: EorH*1
      integer, intent(in) :: is

c==== .        isE -->-- isE
c==== .       <             >
c==== .        isH --<-- isH
c==== .
c==== . A pair of spin (isE, isH) for an electron-hole pair.
c==== . isE = Ispinpair(is,'E'),  isH = Ispinpair(is,'H')
c==== .
c==== .         |  isE = 1         isE = 2
c==== . --------+--------------------------
c==== . isH = 1 |   is = 1          is = 4
c==== .         |
c==== . isH = 2 |   is = 3          is = 2
c==== .

      Ispinpair = 0
      if (EorH == 'E') then
         if(is == 1)then
            Ispinpair = 1
         else if(is == 2)then
            Ispinpair = 2
         else if(is == 3)then
            ispinpair = 1
         else if(is == 4)then
            ispinpair = 2
         else
            write(*,*)"error in Ispinpair: 'is' should be 1 2 3 or 4"
            write(*,*)"                  : is = ",is
            stop
         end if
      else if (EorH == 'H') then
         if(is == 1)then
            Ispinpair = 1
         else if(is == 2)then
            Ispinpair = 2
         else if(is == 3)then
            Ispinpair = 2
         else if(is == 4)then
            Ispinpair = 1
         else
            write(*,*)"error in Ispinpair: 'is' should be 1 2 3 or 4"
            write(*,*)"                  : is = ",is
            stop
         end if
      else
         write(*,*)"error: 'n' should be 1 or 2:"
         stop
      end if

      return
      end
c#########################################################
