      SUBROUTINE ISFILE (DATNAM,CHEK)
      
      CHARACTER*256 DATNAM 
      INTEGER CHEK
      
c     Testet ob einen File vorhanden ist.

      CHEK=0
      
      OPEN (777,FILE=DATNAM,STATUS='OLD',ERR=700) 
      CLOSE (777)
                         
      CHEK=1
      
      RETURN 
      
700   CONTINUE

      CHEK=2

      RETURN
      
      END