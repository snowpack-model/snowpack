! ---------------------------------------------------------------------
! PREVAH (Precipitation-Runoff-Evapotranspiration HRU model)
! Subroutine of Model core for WINDOWS, DOS, LINUX and SUN
! Authors: M.Zappa 2004
! ---------------------------------------------------------------------
      SUBROUTINE CHARNR (icsb,NR)

!     This subroutine transforms numbers between 0 and 99 (ics, icsb) into a character
!     string (NR)
!--------------------------------------------------------------------------


! ---------------------------------------------------------------------
!     Common Variables
! ---------------------------------------------------------------------      
    
      INTEGER ics,icsb
      CHARACTER NR*2 

! ---------------------------------------------------------------------
!     Begin
! ---------------------------------------------------------------------
   
      ics=icsb
      
4     CONTINUE      
      
      IF (ics .GE. 100) THEN
          ics=ics-100 
          GOTO 4
      ENDIF
           
          write (NR,'(I2.2)') ics
           
! --------------------------------------------------------------------           
      RETURN     
      END

! --------------------------------------------------------------------           
      SUBROUTINE ABEXT (IGST,FILEN,JAHR)

      INTEGER JAHR,IGST,I,POS
        CHARACTER*2 JZ

      CHARACTER*12 NAMED,FILEN

      IF (IGST .EQ. 1) THEN

         I=1

         NAMED=FILEN

         POS=INDEX(NAMED,'.')

         IF (POS .EQ. 0 ) POS=9

         IF (POS .GT. 9 ) POS=LEN (NAMED)

         POS=POS-1

           CALL CHARNR (JAHR,JZ)
         FILEN=NAMED(1:POS) // '.S'//JZ         

      END IF


!194   CONTINUE

      IF (IGST .EQ. 3) THEN

         NAMED=FILEN

         POS=INDEX(NAMED,'.')

         IF (POS .EQ. 0 ) POS=9

         IF (POS .GT. 9 ) POS=LEN (NAMED)

         POS=POS-1

           CALL CHARNR (JAHR,JZ)
         FILEN=NAMED(1:POS) // '.H'//JZ  

      END IF

      RETURN

      END

! --------------------------------------------------------------------           
       SUBROUTINE WRITEGRIDS (write2d,dom2d,ncol,nrow,uwas,
     * rcol,rrow,rnodat,rgridsize,rxll,ryll)


       real rcol,rrow,rnodat,rgridsize,rxll,ryll
       real wstats(6)
       integer ncol,nrow
       integer uwas,xcoor,ycoor,rloop
       real write2D(ncol,nrow),dom2D(ncol,nrow)

! ------------------------------------------------------------------------------
!     Writing the grids
! ------------------------------------------------------------------------------

c computes and writes wasim headers and writes grids

        wstats(1)=0
        wstats(2)=1e8
        wstats(3)=-1e8
        wstats(4)=0
        wstats(5)=0
        wstats(6)=0

        do ycoor=1,nrow
        do xcoor=1,ncol
        if (dom2D(xcoor,ycoor) .ne. rnodat) then
                if (int(write2D(xcoor,ycoor)) .eq. int(rnodat)) cycle
                wstats(1)=wstats(1)+1.
                wstats(2)=min(wstats(2),write2D(xcoor,ycoor))
                wstats(3)=max(wstats(3),write2D(xcoor,ycoor))
                wstats(4)=wstats(4)+write2D(xcoor,ycoor)
          else
            write2D(xcoor,ycoor)=rnodat
          endif
        end do
        end do

        wstats(5)=wstats(4)/wstats(1)

        do ycoor=1,nrow
        do xcoor=1,ncol
        if (write2D(xcoor,ycoor) .ne. rnodat) then
                wstats(6)=wstats(6)+(write2D(xcoor,ycoor)-wstats(5))**2
          endif
        end do
        end do

        if (wstats(1) .gt. 1.) then  
                wstats(6)=(wstats(6)/(wstats(1)-1.))**(0.5)
        else
                wstats(6)=0.0
        end if

        write (uwas) rcol
        write (uwas) rrow
        write (uwas) rxll
        write (uwas) ryll
        write (uwas) rgridsize
        write (uwas) rnodat

        do rloop=1,6
            write (uwas) wstats(rloop)
        end do

        do ycoor=nrow,1,-1
                write (uwas) (write2D(xcoor,ycoor),xcoor=1,ncol)
        end do

        return
        end


* DISPLAY A MESSAGE (1)
********************************************************************************
* (1) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        SUBROUTINE MESSAGEOUT(CONTENT)

        CHARACTER*(*) CONTENT

      OPEN (999,FILE='RMODULE.LOG',POSITION="APPEND")
        WRITE (999,'(A)') TRIM(CONTENT)
        print* ,TRIM(CONTENT)
        CLOSE (999)

        END SUBROUTINE

* DISPLAY A MESSAGE (1)
********************************************************************************
* (1) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        SUBROUTINE MESSAGERUN(OPT,CONTENT)

        CHARACTER*(*) CONTENT
        INTEGER OPT

        IF (OPT .ne. 3) THEN
        print* ,TRIM(CONTENT)
      OPEN (999,FILE='RMODULE.LOG',POSITION="APPEND")
        WRITE (999,'(A)') TRIM(CONTENT)
        CLOSE (999)
        END IF
        END SUBROUTINE 
           