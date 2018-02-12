! ------------------------------------------------------------------------------
!
! RUNOFFALPHA
!
! Authors: M. Zappa et al. 2003-2004
! ------------------------------------------------------------------------------
!
! Reviewed literature in international journals
!
! More general documentation
!
! ------------------------------------------------------------------------------
!
! Contacts
! Massimiliano Zappa, Dr. sc. nat. ETH         - massimiliano.zappa@wsl.ch
! Eidgen�sische Forschungsanstalt WSL
! Abteilung Wasser-, Erd- und Felsbewegungen
! Zrcherstr. 111
! CH-8903 Birmensdorf
!
! Tel. ++41 1 739 24 33     Fax ++ 41 1 739 24 88
! http://www.wsl.ch/
!
! ------------------------------------------------------------------------------
!
!      SUBROUTINE INITRUNOFF (NAMEIN,efficnam,ernam,gken,
!     *  diroutput,docali)
      SUBROUTINE initrunoff (anx,any,NAMEIN)

! ------------------------------------------------------------------------------
!
! Features:
! Spatially distributed runoff-generation module
!
! Credits:  H. Jensen (HBV3ETH)  ETH Zuerich
!           J. Gurtz ETH Zuerich
! ------------------------------------------------------------------------------

        use RUNGLOBALS
        character*256 diroutput,ernam,namein
        character*1 docali
        character*4 gken
        integer any,anx,initl,rec
! ------------------------------------------------------------------------------

!     parameters to be initialized are declared in a include file
      include 'init.inc'

      docali='N'
        ERNAM(1:11)='nocalib.out'
! ------------------------------------------------------------------------------
! Reads the control parameters
! ------------------------------------------------------------------------------
!      print* ,anx,any,trim(NAMEIN)
      include 'readalpine.inc'

        IDMIN=int(dtmin(3))

        beg=0


!     *****************************************
!     Read Grid Properties
!     *****************************************

      IDNAM=dirinput(1:pathln(1))//
     *     GRIDNAM(1:LNGRID)//'.'//rdgrid(1)(1:extln(1))

!      OPEN (731, file=IDNAM, status='old',form='binary')
      OPEN (731, file=IDNAM, status='old',
     *      form='unformatted', recl=4, access='direct')
      rec=1

      write (*,*) 'opened file ', trim(IDNAM), ' for reading headers'
      DO CLNAM=1,256
         IDNAM(CLNAM:CLNAM)=' '
      END DO

      uwas=731

      include 'rwasim.inc'

      close (731)

!     *********************
!     Allocation Of Arrays
!     *********************

        if (any .ne. nrow) stop 'initrunoff: bad row number'
        if (anx .ne. ncol) stop 'initrunoff: bad col number'

      ALLOCATE(mzones(ncol,nrow))
      ALLOCATE(dom2D(ncol,nrow))
!      ALLOCATE(muse(ncol,nrow))! 1) not needed by runoff and 2) file not read correctly in fortran

!     ************************
!     Read the Landuse
!     ************************


!      UNIT=731

! 1) not needed by runoff and 2) file not read correctly in fortran
!      IDNAM=dirinput(1:pathln(1))//
!     *      GRIDNAM(1:LNGRID)//'.'//rdgrid(3)(1:extln(3))
!      CALL ISFILE (idnam,LN)

!      IF (LN .EQ. 2) STOP

!!      OPEN (UNIT, file=IDNAM, status='old',form='binary')
!      OPEN (731, file=IDNAM, status='old',
!     *      form='unformatted', recl=4, access='direct') ! don't read ascii file as unformatted file!
!      rec=1

!      write (*,*) 'opened landuse file "', trim(IDNAM), '" for reading'
!      DO CLNAM=1,256
!         IDNAM(CLNAM:CLNAM)=' '
!      END DO

!      uwas=731

!      include 'rwasim.inc'
!      include 'rw2dr.inc'

!      close (731)

!      do ycoor=1,nrow
!      do xcoor=1,ncol
!         muse(xcoor,ycoor)=write2D(xcoor,ycoor)
!      enddo
!      enddo

!      DEALLOCATE(write2D)

!     ************************
!     Read the Domain
!     ************************

      UNIT=731

      IDNAM=dirinput(1:pathln(1))//
     *      GRIDNAM(1:LNGRID)//'.'//rdgrid(2)(1:extln(2))
      CALL ISFILE (idnam,LN)

      IF (LN .EQ. 2) STOP

!      OPEN (UNIT, file=IDNAM, status='old',form='binary')
      OPEN (731, file=IDNAM, status='old',
     *      form='unformatted', recl=4, access='direct')
      rec=1

      write (*,*) 'opened domain file "', trim(IDNAM), '" for reading'
      DO CLNAM=1,256
         IDNAM(CLNAM:CLNAM)=' '
      END DO

      uwas=731

      include 'rwasim.inc'
      include 'rw2dr.inc'

        nvalids=wstats(1)

      close (731)

      do ycoor=1,nrow
      do xcoor=1,ncol
         DOM2D(xcoor,ycoor)=write2D(xcoor,ycoor)
      enddo
      enddo

      DEALLOCATE(write2D)

!     ************************
!     Read the Meteorological zones
!     ************************

      UNIT=731

      IDNAM=dirinput(1:pathln(1))//
     *      GRIDNAM(1:LNGRID)//'.'//rdgrid(1)(1:extln(1))
      CALL ISFILE (idnam,LN)

      IF (LN .EQ. 2) STOP

!      OPEN (UNIT, file=IDNAM, status='old',form='binary')
      OPEN (731, file=IDNAM, status='old',
     *      form='unformatted', recl=4, access='direct')
      rec=1

      write (*,*) 'opened meteo file "', trim(IDNAM), '" for reading'
      DO CLNAM=1,256
         IDNAM(CLNAM:CLNAM)=' '
      END DO

      uwas=731

      include 'rwasim.inc'
      include 'rw2dr.inc'

        nvalids=min(nvalids,wstats(1))

      close (731)

      do ycoor=1,nrow
      do xcoor=1,ncol
         mzones(xcoor,ycoor)=write2D(xcoor,ycoor)
      enddo
      enddo

      DEALLOCATE(write2D)

      ALLOCATE(ZONESLONG(INT(wstats(3))))

      DO I=1,INT(wstats(3))
         ZONESLONG(I)=0.0
      ENDDO

      do ycoor=1,nrow
        do xcoor=1,ncol
           if (dom2D(xcoor,ycoor) .eq. rnodat) cycle
           if (mzones(xcoor,ycoor) .eq. rnodat) cycle
           J=INT(mzones(xcoor,ycoor))
           ZONESLONG(J)=ZONESLONG(J)+1
        enddo
      enddo

      ntg=0

      DO I=1,INT(wstats(3))
         IF (ZONESLONG(I) .ne. 0.) ntg=ntg+1
      END DO


      ALLOCATE (ZONES(ntg,3))

      J=0

      DO I=1,INT(wstats(3))
          IF (ZONESLONG(I) .ne. 0.) then
              J=J+1
              ZONES(J,2)=ZONESLONG(I)
              ZONES(J,1)=Float(I)
              ZONES(J,3)=0.0
          end if
      END DO

      DEALLOCATE (zoneslong)

      do ycoor=1,nrow
         do xcoor=1,ncol
           if (dom2D(xcoor,ycoor) .eq. rnodat) cycle
           DO I=1,ntg
             IF (INT(mzones(xcoor,ycoor)) .eq. INT(ZONES(I,1))) THEN
               mzones(xcoor,ycoor)=float(I)
             end if
           end do
         enddo
      enddo

! ------------------------------------------------------------------------------
! Initialization - allocation of the arrays
! ------------------------------------------------------------------------------

      INCLUDE 'allocate.inc'
      INCLUDE 'nullsetz.inc'

!167   CONTINUE

!132   FORMAT (A72)

      AEGES=AELHG

! ------------------------------------------------------------------------------
!     Decides the number of output tables

        CALL POSNAM (NAMAUS,1,POS)
        CALL POSNAM (diroutput,1,POSI3)
        NAMAUS=diroutput(1:posi3)//NAMAUS(1:POS)
        CALL POSNAM (NAMAUS,1,POS)

      NAMDRU=NAMAUS(1:POS) // '.pri'        ! Monthly results and statistics at end of period for the whole catchment
      NAMSTD=NAMAUS(1:POS) // '.std'        ! Hourly values for the whole catchment
      NAMQ=NAMAUS(1:POS) // '.q'        ! Every time step for the whole catchment
      NAMMIT=NAMAUS(1:POS) // '.mit'         ! Daily results and statistics for the whole catchment
      NAMLOG=NAMAUS(1:POS) // '.log'        ! Logfile

        IF (ERNAM(1:11) .eq. 'nocalib.out') THEN
           ERNAM=NAMMIT
        ELSE
           CALL POSNAM (ERNAM,1,POS)
           ERNAM=diroutput(1:posi3)//ERNAM(1:POS)
        END IF

! ------------------------------------------------------------------------------

      DO I=1,14
        IF (whgrid(I)(9:9) .EQ. 'N') CYCLE
        IF (whgrid(I)(9:9) .EQ. 'n') CYCLE
        tabnam=NAMAUS(1:POS)//labelgrid(I)//'.tab'
        UNIT=1400+I
        OPEN (UNIT,FILE=tabnam,STATUS='UNKNOWN',err=999)
        DO CLNAM=1,256
         tabnam(CLNAM:CLNAM)=' '
        END DO

        LN3=999

        write (UNIT,'(A19,1X,A4)') 'Zonal summary file:',
     *         labelgrid(I)
        write (UNIT,'(A13,999(1X,I9))') 'YYYY MM DD HH',
     *        (INT(zones(J,1)),J=1,ntg),LN3
        write (UNIT,'(A13,999(1X,F9.4))') 'YYYY MM DD HH',
     *        ((zones(J,2)/nvalids),J=1,ntg),nvalids/nvalids

      END DO

! ------------------------------------------------------------------------------
!   Output tables are opened
! ------------------------------------------------------------------------------

        OPEN (MLOG,FILE=NAMLOG,STATUS='UNKNOWN',err=999)

      IF ((ZIERG(1:1) .EQ. 'y') .OR. (ZIERG(1:1) .EQ. 'Y')) THEN
         OPEN (MDQ,FILE=NAMQ,STATUS='UNKNOWN',err=999)
        write (MDQ,'(A19)') 'Summary file runoff'
        write (MDQ,'(A13,3X,99(6X,A4))') 'YYYY MM DD HH',
     *         'qlps',(labelgrid(J),J=2,8)
      END IF

      IF ((ZIERG(3:3) .EQ. 'y') .OR. (ZIERG(3:3) .EQ. 'Y')) THEN
       OPEN (MDSTD,FILE=NAMSTD,STATUS='UNKNOWN',err=999)
        write (MDSTD,'(A26)') 'Summary file all variables'
        write (MDSTD,'(A13,99(6X,A4))') 'YYYY MM DD HH',
     *        (labelgrid(J),J=1,14)
      END IF

        IF ((ZIERG(5:5) .EQ. 'y') .OR. (ZIERG(5:5) .EQ. 'Y'))
     * OPEN (MDMIT,FILE=NAMMIT,STATUS='UNKNOWN',err=999)

      IF ((ZIERG(7:7) .EQ. 'y') .OR. (ZIERG(7:7) .EQ. 'Y'))
     * OPEN (MDR,FILE=NAMDRU,STATUS='UNKNOWN',err=999)

! ------------------------------------------------------------------------------
!      Setting initial conditions
! ------------------------------------------------------------------------------

! Determining the Initial content of the runoff storage of the saturated zone SLZ

      IF (SLZ .LT. 0.0) THEN
         SLZ=0.0
         IF (QGA .LT. 0.0) THEN
            QGA=-QGA
! Transforming from l/s into mm/dt
            QGA=QGA*3600/1000/1000/AELHG
         ENDIF
! Filing of the storage accoring to the current water flow in l/s
         IF (K2 .GT. 0.0) SLZ=(QGA/K2)-QGA
      ENDIF

      do ycoor=1,nrow
      do xcoor=1,ncol
         SUZTF(xcoor,ycoor)= SUZ
         SLZ1TF(xcoor,ycoor)=0.0
         SLZ2TF(xcoor,ycoor)=SLZ
         SLZ3TF(xcoor,ycoor)=0.0
         SLZTF(xcoor,ycoor)=SLZ
      enddo
      enddo

! ------------------------------------------------------------------------------
!      Defining Begin DAY-1
! ------------------------------------------------------------------------------
      TAG=TAGA-1
      MON=MONA
      JAHR=JAHRA

      MONUM(2)=28
      IF (JAHR.EQ.(INT(float(JAHR)/4.)*4)) MONUM(2)=29

      IF (TAG .eq. 0) THEN
          MON=MON-1
          IF (MON .eq. 0) THEN
          JAHR=JAHR-1
          MON=12
          TAG=31
          ELSE
          TAG=MONUM(MON)
          END IF
      END IF

      WRITE (DATUMCHAR,'(2I2.2,I4.4)') TAG,MON,JAHR

999   CONTINUE
      print* ,'Initialization of runoff module successful'
      RETURN
      END

! ------------------------------------------------------------------------------
