

! ------------------------------------------------------------------------------
! Reads the new control file
! Implemented by MZ in September 2003
! ------------------------------------------------------------------------------


      MDAT6=333
       nmetvar=1
       initl=1
!
       allocate(outname(9))
       allocate(missing(9))
       allocate(adjust(9))
       allocate(buildnam(3,3,9))
       allocate(method(9))

      MODUL(1)='[Identification]'
      MODUL(2)='[Paths]     '
      MODUL(3)='[Grids]     '
      MODUL(4)='[Date and save]'
      MODUL(5)='[Snowpack]  '
      MODUL(6)='[Tuneable Parameters]'
      MODUL(7)='[Forecast Parameters]'
      MODUL(8)='[Storage content]'
      MODUL(9)='[Model output]'
      MODUL(10)='[Efficiency]'
      MODUL(11)='[Routing]   '
      MODUL(12)='[Begin Extended Options]'
! ------------------------------------------------------------------------------
c     STANDARDWERTE
! ------------------------------------------------------------------------------

      NTG=1

! ------------------------------------------------------------------------------
c     Eingabe der Name der Inputdatei
! ------------------------------------------------------------------------------
      print* ,trim(namein)
      OPEN (MDAT6,FILE=NAMEIN,STATUS='OLD',ERR=197)
! ------------------------------------------------------------------------------
c    Einlesen der Inputdatei
! ------------------------------------------------------------------------------

33    CONTINUE

      LN=-1

      READ (MDAT6,'(A78)',END=199,err=199) kopf

      DO I=1,12
         IF (kopf(1:12) .EQ. MODUL(I)(1:12)) THEN
            LN=I
            print* ,MODUL(I),'OK'
         ENDIF
      END DO

! ------------------------------------------------------------------------------
!     '[Identification]'
! ------------------------------------------------------------------------------

      IF (LN .EQ. 1) THEN
          READ(MDAT6,'(A32)',err=198) ZN
          ZN=TRIM(ADJUSTL(ZN))
              GKEN=ZN(1:3)//'_'
          READ(MDAT6,*,err=198) GEOBREITE
          READ(MDAT6,*,err=198) AELHG
          READ(MDAT6,*,err=198) AEGRID
      ENDIF


! ------------------------------------------------------------------------------
!     '[Paths]'
! ------------------------------------------------------------------------------

      IF (LN .EQ. 2) THEN
          READ (MDAT6,'(A256)',err=198) dirinput
          READ (MDAT6,'(A256)',err=198) diroutput
          READ (MDAT6,'(A256)',err=198) dirgrids
          READ (MDAT6,'(A256)',err=198) dirstate
          READ (MDAT6,'(A256)',err=198) dirmeteo
          call posnam (dirinput,1,pathln(1))
          call posnam (diroutput,1,pathln(2))
          call posnam (dirgrids,1,pathln(3))
          call posnam (dirstate,1,pathln(4))
          call posnam (dirmeteo,1,pathln(5))
      ENDIF

! ------------------------------------------------------------------------------
!     '[Grids]'
! ------------------------------------------------------------------------------

      IF (LN .EQ. 3) THEN
          READ (MDAT6,'(A256)',err=198) gridnam
          READ (MDAT6,'(A78)',END=199,err=199) kopf
          READ (MDAT6,'(A78)',END=199,err=199) kopf
          READ (MDAT6,'(A256)',err=198) rdgrid(1)
          READ (MDAT6,'(A256)',err=198) rdgrid(2)
          READ (MDAT6,'(A256)',err=198) rdgrid(3)
          call posnam (gridnam,1,LNGRID)
          do J=1,3
          call posnam (rdgrid(J),1,extln(J))
          end do
      ENDIF

! ------------------------------------------------------------------------------
!     '[Date and save]'
! ------------------------------------------------------------------------------

      IF (LN .EQ. 4) THEN
      READ(MDAT6,*,err=198) FIRSTHH
      READ(MDAT6,*,err=198) TAGA,MONA,JAHRA,STDA
      READ (MDAT6,'(A1)',err=198) LOADID
      READ(MDAT6,*,err=198) TAGE,MONE,JAHRE,STDE
      WRITE (LOADCHAR,'(2I2.2,I4.4)') TAGA,MONA,JAHRA
            DATUMCHAR='01011901'

      WRITE (CHHHMM,'(I4.4)') STDA
      READ (CHHHMM,'(2I2.2)') CHRHH,CHRMIN
      MINBEG=CHRHH*60+CHRMIN

 !     print* ,MINBEG,CHRHH,CHRMIN

      WRITE (CHHHMM,'(I4.4)') STDE
      READ (CHHHMM,'(2I2.2)') CHRHH,CHRMIN
      MINEND=CHRHH*60+CHRMIN

!      print* ,MINEND,CHRHH,CHRMIN

      READ (MDAT6,'(A1)',err=198) SAVEID
      READ (MDAT6,'(A1)',err=198) SAVEMMID
      READ (MDAT6,'(A1)',err=198) SAVEMOREID
      READ (MDAT6,*,err=198) SAVEMORENN

      allocate (SAVEDD(SAVEMORENN))
      allocate (SAVEMM(SAVEMORENN))
      allocate (SAVEYY(SAVEMORENN))
      allocate (SAVECHAR(SAVEMORENN))

      DO I=1,SAVEMORENN
      READ (MDAT6,*,err=198) SAVEDD(i),SAVEMM(i), SAVEYY(i)
      WRITE (SAVECHAR(i),'(2I2.2,I4.4)') SAVEDD(i),SAVEMM(i), SAVEYY(i)
      END DO
      ENDIF

! ------------------------------------------------------------------------------
!     '[Snowpack]'
! ------------------------------------------------------------------------------

       IF (LN .EQ. 5) THEN
       READ(MDAT6,*,err=198) nmetvar
          READ(MDAT6,'(A1)',err=198) SPCLEAN

       DO I=1,NMETVAR
       READ (MDAT6,103,END=199,err=198) kopf
       print* ,kopf(1:28)
       READ (MDAT6,*,err=198) method(I)
       READ (MDAT6,*,err=198) dtmin(3)
       READ (MDAT6,'(A256)',err=198) outname(I)

      call posnam (outname(I),1,buildnam(3,1,I))

      READ(MDAT6,*,err=198) missing(I)
      READ(MDAT6,*,err=198) adjust(I)
      READ (MDAT6,103,END=199,err=198) kopf
      print* ,kopf(1:28)
      END DO
      END IF

! ------------------------------------------------------------------------------
!     '[Tuneable Parameters]'
! ------------------------------------------------------------------------------

      IF (LN .EQ. 6) THEN
! ------------------------------------------------------------------------------
! Output Name for Grids and Tables
         READ(MDAT6,'(A256)',err=198) NAMAUS
! ------------------------------------------------------------------------------
! Runoff-generation module
         READ(MDAT6,*,err=198) SGRLUZ
         READ(MDAT6,*,err=198) K0H
         READ(MDAT6,*,err=198) K1H
         READ(MDAT6,*,err=198) K2H
         READ(MDAT6,*,err=198) CG1H
         READ(MDAT6,*,err=198) SLZ1MAX
         READ(MDAT6,*,err=198) CPERC

! ------------------------------------------------------------------------------
! Transformation of storage coefficient from hours to exponential factor
      K0=0.0
      IF (K0H .GT. 0.0) K0=1.0-EXP(-1.0/K0H)
      K1=0.0
      IF (K1H .GT. 0.0) K1=1.0-EXP(-1.0/K1H)
      CG2H=K2H
      K2=0.0
      IF (K2H .GT. 0.0) K2=1.0-EXP(-1.0/K2H)
      K3=0.0
      K3H=0.0
      CG3H=0.11111*CG2H                        ! After Schwarze
      ENDIF

! ------------------------------------------------------------------------------
!     '[Forecast Parameters]'
! ------------------------------------------------------------------------------

      IF (LN .EQ. 7) THEN
         READ(MDAT6,*,err=198) ADJSUZ
         READ(MDAT6,*,err=198) ADJSLZ
      ENDIF

! ------------------------------------------------------------------------------
!     '[Storage content]'
! ------------------------------------------------------------------------------

      IF (LN .EQ. 8) THEN
         READ(MDAT6,*,err=198) SUZ
         READ(MDAT6,*,err=198) SLZ
         READ(MDAT6,*,err=198) QGA
      ENDIF

 ! ------------------------------------------------------------------------------
 !     '[Model output]'
 ! ------------------------------------------------------------------------------

      IF (LN .EQ. 9) THEN
         READ(MDAT6,'(A15)',err=198) ZIERG
         READ(MDAT6,'(A1)',err=198) GRDAKT
         DO I=1,14
            READ (MDAT6,'(A9)') WHGRID(I)
         END DO
      ENDIF

 ! ------------------------------------------------------------------------------
 !     '[Efficiency]'
 ! ------------------------------------------------------------------------------
      IF (LN .EQ. 10) THEN
        READ(MDAT6,*,err=198) neffis
        DO I=1,neffis
           READ (MDAT6,'(A256)',err=198) effinams(i)
        END DO
      ENDIF

 ! ------------------------------------------------------------------------------
 !     '[Routing]'
 ! ------------------------------------------------------------------------------

      IF (LN .EQ. 11) THEN
             READ (MDAT6,'(A1)') dorout
         READ(MDAT6,*,err=198) nrouts
         DO I=1,nrouts
            READ (MDAT6,'(A256)',err=198) routnams(i)
            read (MDAT6,*) translation(I)
            read (MDAT6,*) storage(I)
         END DO
       ENDIF

      goto 33


197   CONTINUE

      STOP 'Runoff control file not found'

198   CONTINUE

      IF (LN .GT. 0) THEN
         print* ,'Fehler in Modul ',MODUL(LN)
      ENDIF
      print*
      print* ,'Fehler in der Inputdatei, moegliche Ursachen:'
      print* ,'- REAL statt INTEGER '
      print* ,'- ZAHL statt BUCHSTAB, oder umgekehrt '
      print* ,'- Dateiname laenger als 12 Buchstaben '
      IF (LN .EQ. 9) print* ,'- Falls SSO(1) Negativ, nur SSO(1)',
     *                       ' muss vorhanden sein '
      print*

199   CONTINUE

      CLOSE(MDAT6)

103   FORMAT (A75)


