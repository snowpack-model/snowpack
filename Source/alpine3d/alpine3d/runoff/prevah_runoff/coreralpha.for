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
! Eidgen�ssische Forschungsanstalt WSL
! Abteilung Wasser-, Erd- und Felsbewegungen
! Z�rcherstr. 111
! CH-8903 Birmensdorf
!
! Tel. ++41 1 739 24 33     Fax ++ 41 1 739 24 88
! http://www.wsl.ch/
!
! ------------------------------------------------------------------------------
!
!      SUBROUTINE RUNOFFALPHA (assspack)
        SUBROUTINE corerunoff (anx,any,adate,rfield)
! ------------------------------------------------------------------------------
!
! Features:
! Spatially distributed runoff-generation module
!
! Credits:  H. Jensen (HBV3ETH)  ETH Zuerich
!           J. Gurtz ETH Zuerich
! ------------------------------------------------------------------------------

      use RUNGLOBALS
        character*1 docali
      character*18 assspack
        INTEGER anx,any,adate,aseq
        character*10 cdate
        DOUBLE PRECISION rfield(anx*any),dblein

! ------------------------------------------------------------------------------

!      print* ,'Running runoff module....'

        docali='N'


! ------------------------------------------------------------------------------
!     Assimilating the output of snowpack
! ------------------------------------------------------------------------------

        write(cdate,'(i10.10)') adate
        read(cdate,'(I4,i2,i2,i2)') JAHR,MON,TAG,ASSHH
      aseq=0
      do ycoor=1,nrow
      do xcoor=1,ncol
           aseq=aseq+1
           dblein=rfield(aseq)
         spackin(xcoor,ycoor)=sngl(dblein)
      enddo
      enddo
!1114  continue

! ------------------------------------------------------------------------------
!     Loading the storages from a prervious run
! ------------------------------------------------------------------------------

      DATUMCHARO=DATUMCHAR
      WRITE (DATUMCHAR,'(2I2.2,I4.4)') TAG,MON,JAHR

      if (asshh .eq. FIRSTHH) then
      if ((loadid .eq. 'y') .or. (loadid .eq. 'Y')) THEN
           if (DATUMCHAR .eq. LOADCHAR) THEN
               loadnam=dirstate(1:pathln(4))//zn(1:5)//
     *         'state'//datumcharo//'.sav'
               include 'readstate.inc'
                 loadid='n'
             else
                 stop 'Model State File not Found'
           end if
      end if
        end if

! ------------------------------------------------------------------------------
!     Computing runoff generation
! ------------------------------------------------------------------------------

      do ycoor=1,nrow
      do xcoor=1,ncol
      if (dom2D(xcoor,ycoor) .eq. rnodat) cycle
! ------------------------------------------------------------------------------
!    Soil moisture and runoff-generation module
!
!    Literature
!    Bergstr�m, S. (1976): Development and Application of a Conceptual Runoff Model
!     for Scandinavian Catchments. Bulletin Series A, No. 52, University of Lund.
!    Bergstr�m, S. (1992): The HBV model - its structure and applications. SMHI, RH
!     No. 4, Norrk�ping, Sd.
!    Gurtz, J., Zappa, M., Jasper, K., Lang, H., Verbunt, M., Badoux, A. and
!     Vitvar, T. (2003): A Comparative Study in Modelling Runoff and its Components
!     in Two Mountainous Catchments. Hydrol. Processes, 17, 297-311.
!    Schwarze, R., Droege, W. and Opherden, K. (1999): Regional analysis and modelling
!     of groundwater runoff components from catchments in hard rock areas.
!     In: Diekkr�ger B, Kirkby MJ, Schr�der U (Editors): Regionalisation in hydrology,
!     IAHS Publication No. 254, 221-232.
!    Zappa, M. and Gurtz, J. (2003): Plot-scale simulation of soil moisture and
!     evapotranspiration during the 1999 MAP-Riviera campaign.
!     Hydrology and Earth System Sciences. Accepted.
! ------------------------------------------------------------------------------

      SUZ=SUZTF(xcoor,ycoor)
      SLZ=SLZTF(xcoor,ycoor)
      SLZ1=SLZ1TF(xcoor,ycoor)
      SLZ2=SLZ2TF(xcoor,ycoor)
      SLZ3=SLZ3TF(xcoor,ycoor)
      PRL=SPACKIN(xcoor,ycoor)

      CALL RUNOFF

      IF (r0 .LT. 0.00001) r0=0.0
      IF (r1 .LT. 0.00001) r1=0.0
      IF (r2 .LT. 0.00001) r1=0.0
      IF (r3 .LT. 0.00001) r3=0.0

      spr0(xcoor,ycoor)=r0
      spr1(xcoor,ycoor)=r1
      spr2(xcoor,ycoor)=r2
      sprgq(xcoor,ycoor)=rg1s
      sprgs(xcoor,ycoor)=rg2s
      sprg3(xcoor,ycoor)=rg3s
      sprges(xcoor,ycoor)=rges
      sprperc(xcoor,ycoor)=addperc
      SUZTF(xcoor,ycoor)=SUZ
      SLZTF(xcoor,ycoor)=SLZ
      SLZ1TF(xcoor,ycoor)=SLZ1
      SLZ2TF(xcoor,ycoor)=SLZ2
      SLZ3TF(xcoor,ycoor)=SLZ3

      END DO
      END DO

! ------------------------------------------------------------------------------
!     Writing the grids
! ------------------------------------------------------------------------------

      IF (GRDAKT .EQ. 'y') GRDAKT='Y'

        IF ((docali .eq. 'y') .or. (docali .eq. 'Y')) GRDAKT='N'

      IF (GRDAKT .EQ. 'Y') THEN

      DO I=1,14
        IF (whgrid(I)(1:1) .EQ. 'N') CYCLE
        IF (whgrid(I)(1:1) .EQ. 'n') CYCLE

      write (assspack,'(A4,3I2.2,A4)') labelgrid(I),MON,
     *       TAG,ASSHH,'.bin'
!      write (assspack,'(A4,3I2.2,A4)') labelgrid(I),cdate,'.bin'
      assnam=dirgrids(1:pathln(3))//assspack

      UNIT=731
      OPEN (UNIT, file=ASSNAM, status='unknown',form='unformatted')

      DO CLNAM=1,256
         ASSNAM(CLNAM:CLNAM)=' '
      END DO

      do ycoor=1,nrow
      do xcoor=1,ncol
      if (dom2D(xcoor,ycoor) .eq. rnodat) cycle
      if (I .eq. 1) exportgrid(xcoor,ycoor)=spackin(xcoor,ycoor)
      if (I .eq. 2) exportgrid(xcoor,ycoor)=sprges(xcoor,ycoor)
      if (I .eq. 3) exportgrid(xcoor,ycoor)=spr0(xcoor,ycoor)
      if (I .eq. 4) exportgrid(xcoor,ycoor)=spr1(xcoor,ycoor)
      if (I .eq. 5) exportgrid(xcoor,ycoor)=sprgq(xcoor,ycoor)
      if (I .eq. 6) exportgrid(xcoor,ycoor)=sprgs(xcoor,ycoor)
      if (I .eq. 7) exportgrid(xcoor,ycoor)=sprg3(xcoor,ycoor)
      if (I .eq. 8) exportgrid(xcoor,ycoor)=spr2(xcoor,ycoor)
      if (I .eq. 9) exportgrid(xcoor,ycoor)=sprperc(xcoor,ycoor)
      if (I .eq. 10) exportgrid(xcoor,ycoor)=SUZTF(xcoor,ycoor)
      if (I .eq. 11) exportgrid(xcoor,ycoor)=SLZTF(xcoor,ycoor)
      if (I .eq. 12) exportgrid(xcoor,ycoor)=SLZ1TF(xcoor,ycoor)
      if (I .eq. 13) exportgrid(xcoor,ycoor)=SLZ2TF(xcoor,ycoor)
      if (I .eq. 14) exportgrid(xcoor,ycoor)=SLZ3TF(xcoor,ycoor)
      END DO
      END DO

       call WRITEGRIDS (exportgrid,dom2d,ncol,nrow,unit,
     * rcol,rrow,rnodat,rgridsize,rxll,ryll)

       close (unit)

       END DO
       END IF

! ------------------------------------------------------------------------------
!     Writing tables
! ------------------------------------------------------------------------------

      DO I=1,14
        allmeans(I)=-999.

      do ycoor=1,nrow
      do xcoor=1,ncol
      if (dom2D(xcoor,ycoor) .eq. rnodat) cycle
      if (I .eq. 1) exportgrid(xcoor,ycoor)=spackin(xcoor,ycoor)
      if (I .eq. 2) exportgrid(xcoor,ycoor)=sprges(xcoor,ycoor)
      if (I .eq. 3) exportgrid(xcoor,ycoor)=spr0(xcoor,ycoor)
      if (I .eq. 4) exportgrid(xcoor,ycoor)=spr1(xcoor,ycoor)
      if (I .eq. 5) exportgrid(xcoor,ycoor)=sprgq(xcoor,ycoor)
      if (I .eq. 6) exportgrid(xcoor,ycoor)=sprgs(xcoor,ycoor)
      if (I .eq. 7) exportgrid(xcoor,ycoor)=sprg3(xcoor,ycoor)
      if (I .eq. 8) exportgrid(xcoor,ycoor)=spr2(xcoor,ycoor)
      if (I .eq. 9) exportgrid(xcoor,ycoor)=sprperc(xcoor,ycoor)
      if (I .eq. 10) exportgrid(xcoor,ycoor)=SUZTF(xcoor,ycoor)
      if (I .eq. 11) exportgrid(xcoor,ycoor)=SLZTF(xcoor,ycoor)
      if (I .eq. 12) exportgrid(xcoor,ycoor)=SLZ1TF(xcoor,ycoor)
      if (I .eq. 13) exportgrid(xcoor,ycoor)=SLZ2TF(xcoor,ycoor)
      if (I .eq. 14) exportgrid(xcoor,ycoor)=SLZ3TF(xcoor,ycoor)
      END DO
      END DO

      gebmean=0.0

      do ycoor=1,nrow
      do xcoor=1,ncol
         if (dom2D(xcoor,ycoor) .eq. rnodat) cycle
         J=INT(mzones(xcoor,ycoor))
         zones(J,3)=zones(J,3)+exportgrid(xcoor,ycoor)/zones(J,2)
         gebmean=gebmean+exportgrid(xcoor,ycoor)/nvalids
      enddo
      enddo
      allmeans(I)=gebmean

      IF ((whgrid(I)(9:9) .EQ. 'Y') .or.
     *   (whgrid(I)(9:9) .EQ. 'y')) THEN
        UNIT=1400+I

      write (unit,'(I4,3(1X,I2),999(1X,F9.3))')
     *        JAHR,MON,TAG,ASSHH,
     *       (zones(J,3),J=1,ntg),gebmean
      END IF

!      print*,'Now loop over',ntg
      do j=1,ntg
        zones(J,3)=0.0
      end do

      end do

      qlps=allmeans(2)*aelhg*1000*1000/60/idmin

      IF ((ZIERG(1:1) .EQ. 'y') .OR. (ZIERG(3:3) .EQ. 'Y'))
     *       write (MDQ,'(I4,3(1X,I2),1X,F12.2,99(1X,F9.5))')
     *       JAHR,MON,TAG,ASSHH,qlps,(allmeans(J),J=2,8)

      IF ((ZIERG(3:3) .EQ. 'y') .OR. (ZIERG(3:3) .EQ. 'Y'))
     *       write (MDSTD,'(I4,3(1X,I2),9(1X,F9.5),5(1X,F9.3))')
     *       JAHR,MON,TAG,ASSHH,(allmeans(J),J=1,14)

! ------------------------------------------------------------------------------
!     Saving the model state
! ------------------------------------------------------------------------------
      if (asshh .eq. (firsthh+23)) then
      if ((saveid .eq. 'y') .or. (saveid .eq. 'Y')) THEN

      IF (TAG .EQ. MONUM(MON)) THEN
      if ((savemmid .eq. 'y') .or. (savemmid .eq. 'Y')) THEN
         savenam=dirstate(1:pathln(4))//zn(1:5)//
     *   'state'//datumchar//'.sav'
         include 'writestate.inc'
      end if
      end if

      if ((savemoreid .eq. 'y') .or. (savemoreid .eq. 'Y')) THEN
      DO I=1,SAVEMORENN
      if (DATUMCHAR .eq. SAVECHAR(i)) THEN
         savenam=dirstate(1:pathln(4))//zn(1:5)//
     *          'state'//datumchar//'.sav'
         include 'writestate.inc'
      end if
      END DO
      end if

      end if
        end if

!999   CONTINUE
      RETURN
      END

! ------------------------------------------------------------------------------
