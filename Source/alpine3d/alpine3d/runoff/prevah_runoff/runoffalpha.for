! ------------------------------------------------------------------------------
!      SUBROUTINE RUNOFF (PRL,CPERC,SGRLUZ,K0,K1,K3,CG1H,CG2H,CG3H,
!     *            SLZ1MAX,IDMIN,SUZ,SLZ,SLZ1,SLZ2,SLZ3,RG1S,RG2S,
!     *            RG3S,R0,R1,R2,R3,RGES,addperc)
       SUBROUTINE RUNOFF
!      
!      PROGRAM RUNOFF
!    
!    Runoff-generation module for GridAlpHA
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
!
!   
!    M. Zappa - Sept.-Dec. 2003
! ------------------------------------------------------------------------------

      use RUNGLOBALS

      INTEGER IM                                ! Integration steps within the current time step

      REAL Q0,                                        ! Direct runoff in portion of time step
     *     Q1,                                        ! Interflow in portion of time step
     *     Q2,                                        ! Total baseflow in portion of time step
     *     Q3,                                        ! Karst flow in portion of time step (not active) 
     *     RG1,                                        ! Quick response baseflow in portion of time step
     *     RG2,                                        ! Delayed Baseflow  in portion of time step
     *     RG3                                        ! Additional baseflow (Schwarze)  in portion of time step
     
      REAL TS,                                        ! portion of time step (0..1)
     *     TN,                                        ! portion of time step in hours
     *     TH,                                        ! Time step in hours
     *     SUZE,                                ! Storage reservoirs in the saturated zone (from previous and to next time step)
     *     SLZE,SLZ1E,SLZ2E,SLZ3E,                ! Storage reservoirs in the saturated zone (from previous and to next time step)
     *     QGES,                                ! Storage coefficient in hours for karst flow (not active)
     *     DIFF,                                ! Help variable for baseflow
     *     RIM,                                        ! Integration steps within the current time step (as float number)    
     *     perc,                                 ! Help variable for percolation    
     *     dsuz,dmin                                ! Help variable for inflow
     
     
! ------------------------------------------------------------------------------
! Begin
! ------------------------------------------------------------------------------
!    Initialization of the storages and of the variables

      DMIN=FLOAT(IDMIN)

!5     continue



     
      SUZE=SUZ
      SLZE=SLZ
      SLZ1E=SLZ1
      SLZ2E=SLZ2
      SLZ3E=SLZ3
!      print*, '1.5',dmin,idmin,slz1e,slz2e,slz3e,SLZ1,SLZ2,SLZ3
      QGES=0.0
      Q0=0.0
      Q1=0.0
      Q2=0.0
      RGES=0.0
      R0=0.0
      R1=0.0
      R2=0.0
      RG1s=0.0
      RG2s=0.0
      RG3s=0.0
      RG1=0.0
      RG2=0.0
      RG3=0.0
      R3=0.0
      PERC=0.0
      addperc=0.0

! ------------------------------------------------------------------------------
!    PRL determines the number of integration steps (at least 3)

      IM=INT(((PRL*10.0)/30.0)+1.5)
      IM=MAX(IM,3)
      TS=1.0/IM
      TN=TS*DMIN/60
      TH=DMIN/60
      RIM=IM*1.0

!    After Schwarze CG3H is 1/9 of CG2H
      CG3H=0.11111*CG2H

      if (im .gt. 5) print*,'Now loop over',IM,PRL
! ------------------------------------------------------------------------------
!    Core of the runoff-generation module - Integration within IDMIN
      DO I=1,IM

! ------------------------------------------------------------------------------
!     Available water for runoff generation DSUZ
       DSUZ=PRL*TS

! ------------------------------------------------------------------------------
!         Computation of surface runoff r0 basing on the threshold parameter SGRLUZ
          IF (SUZE .GT. SGRLUZ) THEN
             Q0=K0*(SUZE-SGRLUZ)
          ELSE
             Q0=0.0
          ENDIF
          IF (Q0 .LE. 0.0) Q0=0.0

! ------------------------------------------------------------------------------
!       Computation of interflow Q1 and, eventually karstflow Q3.
!       Q1 and Q3 generates from the storage of gravitative water in the
!       unsaturated zone SUZ
        IF (SUZE .GT. 0.0) THEN
          Q1=K1*SUZE
          Q3=0.0
          IF (K3 .GT. 0.0) Q3=K3*SUZE
        ELSE
          Q1=0.0
          Q3=0.0
        ENDIF
! ------------------------------------------------------------------------------
       PERC=CPERC*TN
       IF (PERC .GT. SUZE) PERC=SUZE
       addperc=addperc+perc

!      Actualization of the storage of gravitative water in the
!      unsaturated zone SUZ
       SUZE=SUZE+DSUZ-(PERC+(Q0*TN)+(Q1*TN)+(Q3*TN))

       IF (SUZE .LT. 0.0) THEN
            IF (Q3 .LT. 0.0) THEN
               Q1=Q1+Q3
               Q3=0.0
            ENDIF
            IF (Q1 .LT. 0.0) THEN
               Q0=Q0+Q1
               Q1=0.0
            ENDIF
            SUZE=0.0
       ENDIF

! ------------------------------------------------------------------------------
!     Slowcomp

!12    CONTINUE
      IF ((SLZ1MAX .GT. 0.0) .and. (CG1H .GT. 0.0)) THEN
         SLZ1=SLZ1E*EXP(-TS/CG1H/TH)
         IF (SLZ1 .GE. SLZ1MAX) THEN
            DIFF=PERC
            GOTO 10
         ENDIF
         DIFF=(SLZ1MAX-SLZ1)/(CG1H*RIM*TH)
         IF (DIFF .GT. PERC) THEN
            SLZ1=SLZ1+((1.0-EXP(-TS/CG1H/TH))*PERC*(CG1H*RIM*TH))
            DIFF=0.0
           ELSE
            SLZ1=SLZ1+((1.0-EXP(-TS/CG1H/TH))*(DIFF*CG1H*RIM*TH))
            DIFF=PERC-DIFF
         ENDIF
       ELSE
         DIFF=PERC
         SLZ1=0.0
         RG1=0.0
      ENDIF
10    CONTINUE
      IF (CG2H .GT. 0.0) THEN
         IF (DIFF .LE. 0.0) THEN
!          no infiltration in the slow groundwater component SLZ2+SLZ3
            DIFF=0.0
!              print*,'1.66',SLZ2,SLZ2E,EXP(-TS/CG2H/TH),TS,TH
            SLZ2=SLZ2E*EXP(-TS/CG2H/TH)
            SLZ3=SLZ3E*EXP(-TS/CG3H/TH)
!              print*,SLZ1,SLZ2,SLZ3
           ELSE
!    infiltration in the slow groundwater component SLZ2+SLZ3
            SLZ2=(SLZ2E*EXP(-TS/CG2H/TH))+((1.0-EXP(-TS/CG2H/TH))*
     *           (DIFF*0.8889*CG2H*RIM*TH))
            SLZ3=(SLZ3E*EXP(-TS/CG3H/TH))+((1.0-EXP(-TS/CG3H/TH))*
     *           (DIFF*0.1111*CG3H*RIM*TH))
         ENDIF
       ELSE
         SLZ2=0.0
         RG2=0.0
         SLZ3=0.0
         RG3=0.0
      ENDIF
      
      IF (CG1H .GT. 0.0) RG1=SLZ1/(CG1H*RIM*TH)
      IF (CG2H .GT. 0.0) THEN
         RG2=SLZ2/(CG2H*RIM*TH)
         RG3=SLZ3/(CG3H*RIM*TH)
      ENDIF

! ------------------------------------------------------------------------------
!     Storage change in the saturated zone SLZ*.
       SLZ1E=SLZ1
       SLZ2E=SLZ2
       SLZ3E=SLZ3
! ------------------------------------------------------------------------------
!    Determination of the total runoff
       Q0=Q0*TN
       Q1=Q1*TN
!11     CONTINUE
       RG1S=Rg1S+RG1
       RG2S=Rg2S+RG2
       RG3S=Rg3S+RG3
       Q2=RG1+RG2+RG3
       Q3=Q3*TN
       QGES=Q0+Q1+Q2
       R0=R0+Q0
       R1=R1+Q1
       R2=R2+Q2
       R3=R3+Q3
       RGES=RGES+QGES

      END DO

! ------------------------------------------------------------------------------
!     Actualization of the state variables
      SUZ=SUZE
      SLZ=SLZ1E+SLZ2E+SLZ3E
!        print*,'1.75',slz
      SLZ1=SLZ1E
      SLZ2=SLZ2E
      SLZ3=SLZ3E


!300   CONTINUE


! ------------------------------------------------------------------------------
      RETURN
      END
! ------------------------------------------------------------------------------

