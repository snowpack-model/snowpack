      MODULE RUNGLOBALS


      integer i,j,POS,clnam,unit,MINS,ln3,lngrid,pos3
      integer tag,taga,tage,mon,mona,mone
      integer jahr,jahra,jahre
      integer tagbeg,tagend,monbeg,monend
      integer mta,mbi,mdr,mdtag,mdmon,mdstd,mdtmp,mdmit,midat,mdgkw,mlog
      integer mdat6,mdhou,mdq
      integer mdatn(6),idmin
      integer mpri,rmcycle
      integer posi3,ntg,monum(12)

! ------------------------------------------------------------------------------

      real aeges,nvalids
      real geobreite,aelhg
      real cln,qga,qlps,translation(99),storage(99)
      real aegrid,allmeans(14),gebmean,addperc

      real, allocatable :: r0js(:,:)
      real, allocatable :: r0s(:,:),r1s(:,:)
      real, allocatable :: r1ms(:,:),r1js(:,:),r2s(:,:),r2ms(:,:)
      real, allocatable :: r2js(:,:),bil(:,:)
      real, allocatable :: exportgrid(:,:)
      real, allocatable :: ssotf(:,:),ssnowtf(:,:)
      real, allocatable :: sliqtf(:,:),sitf(:,:),ssztf(:,:),slz2tf(:,:)
      real, allocatable :: ssmtf(:,:),suztf(:,:),slztf(:,:),slz1tf(:,:)
      real, allocatable :: mzones(:,:),muse(:,:),spackin(:,:)
      real, allocatable :: spr0(:,:),spr1(:,:),spr2(:,:),sprges(:,:)
      real, allocatable :: sprgq(:,:),sprgs(:,:),sprg3(:,:),sprperc(:,:)
      real, allocatable :: ZONESLONG(:),ZONES(:,:),slz3tf(:,:)
 
! ------------------------------------------------------------------------------      

      character*1  dorout,spclean
      character*4  dum
      character*15 zierg
      character*256 namaus,filen(6),namlog,
     *              namdru,namstd,nammit,routeffnam,
     *              namq,gridnam,assnam,tabnam
      character*32 zn
      CHARACTER*78 kopf
      CHARACTER*30 MODUL(20)
      
! ---------------------------------------------------------------------
!     Tuneable parameters
! ---------------------------------------------------------------------

      REAL
!for runoff-generation module
     *     SGRLUZ, ! treshold coefficient for surface runoff 
     *     CPERC, ! Percolation rate
     *     K0H, ! Storage coefficient in hours for surface runoff
     *     K0, ! Storage coefficient for surface runoff as exponent
     *     K1H, ! Storage coefficient in hours for interflow
     *     K1, ! Storage coefficient for interflow as exponent
     *     K2H,K3H, ! Storage coefficient in hours for baseflow
     *     K2,K3, ! Storage coefficient for delayed baseflow as exponent
     *     CG1H,CG2H, ! Storage coefficient in hours for quick basflow
     *     SLZ1MAX, ! Maximal content of the quick baseflow storages SLZ1
     *     CG3H, ! 1/9 of CG1H after Schwarze
     *     ADJSUZ,ADJSLZ          ! Adjustment factors for forecast runs

      REAL PRL ! Inflow from SNOWPACK in the runoff generation module

      REAL SUZ, ! Storage reservoir in the unsaturated zone
     *     SLZ,SLZ1,SLZ2,SLZ3, ! Storage reservoirs in the saturated zone
     *     RGES, ! Total runoff in portion of time step
     *     R1, ! Quick response baseflow in portion of time step
     *     R2, ! Delayed Baseflow  in portion of time step
     *     R0, ! Additional baseflow (Schwarze)  in portion of time step
     *     RG1s, ! Quick response baseflow for this time step
     *     RG2s, ! Delayed Baseflow for this time step
     *     RG3s, ! Additional baseflow (Schwarze) for this time step
     *     R3 ! Karst flow in this time step (not active) for this time step

!    Grid Ausgabe

      INCLUDE 'wasim.h'
      INTEGER LN,opnssm,beg
      REAL, allocatable :: wrgrd(:,:)  
      character*9 whgrid(14)          
      character*4 ergrid(14),labelgrid(14)
      character*256 idnam,rdgrid(3),gkwfile
      CHARACTER dpunkt*1,grdakt*1     
! ------------------------------------------------------------------------------
!     Load and Save state
      character*256 dirstate,dirmeteo,dirgrids,dirinput,rknam
      character*1 LOADID,SAVEID,SAVEMMID,SAVEMOREID,ZK
      character*8 LOADCHAR,DATUMCHAR,DATUMCHARO
      character*8, allocatable :: savechar(:)
      character*256 savenam,loadnam,effinams(9),routnams(99)
      character*4 chhhmm
      integer CHRMIN,CHRHH,FIRSTHH,nrouts
      integer SAVEMORENN,mdstate,pathln(5),extln(3)
      integer, allocatable ::SAVEDD(:),SAVEMM(:), SAVEYY(:)
      integer asshh,STDA,STDE,MINBEG,MINEND
      INTEGER dtmin(6),neffis,nmetvar   

      CHARACTER*256, allocatable :: outname(:)
      REAL, allocatable :: missing(:)
      REAL, allocatable :: adjust(:)
      INTEGER, allocatable :: buildnam(:,:,:)
      INTEGER, allocatable :: method(:)


      END MODULE RUNGLOBALS