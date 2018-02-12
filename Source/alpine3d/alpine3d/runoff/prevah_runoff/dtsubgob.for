! ------------------------------------------------------------------------------      
      MODULE DTGLOBALS
!     Meteorological interpolation tool
!     From DTIDW.EXE
!     Version 2.0
!     Massimiliano Zappa. massimiliano.zappa@wsl.ch
!
!     Methods:
!     1. Detrended IDW
!     2. Elevation Dependent Regression
!     3. IDW
!     4. 2D Kriging (basic routines www.gslib.com)
!     5. Detrended 2D Kriging
! ------------------------------------------------------------------------------

      CHARACTER*256 IDNAM,ODNAM,METNAM,IDWNAM,REGRNAM
      CHARACTER*256 filop,filopi
      CHARACTER*4 kogrid(3)
      CHARACTER*2 YZ,MZ,DZ,HZ
      CHARACTER message*200,MASKID*1,COMPR*1,KOZZ*1000
      REAL COLUMNS,ROWS,xcoorOR,ycoorOR,SIZE,STATS(6,50),MASK
      REAL dxx,dyy,dtt
      REAL MAXDIST,IDWW,gebmean
      REAL, allocatable :: TDARR(:,:,:)
      REAL, allocatable :: STIDW(:,:,:,:)
      REAL, allocatable :: FIELDS(:,:,:)
      REAL, allocatable :: RESIDW(:,:,:)
      REAL, allocatable :: REGFIELD(:,:)
      REAL, allocatable :: ZONESLONG(:),ZONES(:,:)
      INTEGER LN,LN2,INCR,LN3,ertyp(25),GRIDID(50)
      INTEGER UNIT,UID,varm,ntg,clnam,endinta,endintb
      LOGICAL(4) OUCOMM
      INCLUDE 'wasim.h'

! ------------------------------------------------------------------------------
! Input file

      INTEGER TAGA,MONA,JAHRA,TAGE,MONE,JAHRE
      INTEGER JJOLD,DDOLD,MMOLD,MMMAX
      INTEGER I,J,MDAT6,LNGRID,k,pos
      INTEGER nmetvar,dtmin,pathln(3),extln(3)
      CHARACTER*256 gridnam
      CHARACTER*30 MODUL(6)
      CHARACTER ZN*32,kopf*75
      CHARACTER*9, allocatable :: outdt(:)
      CHARACTER*256, allocatable :: tables(:)
      CHARACTER*256, allocatable :: regress(:)
      CHARACTER*256, allocatable :: outname(:)
      REAL, allocatable :: LIMITS(:,:)
      REAL, allocatable :: missing(:)
      REAL, allocatable :: adjust(:)
      INTEGER, allocatable :: buildnam(:,:,:)
      INTEGER, allocatable :: method(:),ndigits(:)
      INTEGER, allocatable :: yopt(:,:)       
      character*256 ergrid(3)
      character yychar*4,statn*3,regform*1,gstack*1
	character datgr*6
	character*8 dateact(9),dateold(9)
      character*256 dirinput,dirgrids,diroutput   
      
! ------------------------------------------------------------------------------      
!     Stations
      INTEGER XK(200,9),YK(200,9),MUM(200,9)
      INTEGER allostat,nstat(9)
      REAL areg(3),breg(3),igu,igo,af,bf(2)

! ------------------------------------------------------------------------------      
!     Kriging
      integer nx,ny,rnst,rktype,rstmin,rstmax
      REAL, allocatable :: KRIGFIELD(:,:)
	REAL ksum,nugget,aak,aa2k,cck,maxdistk,rkmean,rkazm,rrsts
      
! ------------------------------------------------------------------------------
!     Data and date
      INTEGER JJ,MM,DD,HH
      INTEGER JJR,MMR,DDR,HHR,JJIN
      REAL DATEN(900),DATAREG(900),RESID(900)

	END MODULE DTGLOBALS

