!     Initialization of gridded output

      gkwfile='remap.tmp'
      dpunkt=':'       
      grdakt='Y'
      opnssm=0

!812   FORMAT (2(F10.2,1X),A1,1X,F15.5)
      
      ergrid(1)='.SPK'
      ergrid(2)='.RGS'
      ergrid(3)='.R0' 
      ergrid(4)='.R1'
      ergrid(5)='.RGQ'
      ergrid(6)='.RGS'
      ergrid(7)='.RG3'
      ergrid(8)='.R2'
      ergrid(9)='.GWN'
      ergrid(10)='.SUZ'
      ergrid(11)='.SLZ'
      ergrid(12)='.SL1'
      ergrid(13)='.SL2'
      ergrid(14)='.SL3'

      labelgrid(1)='snpk'
      labelgrid(2)='rges'
      labelgrid(3)='rsur' 
      labelgrid(4)='rint'
      labelgrid(5)='rb_1'
      labelgrid(6)='rb_2'
      labelgrid(7)='rb_3'
      labelgrid(8)='rbas'
      labelgrid(9)='rech'
      labelgrid(10)='suzs'
      labelgrid(11)='slzs'
      labelgrid(12)='slz1'
      labelgrid(13)='slz2'
      labelgrid(14)='slz3'
      
      MONUM(1)=31
      MONUM(2)=28
      MONUM(3)=31
      MONUM(4)=30
      MONUM(5)=31
      MONUM(6)=30
      MONUM(7)=31
      MONUM(8)=31
      MONUM(9)=30
      MONUM(10)=31
      MONUM(11)=30
      MONUM(12)=31
      
      DO I=1,14
        whgrid(I)='N-N-N-N-N'
      END DO

!     First initialization

      MTA=5
      MBI=6
      MDTAG=8
      MDMON=9
      MDR=42
      MLOG=48
      MPRI=MDR
      MIDAT=10
      MDSTD=11
      MDGKW=12
      MDMIT=41
      MDTMP=40
          MDHOU=13
          MDQ=14

      AEGRID=0.01
      SUZ=0.0
      SLZ1=0.0
      SLZ2=0.0
      SLZ3=0.0
 
       
              DO i=1,32
                  zn(i:i)=' '
            ENDDO
      
            DO i=1,6
                  filen(i)='            '
                  mdatn(i)=1111+i
      ENDDO