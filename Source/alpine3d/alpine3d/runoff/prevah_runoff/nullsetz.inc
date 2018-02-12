! set variables to zero at model initialization

      do ycoor=1,nrow
      do xcoor=1,ncol
          sitf(xcoor,ycoor)=0.0
          ssotf(xcoor,ycoor)=0.0
          ssnowtf(xcoor,ycoor)=0.0
          sliqtf(xcoor,ycoor)=0.0
          ssztf(xcoor,ycoor)=0.0
          suztf(xcoor,ycoor)=0.0
          slztf(xcoor,ycoor)=0.0
          slz1tf(xcoor,ycoor)=0.0
          slz2tf(xcoor,ycoor)=0.0
          slz3tf(xcoor,ycoor)=0.0
          ssmtf(xcoor,ycoor)=0.0
          bil(xcoor,ycoor)=0.0
          ENDDO   
      ENDDO
      