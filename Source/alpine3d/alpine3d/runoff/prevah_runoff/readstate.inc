! -----------------------------------------------------
! Loads the model state
! -----------------------------------------------------

      open (9876, file=loadnam, status='old', form='unformatted')

      MDSTATE=9876
      
      do ycoor=1,nrow
      do xcoor=1,ncol
        read (MDSTATE) suztf(xcoor,ycoor)
        read (MDSTATE) slz1tf(xcoor,ycoor)
        read (MDSTATE) slz2tf(xcoor,ycoor)
        read (MDSTATE) slz3tf(xcoor,ycoor)
        read (MDSTATE) slztf(xcoor,ycoor)
      end do
      end do
      
      close (9876)
     
 
! -----------------------------------------------------