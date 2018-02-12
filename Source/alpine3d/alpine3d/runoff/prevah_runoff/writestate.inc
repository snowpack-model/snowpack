! -----------------------------------------------------
! Saves the model state
! -----------------------------------------------------

      open (9876, file=savenam, status='unknown', form='unformatted')

      MDSTATE=9876
    
      do ycoor=1,nrow
      do xcoor=1,ncol
        write (MDSTATE) suztf(xcoor,ycoor)
        write (MDSTATE) slz1tf(xcoor,ycoor)
        write (MDSTATE) slz2tf(xcoor,ycoor)
        write (MDSTATE) slz3tf(xcoor,ycoor)
        write (MDSTATE) slztf(xcoor,ycoor)
      end do
      end do
      
      close (9876)
      
! -----------------------------------------------------