c reads wasim grid
       ALLOCATE (WRITE2D(ncol,nrow))
       do ycoor=nrow,1,-1
!          read (uwas, rec=rec) (write2D(xcoor,ycoor),xcoor=1,ncol)
          do xcoor=1,ncol
             read (uwas, rec=rec) rtmp
             write2D(xcoor,ycoor)=rtmp
             rec=rec+1
          end do
       end do
!       close (uwas)
