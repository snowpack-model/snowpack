      SUBROUTINE POSNAM (namfile,opt,ext)
      
      integer ext,opt,i
      character namfile*256
      
      if (opt .eq. 1) then
      do i=2,256
              if (namfile((i-1):i) .eq. '  ') THEN
              ext=i-2
              if (namfile(ext:ext) .eq. '.') ext=ext-1
              exit
              end if
      end do      
      end if
      
      if (opt .eq. 2) then
      do i=2,256
              if (namfile(i:i) .eq. '.') THEN
              ext=i-1
              exit
              end if
      end do      
      end if
      
      return
      end
      
      
      
      
