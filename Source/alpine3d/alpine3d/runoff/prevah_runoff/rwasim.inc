c reads wasim headers

       read (uwas, rec=rec) rcol ; rec=rec+1
       read (uwas, rec=rec) rrow ; rec=rec+1
       read (uwas, rec=rec) rxll ; rec=rec+1
       read (uwas, rec=rec) ryll ; rec=rec+1
       read (uwas, rec=rec) rgridsize ; rec=rec+1
       read (uwas, rec=rec) rnodat ; rec=rec+1

       ncol=int(rcol)
       nrow=int(rrow)
       nxll=int(rxll)
       nyll=int(ryll)

       ngridsize=int(rgridsize)
       nodat=int(rnodat)

       do rloop=1,6
         read (uwas, rec=rec) wstats(rloop) ; rec=rec+1
       end do
