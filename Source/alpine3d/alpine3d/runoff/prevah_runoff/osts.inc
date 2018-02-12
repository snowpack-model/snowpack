c     **************************************
c     COMPUTATION OF THE STATISTICAL
c     PARAMETERS
c     **************************************

	wstats(1)=0
	wstats(2)=1e8
	wstats(3)=-1e8
	wstats(4)=0
	wstats(5)=0
	wstats(6)=0

	do I=1,stvar
		wstats(1)=wstats(1)+1
		wstats(2)=min(wstats(2),value(I))		
		wstats(3)=max(wstats(3),value(I))
		wstats(4)=wstats(4)+value(I)
	enddo

	wstats(5)=wstats(4)/wstats(1)

	do I=1,stvar
 		wstats(6)=wstats(6)+(value(I)-wstats(5))**2
	enddo

      	If (stvar .LE. 1) wstats(6)=0.0

	If (stvar .GT. 1) wstats(6)=(wstats(6)/(wstats(1)-1))**(0.5)
	
