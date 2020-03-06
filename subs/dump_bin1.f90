

       WRITE(which,'(I5)') 10000 + icount
    
       string1 = 'data/w' // '_' // which(1:5)       
       string2 = 'data/w_gnu' // '_' // which(1:5)       
  

       open(unit=14,file=string1,access='DIRECT',&
           & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2))
       write(14,REC=1) div_ek(:,:)
       close(14)

!      open(unit=91, file = string2)
!      write(91,*) 'set pm3d map'
!      write(91,*) 'splot "-"'
!      do j = 1,ny
!      do i = 1,nx
!      write(91,*) i, j, div_ek(i,j)
!      enddo
!      write(91,*)
!      enddo
!      close(91)

