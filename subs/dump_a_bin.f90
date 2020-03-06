
  eta(:,:,1,3) = Psurf(:,:)

  WRITE(which,'(I5)') 10000 + icount
  
   string1 = 'data/u_o'  // '_' // which(1:5)
   string2 = 'data/v_o'  // '_' // which(1:5)
   string3 = 'data/eta_o'  // '_' // which(1:5)
   string4 = 'data/w_o'  // '_' // which(1:5)
   string5 = 'data/Uek'  // '_' // which(1:5)
   string6 = 'data/Vek'  // '_' // which(1:5)
!  string5 = 'data/Vek'  // '_' // which(1:5)
!  string6 = 'data/forci_qg'  // '_' // which(1:5)
!  string7 = 'data/forci_ag'  // '_' // which(1:5)
!  string8 = 'data/forci_to'  // '_' // which(1:5)
!  string9 = 'data/dissi_u'  // '_' // which(1:5)
!  string10 = 'data/dissi_v'  // '_' // which(1:5)
!  string11 = 'data/q'  // '_' // which(1:5)
!  string12 = 'data/taux'  // '_' // which(1:5)
!  string13 = 'data/u_qg'  // '_' // which(1:5)
!  string14 = 'data/u_ag'  // '_' // which(1:5)
!  string15 = 'data/v_qg'  // '_' // which(1:5)
!  string16 = 'data/v_ag'  // '_' // which(1:5)
!  string17 = 'data/eta_qg'  // '_' // which(1:5)
!  string18 = 'data/eta_ag'  // '_' // which(1:5)
!  string19 = 'data/div_ek'  // '_' // which(1:5)


   do j = 1,ny,2
   jj = 1 + (j-1)/2
   do i = 1,nx,2
      ii = 1 + (i-1)/2
      div_ek_out(ii,jj) = div_ek(i,j)
   enddo
   enddo

   do k = 1,nz
   do j = 1,ny,2
   jj = 1 + (j-1)/2
   do i = 1,nx,2
      ii = 1 + (i-1)/2
      u_out(ii,jj,k) = u(i,j,k,3)
      v_out(ii,jj,k) = v(i,j,k,3)
      eta_out(ii,jj,k) = eta(i,j,k,3)
      Uek_out(ii,jj) = Uek(i,j,2)
      Vek_out(ii,jj) = Vek(i,j,2)
   enddo
   enddo
   enddo
   open(unit=14,file=string1,access='DIRECT',&
        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx/2+2)*(ny/2+2)*2)
!       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!  write(14,REC=1) u(:,:,:,3)
   write(14,REC=1) u_out(:,:,:)
   close(14)
  
   open(unit=15,file=string2,access='DIRECT',&
        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx/2+2)*(ny/2+2)*2)
!       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!  write(15,REC=1) v(:,:,:,3)
   write(15,REC=1) v_out(:,:,:)
   close(15)
  
  open(unit=16,file=string3,access='DIRECT',&
        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx/2+2)*(ny/2+2)*2)
!      & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
! write(16,REC=1) eta(:,:,:,3)
  write(16,REC=1) eta_out(:,:,:)
  close(16)
  
  open(unit=17,file=string4,access='DIRECT',&
        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx/2+2)*(ny/2+2))
!      & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
! write(17,REC=1) div_ek(:,:)
  write(17,REC=1) div_ek_out(:,:)
  close(17)
  
  open(unit=18,file=string5,access='DIRECT',&
        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx/2+2)*(ny/2+2))
  write(18,REC=1) Uek_out(:,:)
  close(18)
  
  open(unit=19,file=string6,access='DIRECT',&
        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx/2+2)*(ny/2+2))
  write(19,REC=1) Vek_out(:,:)
  close(19)
  
!  open(unit=18,file=string5,access='DIRECT',&
!       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!  write(18,REC=1) Vek(:,:,2)
!  close(18)

!  open(unit=19,file=string6,access='DIRECT',&
!       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!  write(19,REC=1) forcing_qg(:,:)
!  close(19)
  
!  open(unit=20,file=string7,access='DIRECT',&
!       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!  write(20,REC=1) forcing_ag(:,:)
!  close(20)
 
!  open(unit=21,file=string8,access='DIRECT',&
!       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!  write(21,REC=1) forcing_total(:,:)
!  close(21)
 
!  open(unit=22,file=string9,access='DIRECT',&
!       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!  write(22,REC=1) dissi_u(:,:)
!  close(22)

!  open(unit=23,file=string10,access='DIRECT',&
!       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!  write(23,REC=1) dissi_v(:,:)
!  close(23)

!  open(unit=24,file=string11,access='DIRECT',&
!       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!  write(24,REC=1) q(:,:,:)
!  close(24)

!  open(unit=25,file=string12,access='DIRECT',&
!       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!  write(25,REC=1) taux(:,:)
!  close(25)

!   open(unit=26,file=string13,access='DIRECT',&
!        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!   write(26,REC=1) u_qg(:,:,:)
!   close(26)

!   open(unit=27,file=string14,access='DIRECT',&
!        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!   write(27,REC=1) u_ag(:,:,:)
!   close(27)

!   open(unit=28,file=string15,access='DIRECT',&
!        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!   write(28,REC=1) v_qg(:,:,:)
!   close(28)

!   open(unit=29,file=string16,access='DIRECT',&
!        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!   write(29,REC=1) v_ag(:,:,:)
!   close(29)

!   open(unit=30,file=string17,access='DIRECT',&
!        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!   write(30,REC=1) eta_qg(:,:)
!   close(30)
   
!   open(unit=31,file=string18,access='DIRECT',&
!        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!   write(31,REC=1) eta_ag(:,:)
!   close(31)

!   open(unit=32,file=string19,access='DIRECT',&
!        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!   write(32,REC=1) div_ek(:,:)
!   close(32)

