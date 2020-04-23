
  eta(:,:,1,3) = Psurf(:,:)
  
  WRITE(which,'(I6)') 100000 + icount
  
   string1 = 'data/u_o'  // '_' // trim(which)
   string2 = 'data/v_o'  // '_' // trim(which)
   string3 = 'data/eta_o'  // '_' // trim(which)
   string4 = 'data/div_o'  // '_' // trim(which)
!  string5 = 'data/Uek'  // '_' // trim(which)
!  string6 = 'data/Vek'  // '_' // trim(which)


!  string5 = 'data/Vek'  // '_' // trim(which)
!  string6 = 'data/forci_qg'  // '_' // trim(which)
!  string7 = 'data/forci_ag'  // '_' // trim(which)
!  string8 = 'data/forci_to'  // '_' // trim(which)
!  string9 = 'data/dissi_u'  // '_' // trim(which)
!  string10 = 'data/dissi_v'  // '_' // trim(which)
!  string11 = 'data/q'  // '_' // trim(which)
!  string12 = 'data/taux'  // '_' // trim(which)
!  string13 = 'data/u_qg'  // '_' // trim(which)
!  string14 = 'data/u_ag'  // '_' // trim(which)
!  string15 = 'data/v_qg'  // '_' // trim(which)
!  string16 = 'data/v_ag'  // '_' // trim(which)
!  string17 = 'data/eta_qg'  // '_' // trim(which)
!  string18 = 'data/eta_ag'  // '_' // trim(which)
!  string19 = 'data/div_ek'  // '_' // trim(which)


! Note indices for (u,v,eta ...) starting with 0, useful part is 1:256
  !  real u_out(0:nx/subsmprto+1,0:ny/subsmprto+1,nz), v_out(0:nx/subsmprto+1,0:ny/subsmprto+1,nz)
  !  real div_ek_out(0:nx/subsmprto+1,0:ny/subsmprto+1), eta_out(0:nx/subsmprto+1,0:ny/subsmprto+1,nz)
 
   div_ek_out=div_ek(isubx,isuby)

   do k = 1,nz
   u_out(:,:,k)=u(isubx,isuby,k,3)
   v_out(:,:,k)=v(isubx,isuby,k,3)
   eta_out(:,:,k)=eta(isubx,isuby,k,3)
   enddo

   open(unit=14,file=string1,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isubx)*nz))
   write(14,REC=1) (((u_out(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto),k=1,nz)

   close(14)
  
   open(unit=15,file=string2,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
        write(15,REC=1) (((v_out(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto),k=1,nz)
   close(15)
  
  open(unit=16,file=string3,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
        write(16,REC=1) (((eta_out(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto),k=1,nz)
  close(16)
  
  open(unit=17,file=string4,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(17,REC=1) ((div_ek_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
  close(17)
  
!  open(unit=18,file=string5,access='DIRECT',&
!        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx/2+2)*(ny/2+2))
!  write(18,REC=1) Uek_out(:,:)
!  close(18)
  
!  open(unit=19,file=string6,access='DIRECT',&
!        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx/2+2)*(ny/2+2))
!  write(19,REC=1) Vek_out(:,:)
!  close(19)
  

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

