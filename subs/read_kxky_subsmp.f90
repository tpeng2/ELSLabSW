! FFT subsample setting
if (ftsubsmprto.ne.1) then
    write(*,*) 'ftsubsmprtio',ftsubsmprto,'==> Fourier spaces is partialy sampled'
    if(nx==256) then
        szftrdrow=4812
        open (unit=99, file='kxky_subsmp_nx_256.txt', status='old', action='read')
    elseif(nx==512) then
        szftrdrow=6207
        open (unit=99, file='kxky_subsmp_nx_512.txt', status='old', action='read')
    end if
    
    szftrdcol=2
    allocate(iftsubkl(szftrdrow,szftrdcol))

    do inkrow = 1,szftrdrow
        read(99, *) rdsubk,rdsubl
        write(*,*) '(n_k,n_l)',rdsubk,rdsubl
        iftsubkl(inkrow,1)=rdsubk+1 !nk=0 ==> index =1
        iftsubkl(inkrow,2)=rdsubl+1 !nl
    end do
else !if (fftsumbsmp==1)
    write(*,*) 'ftsubsmprtio',ftsubsmprto,'==> Fourier spaces is fully sampled'
    if(nx==256) then
        szftrdrow=33024
        open (unit=99, file='kxky_full_nx_256.txt', status='old', action='read')
    elseif(nx==512) then
        szftrdrow=131584
        open (unit=99, file='kxky_full_nx_512.txt', status='old', action='read')
    end if
    
    szftrdcol=2
    allocate(iftsubkl(szftrdrow,szftrdcol))

    do inkrow = 1,szftrdrow
        read(99, *) rdsubk,rdsubl
        write(*,*) '(n_k,n_l)',rdsubk,rdsubl
        iftsubkl(inkrow,1)=rdsubk+1 !nk=0 ==> index =1
        iftsubkl(inkrow,2)=rdsubl+1 !nl
    end do
end if