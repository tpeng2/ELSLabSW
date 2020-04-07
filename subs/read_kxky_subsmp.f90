

    if(nx==256) then
        szftrdrow=4633
        open (unit=99, file='kxky_subsmp_256.txt', status='old', action='read')
    elseif(nx==512) then
        szftrdrow=6143
        open (unit=99, file='kxky_subsmp_512.txt', status='old', action='read')
    end if
    szftrdcol=2
    allocate(iftsubkl(szftrdrow,szftrdcol))

    do inkrow = 1,szftrdrow
    read(99, *) rdsubk,rdsubl
    write(*,*) n,m
    iftsubkl(inkrow,1)=rdsubk+1 !nk=0 ==> index =1
    iftsubkl(inkrow,2)=rdsubl+1 !nl
    end do
