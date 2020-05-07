
! Tianze 5/5/2020
! Now separating AG+ and AG-
! We have eta_ag and div_ag=div, for each field apply the thrid leap-frog stage

! wave number k,l
     real k,l,M,c
     complex eye,M,eta_ag_pl,eta_ag_mi
     eye=(0.0,1.0)
     c=g_rd*He
     M=eye*sqrt(c**2*kappa**2+f0**2)*g/c**2 !dim: same as eta_A
     etaagfft_pls=1/(2*M)*(M)
     