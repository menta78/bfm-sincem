 MODULE  Service
!
 use POM,ONLY: KB
 use global_mem, ONLY:RLEN
 implicit none
 private
 real(RLEN),public                 :: pon1p,pon3n,pon4n,pon5s,smothass
 real(RLEN),public                 :: o2surf,n1psurf,n3nsurf,n5ssurf,n4nsurf 
 real(RLEN),public,dimension(KB-1) :: GPP, NPP,bacP
 real(RLEN),public,dimension(KB-1) :: GP1, GP2,GP3,NP1, NP2,NP3
 real(RLEN),public,dimension(KB-1) :: ISM
 real(rlen),public                 :: deltat
 integer,public                    :: savef,nitend
 character(200),public             :: wind_input,surfaceS_input,radiance_input,ism_input, &
                                      Sal_input,Temp_input,Sprofile_input,Tprofile_input, &
                                      heat_input,surfNut_input

 END  
