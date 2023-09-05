subroutine atmdisrom(N_in, N_receptor, Nstep, param, x_in, y_in, x_re, y_re, rate_in, out_flag, & 
           num_source, x_new, y_new, critical_distance, log_message) &
           bind(C, name='atmdisrom')
!C***********************************************************************
!CPS
! Note: The third equation in B&M solution in Table 2.15 should be: beta=-0.5alpha+1.78
! To compile: 1. gfortran -c AtmDisROM.f90
!             2. gfortran -shared -o atmdisrom.dll AtmDisRom.o  (windows)
!                or gfortran -dynamiclib -o atmdisrom.dylib AtmDisROM.o (MacOS)
! N_in: number of leakage points;   N_receptor: number of receptors; Nstep: number of time steps
! x_in, Y_in: leakage locations; x_re,y_re: receptor
!
!C***********************************************************************


     implicit none

     integer, intent(in):: N_in,  N_receptor, Nstep
     integer, intent(out):: log_message
     real*8, dimension(5), intent(in):: param
     real*8, dimension(N_in), intent(in):: x_in, y_in, rate_in
     real*8, dimension(N_receptor), intent(in):: x_re, y_re
     integer, dimension(N_receptor), intent(out):: out_flag   

     integer nc, nloc, i, j, mtotal, k, mloc, error, num_check, m_check, i_source, num_count,ntotal_above
     integer i1, i2, i3
     real*8, dimension(N_in):: xx_coor, yy_coor, rate_final                   
     real*8, dimension (N_in, N_in):: dis
 
     real*8 T_amb, P_amb, V_wind, C0_critical,C_critical, AirDensity_ambient, CO2_density, Initial_buoyancy_factor
     real*8 Cloud_Rep_criteria, CO2_leakageRate_Vol, D_characteristic, alpha, Beta_intepolated, T_source, gas_density
     real*8 Cupper_select, Clower_select, beta_BM_upper, beta_BM_lower, distance
     real*8 t1, t2, t3, rate1, rate2, dtime
     real*8, allocatable, dimension(:)::  CriticalX_downwind
     integer, allocatable, dimension(:):: mflag, kflag, m_rep
     real*8, dimension(N_in), intent(out)::  x_new, y_new, critical_distance
     integer, intent(out) :: num_source
    
     integer iflag_source


! Body 

     log_message = 0        ! valid B&M

! Parameters of the model   
     T_amb = param(1)       ! ambient temperature for atmospheric dispersion model
     P_amb = param(2)       ! ambient pressure
     V_wind = param(3)      ! wind speed at 10m height
     C0_critical = param(4) ! critial concentration
     T_source = param(5)    ! source temperature

          
! Receptor locations - locations where concentration of CO2 is a concern
     mloc = N_receptor  ! number of receptors - locations that CO2 concentration is a concern, these locations are read in through a file "CheckPoints_CO2concentration.txt", max=100
                  
          
! read in the well locations where there is a CO2 leakage, as well as the leakage rate         
     nloc = N_in   ! number of leakage locations
          
     rate_final=0
 
     i=1
     j=1
          
     rate1=0.
          
     Do i = 1, nloc
        xx_coor(i) =  x_in(i)
        yy_coor(i) =  y_in(i) 
        rate_final(i) =  rate_in(i)      
        rate1 = rate1 +rate_final(i)
     enddo
     mtotal = nloc
                   
     allocate (CriticalX_downwind(mtotal), stat=error)   
     allocate (mflag(mloc), stat=error)   
     allocate (kflag(mtotal), stat=error)        
     allocate (m_rep(mtotal), stat=error)             
          
     AirDensity_ambient = P_amb*29.92/((T_amb+273.15)*0.082057) ! unit: kg/m^3
     CO2_density = 44.0 *P_amb/(0.082057*(T_source+273.15))        ! unit: g/L
     Initial_buoyancy_factor = 9.8*(CO2_density-AirDensity_ambient)/AirDensity_ambient
     C_critical=C0_critical/(C0_critical+(1-C0_critical)*((T_amb+273.15)/(T_source+273.15)))    
          
          if (C_critical>=0.1) then
              Cupper_select = 0.1
              Clower_select = 0.1              
          elseif ((C_critical<0.1).and.(C_critical>=0.05)) then
              Cupper_select = 0.1
              Clower_select = 0.05              
          elseif ((C_critical<0.05).and.(C_critical>=0.02)) then   
              Cupper_select = 0.05
              Clower_select = 0.02             
          elseif ((C_critical<0.02).and.(C_critical>=0.01)) then   
              Cupper_select = 0.02  
              Clower_select = 0.01             
          elseif ((C_critical<0.01).and.(C_critical>=0.005)) then
              Cupper_select = 0.01
              Clower_select = 0.005              
          elseif ((C_critical<0.005).and.(C_critical>=0.002)) then                 
              Cupper_select = 0.005
              Clower_select = 0.002             
          else
              Cupper_select = 0.002
              Clower_select = 0.002
          endif
 
 ! iterations start here: multiple sources, if their critical radius overlap, add the sources together. iflag_source is the loop flag.
 ! kflag records the source terms that are added to others, so they become inactive next time in the loop.
          iflag_source = 0
          kflag = 0
          m_rep = 0
          Do while (iflag_source<2)         
          num_check = 0
          Do i = 1, mtotal      ! calculate the critical distance for each leakage       
              
            if (kflag(i)<1) then     ! Only do the calculation of radius for active source                   
               CO2_leakageRate_Vol = rate_final(i)/CO2_density      ! unit: m^3/s
               D_characteristic = sqrt(CO2_leakageRate_Vol/V_wind ) ! unit: m
               Cloud_Rep_criteria = (Initial_buoyancy_factor*CO2_leakageRate_Vol/(V_wind**3.*D_characteristic))**(1./3.)
               if (Cloud_Rep_criteria<0.15) then
 !                write (*,*) "Not qualified as dense gas release, please consider alternative methods"
 !                  stop
                  log_message=1
               endif
 ! Notice there is a typo in the Table 1 of user's guide. Based on Figure 1, there is a coefficient of 1/5 in alpha            
               alpha = 0.2*log10(Initial_buoyancy_factor**2.*CO2_leakageRate_Vol/V_wind**5.)    
               If (alpha>1.) then
                   write (*,*) "Alpha in Table 1 (Users Guide) is out of range. B&M is not a valid approach"
				   write (*,*) "alpha=", alpha
 !                   stop
                   log_message=log_message+2
               endif

 ! Calculate beta upper limit for interpolation            
                 if ((Cupper_select==0.1) .and. (alpha<=-0.55)) then
                    beta_BM_upper = 1.75
                 elseif ((Cupper_select==0.1) .and. (alpha>-0.55) .and. (alpha<=-0.14)) then
                    beta_BM_upper = 1.88+0.24*alpha
                 elseif ((Cupper_select==0.1) .and. (alpha>-0.14) .and. (alpha<=1)) then
                    beta_BM_upper = 1.78-0.5*alpha     ! another typo in table 1, should be -0.5 instead of 0.5
                 elseif ((Cupper_select==0.05) .and. (alpha<=-0.68)) then
                    beta_BM_upper = 1.92
                 elseif ((Cupper_select==0.05) .and. (alpha>-0.68) .and. (alpha<=-0.29)) then
                    beta_BM_upper = 2.16+0.36*alpha
                 elseif ((Cupper_select==0.05) .and. (alpha>-0.29) .and. (alpha<=-0.18)) then                 
                    beta_BM_upper = 2.06
                 elseif ((Cupper_select==0.05) .and. (alpha>-0.18) .and. (alpha<=1)) then                
                    beta_BM_upper = 1.96-0.56*alpha
                 elseif ((Cupper_select==0.02) .and.  (alpha<=-0.69)) then
                    beta_BM_upper = 2.08
                 elseif ((Cupper_select==0.02) .and. (alpha>-0.69) .and. (alpha<=-0.31)) then                 
                    beta_BM_upper = 2.39+0.45*alpha
                 elseif ((Cupper_select==0.02) .and. (alpha>-0.31) .and. (alpha<=-0.16)) then                                 
                    beta_BM_upper = 2.25
                 elseif ((Cupper_select==0.02) .and. (alpha>-0.16) .and. (alpha<=1)) then                 
                    beta_BM_upper = 2.16-0.54*alpha
                 elseif ((Cupper_select==0.01) .and.  (alpha<=-0.70)) then                    
                    beta_BM_upper = 2.25
                 elseif ((Cupper_select==0.01) .and. (alpha>-0.70) .and. (alpha<=-0.29)) then   
                    beta_BM_upper = 2.59+0.49*alpha
                 elseif ((Cupper_select==0.01) .and. (alpha>-0.29) .and. (alpha<=-0.20)) then                    
                    beta_BM_upper = 2.45
                 elseif ((Cupper_select==0.01) .and. (alpha>-0.20) .and. (alpha<=1)) then                    
                    beta_BM_upper = 2.35-0.52*alpha
                 elseif ((Cupper_select==0.005) .and. (alpha<=-0.67)) then                    
                    beta_BM_upper = 2.40
                 elseif ((Cupper_select==0.005) .and. (alpha>-0.67) .and. (alpha<=-0.28)) then                    
                    beta_BM_upper = 2.80+0.59*alpha
                 elseif ((Cupper_select==0.005) .and. (alpha>-0.28) .and. (alpha<=-0.15)) then                    
                    beta_BM_upper = 2.63
                 elseif ((Cupper_select==0.005) .and. (alpha>-0.15) .and. (alpha<=1)) then                           
                    beta_BM_upper = 2.56-0.49*alpha
                 elseif ((Cupper_select==0.002) .and.  (alpha<=-0.69)) then                    
                    beta_BM_upper = 2.60
                 elseif ((Cupper_select==0.002) .and. (alpha>-0.69) .and. (alpha<=-0.25)) then                    
                    beta_BM_upper = 2.87+0.39*alpha
                 elseif ((Cupper_select==0.002) .and. (alpha>-0.25) .and. (alpha<=-0.13)) then                           
                    beta_BM_upper = 2.77
                 elseif ((Cupper_select==0.002) .and. (alpha>-0.13) .and. (alpha<=1)) then                    
                    beta_BM_upper = 2.71-0.5*alpha
                 else                            
                    write (*,*) "Use other models. beta_BM_upper out of range:", beta_BM_upper
                 endif                     

  ! Calculate beta lower limit for interpolation                  
                 if ((Clower_select==0.1) .and. (alpha<=-0.55)) then
                    beta_BM_lower =  1.75         
                 elseif ((Clower_select==0.1) .and. (alpha>-0.55) .and. (alpha<=-0.14)) then
                    beta_BM_lower =  1.88+0.24*alpha         
                 elseif ((Clower_select==0.1) .and. (alpha>-0.14) .and. (alpha<=1)) then
                    beta_BM_lower =  1.78-0.5*alpha         
                 elseif ((Clower_select==0.05) .and. (alpha<=-0.68)) then    
                    beta_BM_lower =  1.92         
                 elseif ((Clower_select==0.05) .and. (alpha>-0.68) .and. (alpha<=-0.29)) then
                    beta_BM_lower = 2.16+0.36*alpha        
                 elseif ((Clower_select==0.05) .and. (alpha>-0.29) .and. (alpha<=-0.18)) then
                    beta_BM_lower = 2.06         
                 elseif ((Clower_select==0.05) .and. (alpha>-0.18) .and. (alpha<=1)) then                 
                    beta_BM_lower = 1.96-0.56*alpha         
                 elseif ((Clower_select==0.02) .and.  (alpha<=-0.69)) then
                    beta_BM_lower = 2.08          
                 elseif ((Clower_select==0.02) .and. (alpha>-0.69) .and. (alpha<=-0.31)) then
                    beta_BM_lower = 2.39+0.45*alpha          
                 elseif ((Clower_select==0.02) .and. (alpha>-0.31) .and. (alpha<=-0.16)) then                 
                    beta_BM_lower = 2.25          
                 elseif ((Clower_select==0.02) .and. (alpha>-0.16) .and. (alpha<=1)) then
                    beta_BM_lower = 2.16-0.54*alpha          
                 elseif ((Clower_select==0.01) .and.  (alpha<=-0.70)) then
                    beta_BM_lower = 2.25         
                 elseif ((Clower_select==0.01) .and. (alpha>-0.70) .and. (alpha<=-0.29)) then                 
                    beta_BM_lower = 2.59+0.49*alpha          
                 elseif ((Clower_select==0.01) .and. (alpha>-0.29) .and. (alpha<=-0.20)) then
                    beta_BM_lower = 2.45          
                 elseif ((Clower_select==0.01) .and. (alpha>-0.20) .and. (alpha<=1)) then
                    beta_BM_lower = 2.35-0.52*alpha          
                 elseif ((Clower_select==0.005) .and. (alpha<=-0.67)) then    
                    beta_BM_lower = 2.40          
                 elseif ((Clower_select==0.005) .and. (alpha>-0.67) .and. (alpha<=-0.28)) then
                    beta_BM_lower = 2.80+0.59*alpha          
                 elseif ((Clower_select==0.005) .and. (alpha>-0.28) .and. (alpha<=-0.15)) then
                    beta_BM_lower = 2.63          
                 elseif ((Clower_select==0.005) .and. (alpha>-0.15) .and. (alpha<=1)) then                 
                    beta_BM_lower = 2.56-0.49*alpha           
                 elseif ((Clower_select==0.002) .and. (alpha<=-0.69)) then
                    beta_BM_lower = 2.60         
                 elseif ((Clower_select==0.002) .and. (alpha>-0.69) .and. (alpha<=-0.25)) then
                    beta_BM_lower = 2.87+0.39*alpha           
                 elseif ((Clower_select==0.002) .and. (alpha>-0.25) .and. (alpha<=-0.13)) then                 
                    beta_BM_lower = 2.77          
                 elseif ((Clower_select==0.002) .and. (alpha>-0.13) .and. (alpha<=1)) then
                    beta_BM_lower = 2.71-0.5*alpha          
                 else
                    write (*,*) "Use other models. beta_BM_lower out of range:", beta_BM_lower           
                 endif     
             
                 Beta_intepolated=beta_BM_lower+(beta_BM_upper-beta_BM_lower)*log10(C_critical/Clower_select)  &
				                   /log10(Cupper_select/Clower_select)
                 if (abs(beta_BM_upper-beta_BM_lower)<=1.e-8)  Beta_intepolated = beta_BM_upper
                 CriticalX_downwind(i) = 10**(Beta_intepolated)*D_characteristic
            else 
                num_check = num_check + 1              
            endif
          Enddo  
          
          
          if (iflag_source<1) then
                 
             i_source = 0         
             Do i = 1, mtotal - 1                 
                if (kflag(i)<1) then

                   Do j=i+1, mtotal
                   if (kflag(j)<1) then    
                      dis(i, j)= sqrt((xx_coor(i)-xx_coor(j))**2.+(yy_coor(i)-yy_coor(j))**2.)
                      if ((dis(i, j).le.CriticalX_downwind(i)).or.(dis(i, j).le.CriticalX_downwind(j))) then
                         xx_coor(i) = CriticalX_downwind(j)/(CriticalX_downwind(i)+CriticalX_downwind(j))      &
                                      *(xx_coor(j)-xx_coor(i))+xx_coor(i)
                         yy_coor(i) = CriticalX_downwind(j)/(CriticalX_downwind(i)+CriticalX_downwind(j))      &
                                      *(yy_coor(j)-yy_coor(i))+yy_coor(i)
                         rate_final(i) = rate_final(i) + rate_final(j)
                         rate_final(j) = 0.
                         kflag(j)=1
                         m_rep(i) = m_rep(i) + 1     
                         i_source = i_source + 1
                      endif  
                   endif                    
                   End do
 
                endif
                
             enddo   
          endif   
          iflag_source = iflag_source + 1       
          
          Enddo
          
          m_check = 0
          Do i = 1, mtotal
             m_check = m_check + m_rep(i)
          enddo
          
          if (m_check.ne.num_check) write (*,*) "mistake in source adding"
 ! Output the critical distance to a file for debugging or independent use.
 

!          write (21,*)  "LocationNumber    x_cor   y_cor   leakage_rate   Critical_downwind_distance"

          num_source=0                    
          Do i = 1,  mtotal
              if (kflag(i)<1) then
                 num_source=num_source+1
                 x_new(num_source)=xx_coor(i)
                 y_new(num_source)=yy_coor(i)
                 critical_distance(num_source)=CriticalX_downwind(i)
              endif
          enddo          
          

  ! Now check all the locations if they are within the critical distance radius of leakage wells.  
  ! If mflag(i)>0, then location i has a risk that CO2 concentration is above critical value, The higher the value, the more likely it happens. 
 
 
 !      open (22, file=file_receptor)
    
 !          write (22,*)  "x_cor   y_cor   Flag"
          out_flag(1:mloc) = 0 
          mflag(1:mloc) = 0     ! assuming the location is not within the critical distance radius of a leakage w    
          Do i = 1, mloc           
             Do j = 1, mtotal
                distance = sqrt((x_re(i)-xx_coor(j))**2+(y_re(i)-yy_coor(j))**2)
                if (distance.le. CriticalX_downwind(j).and.(CriticalX_downwind(j).gt.1.e-20)) then
                   mflag(i) = mflag(i) + 1
                endif               
             enddo
                        
          Enddo       
          out_flag(1:mloc) = mflag(1:mloc)                   

       close (22)
 
 
 
        deallocate (CriticalX_downwind, mflag)   

     
       
    return

    end subroutine atmdisrom