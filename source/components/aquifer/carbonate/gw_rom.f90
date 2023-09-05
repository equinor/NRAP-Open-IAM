    subroutine gw_plumes(ithresh,rmin, aquifer,time,    &
        n_leak,x,y,cl,co2_rate,co2_mass,brine_rate,brine_mass,   &
        vol,chem,logf) bind(C, name='gw_plumes')

    use rs_param
    use mars_coefs
    ! method specifies linear (1) quadratic (2) cubic (3) mars (4) response function
    ! ithresh specifies metric:  conc > MCL (1)  or conc > background (2)
    ! min_disance specifies how close 2 leaky wells can be before they are combined into 1 source
    ! this routine requires a call of time=0 for each leakage scenario to properly initialize   
    ! if chem(1)<0 chemical scaling is turned off

    real*8 x(n_leak),y(n_leak),aquifer(6),cl(n_leak),co2_rate(n_leak),co2_mass(n_leak),rmin,r,outArgs(24),inArgs(24)    
    real*8 brine_rate(n_leak),brine_mass(n_leak),vol(12),xplume,yplume,time,flux
    real*8 input1(13)
    integer n_leak,iw,groups(250),ng,k,ithresh,nin,im,rs_method(12),nrs
    integer quad_inputs(12,12),nhout,logf(12)
    real*8 time_stop(1000),outlink(12),xp(24),yout,chem(8),y_result
    real*8 plume_final(1000,12)
    logical debug,scale,quad,mars,leak(1000),stop_leak(1000),dis(1000,12),app(1000,12)
    save time_stop,groups,plume_final,leak,stop_leak,scale
    character(len=80) f,dir
    character(len=10) dum,dum1,fname(2),type(2)
    CHARACTER(len=255) :: cwd
    data rs_method/4,4,4,4,4,4,4,4,4,4,4,4/
    ! for now if rs_method<4 all integers have to be the same    
    data fname/'MCL/MCL','THR/THR'/
    data type/'pH','TDS'/
    ! it =  1:MCL  2:THR
    ! nrs = 1:pH-related variables  2:TDS-related variables 
    vol=-99;flux=0;xplume=0;yplume=0
    scale=.true.
    debug=.false.
    if(debug) write(*,*) 'ithresh', ithresh
    if(debug) write(*,*) 'rmin', rmin
    if(debug) write(*,*) 'aquifer', aquifer
    if(debug) write(*,*) 'time', time
    if(debug) write(*,*) 'n_leak', n_leak
    if(debug) write(*,*) 'x,y', x,y
    if(debug) write(*,*) 'cl',cl
    if(debug) write(*,*) 'co2_rate,co2_mass,brine_rate,brine_mass',co2_rate,co2_mass,brine_rate,brine_mass
    if(debug) write(*,*) 'chem',chem
    if(debug) write(*,*) 'logf',logf
    CALL getcwd(cwd)
    if(debug) WRITE(*,*) 'cwd ',TRIM(cwd)

    if (cwd(1:1) == '/') then
        dir='./carbonate/RScoeff/hydro_'
    else
        dir='.\carbonate\RScoeff\hydro_' 
    endif

    if(ithresh.eq.1) then
        nhout=10
    else
        nhout=12
    endif
    quad=.false.;mars=.false.
    do i=1,nhout
        if(rs_method(i)<4) quad=.true.
        if(rs_method(i)==4) mars=.true.
    end do
    if(chem(1)<0) scale=.false.
    if(iread==0) then
        !     write(*,*) 'reading coefficient files'
        do i=1,nhout
            if(rs_method(i)<4) then
                if(i>4) then
                    nrs=2
                else
                    nrs=1
                endif
                f=trim(dir)//trim(fname(ithresh))//'_RScoefficient'//trim(type(nrs))//'.dat' 
                ! for now, this gets called more times than necessary, but it doesn't hurt anything          
                call read_param(f,nrs,rs_method(i)) 
            else      
                if(logf(i).eq.1) then
                    f=trim(dir)//trim(fname(ithresh))//'_rsmlog'//trim(out_names(i))//'_pre'//'.h' 
                else 
                    f=trim(dir)//trim(fname(ithresh))//'_rsm'//trim(out_names(i))//'_pre'//'.h' 
                endif
                !
                !           if(debug) write(*,*) ithresh,f            
                call read_rsm( f, i )   ! need to add loop over pre and post so we read in 24 rows
                ! Flux is handled differently
                if(i.ne.2) then
                    f=trim(dir)//trim(fname(ithresh))//'_rsm'//trim(out_names(i))//'f_post'//'.h' 
                else
                    f=trim(dir)//trim(fname(ithresh))//'_rsm'//trim(out_names(i))//'_post'//'.h' 
                end if
                !           if(debug) write(*,*) ithresh,f            
                call read_rsm( f, i+12 )   ! need to add loop over pre and post so we read in 24 rows
            end if 
        end do 
        if(scale) then
            call read_geochem_param(ithresh)
        endif  
        iread=1
    end if
    ! load geochemical parameters that do not vary from well-to-well
    do k=13,20
       xp(k)=chem(k-12)
   end do
   ! read input variable matrices
   if(quad) then
    open(1,file='quad_input_variables')
    read(1,'(a)')       
    do i=1,nhout
        read(1,*) dum,dum1,(quad_inputs(i,nin),nin=1,12)
    end do
    close(1)
endif
! intialize at beginning of each 200 year run   
!    leak=.false.
if(time<=1) then
    time=abs(time)
    dis=.false.   ! this flag watches for pH plume disappearance during recovery
    leak=.false.
    stop_leak=.false.
    time_stop=0.
    app=.false.

    time_zero=200.
    ! loop to see if any wells need to be combined
    groups=0
    ng=0
    do i=1,n_leak
        groups(i)=i
    end do
    do i=1,n_leak
        if(groups(i)>0) then
            do m=i+1,n_leak
               r=sqrt((x(m)-x(i))**2+(y(m)-y(i))**2)
               if(groups(m)>0.and.r.lt.rmin) then
                groups(m)=-1*i
                groups(i)=i 
                ng=ng+1 
            endif
        end do
    endif
end do
!!  if(debug) write(*,*) ng,' wells will disappear ',(groups(i),i=1,n_leak)
endif

! combine leakage from closely spaced wells if necessary
! how are we combining chloride concentration?
do iw=1,n_leak
    if(groups(iw)>0) then       
        do k=iw+1,n_leak
            if(groups(k).eq.-1*iw) then
               co2_rate(iw)=co2_rate(iw)+co2_rate(k)
               co2_mass(iw)=co2_mass(iw)+co2_mass(k)
               brine_rate(iw)=brine_rate(iw)+brine_rate(k)
               brine_mass(iw)=brine_mass(iw)+brine_mass(k)
           endif
       end do
   endif
end do       

vol=0.;flux=0.;xplume=0;yplume=0    
! load aquifer input variables
input1(1)=aquifer(5);input1(2)=aquifer(6);input1(3)=aquifer(1);input1(4)=aquifer(2);input1(5)=aquifer(3)
input1(6)=aquifer(4)
! aquifer(5) is porosity - not used  
do iw=1,n_leak
    ! check for leak start/stop  
    if(.not.leak(iw).and.co2_rate(iw)>0) leak(iw)=.true.
    if(leak(iw).and.co2_rate(iw)<=0.) then
      stop_leak(iw)=.true. 
      if(debug) write(*,*) 'leak has stopped'
      time_stop(iw)=time
      leak(iw)=.false.
  endif   
  ! load leak variables
  input1(7)=cl(iw)
  input1(8)=time-time_stop(iw)
  input1(9)=co2_rate(iw)
  input1(10)=co2_mass(iw)
  input1(11)=brine_rate(iw)
  input1(12)=brine_mass(iw)   
  !        write(355,'(14f15.4)') (input1(k),k=1,12),time,time_stop(iw)

  ! only calculated for wells that have not been combined
  if(groups(iw)>=0) then
    if(scale) then
        xp(5)=cl(iw)
        ! scale trace metals and organics for geochem model         
        do k=6,12
            xp(k)=xp(5)*brine_ratios(k-5)
            xp(k)=log10(xp(k))
        end do  
        xp(5)=log10(xp(5)) 
        xp(1)=co2_mass(iw);xp(2)=co2_rate(iw)*1E-3;xp(3)=brine_mass(iw)
        xp(4)=brine_rate(iw)*1E-3 !! geochem model expects leak rates in kg/s - hydro ROM receiveds them in g/s
        xp(21)=time         
    endif     
    455        format(i3,19e10.3)
    ! do we want to filter for CO2_mass)>0 ??
    vol=0.0; outArgs=0.0; outlink=0.0; yout=0.0   ! not sure if we need this now
    do im=1,nhout  ! loop over outputs
        iuse=0
        if(stop_leak(iw)) then
            if(debug) write(*,*) 'in recovery'
            ! we are in recovery
            input1(13) = plume_final(iw,im)
            call mars_rsm(input1,y_result,im+nhout)   
            ! Flux (im==2) is handled differently than other post leak
            ! outputs
            if(im.eq.1) then
            endif
            if(im.ne.2) then  ! dont apply the correction to flux
                y_result=y_result*plume_final(iw,im)
            endif
        else

            ! we are pre-leak or during leak
            if (logf(im).eq.1) then
                call mars_rsm(inArgs,y_result,im)                           
                yresult = 10**yresult - 1
            else
                call mars_rsm(input1,y_result,im)  
                if(im.eq.12) then
                endif
            endif
        endif   
        if (y_result<0) y_result = 0   !  recovery ROM is sending negative #'s


        ! make sure the plume is not reappearing

        ! the following 2 blocks make sure once a ph or tds plume disappears, it stays away                        

        if(.not.app(iw,im).and.y_result>0) app(iw,im)=.true.

        if (app(iw,im).and.y_result<=0.) dis(iw,im)=.true.

        if(dis(iw,im)) y_result=0

        ! scale for chemistry
        if(y_result>0) then  
            xp(22)=y_result 
            if(scale) then
               call rs(xp,y_result,im+12)
               308            format(a15,1x,i4,25e12.4)                
           endif
       endif





       outlink(im)=y_result
   end do
   ! apply superposition
   if(input1(10).le.0.d0) then
    outlink(1) = 0.d0
    outlink(2) = 0.d0
    outlink(3) = 0.d0
    outlink(4) = 0.d0
endif
if(input1(12).le.0.d0) then
    outlink(5) = 0.d0
    outlink(6) = 0.d0
    outlink(7) = 0.d0
    outlink(8) = 0.d0
    outlink(9) = 0.d0
    outlink(10) = 0.d0
    outlink(11) = 0.d0
    outlink(12) = 0.d0
endif
do k=1,12                
    vol(k)=vol(k)+outlink(k)
end do

! save plume volumes if leak is still active   

if(leak(iw)) then

    do k=1,nhout

        plume_final(iw,k)=outlink(k)

    end do

endif
endif
344        format(a6,f5.0,i5,5f6.1,e10.2)
! next well
end do
if(debug) write(*,*) 'vol',vol
return
end
