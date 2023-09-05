        subroutine read_param(ff,nrs,model)
        use rs_param    
        implicit none
        integer i,ii,k,k2,k3,jj,nrs,np,nout2,nn,model,n
        character skip*10,ff*80
        logical debug
! nrs = 1:pH-related variables  2:TDS-related variables 
        debug=.true.
        ymin=0;ymax=1e10
            if(nrs==1) then
                nn=1;nout2=4
                np=5
            else
                nn=5;nout2=12   
                np=6
            endif
            do n=nn,nout2
                if ( model == 1 ) then
                    ivar(n) = np
                    jvar(n) = 0
                    kvar(n) = 0
                    lvar(n) = 0
                elseif ( model == 2 ) then
                    ivar(n) = np
                    jvar(n) = np
                    kvar(n) = 0
                    lvar(n) = 0
                elseif ( model == 3 ) then
                     ivar(n) = np
                     jvar(n) = np
                     kvar(n) = np
                     lvar(n) = 0
                end if  
             end do             
        open(3,file=ff) 
            call progress(3,2)  
! read linear parameters
            read(3,*) skip,(const(n),n=nn,Nout2)

            !write(*,*) const(2)
            do i=1,Np
                 read(3,*) skip,k, (coef(i,n),n=nn,Nout2)
            end do
         if(model==1) return
! read quadratic parameters         
        call progress(3,4)  
         read(3,*)skip,(const(n),n=nn,Nout2)
            do i=1,Np
                read(3,*) skip,k, (coef(i,n),n=nn,Nout2)
             end do        
            do i=1,Np                    
            do ii=1,Np
            read(3,*)skip, k,k2,(cross(i,ii,n),n=nn,Nout2)
            end do
            end do
        if(model==2) return      
        call progress(3,4)
! read cubic parameters        
         read(3,*)skip,(const(n),n=nn,Nout2)
            do i=1,Np
                read(3,*) skip,k, (coef(i,n),n=nn,Nout2)
             end do        
            do i=1,Np                    
                do ii=1,Np
                    read(3,*)skip, k,k2,(cross(i,ii,n),n=nn,Nout2)
                end do
            end do
             do k=1,Np
            do k2=1,Np
                do k3=1,Np
                    read(3,*)skip,i,ii,jj,(cubic(k,k2,k3,n),n=nn,Nout2)
                 end do
             end do
             end do  
    close(3)              
    return
    end                    
    subroutine progress(iu,n)
    do i=1,n
    read(iu,'(a)')
    end do
    return
    end
!


