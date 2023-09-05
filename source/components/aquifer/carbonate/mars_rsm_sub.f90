        MODULE mars_coefs
            integer, parameter :: nrmax = 200 ! Max number of rsms
            integer, parameter :: ntmax = 200 ! Max number of terms
            integer, parameter :: npmax = 50  ! Max number of inputs
            integer, parameter :: nsmax = 200 ! Max number of selected terms
            integer :: ntm(nrmax), npm(nrmax), nsm(nrmax)
            real(kind=8) :: cuts(nrmax,ntmax,npmax), c(nrmax,nsmax) 
            real(kind=8) :: mins(nrmax,npmax), maxs(nrmax,npmax)
            real(kind=8) :: trans(nrmax,3)
            integer :: dirs(nrmax,ntmax,npmax), terms(nrmax,nsmax)
            integer :: mars_inputs(nrmax,13)
            character :: in_names(13)*15,out_names(12)*15
            integer :: inputs(nrmax,npmax)
       data in_names/'thick','grad','var','corr','anis','perm','cl','time','co2','co2_c','br','br_c','plume'/
      data out_names/'pH','flux','dx','dy','TDS','as','pb','cd','ba','benz','nap','phenol'/

         END MODULE mars_coefs

      subroutine mars_rsm(x_all, y, nr)
        use mars_coefs
        implicit none
        integer i, j, nr
        real*8, dimension(*) :: x_all
        real*8 y,x(npmax)
        real*8, allocatable :: bx(:)
        real*8 temp 
        character outn*15
 
        ! Allocate basis function array
        allocate( bx(nsm(nr)) )
		! load input parameter array
		do i=1,npm(nr)
			x(i)=x_all(inputs(nr,i))
!			if(nr<=12) then
!			write(*,'(2i5,2(a15,1x),e15.3)') nr,i,out_names(nr),in_names(inputs(nr,i)),x(i)
!			else
!			write(*,'(2i5,3(a15,1x),e15.3)') nr,i,'post',out_names(nr-12),in_names(inputs(nr,i)),x(i)	
!			endif		
		end do		

        ! Constrain inputs to be within bounds
        do i =1,npm(nr)
            if( x(i) < mins(nr,i) ) then
            	if(nr<13) then
            		outn=out_names(nr)
            	else
            		outn=trim(out_names(nr-12))//'.post '
            	endif	
!             	write(222,333) outn,in_names(inputs(nr,i)),x(i),' truncated at min to ',mins(nr,i)
               x(i) = mins(nr,i)

            endif
            if( x(i) > maxs(nr,i) ) then
            	if(nr<13) then
            		outn=out_names(nr)
            	else
            		outn=trim(out_names(nr-12))//'.post '
            	endif	
!             	if(i.eq.7) write(222,333)outn,in_names(inputs(nr,i)),x(i),' truncated at max to ',maxs(nr,i)
              x(i) = maxs(nr,i)
            endif
             
        end do
 333	format(2a15,e15.2,a20,e15.2)
        ! Compute
        y = 0.0
        do i = 1,nsm(nr)
            bx(i) = 1
            do j = 1,npm(nr)
                !print *, x(j)
                if(dirs(nr,terms(nr,i),j) == 2) then 
                    bx(i) = bx(i) * x(j)
                else if(dirs(nr,terms(nr,i),j)==-1 .or. &
                                 dirs(nr,terms(nr,i),j)==1) then
                    temp = dirs(nr,terms(nr,i),j)*(x(j)- &
                                 cuts(nr,terms(nr,i),j))
                    if(temp>0) then 
                        bx(i) = bx(i) * temp
                    else 
                        bx(i) = 0
                    end if
               end if
           end do
           y = y + bx(i) * c(nr,i)
       end do
       ! Apply back transformations
       y = y + trans(nr,1)
       if(trans(nr,2).ne.0) y = 10**y
       y = y + trans(nr,3)

      end subroutine

      subroutine read_rsm( f, ri )
        use mars_coefs
        character(len=132) :: line
        character(len=80) :: f
        character(len=32) :: dum
        integer           :: success, ri, indx, prev, beginning
        integer           :: i, j, k, intval
        real(kind=8)      :: value, scl
		logical           :: debug
		debug=.false.
        open(1,file=trim(f),status='OLD')
        read(1,*,err=666)
        read(1,*,err=666)
        read(1,*,err=666)
        read(1,*,err=666) dum, ntm(ri)
        read(1,'(A)',err=666)
        read(1,*,err=666) dum, npm(ri)
        read(1,'(A)',err=666)
        read(1,*,err=666) dum, nsm(ri)
        ! Read in cuts
        read(1,'(A)',err=666)
        do i=1,npm(ri)
            read(1,*,err=666) dum, scl
            k = 1
            do
                read(1,'(A)',iostat=success) line
                prev      = 1 
                beginning = 1 

                do j=1,len(line)

                    indx = index('-0123456789.eE', line(j:j))

                    ! store value when you have reached a blank (or any other 
                    ! non-real number character)
                    if (indx.eq.0 .and. prev.gt.0) then
                       read(line(beginning:j-1), *) value
                       cuts(ri,k,i) = value
                       k = k + 1
                    else if (indx.gt.0 .and. prev.eq.0) then
                       beginning = j 
                    end if

                    prev = indx
                end do
                if (k .gt. ntm(ri)) EXIT
            end do
            if (scl .ne. 0.) then
                k = 1
                do
                    read(1,'(A)',iostat=success) line
                    prev      = 1 
                    beginning = 1 

                    do j=1,len(line)

                        indx = index('-0123456789.eE', line(j:j))

                        ! store value when you have reached a blank (or any other 
                        ! non-real number character)
                        if (indx.eq.0 .and. prev.gt.0) then
                           read(line(beginning:j-1), *) value
                           cuts(ri,k,i) = (cuts(ri,k,i) + value) / scl
                           k = k + 1
                           !print *, cuts(ri,k,i)
                        else if (indx.gt.0 .and. prev.eq.0) then
                           beginning = j 
                        end if

                        prev = indx
                    end do
                    if (k .gt. ntm(ri)) EXIT
                end do
            end if
        end do

        ! Read in dirs
        read(1,'(A)',err=666)
        do i=1,npm(ri)
            k = 1
            do
                read(1,'(A)',iostat=success,err=666) line
                prev      = 1 
                beginning = 1 
                do j=1,len(line)
                    indx = index('-0123456789', line(j:j))
                    if (indx.eq.0 .and. prev.gt.0) then
                       read(line(beginning:j-1), *) intval
                       dirs(ri,k,i) = intval
                       k = k + 1
                       !print *, dirs(ri,k,i)
                    else if (indx.gt.0 .and. prev.eq.0) then
                       beginning = j 
                    end if
                    prev = indx
                end do
                if (k .gt. ntm(ri)) EXIT
            end do
        end do

        ! Read in coefficients
        read(1,'(A)')
        k = 1
        do
            read(1,'(A)',iostat=success,err=666) line
            prev      = 1 
            beginning = 1 
            do j=1,len(line)
                indx = index('-0123456789.eE', line(j:j))
                if (indx.eq.0 .and. prev.gt.0) then
                   read(line(beginning:j-1), *) value
                   c(ri,k) = value
                   !print *, c(ri,k)
                   k = k + 1
                else if (indx.gt.0 .and. prev.eq.0) then
                   beginning = j 
                end if
                prev = indx
            end do
            if (k .gt. nsm(ri)) EXIT
        end do

        ! Read in selected terms
        read(1,'(A)')
        k = 1
        do
            read(1,'(A)',iostat=success,err=666) line
            prev      = 1 
            beginning = 1 
            do j=1,len(line)
                indx = index('-0123456789', line(j:j))
                if (indx.eq.0 .and. prev.gt.0) then
                   read(line(beginning:j-1), *) intval
                   terms(ri,k) = intval
                   k = k + 1
                   !print *, terms(ri,k)
                else if (indx.gt.0 .and. prev.eq.0) then
                   beginning = j 
                end if
                prev = indx
            end do
            if (k .gt. nsm(ri)) EXIT
        end do

        ! Read in mins
        read(1,'(A)')
        k = 1
        do
            read(1,'(A)',iostat=success,err=666) line
            prev      = 1 
            beginning = 1 
            do j=1,len(line)
                indx = index('-0123456789.eE', line(j:j))
                if (indx.eq.0 .and. prev.gt.0) then
                   read(line(beginning:j-1), *) value
                   mins(ri,k) = value
                   !print *, value, mins(ri,k)
                   k = k + 1
                else if (indx.gt.0 .and. prev.eq.0) then
                   beginning = j 
                end if
                prev = indx
            end do
            if (k .gt. npm(ri)) EXIT
        end do

        ! Read in maxs
        read(1,'(A)')
        k = 1
        do
            read(1,'(A)',iostat=success,err=666) line
            prev      = 1 
            beginning = 1 
            do j=1,len(line)
                indx = index('-0123456789.eE', line(j:j))
                if (indx.eq.0 .and. prev.gt.0) then
                   read(line(beginning:j-1), *) value
                   maxs(ri,k) = value
                   !print *, value, maxs(ri,k), ri, k
                   k = k + 1
                else if (indx.gt.0 .and. prev.eq.0) then
                   beginning = j 
                end if
                prev = indx
            end do
            if (k .gt. npm(ri)) EXIT
        end do

        ! Read in input flags (flags indicating which inputs are used)
        read(1,'(A)',err=666)
        k = 1
        do
            read(1,'(A)',iostat=success,err=666) line
             prev      = 1 
            beginning = 1 
            do j=1,len(line)
                indx = index('-0123456789', line(j:j))
                if (indx.eq.0 .and. prev.gt.0) then
                   read(line(beginning:j-1), *) intval
                   mars_inputs(ri,k) = intval
                   !print *, ri,k,mars_inputs(ri,k)
                   k = k + 1
                else if (indx.gt.0 .and. prev.eq.0) then
                   beginning = j 
                end if
                prev = indx
            end do
            if (k .gt. 13) EXIT
        end do
        k=0
		do j=1,13
		if(mars_inputs(ri,j)>0) then
			k=k+1
			inputs(ri,k)=j
		endif
		end do
		if(debug) then
		if(ri.le.12) then
		write(*,'(a15,i3,10a15)') out_names(ri),k,(in_names(inputs(ri,m)),m=1,npm(ri))	
		else
		write(*,'(a15,i3,10a15)') 'post.'//out_names(ri-12),k,(in_names(inputs(ri,m)),m=1,npm(ri))	
		endif		
		endif	
        ! Read in output tranformation flags/values
        read(1,'(A)',err=666)
        k = 1
        do
            read(1,'(A)',iostat=success,err=666) line
            prev      = 1 
            beginning = 1 
            do j=1,len(line)
                indx = index('-0123456789', line(j:j))
                if (indx.eq.0 .and. prev.gt.0) then
                   read(line(beginning:j-1), *) intval
                   trans(ri,k) = intval
                   !print *, trans(ri,k)
                   k = k + 1
                else if (indx.gt.0 .and. prev.eq.0) then
                   beginning = j 
                end if
                prev = indx
            end do
            if (k .gt. 3) EXIT
        end do
        close(1)
		return
 666	write(*,*) 'error reading hydrologic ROM coefficient files'
 		write(*,*) trim(f)
 		stop
      end subroutine
