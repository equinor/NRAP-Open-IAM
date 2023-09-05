    module rs_param
      parameter(nvar=24,nout=24)
      real(kind=8) :: coef(nvar,nout),cross(nvar,nvar,nout), &
      cubic(nvar,nvar,nvar,nout),quartic(nvar,nvar,nvar,nvar,nout)  
      real(kind=8) :: ymin(nout),ymax(nout),const(nout)
      integer :: ivar(nout),jvar(nout),kvar(nout),lvar(nout), iread
      real(kind=8) :: brine_ratios(7)
      data brine_ratios/1.71863E-07,9.66456E-08,1.71863E-08,1.216E-4,1.531E-7,5.433E-8,8.61E-8/      
    end module    
    subroutine read_geochem_param(it)
     use rs_param
     ! it =  1:MCL  2:THR
     ! nout is number of outputs
     ! nvar is number in input variables
     CHARACTER(len=255) :: cwd
     character outnam(nout/2)*6, linkmodel*48, linkfile*80, line*80
     data outnam/'phvol','flux','phlen','phwid','tdsvol','asvol','pbvol',&
     'cdvol','bavol','bzvol','npvol','pnvol'/
     integer :: model,i,j,k,l,indx,n,nhalf
     ! Zero coefficients
     !  const = 0.d+0
     nhalf=nout/2
     do n=nhalf+1,nout

      const(n)=0.
      do i = 1,nvar
        coef(i,n) = 0.d+0
        do j = 1,nvar
          cross(i,j,n) = 0.d+0
          do k = 1,nvar
            cubic(i,j,k,n) = 0.d+0
            do l = 1,nvar
              quartic(i,j,k,l,n) = 0.d+0
            end do
          end do
        end do
      end do
    end do
    !!    write(*,*) 'Diana coefficient file:', linkfile
    CALL getcwd(cwd)
    if (cwd(1:1) == '/') then
      if(it==1) linkmodel='./carbonate/RScoeff/chem_MCL'
      if(it==2) linkmodel='./carbonate/RScoeff/chem_THR'
    else
      if(it==1) linkmodel='.\carbonate\RScoeff\chem_MCL'
      if(it==2) linkmodel='.\carbonate\RScoeff\chem_THR'
    endif
    do n=nhalf+1,nout
      linkfile=trim(linkmodel)//'/quadratic_model_'//trim(outnam(n-nhalf))//'.txt'
      open(1,file=linkfile,status='OLD')
      ! Find coefficients in input file
      1 read(1,'(a)',err=666) line
      if (index(line,'Linear regression model') /= 0) model = 1
      if (index(line,'Quadratic regression model') /= 0) model = 2
      if (index(line,'Cubic regression model') /= 0) model = 3
      if (index(line,'Quartic regression model') /= 0) model = 4

      ! Read max/min values
      if (index(line,'Maximum') /= 0) then
        indx = index(line,'=')
        read(line(indx+1:),*,err=666) ymax(n),ymin(n)
      end if

      if (index(line,'coefficient') == 0) go to 1

      if ( model == 1 ) then
        ivar(n) = nvar
        jvar(n) = 0
        kvar(n) = 0
        lvar(n) = 0
      elseif ( model == 2 ) then
        ivar(n) = nvar
        jvar(n) = nvar
        kvar(n) = 0
        lvar(n) = 0
      elseif ( model == 3 ) then
        ivar(n) = nvar
        jvar(n) = nvar
        kvar(n) = nvar
        lvar(n) = 0
      elseif ( model == 4 ) then
        ivar(n) = nvar
        jvar(n) = nvar
        kvar(n) = nvar
        lvar(n) = nvar
      end if

      ! Read constant, if specified
      4 read(1,'(a)',err=666) line
      if (index(line,'Constant') /= 0) then
        indx = index(line,'=')
        read(line(indx+1:),*) const(n)
      else
        backspace(1)
      end if

      ! Read coefficients
      2 read(1,'(a)',err=666) line
      if (index(line,'Input') /= 0) then
        indx = index(line,'=')
        if ( model == 1 ) then
          read(line(8:indx-1),'(i4)') i
          j = 0
          k = 0
          l = 0
        elseif ( model == 2 ) then
          read(line(8:indx-1),'(2i4)') i,j
          k = 0
          l = 0
        elseif ( model == 3 ) then
          read(line(8:indx-1),'(3i4)') i,j,k
          l = 0
        elseif ( model == 4 ) then
          read(line(8:indx-1),'(4i4)') i,j,k,l
        end if
        if ( l == 0. ) then
          if ( k == 0. ) then
            if ( j == 0. ) then
              !         Read linear coefficients
              read(line(indx+1:),*) coef(i,n)
            else
              !         Read quadratic coefficients
              read(line(indx+1:),*) cross(i,j,n)
            end if
          else
            !       Read qubic coefficients
            read(line(indx+1:),*) cubic(i,j,k,n)
          end if
        else
          !     Read quartic coefficients
          read(line(indx+1:),*) quartic(i,j,k,l,n)
        end if
        go to 2
      end if                  

      close(1)
      ! read next coefficient file    
    end do
    return
    666  write(*,*) 'Error during read of geochemical ROM coefficient file'    
    write(*,*) linkfile
    stop
  end

