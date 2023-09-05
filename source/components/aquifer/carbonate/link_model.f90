subroutine rs (x,y,n)
! Calculates model values using polynomial fit to PSUADE response surface
! Supports linear, quadratic, cubic or quartic curve fits
        use rs_param
  implicit none
  real(kind=8) :: y,x(nvar),delta_max,rmax,delta,ratio
  integer :: i,j,k,l,n

        delta_max=1.e7;rmax=0.5
! n=output flag.  1: ph plume vol, 2: flux 3: x-dim ph plume  4: y-dimen ph plume
! 5: tds , 6: As 7: Pb 8: Cd  9: Ba  10: Org1  11: org2  12:org3
! n=1,12 hydro model
! n=13,24 geochem scaling factors
! geochem model is truncated to control behavior outside intended input param value ranges
! flux, xy dimensions will only be controlled by rmax
!
! first make sure geochem model returns 0 output for 0 input
! geochem model requires 22 parameters:  co2m,co2f,brine_mass,brinef, 7 solute concentrations, 8 chem parameters, time, plume vol
! premise is that we only use this for chem model now?
!	x(1)=x(1)*100
!	if(x(4)<=0.or.x(5)<=0.) then
!		y=0.
!		return
!	endif
    if(n>12.and.x(22).eq.0.) then
        y=0
        return
    endif
    y = const(n)
    do i = 1,ivar(n)
      y = y + x(i)*coef(i,n)
      do j = i,jvar(n)
        y = y + x(i)*x(j)*cross(i,j,n)
        do k = j,kvar(n)
          y = y + x(i)*x(j)*x(k)*cubic(i,j,k,n)
          do l = k,lvar(n)
            y = y + x(i)*x(j)*x(k)*x(l)*quartic(i,j,k,l,n)
          end do
        end do
      end do
    end do
    if ( y < ymin(n) ) y = ymin(n)
! truncate geochemistry ROM here; ideally it should be output variable specific criterion
    if(n.gt.12) then
        delta=abs(y-x(22))
        ratio=delta/x(22)
        if(ratio>rmax) then
            if(y>x(22)) then
                y=x(22)*(1.+rmax)
            else
                y=x(22)*(1.-rmax)
            endif
            return
        endif
        if(delta>delta_max) then
            if(y>x(22)) then
                y=x(22)+delta_max
            else
                y=x(22)-delta_max
            endif
            return          
        endif
    endif
!    write(*,*) 'in rs,',(x(k),k=1,5),y
!    if ( y > ymax(n) ) y = ymax(n)
!    if ( x(nvar) == 0. ) y = 0.d+0
return
end subroutine rs

