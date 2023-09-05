! to compile I use gfortran 6.4.0
! 1st step:
! gfortran -c iam_quad_eq_fortran_model.f90     or
! x86_64-w64-mingw32-gfortran -c iam_quad_eq_fortran_model.f90
! 2nd step: in Windows
! gfortran -shared -o quad_eq_fun.dll iam_quad_eq_fortran_model.o   or
! x86_64-w64-mingw32-gfortran -shared -o quad_eq_fun.dll iam_quad_eq_fortran_model.o
! or in Mac
! gfortran -dynamiclib -o quad_eq_fun.dylib iam_quad_eq_fortran_model.o
subroutine quad_eq_fun(a, b, c, x, N, x1, x2, y, flag) bind(C, name='quad_eq_fun')
implicit none
real*8, intent(in) :: a
real*8, intent(in) :: b
real*8, intent(in) :: c
real*8 :: D
real*8, intent(out) :: x1
real*8, intent(out) :: x2
integer, intent(in) :: N
integer, intent(out) :: flag
integer :: i
real*8, parameter :: epsil = 1d-20
real*8, dimension(1:N), intent(in) :: x
real*8, dimension(1:N), intent(out) :: y

D = b*b - 4.0*a*c

if (abs(D) > epsil) then
    if (D .GT. 0.0) then
        x1 = (-b + sqrt(D))/(2.0*a)
        x2 = (-b - sqrt(D))/(2.0*a)
        flag = 1                 ! indicates two real distinct roots
    else
        x1 = -b/(2.0*a)          ! returns real part of the roots
        x2 = sqrt(-D)/(2.0*a)    ! returns complex part of the roots
        flag = 3
    endif
else
    x1 = -b/(2.0*a)
    x2 = x1
    flag = 2
endif

do i = 1, N
    y(i) = a*x(i)*x(i) + b*x(i) + c   ! calculates parabola y-values for the entered parameters a, b, c and x-values
end do

return
end subroutine quad_eq_fun
