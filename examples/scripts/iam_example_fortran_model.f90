! to compile I use gfortran 4.9.2
! 1st step:
! gfortran -c iam_example_fortran_model.f90             or
! x86_64-w64-mingw32-gfortran -c iam_example_fortran_model.f90
! 2nd step: in Windows
! gfortran -shared -o example_rom.dll iam_example_fortran_model.o   or
! x86_64-w64-mingw32-gfortran -shared -o example_rom.dll iam_example_fortran_model.o
! or in Mac
! gfortran -dynamiclib -o example_rom.dylib iam_example_fortran_model.o (MacOS)
subroutine example_rom(var1, array1, N, out_var1, out_array1) bind(C, name='example_rom')
implicit none
real*8, intent(in) :: var1
real*8, intent(out) :: out_var1
integer, intent(in) :: N
integer :: i
real*8, dimension(1:N), intent(in) :: array1
real*8, dimension(1:N), intent(out) :: out_array1

if (var1>5) then
   out_var1 = var1 + 0.5
else
   out_var1 = var1 - 0.5
endif

do i = 1, N
  out_array1(i) = array1(i)*1000.0
end do

return
end subroutine example_rom
