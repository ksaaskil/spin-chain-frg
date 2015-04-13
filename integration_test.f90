PROGRAM test
  USE nrtype
  USE qromomod, only: qromo,func_to_int
  USE integrators, only: midinf, midpnt, polint
  IMPLICIT NONE
  ! gfortran -o test nrtype.f90 nrutil.f90 nr.f90 integrators.f90 qromo.f90 integration_test.f90
  real(sp) :: a,b

  a=pi
  b=3*pi
  !print*, qromo(func_to_int,a,b,midinf)
  print*, qromo(func_to_int,a,b,midpnt)

contains
  function func_to_int2(x) result(value)
    use nrtype
    implicit none
    real(sp)::x(:)
    real(sp)::value(size(x))
    value=sin(x)
  end function func_to_int2

END PROGRAM test
