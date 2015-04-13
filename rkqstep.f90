MODULE rkqstep

CONTAINS
  SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
    USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
    USE rkckmod, ONLY : rkck
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
    REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
    REAL(SP), INTENT(INOUT) :: x
    REAL(SP), INTENT(IN) :: htry,eps
    REAL(SP), INTENT(OUT) :: hdid,hnext
    INTERFACE
       SUBROUTINE derivs(x,y,dydx)
         USE nrtype
         IMPLICIT NONE
         REAL(SP), INTENT(IN) :: x
         REAL(SP), DIMENSION(:), INTENT(IN) :: y
         REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
       END SUBROUTINE derivs
    END INTERFACE
    !Fifth order Runge-Kutta step with monitoring of local truncation error to ensure accuracy
    !and adjust stepsize. Input are the dependent variable vector y and its derivative dydx at
    !the starting value of the independent variable x. Also input are the stepsize to be attempted
    !htry, the required accuracy eps, and the vector yscal against which the error is scaled. y,
    !dydx, and yscal are all of the same length. On output, y and x are replaced by their new
    !values, hdid is the stepsize that was actually accomplished, and hnext is the estimated
    !next stepsize. derivs is the user-supplied subroutine that computes the right-hand-side
    !derivatives.
    INTEGER(I4B) :: ndum, errloc(1)
    REAL(SP) :: errmax,h,htemp,xnew
    REAL(SP), DIMENSION(size(y)) :: yerr,ytemp
    REAL(SP), PARAMETER :: SAFETY=0.9_sp,PGROW=-0.2_sp,PSHRNK=-0.25_sp,&
         ERRCON=1.89e-4
    REAL(SP), PARAMETER :: TINY=1e-15_sp
    ndum=assert_eq(size(y),size(dydx),size(yscal),'rkqs')
    h=htry
    !write(*,*) h
    do
       call rkck(y,dydx,x,h,ytemp,yerr,derivs)
       !write(*,*) maxval(abs(dydx))
       errmax=maxval(abs(yerr(:)/yscal(:)))/eps
       !errmax=sqrt(sum(yerr(:)**2))/sqrt(sum(yscal(:)**2))/eps
       write(*,*) h,maxval(yscal),errmax
       !errloc=maxloc(abs(yerr(:)/yscal(:)))
       !write(*,*) h,maxval(yscal),errmax,errloc,yscal(errloc)
       
       if (errmax <= 1.0) exit
       !write(*,*) 'here'
       htemp=SAFETY*h*(errmax**PSHRNK)
       h=sign(max(abs(htemp),0.1_sp*abs(h)),h)
       xnew=x+h
       !write(*,*) h
       !if (xnew == x) call nrerror('stepsize underflow in rkqs')
       if (abs(xnew-x)<TINY) call nrerror('stepsize underflow in rkqs')
    end do
    if (errmax > ERRCON) then
       hnext=SAFETY*h*(errmax**PGROW)
    else
       hnext=5.0_sp*h
    end if
    hdid=h
    x=x+h
    y(:)=ytemp(:)
    hnext=min(hnext,1e-1_sp)
    
    write(*,*) 'hdid=',hdid,'hnext=',hnext
  END SUBROUTINE rkqs
END MODULE RKQSTEP
