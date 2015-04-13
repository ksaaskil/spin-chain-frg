MODULE rk4step
PRIVATE
PUBLIC :: rk4
CONTAINS
  SUBROUTINE rk4(y,dydx,x,nvar,h,yout,derivs)
    ! Global variables not needed
    IMPLICIT NONE
    !EXTERNAL derivs
    INTERFACE
       SUBROUTINE derivs(x,y,dydx,nvar)
         !USE global_variables, ONLY: nvar
         !INTEGER :: nf,sites,freqs,nvar
         !COMMON /indices/ nf,sites,freqs,nvar
         INTEGER, INTENT(IN) :: nvar
         REAL(KIND(1.d0)),INTENT(IN) :: y(nvar),x
         REAL(KIND(1.d0)),INTENT(OUT):: dydx(nvar)
       END SUBROUTINE derivs
    END INTERFACE
    INTEGER, PARAMETER:: NMAX=1e5
    INTEGER,INTENT(IN) :: nvar
    !,nf,sites,freqs
    !COMMON /indices/ nf,sites,freqs,nvar
    REAL(KIND(1.d0)), INTENT(IN):: h,x,dydx(nvar),y(nvar)
    REAL(KIND(1.d0)), INTENT(OUT):: yout(nvar)
    
    !Given values for the variables y(l:n) and their derivatives dydx(l:n) known at x, use
    !the fourth-order Runge-Kutta method to advance the solution over an interval h and return
    !the incremented variables as yout(l:n), which need not be a distinct array from y. The
    !user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
    INTEGER :: i,n
    REAL(KIND(1.d0)):: h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
    IF (nvar>NMAX) THEN
       write(*,*) 'error, NMAX too small!',nvar,NMAX
       stop
    END IF

    n=nvar
    hh=h*0.5
    h6=h/6.
    xh=x+hh
    do i=1,n !First step.
       yt(i)=y(i)+hh*dydx(i)
    enddo
    call derivs(xh,yt,dyt,nvar) !Second step,
    do i=1,n
       yt(i)=y(i)+hh*dyt(i)
    enddo
    call derivs(xh,yt,dym,nvar) !Third step,
    do i=1,n
       yt(i)=y(i)+h*dym(i)
       dym(i)=dyt(i)+dym(i)
    enddo
    call derivs(x+h,yt,dyt,nvar) !Fourth step,
    do i=1,n !Accumulate increments with proper weights.
       yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
    enddo

  END SUBROUTINE RK4
END MODULE rk4step
