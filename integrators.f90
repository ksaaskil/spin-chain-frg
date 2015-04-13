MODULE integrators

CONTAINS

  SUBROUTINE midinf(funk,aa,bb,s,n)
    USE nrtype; USE nrutil, ONLY : arth,assert
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: aa,bb
    REAL(SP), INTENT(INOUT) :: s
    INTEGER(I4B), INTENT(IN) :: n
    INTERFACE
       FUNCTION funk(x)
         USE nrtype
         REAL(SP), DIMENSION(:), INTENT(IN) :: x
         REAL(SP), DIMENSION(size(x)) :: funk
       END FUNCTION funk
    END INTERFACE
    REAL(SP) :: a,b,del
    INTEGER(I4B) :: it
    REAL(SP), DIMENSION(2*3**(n-2)) :: x
    call assert(aa*bb > 0.0, 'midinf args')
    b=1.0_sp/aa
    a=1.0_sp/bb
    if (n == 1) then
       s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
    else
       it=3**(n-2)
       del=(b-a)/(3.0_sp*it)
       x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
       x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
       s=s/3.0_sp+del*sum(func(x))
    end if
  CONTAINS
    !BL
    FUNCTION func(x)
      REAL(SP), DIMENSION(:), INTENT(IN) :: x
      REAL(SP), DIMENSION(size(x)) :: func
      func=funk(1.0_sp/x)/x**2
    END FUNCTION func
  END SUBROUTINE midinf
  
  SUBROUTINE midpnt(func,a,b,s,n)
    USE nrtype; USE nrutil, ONLY : arth
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: a,b
    REAL(SP), INTENT(INOUT) :: s
    INTEGER(I4B), INTENT(IN) :: n
    INTERFACE
       FUNCTION func(x)
         USE nrtype
         REAL(SP), DIMENSION(:), INTENT(IN) :: x
         REAL(SP), DIMENSION(size(x)) :: func
       END FUNCTION func
    END INTERFACE
    REAL(SP) :: del
    INTEGER(I4B) :: it
    REAL(SP), DIMENSION(2*3**(n-2)) :: x
    if (n == 1) then
       s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
    else
       it=3**(n-2)
       del=(b-a)/(3.0_sp*it)
       x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
       x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
       s=s/3.0_sp+del*sum(func(x))
    end if
  END SUBROUTINE midpnt
  
  
  SUBROUTINE midsql(funk,aa,bb,s,n)
    USE nrtype; USE nrutil, ONLY : arth
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: aa,bb
    REAL(SP), INTENT(INOUT) :: s
    INTEGER(I4B), INTENT(IN) :: n
    INTERFACE
       FUNCTION funk(x)
         USE nrtype
         REAL(SP), DIMENSION(:), INTENT(IN) :: x
         REAL(SP), DIMENSION(size(x)) :: funk
       END FUNCTION funk
    END INTERFACE
    REAL(SP) :: a,b,del
    INTEGER(I4B) :: it
    REAL(SP), DIMENSION(2*3**(n-2)) :: x
    b=sqrt(bb-aa)
    a=0.0
    if (n == 1) then
       s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
    else
       it=3**(n-2)
       del=(b-a)/(3.0_sp*it)
       x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
       x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
       s=s/3.0_sp+del*sum(func(x))
    end if
  CONTAINS
  !  !BL
    FUNCTION func(x)
      REAL(SP), DIMENSION(:), INTENT(IN) :: x
      REAL(SP), DIMENSION(size(x)) :: func
      func=2.0_sp*x*funk(aa+x**2)
    END FUNCTION func
  END SUBROUTINE midsql
  
  SUBROUTINE midsqu(funk,aa,bb,s,n)
    USE nrtype; USE nrutil, ONLY : arth
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: aa,bb
    REAL(SP), INTENT(INOUT) :: s
    INTEGER(I4B), INTENT(IN) :: n
    INTERFACE
       FUNCTION funk(x)
         USE nrtype
         REAL(SP), DIMENSION(:), INTENT(IN) :: x
         REAL(SP), DIMENSION(size(x)) :: funk
       END FUNCTION funk
    END INTERFACE
    REAL(SP) :: a,b,del
    INTEGER(I4B) :: it
    REAL(SP), DIMENSION(2*3**(n-2)) :: x
    b=sqrt(bb-aa)
    a=0.0
    if (n == 1) then
       s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
    else
       it=3**(n-2)
       del=(b-a)/(3.0_sp*it)
       x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
       x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
       s=s/3.0_sp+del*sum(func(x))
    end if
  CONTAINS
    !BL
   FUNCTION func(x)
      REAL(SP), DIMENSION(:), INTENT(IN) :: x
      REAL(SP), DIMENSION(size(x)) :: func
      func=2.0_sp*x*funk(bb-x**2)
    END FUNCTION func
  END SUBROUTINE midsqu
  
  SUBROUTINE midexp(funk,aa,bb,s,n)
    USE nrtype; USE nrutil, ONLY : arth
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: aa,bb
    REAL(SP), INTENT(INOUT) :: s
    INTEGER(I4B), INTENT(IN) :: n
    INTERFACE
       FUNCTION funk(x)
         USE nrtype
         REAL(SP), DIMENSION(:), INTENT(IN) :: x
         REAL(SP), DIMENSION(size(x)) :: funk
       END FUNCTION funk
    END INTERFACE
    REAL(SP) :: a,b,del
    INTEGER(I4B) :: it
    REAL(SP), DIMENSION(2*3**(n-2)) :: x
    b=exp(-aa)
    a=0.0
    if (n == 1) then
       s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
    else
       it=3**(n-2)
       del=(b-a)/(3.0_sp*it)
       x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
       x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
       s=s/3.0_sp+del*sum(func(x))
    end if
  CONTAINS
    !BL
    FUNCTION func(x)
      REAL(SP), DIMENSION(:), INTENT(IN) :: x
      REAL(SP), DIMENSION(size(x)) :: func
      func=funk(-log(x))/x
    END FUNCTION func
  END SUBROUTINE midexp

  SUBROUTINE polint(xa,ya,x,y,dy)
    USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
    REAL(SP), INTENT(IN) :: x
    REAL(SP), INTENT(OUT) :: y,dy
    INTEGER(I4B) :: m,n,ns
    REAL(SP), DIMENSION(size(xa)) :: c,d,den,ho
    n=assert_eq(size(xa),size(ya),'polint')
    c=ya
    d=ya
    ho=xa-x
    ns=iminloc(abs(x-xa))
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(den(1:n-m) == 0.0)) &
            call nrerror('polint: calculation failure')
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
  END SUBROUTINE polint
  
  
END MODULE integrators
