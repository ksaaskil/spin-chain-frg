MODULE odemod
  USE nrtype

CONTAINS 

  SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
  USE nrtype 
  USE nrutil, ONLY : nrerror,reallocate
  USE ode_path
  USE globalvariables, ONLY: nvar, om0, aa
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: ystart
  REAL(SP), INTENT(IN) :: x1,x2,eps,h1,hmin
  INTERFACE
     SUBROUTINE derivs(x,y,dydx)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: x
       REAL(SP), DIMENSION(:), INTENT(IN) :: y
       REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
     END SUBROUTINE derivs
     !BL
     SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
       USE nrtype
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
     END SUBROUTINE rkqs
  END INTERFACE
  REAL(SP), PARAMETER :: TINY=1.0e-30_sp
  INTEGER(I4B), PARAMETER :: MAXSTP=1023
  ! Runge-Kutta driver with adaptive stepsize control. Integrate the array of starting values
  ! ystart from x1 to x2 with accuracy eps, storing intermediate results in the module
  ! variables in ode path. h1 should be set as a guessed ﬁrst stepsize, hmin as the minimum
  ! allowed stepsize (can be zero). On output ystart is replaced by values at the end of the
  ! integration interval. derivs is the user-supplied subroutine for calculating the right-hand-
  ! side derivative, while rkqs is the name of the stepper routine to be used.
  INTEGER(I4B) :: nstp,i
  REAL(SP) :: h,hdid,hnext,x,xsav
  REAL(SP), DIMENSION(size(ystart)) :: dydx,y,yscal
  x=x1
  h=sign(h1,x2-x1)
  nok=0
  nbad=0
  kount=0
  y(:)=ystart(:)
  nullify(xp,yp)
  if (save_steps) then
     xsav=x-2.0_sp*dxsav
     allocate(xp(256))
     allocate(yp(size(ystart),size(xp)))
  end if
  do nstp=1,MAXSTP+1
     write(*,*) 'step=', nstp,'x=',x
     if (nstp==10) then
        do i=1,size(y)
           write(456,*) y(i)
        end do
        print*,i,size(y),nvar
        print*, 'Wrote data to file fort.456.'
     end if
     call derivs(x,y,dydx)
     !write(*,*) 'max(abs(dv))=',maxval(abs(dydx))
     !yscal(:)=abs(y(:))+abs(h*dydx(:))+1e-2_sp
     !yscal(:)=abs(y(:))
     !yscal=maxval(abs(y))
     ! yscal(:)=abs(h*dydx(:))+TINY
     do i=1,size(y)
        !yscal(i)=max(abs(y(i))+abs(h*dydx(i)),1e-3_sp)
        !yscal(i)=min(abs(y(i)),1._sp)
        yscal(i)=aa**(i-1)*om0*min(abs(y(i)),1._sp)
     end do
     if (save_steps .and. (abs(x-xsav) > abs(dxsav))) &
          call save_a_step
     if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x
     call rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)
     if (hdid == h) then
        nok=nok+1
     else
        nbad=nbad+1
     end if
     if ((x-x2)*(x2-x1) >= 0.0) then
        ystart(:)=y(:)
        if (save_steps) call save_a_step
        !STOP
        RETURN
     end if
     if (abs(hnext) < hmin) THEN
        write(*,*) 'stepsize smaller than minimum in odeint, exiting...'
        return
        !call nrerror('stepsize smaller than minimum in odeint')
     end if
     if (nstp==maxstp) then
        write(*,*) 'Maximum number of steps reached, exiting...'
        return
     end if
       h=hnext
       
    end do
    call nrerror('too many steps in odeint')
    
    CONTAINS
      !BL
      SUBROUTINE save_a_step
        kount=kount+1
        if (kount > size(xp)) then
           xp=>reallocate(xp,2*size(xp))
           yp=>reallocate(yp,size(yp,1),size(xp))
        end if
        xp(kount)=x
        yp(:,kount)=y(:)
        xsav=x
      END SUBROUTINE save_a_step
      
    END SUBROUTINE odeint
  END MODULE ODEMOD
