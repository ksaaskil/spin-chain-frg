PROGRAM main
  USE equation, ONLY: derivs
  USE rk4step, ONLY: rk4
  IMPLICIT NONE
  INTEGER, PARAMETER:: nsteps=20, nvar=2
  INTEGER, PARAMETER:: NMAX=50, NSTPMX=200 ! Maximum number of functions, max number of variables
  REAL(KIND(1.d0)):: x1=0.d0,x2=1.d0,vstart(nvar),xx(NSTPMX),y(NMAX,NSTPMX) 
  COMMON /path/ xx,y !Storage of results.
  INTEGER i,k
  REAL(KIND(1.d0)):: h,x,dv(NMAX),v(NMAX)
  OPEN(55,file='testi.dat')
  
  vstart=1.d0
  do i=1,nvar !Load starting values.
     v(i)=vstart(i)
     y(i,1)=v(i)
  enddo
  xx(1)=x1
  write(55,*) xx(1),y(2,1)
  x=x1
  h=(x2-x1)/nsteps

  do k=1,nsteps !Take nstep steps
     call derivs(x,v,dv,nvar) ! Calculate the derivatives dv, used in the first step of RK
     call rk4(v,dv,x,nvar,h,v,derivs) ! The Runge-Kutta step
     !if(x+h.eq.x)pause 'stepsizenot significant in rkdumb'
     x=x+h
     xx(k+1)=x !Store intermediate steps.
     do i=1,nvar ! Loop over the components of y
        y(i,k+1)=v(i)
     enddo
     write(55,*) xx(k+1),y(2,k+1)
  enddo

  close(55)


END PROGRAM main
