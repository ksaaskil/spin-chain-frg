PROGRAM init
USE global_variables, ONLY : nf,sites,freqs,hf
USE nrutil, ONLY: assert_eq
USE nrtype
IMPLICIT NONE
INTEGER ::nsteps=200
INTEGER,PARAMETER :: nvar=nf+2*sites*freqs
REAL(KIND(1.d0)) :: omv(nf),x1=0.d0,x2=1.d0
INTEGER:: i,j,k,lask
REAL(KIND(1.d0)) :: gamma(nf)=0.d0,vfs(freqs,sites)=0.d0,vfd(freqs,sites)=0.d0 ! Self-energy + vertex functions in matrix form
REAL(KIND(1.d0)) :: vfsvec(sites*freqs), vfdvec(sites*freqs) ! vertex functions in vector form
REAL(KIND(1.d0)) :: y(nvar)=0.d0,y0(nvar)=0.d0,yout(nvar)=0.d0
INTEGER :: ind(freqs,3),ind2(nf,nf,nf)

lask=0
DO i=1,nf
   DO j=1,nf
      DO k=1,nf
         ind(lask,1)=i  ! indexing of frequencies
         ind(lask,2)=j
         ind(lask,3)=k
         ind2(i,j,k)=lask ! inverse indexing
      END DO
   END DO
END DO

write(*,*) 'Frequency step hf is ', hf

! Initial values
do i=1,freqs
   vfs(i,2)=.25d0 ! The nearest neighbor spin interaction
END DO

CALL compiley(y0,gamma,vfs,vfd,nvar) ! y contains the initial values (in vector form)
write(*,*) maxval(y0)


CALL rksolver(yout,y0,nvar,x1,x2,nsteps) ! Call the RK solver

!write(*,*) maxval(abs(yout))

CONTAINS 

  SUBROUTINE rksolver(yout,vstart,nvar,x1,x2,nsteps)
    USE rk4step, ONLY: rk4
    USE equation, ONLY: derivs
    USE global_variables, ONLY: nf
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nsteps, nvar
    REAL(KIND(1.d0)), INTENT(OUT) :: yout(nvar)
    REAL(KIND(1.d0)), INTENT(IN) :: vstart(nvar),x1,x2
    INTEGER, PARAMETER:: NMAX=1e6, NSTPMX=200 ! Max number of variables, max number of steps
    REAL(KIND(1.d0)):: xx(NSTPMX),y(NMAX,NSTPMX) ! Contain the values of x and the corresponding values of y
    INTEGER i,k
    REAL(KIND(1.d0)):: h,x,dv(NMAX),v(NMAX)

    write(*,*) nvar
    write(*,*) NMAX
    IF (nvar>NMAX) THEN
       write(*,*) 'error, NMAX too small!',nvar,NMAX
       stop
    END IF

    OPEN(55,file='testi.dat')

    do i=1,nvar !Load starting values.
       v(i)=vstart(i)
       y(i,1)=v(i)
    enddo
    !

    xx(1)=x1
    write(55,*) xx(1),y(nf+freqs+freqs+1,1)
    x=x1
    h=(x2-x1)/nsteps

    do k=1,nsteps !Take nsteps steps
       !write(*,*) maxval(abs(v))
       call derivs(x,v,dv,nvar) ! Calculate the derivatives dv, used in the first step of RK
       write(*,*) 'x=',x,'max(y)=',maxval(abs(v)),',max(abs(dv))=', maxval(abs(dv))
       write(*,*) maxloc(v)
       call rk4(v,dv,x,nvar,h,v,derivs) ! The Runge-Kutta step, changes v
       !if(x+h.eq.x)pause 'stepsizenot significant in rkdumb'
       x=x+h
       xx(k+1)=x !Store intermediate steps.
       do i=1,nvar ! Loop over the components to store into y at this x
          y(i,k+1)=v(i)
       enddo
       write(55,*) xx(k+1),y(nf+freqs+freqs+2,k+1)
    enddo
    yout=y(:,k+1)

  close(55)
 
  END SUBROUTINE RKSOLVER

  SUBROUTINE vec2mat(vec,mat,m,n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m,n
    REAL(KIND(1.d0)),INTENT(IN) :: vec(m*n)
    REAL(KIND(1.d0)),INTENT(OUT) :: mat(m,n)
    INTEGER :: lask, i,j

    lask=0
    DO j=1,n
       DO i=1,m
          lask=lask+1
          mat(i,j)=vec(lask)
       END DO
    END DO
  END SUBROUTINE VEC2MAT

  SUBROUTINE mat2vec(vec,mat,m,n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m,n
    REAL(KIND(1.d0)),INTENT(OUT) :: vec(m*n)
    REAL(KIND(1.d0)),INTENT(IN) :: mat(m,n)
    INTEGER :: lask, i,j

    lask=0
    DO j=1,n
       DO i=1,m
          lask=lask+1
          vec(lask)=mat(i,j)
       END DO
    END DO
  END SUBROUTINE MAT2VEC

  SUBROUTINE compiley(y,gamma,vfs,vfd,nvar)
    USE global_variables,ONLY: nf,sites,freqs
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nvar
    REAL(KIND(1.d0)),INTENT(OUT) :: y(nvar)
    REAL(KIND(1.d0)), INTENT(IN) :: gamma(nf),vfs(freqs,sites),vfd(freqs,sites)
    INTEGER :: i,j,lask
    !COMMON /indices/ nf,sites,freqs,nvar

    
    lask=0
    DO j=1,nf
       lask=lask+1
       y(j)=gamma(j)
    end DO

    DO j=1,sites
       do i=1,freqs
          lask=lask+1
          y(lask)=vfs(i,j)
       end do
    end do
    
    DO j=1,sites
       do i=1,freqs
          lask=lask+1
          y(lask)=vfd(i,j)
       end do
    end do

  END SUBROUTINE COMPILEY
    

END PROGRAM init
