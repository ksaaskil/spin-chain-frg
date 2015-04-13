PROGRAM main
USE derivatives, ONLY: derivs
USE rkqstep, ONLY: rkqs
USE odemod, ONLY: odeint
USE globalvariables, ONLY: nf, sites,freqs, hf
USE nrtype
IMPLICIT NONE
!INTEGER ::nsteps=200
INTEGER,PARAMETER :: nvar=nf+2*sites*freqs
REAL(SP) :: omv(nf),x1=0._sp,x2=1._sp,eps=0.01_sp,h1=0.01_sp,hmin=0._sp
INTEGER:: i,j,k,lask
REAL(SP) :: gamma(nf)=0.d0,vfs(freqs,sites)=0.d0,vfd(freqs,sites)=0.d0 ! Self-energy + vertex functions in matrix form
REAL(SP) :: vfsvec(sites*freqs), vfdvec(sites*freqs) ! vertex functions in vector form
REAL(SP) :: y(nvar)=0._sp,y0(nvar)=0._sp,yout(nvar)=0._sp
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
   vfs(i,2)=.25_sp ! The nearest neighbor spin interaction
END DO

!CALL compiley(y0,gamma,vfs,vfd,nvar) ! y contains the initial values (in vector form)
!write(*,*) maxval(y0)

!CALL rksolver(yout,y0,nvar,x1,x2,nsteps) ! Call the RK solver
CALL odeint(y0,x1,x2,eps,h1,hmin,derivs,rkqs)
  
END PROGRAM MAIN


