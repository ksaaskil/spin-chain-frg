PROGRAM main
USE nrtype ! Various definitions, such as SP, PI,...
!gfortran -O3 -fopenmp -o main nrtype.f90 nrutil.f90 ode_path.f90 globalvariables.f90 rkck.f90 rkqstep.f90 odemod.f90 integrators_par.f90 interpmod.f90 derivs.f90 main.f90   
USE globalvariables, ONLY: nf, sites_y,sites_x,freqs, lambda0=>LO,svar,ind,&
ind2,nvar, aa, logaa, om0, omv, gaa,log_gaa,symmetry_flag,site_ind, site_ind2, &
vind, vind2, nind, nind2,vomv,vfreqs, gnf, gom0, gomv ! Initial values for the specific problem
USE derivatives, ONLY: derivs, funktio ! The derivatives
USE rkqstep, ONLY: rkqs ! The Runge-Kutta step with adaptive step size
USE odemod, ONLY: odeint ! The driver for the integration
USE ode_path ! Pointers for the storing of results
USE nrutil, ONLY: assert
IMPLICIT NONE
! The start and end points of the flow
REAL(SP) :: x1=.0_sp,x2=LOG(2*lambda0/om0)
!REAL(SP) :: x1=.0_sp,x2=1_sp
!REAL(SP) :: x1=LOG(lambda0/om0),x2=LOG(10*lambda0/om0)
! Parameters for Runge-Kutta
REAL(SP), PARAMETER :: eps=1e-2_sp, hmin=1e-4_sp, h1=0.1_sp !, abstol=1e-3_sp
INTEGER:: i,j,k,m,lask
! The initial conditions for the self-energy and spin and density interactions
REAL(SP) :: gamma(gnf)=0._sp,vfs(freqs,svar)=0._sp,vfd(freqs,svar)=0._sp
REAL(SP) :: chi0(nf,svar)=0._sp,v0(vfreqs,svar)=0._sp
REAL(SP):: rho0(nf,svar)=0._sp, vrho0(vfreqs,svar)=0._sp
REAL(SP) :: y(nvar)=0._sp,y0(nvar)=0._sp,yout(nvar)=0._sp
REAL(SP) :: LO
! The bare interactions
REAL(SP) :: bareint(svar)
CHARACTER (LEN=*), PARAMETER :: wrf='dat_s2_noint2.dat' ! Where the data will be written
CHARACTER (LEN=*), PARAMETER :: wrf2=wrf//'.restart' ! Where the restart data will be written
CHARACTER (LEN=*), PARAMETER :: wrf3=wrf//'.g0' ! Where the flow of gamma will be written
CHARACTER (LEN=*), PARAMETER :: wrf4=wrf//'.chi' ! Where the flow of the static chi will be written
CHARACTER (LEN=*), PARAMETER :: wrf5=wrf//'.freqs' 
CHARACTER (LEN=*), PARAMETER :: wrf55=wrf//'.gfreqs' 
CHARACTER (LEN=*), PARAMETER :: wrf6=wrf//'.chiall' ! Where the flow of the full chi will be written
CHARACTER (LEN=*), PARAMETER :: wrf65=wrf//'.rhoall'
CHARACTER (LEN=*), PARAMETER :: wrf7=wrf//'.all' ! Where the flow of everything will be written
CHARACTER (LEN=*), PARAMETER :: wrf8=wrf//'.gall'
CHARACTER (LEN=*), PARAMETER :: wrf9=wrf//'.en'
CHARACTER (LEN=*), PARAMETER :: restartfile='dat_s2.dat.restart' ! Where the initial values are read from
LOGICAL, PARAMETER :: restart_boolean=.false.
! Variables for the time counter
REAL(SP) :: start1, finish1
INTEGER :: count,count_rate, count_max
INTEGER :: count2
INTRINSIC LOG

CALL SYSTEM_CLOCK(count,count_rate,count_max)
CALL CPU_TIME(start1)

LO=lambda0*exp(-x1)
print*, 'Flowing to x=',x2,', lambda=',lambda0*exp(-x2)

lask=0
do i=1,sites_y/2+1
   do j=1,sites_x/2+1
      lask=lask+1
      site_ind(lask,1)=i ! The y-coordinate
      site_ind(lask,2)=j ! The x-coordinate
      site_ind2(i,j)=lask
   end do
end do

lask=0
DO i=1,nf
   DO j=1,nf
      DO k=i,nf
         lask=lask+1
         nind(lask,1)=i  ! indexing of frequencies
         nind(lask,2)=j
         nind(lask,3)=k
         nind2(i,j,k)=lask ! inverse indexing
      END DO
   END DO
END DO

call assert(lask==freqs,'number of frequencies')


lask=0
DO i=1,nf
   DO j=1,nf
      DO k=1,nf
         lask=lask+1
         ind(lask,1)=i  ! indexing of frequencies
         ind(lask,2)=j
         ind(lask,3)=k
         ind2(i,j,k)=lask ! inverse indexing
      END DO
   END DO
END DO

lask=0
DO i=1,size(vomv)
   DO j=1,size(vomv)
      lask=lask+1
      vind(lask,1)=i  ! indexing of frequencies
      vind(lask,2)=j
      vind2(i,j)=lask ! inverse indexing
   END DO
END DO
!print*, funktio(1.1_sp),funktio(1.2_sp),funktio(1.3_sp)
aa=rtbis(funktio,1.01_sp,3._sp,1e-7_sp)
print*, 'a=', aa
call assert(aa>1,'invalid a!')
logaa=log(aa)
gaa=aa
log_gaa=log(gaa)

! Construct the frequency mesh
omv=0._sp
forall(i=1:nf-1)
   omv(i+1)=om0*(aa**(i)-1._sp)/(aa-1._sp)
end forall

gomv=0._sp
forall(i=1:gnf-1)
   gomv(i+1)=gom0*(gaa**(i)-1._sp)/(gaa-1._sp)
end forall

print*, 'The smallest non-zero frequency is ', omv(2)
print*, omv

! Construct the frequency mesh for the susceptibility part
vomv=omv

open(99,file=wrf5)
do i=1,nf
   write(99,*) i,omv(i)
end do
write(*,*) 'Wrote frequency mesh to file '//wrf5//'.'
close(99)

open(99,file=wrf55)
do i=1,gnf
   write(99,*) i,gomv(i)
end do
write(*,*) 'Wrote g-frequency mesh to file '//wrf55//'.'
close(99)
!stop
! The bare interaction
bareint=0._sp
!bareint(2)=0.25_sp
!print*, site_ind2(1,2), site_ind2(2,1)
!bareint(site_ind2(1,2))=0.25_sp
!bareint(site_ind2(2,1))=0.25_sp
symmetry_flag=.false.
if (symmetry_flag) print*, 'Assuming xy-symmetry!'
!bareint(3)=0.25_sp*g

!do i=1,svar-1 ! Runs from 1 to 8 for 16 sites
!   angle=2*PI*i/sites
!   !bareint(i+1)=1/((cos(angle)-1)**2+sin(angle)**2)
!   bareint(i+1)=1/sin(angle/2)**2
!end do
!bareint=bareint/bareint(2)

call assert(bareint(1)==0._sp, 'local interactions!')
   
! Initial spin interactions
vfs=0._sp
do i=1,svar
   !forall(j=1:freqs)
      vfs(:,i)=bareint(i) ! The nearest neighbor spin interaction
   !end forall
end do
! Density part
vfd=0._sp
! sz-vertex part
v0=0._sp
v0(:,1)=.5_sp 
! ninja part
chi0=0._sp

vrho0(:,1)=1._sp
rho0=0._sp

CALL compiley(y0,gamma,vfs,vfd,chi0,v0,rho0,vrho0,0._sp,nvar) ! y contains the initial values (in vector form)

!if (.true.) then
if (restart_boolean) then
   write(*,*) 'Restarting from file '//restartfile//'.'
   y0=0._sp
   open(25,file=restartfile)
   read(25,*) y0
   close(25)
end if


CALL odeint(y0,x1,x2,eps,h1,hmin,derivs,rkqs)

OPEN(124,file=wrf)
write(124,*) 'om0', om0
write(124,*) 'om_max', omv(nf)
write(124,*) 'aa', aa
write(124,*) 'nf', nf
write(124,*) 'gnf', gnf
write(124,*) 'sites_x', sites_x
write(124,*) 'sites_y', sites_y
write(124,*) 'freqs', freqs
write(124,*) 'LO', lambda0
write(124,*) 'x1',x1
write(124,*) 'x2', min(x2,maxval(xp))
write(124,*) 'eps', eps
write(124,*) 'flow ', 'log'
do i=1,nvar
   write(124,*) yp(i,kount)
end do
write(*,*) 'Wrote data to file '//wrf//'.'
CLOSE(124)

open(1234,file=wrf2)
do i=1,nvar
   write(1234,*) yp(i,kount)
end do
write(*,*) 'Wrote to file '//wrf2//' for restart.'
close(1234)

if (.true.) then
   open(1235,file=wrf3)
   write(1235,*) 'kount', kount
   do j=1,kount
      do i=1,gnf
         write(1235,*) xp(j),yp(i,j)
      end do
   end do
   write(*,*) 'Flow of gamma stored into '//wrf3//'.'
   close(1235)
end if

if (.true.) then
   open(1236,file=wrf4)
   write(1236,*) 'kount', kount
   do j=1,kount
      do i=1,svar
         lask=gnf+2*freqs*svar+1+(i-1)*nf ! The static ninja functions
         write(1236,*) xp(j),yp(lask,j)
      end do
   end do
   write(*,*) 'Flow of chi stored into '//wrf4//'.'
   close(1236)
end if

if (.true.) then
   open(1236,file=wrf6)
   write(1236,*) 'kount', kount
   do j=1,kount
      do i=1,svar
         do m=1,nf
            lask=gnf+2*freqs*svar+(i-1)*nf+m ! The full ninja functions
            write(1236,*) xp(j),yp(lask,j)
         end do
      end do
   end do
   write(*,*) 'Flow of chiall stored into '//wrf6//'.'
   close(1236)
end if

if (.true.) then
   open(1236,file=wrf65)
   write(1236,*) 'kount', kount
   do j=1,kount
      do i=1,svar
         do m=1,nf
            lask=gnf+2*freqs*svar+nf*svar+vfreqs*svar+(i-1)*nf+m ! The full density-density functions
            write(1236,*) xp(j),yp(lask,j)
         end do
      end do
   end do
   write(*,*) 'Flow of rhoall stored into '//wrf65//'.'
   close(1236)
end if

if (.false.) then
   open(1235,file=wrf7)
   write(1235,*) 'kount', kount
   do j=1,kount
      do i=1,nvar
         write(1235,*) xp(j),yp(i,j)
      end do
   end do
   write(*,*) 'Flow of everything stored into '//wrf7//'.'
   close(1235)
end if

if (.true.) then
   open(1235,file=wrf8)
   write(1235,*) 'kount', kount
   do j=1,kount
      do i=1,gnf
         write(1235,*) xp(j),yp(i,j)
      end do
   end do
   write(*,*) 'Flow of self-energies stored into '//wrf8//'.'
   close(1235)
end if

if (.true.) then
   open(1235,file=wrf9)
   write(1235,*) 'kount', kount
   do j=1,kount
      !do i=1,nf
      write(1235,*) xp(j),yp(nvar,j)
      !end do
   end do
   write(*,*) 'Flow of energies stored into '//wrf9//'.'
   close(1235)
end if

CALL CPU_TIME(finish1)
CALL SYSTEM_CLOCK(count2,count_rate,count_max)

print*, 'Elapsed CPU time ', finish1-start1, ' seconds.'
print*, 'Elapsed wall clock time ',(count2-count)/1000._dp, ' seconds.'

CONTAINS 

  SUBROUTINE compiley(y,gamma,vfs,vfd,chi,v,rho,vrho,pot,nvar)
    USE globalvariables,ONLY: nf,svar,freqs,vfreqs, gnf
    USE nrtype
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nvar
    REAL(SP),INTENT(OUT) :: y(nvar)
    REAL(SP), INTENT(IN) :: gamma(gnf),vfs(freqs,svar),vfd(freqs,svar)
    REAL(SP), INTENT(IN) :: chi(nf,svar),v(vfreqs,svar), pot, rho(nf,svar),vrho(vfreqs,svar)
    INTEGER :: i,j,lask
    
    lask=0
    DO j=1,gnf
       lask=lask+1
       y(j)=gamma(j)
    end DO

    DO j=1,svar
       do i=1,freqs
          lask=lask+1
          y(lask)=vfs(i,j)
       end do
    end do
    
    DO j=1,svar
       do i=1,freqs
          lask=lask+1
          y(lask)=vfd(i,j)
       end do
    end do

    DO j=1,svar
       do i=1,nf
          lask=lask+1
          y(lask)=chi(i,j)
       end do
    end do

    DO j=1,svar
       do i=1,vfreqs
          lask=lask+1
          y(lask)=v(i,j)
       end do
    end do   

     DO j=1,svar
       do i=1,nf
          lask=lask+1
          y(lask)=rho(i,j)
       end do
    end do

    DO j=1,svar
       do i=1,vfreqs
          lask=lask+1
          y(lask)=vrho(i,j)
       end do
    end do   

    lask=lask+1
    y(lask)=pot

  END SUBROUTINE COMPILEY
  
  FUNCTION func_i(x) result(value)
    use nrtype
    use nrutil
    real(sp), intent(in) :: x
    real(sp) :: value
    if (abs(x)>1e-10) then
       value=1._sp/(x*PI)*LOG(1+x)
    else
       value=1._sp/PI
    end if
  end FUNCTION func_i

  FUNCTION rtbis(func,x1,x2,xacc)
    use nrtype
    use nrutil, only: assert
    implicit none
    !INTEGER:: JMAX
    REAL(SP) :: rtbis,x1,x2,xacc,func
    EXTERNAL func
    !INTERFACE
    !     FUNCTION func(x)
    !       USE nrtype
    !       IMPLICIT NONE
    !       REAL(SP), INTENT(IN) :: x
    !     END FUNCTION func
    !  end INTERFACE
    INTEGER, PARAMETER:: JMAX=80
    INTEGER:: j
    REAL(SP):: dx,f,fmid,xmid
    fmid=func(x2)
    f=func(x1)
    call assert(f*fmid<0.,'root must be bracketed in rtbis')
    if(f.lt.0.)then
       rtbis=x1
       dx=x2-x1
    else
       rtbis=x2
       dx=x1-x2
    endif

    do j=1,JMAX
       dx=dx*.5
       xmid=rtbis+dx
       fmid=func(xmid)
       if(fmid.le.0.)rtbis=xmid
       if(abs(dx).lt.xacc .or. fmid.eq.0.) return   
    END do

  end FUNCTION rtbis




END PROGRAM MAIN


