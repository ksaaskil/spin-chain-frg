MODULE derivatives
USE nrtype
USE INTERPMOD,ONLY: BLEND_102, BLEND_101
!USE INTERPMOD,ONLY: BLEND_102=>BLEND_102_BF, BLEND_101=>BLEND_101_BF
IMPLICIT NONE
INTEGER :: flag_stu, flag_sd
REAL(SP) :: stuv
PUBLIC:: derivs, funktio
CONTAINS

  SUBROUTINE derivs(x,y,dydx)
    ! Derivatives for a logarihmic frequency mesh

    USE nrtype
    USE globalvariables, ONLY: nf,svar,freqs,LO,ind,ind2,&
         omv, vind, vind2, site2, freq2, gnf, gomv, site_ind,site_ind2,symmetry_flag,&
         nind, nind2, v_in2, vrho_in2, s_in2, d_in2, gout2, g_in2, smat2, dmat2,vfreqs,vomv !,ind,ind2
    USE INTEGRATORS
    IMPLICIT NONE
    REAL(SP),INTENT(IN) :: x
    REAL(SP),DIMENSION(:),INTENT(IN) :: y
    REAL(SP),DIMENSION(:),INTENT(OUT):: dydx
    INTEGER :: nvar
    INTEGER :: indg(gnf),inds(svar*freqs),indd(svar*freqs)
    INTEGER :: indchi(nf*svar),indv(vfreqs*svar), indrho(nf*svar),indvrho(vfreqs*svar)
    INTEGER :: i,j,k,lask,m, site ,freq
    REAL(SP) :: gout(gnf), sout(freqs,svar),dout(freqs,svar) ! The outgoing derivatives of selfenergies and vertex functions
    REAL(SP) :: chiout(nf,svar), vout(vfreqs,svar)
    REAL(SP) :: rhoout(nf,svar), vrhoout(vfreqs,svar)
    REAL(SP) :: stemp, dtemp
    REAL(SP) :: g_in(gnf),s_in_orig(freqs,svar),d_in_orig(freqs,svar) ! vertex function matrices
    REAL(SP) :: chi_in(nf,svar),rho_in(nf,svar)
    INTEGER :: i1,i2
    REAL(SP) :: oms,pf,omega2
    INTEGER :: lambda_ind, om1i, om2i, s1i, u1i, s2i, u2i
    INTEGER :: abstand, nui, nui2
    INTRINSIC EXP, ABS
    REAL(SP) :: lambda
    REAL(SP) :: chitemp, nu, omega, omegaw, omegaf, rhotemp
    REAL(SP) :: vtemp
    INTEGER :: omegai, omegawi, omegafi
    REAL(SP), PARAMETER :: TINY=1.0e-15_sp
    REAL(SP) :: UPLIM, UPLIM2
    REAL(SP) :: s1, s2 
!    INTEGER :: count,count2, count_rate, count_max
    
    gout=0._sp
    gout2=0._sp
    sout=0._sp
    dout=0._sp

    g_in2=0._sp
    s_in2=0._sp
    d_in2=0._sp

    chiout=0._sp
    vout=0._sp

    chi_in=0._sp
    v_in2=0._sp
    !v_in=0._sp
    
    nvar=size(y)
    dydx=0.0_sp
          
    lambda=LO*exp(-x)

    if (abs(x)<1e-10_sp) then
       open(99,file='freqs_test.dat')
       do i=1,10000
         write(99,*) i*10._sp/10000._sp, find_ind_down(i*10._sp/10000._sp)
       end do
       write(*,*) 'Wrote frequency test to file freqs_test.dat.'
       close(99)
       !stop
    end if

    if (abs(x)<1e-10_sp) then
       open(99,file='freqs_test2.dat')
       do i=1,10000
         omega=-10._sp+i*20._sp/10000._sp
         omegai=find_ind_down(omega)
         if (omega>0) then
            omegai=omegai+nf ! The index >= 1
         else
            omegai=nf-omegai ! this can also be zero
         end if
         write(99,*) omega, omegai
       end do
       write(*,*) 'Wrote frequency test to file freqs_test2.dat.'
       close(99)
       !stop
    end if


    !stop

    DO i=1,gnf
       indg(i)=i
    END DO

    g_in=y(indg)

    j=0
    DO i=gnf+1,gnf+svar*freqs
       j=j+1
       inds(j)=i
    END DO

    ! Construct sin, the incoming matrix for spin interactions
    CALL VEC2MAT(y(inds),s_in_orig,freqs,svar)

    j=0
    DO i=gnf+svar*freqs+1,gnf+2*svar*freqs
       j=j+1
       indd(j)=i
    END DO
    
    ! Construct din, the incoming matrix for density interactions
    CALL VEC2MAT(y(indd),d_in_orig,freqs,svar)

    !

    do site=1,svar
       lask=0
       do i=1,nf
          do j=1,nf
             do k=1,nf
                lask=lask+1
                s_in2(lask,site)=s_in_orig(nind2(min(i,k),j,max(i,k)),site)
                d_in2(lask,site)=d_in_orig(nind2(min(i,k),j,max(i,k)),site)
                if (k<i) then ! Antisymmetric
                   d_in2(lask,site)=-d_in2(lask,site)
                elseif (i==k) then
                   d_in2(lask,site)=0._sp
                end if
             end do
          end do
       end do
    end do
    !CALL SYSTEM_CLOCK(count2,count_rate,count_max)
 
    !print*, 'Elapsed wall clock time ',(count2-count)/1000._dp, ' seconds.'
    !stop

    j=0
    DO i=gnf+2*svar*freqs+1,gnf+2*svar*freqs+nf*svar
       j=j+1
       indchi(j)=i
    end do
    
    CALL VEC2MAT(y(indchi),chi_in,nf,svar)

    j=0
    DO i=gnf+2*svar*freqs+nf*svar+1,gnf+2*svar*freqs+nf*svar+vfreqs*svar
       j=j+1
       indv(j)=i
    end do

    CALL VEC2MAT(y(indv),v_in2,vfreqs,svar)   

    j=0
    DO i=gnf+2*svar*freqs+(nf+vfreqs)*svar,gnf+2*svar*freqs+(nf+vfreqs)*svar+nf*svar
       j=j+1
       indrho(j)=i
    end do
    
    CALL VEC2MAT(y(indrho),rho_in,nf,svar)

    j=0
    DO i=gnf+2*svar*freqs+(nf+vfreqs)*svar+nf*svar+1,gnf+2*svar*freqs+(nf+vfreqs)*svar+nf*svar+vfreqs*svar
       j=j+1
       indvrho(j)=i
    end do

    CALL VEC2MAT(y(indvrho),vrho_in2,vfreqs,svar)   
    

    g_in2=g_in
    
    !---------------------------------------------------------------------


    CALL calc_g_i(gout,g_in2,s_in2,d_in2,lambda)

    ! For the Katanin integration
    g_in2=g_in
    gout2=gout/(2*PI)


    do site=1,svar
       smat2=0._sp
       dmat2=0._sp
       forall(i=1:nf,j=1:nf,k=1:nf)
          smat2(i,j,k)=s_in2(ind2(i,j,k),site)
          dmat2(i,j,k)=d_in2(ind2(i,j,k),site)
       end forall
       if (symmetry_flag.and.(site_ind(site,1)>site_ind(site,2))) then ! The y-direction larger
             sout(:,site)=sout(:,site_ind2(site_ind(site,2),site_ind(site,1)))
             dout(:,site)=dout(:,site_ind2(site_ind(site,2),site_ind(site,1)))
       else
          !$OMP parallel do &
          !$OMP DEFAULT(NONE) &
          !$OMP PRIVATE(stuv,stemp,dtemp,oms,uplim,uplim2,s1,s2,flag_sd,flag_stu,freq,pf,omega2) &
          !$OMP shared(omv,nind,lambda,s_in2,d_in2,g_in2,gout2,ind2,sout,dout,site,smat2,dmat2)
       do freq=1,freqs
          do flag_stu=1,3 ! Over the channels
             stuv=omv(nind(freq,flag_stu))
            
             do k=1,4 ! Over the four possibilities for oms
                stemp=0._sp
                dtemp=0._sp

                IF (k==1) THEN
                   oms=lambda 
                   omega2=oms+stuv
                ELSEIF (k==2) THEN
                   oms=-lambda
                   omega2=oms+stuv
                ELSEIF (k==3) THEN
                   oms=lambda-stuv
                   omega2=oms
                ELSEIF (k==4) THEN
                   oms=-lambda-stuv
                   omega2=oms
                END IF

                !if (((k==1.or.k==2).and.abs(stuv+oms)>lambda).or.((k==3.or.k==4).and.abs(oms)>lambda)) then
                if (abs(omega2)+1e-3>=lambda) then

                   if (flag_stu==1) then
                      CALL CALC_SCH_I(oms,stemp,dtemp,smat2,dmat2,freq)
                   elseif (flag_stu==2) then
                      CALL CALC_TCH_I(oms,stemp,dtemp,s_in2,d_in2,site,freq)
                   elseif (flag_stu==3) then
                      CALL CALC_UCH_I(oms,stemp,dtemp,smat2,dmat2,freq)
                   else
                      stop
                   end if
                   
                   IF (k==1.or.k==2) THEN ! The same pf for k=1 and k=2
                      pf=pfunc_i(oms,stuv+oms,lambda,g_in2)
                   ELSEif (k==3.or.k==4) then
                      pf=pfunc_i(stuv+oms,oms,lambda,g_in2)
                   END IF
                   stemp=stemp*pf
                   dtemp=dtemp*pf

                   sout(freq,site)=sout(freq,site)+stemp
                   dout(freq,site)=dout(freq,site)+dtemp
                end if

             END do
                           
             IF (.true.) THEN ! Katanin integration
   
                uplim=max(lambda,stuv-lambda)+20*omv(nf)
                uplim2=lambda+stuv+20*omv(nf)
                
                s1=0._sp
                s2=0._sp
         
                flag_sd=1 ! Spin interactions                
                s1=qromo(func_to_int,max(lambda,stuv-lambda),uplim,midinf,flag_stu,flag_sd,freq,stuv,site)
                sout(freq,site)=sout(freq,site)+s1                 
                s2=qromo(func_to_int,-uplim2,-lambda-stuv,midinf,flag_stu,flag_sd,freq,stuv,site)                 
                sout(freq,site)=sout(freq,site)+s2
                if (stuv>2*lambda) then                 
                   s1=qromo(func_to_int,lambda-stuv,-lambda,midpnt,flag_stu,flag_sd,freq,stuv,site)
                   s2=qromo(func_to_int,lambda,stuv-lambda,midpnt,flag_stu,flag_sd,freq,stuv,site)           
                   sout(freq,site)=sout(freq,site)+s1+s2
                end if
                flag_sd=2 ! Density interactions
                s1=qromo(func_to_int,max(lambda,stuv-lambda),uplim,midinf,flag_stu,flag_sd,freq,stuv,site)
                dout(freq,site)=dout(freq,site)+s1
                s2=qromo(func_to_int,-uplim2,-lambda-stuv,midinf,flag_stu,flag_sd,freq,stuv,site)
                dout(freq,site)=dout(freq,site)+s2
                if (stuv>2*lambda) then
                   s1=qromo(func_to_int,lambda-stuv,-lambda,midpnt,flag_stu,flag_sd,freq,stuv,site)
                   s2=qromo(func_to_int,lambda,stuv-lambda,midpnt,flag_stu,flag_sd,freq,stuv,site)
                   dout(freq,site)=dout(freq,site)+s1+s2
                end if
             END IF
          END do
       END do
!$OMP END PARALLEL DO
    end if
    END do

    ! The flow eq. for the ninja function

    chiout=0._sp

    do site=1,svar
       do freq=1,nf
          omega=omv(freq)
          do k=1,4
             chitemp=0._sp  
             IF (k==1) THEN
                nu=lambda !-> |nu+omega|>lambda
             ELSEIF (k==2) THEN
                nu=-lambda
             ELSEIF (k==3) THEN
                nu=lambda-omega
             ELSEIF (k==4) THEN
                nu=-lambda-omega !-> |nu|>lambda
             END IF
             
             call calc_chi_i(nu,chitemp,v_in2,freq,site)

             IF (k==1.OR.k==2) THEN
                chitemp=chitemp*pfunc_i(nu,nu+omega,lambda,g_in)
             ELSEif (k==3.or.k==4) then
                chitemp=chitemp*pfunc_i(nu+omega,nu,lambda,g_in)
             END IF
             
             chiout(freq,site)=chiout(freq,site)+chitemp

             if (.false.) then
             uplim=max(lambda,omega-lambda)+20*omv(nf)        
             uplim2=lambda+omega+20*omv(nf)
             chiout(freq,site)=chiout(freq,site)+qromo(func_to_int_chi,max(lambda,omega-lambda),uplim,midinf&
                  ,0,0,freq,omega,site)
             chiout(freq,site)=chiout(freq,site)+qromo(func_to_int_chi,-uplim2,-lambda-omega,midinf&
                  ,0,0,freq,omega,site)
             if (omega>2*lambda) then
                chiout(freq,site)=chiout(freq,site)+qromo(func_to_int_chi,lambda-omega,-lambda,midpnt&
                     ,0,0,freq,omega,site)              
                chiout(freq,site)=chiout(freq,site)+qromo(func_to_int_chi,lambda,omega-lambda,midpnt&
                     ,0,0,freq,omega,site)
                !vout(freq,site)=vout(freq,site)+qromo(func_to_int_v,lambda-stuv,lambda-stuv+(stuv-2*lambda)/2._sp,midpnt&
                !     ,0,0,freq,stuv,site)              
                !vout(freq,site)=vout(freq,site)+qromo(func_to_int_v,lambda-stuv+(stuv-2*lambda)/2._sp,-lambda,midpnt&
                !     ,0,0,freq,stuv,site)  
                !vout(freq,site)=vout(freq,site)+qromo(func_to_int_v,lambda,lambda+(stuv-2*lambda)/2._sp,midpnt&
                !     ,0,0,freq,stuv,site)           
                !vout(freq,site)=vout(freq,site)+qromo(func_to_int_v,lambda+(stuv-2*lambda)/2._sp,stuv-lambda,midpnt&
                !     ,0,0,freq,stuv,site)
             end if
          end if
          end do
       end do
    end do

    vout=0._sp
    
    DO site=1,svar ! The flow eqs. for the fermion-boson vertex functions
       do freq=1,vfreqs
          omegai=vind(freq,1)
          omega=vomv(omegai)
          stuv=vomv(vind(freq,1))

          do k=1,4
             vtemp=0._sp
             IF (k==1) THEN
                nu=lambda
             ELSEIF (k==2) THEN
                nu=-lambda 
             ELSEIF (k==3) THEN
                nu=lambda-omega 
             ELSEIF (k==4) THEN
                nu=-lambda-omega
             END IF

             CALL calc_vch_i(nu,vtemp,s_in2,d_in2,v_in2,site,freq)

             IF (k==1.OR.k==2) THEN
                vtemp=vtemp*pfunc_i(nu,nu+omega,lambda,g_in)
             ELSE
                vtemp=vtemp*pfunc_i(nu+omega,nu,lambda,g_in)
             END IF
             vout(freq,site)=vout(freq,site)+vtemp
          end do


          if (.false.) then
             uplim=max(lambda,stuv-lambda)+20*omv(nf)        
             uplim2=lambda+stuv+20*omv(nf)
             vout(freq,site)=vout(freq,site)+qromo(func_to_int_v,max(lambda,stuv-lambda),uplim,midinf&
                  ,flag_stu,flag_sd,freq,stuv,site)
             vout(freq,site)=vout(freq,site)+qromo(func_to_int_v,-uplim2,-lambda-stuv,midinf&
                  ,flag_stu,flag_sd,freq,stuv,site)
             if (stuv>2*lambda) then
                vout(freq,site)=vout(freq,site)+qromo(func_to_int_v,lambda-stuv,-lambda,midpnt&
                     ,0,0,freq,stuv,site)              
                vout(freq,site)=vout(freq,site)+qromo(func_to_int_v,lambda,stuv-lambda,midpnt&
                     ,0,0,freq,stuv,site)
                !vout(freq,site)=vout(freq,site)+qromo(func_to_int_v,lambda-stuv,lambda-stuv+(stuv-2*lambda)/2._sp,midpnt&
                !     ,0,0,freq,stuv,site)              
                !vout(freq,site)=vout(freq,site)+qromo(func_to_int_v,lambda-stuv+(stuv-2*lambda)/2._sp,-lambda,midpnt&
                !     ,0,0,freq,stuv,site)  
                !vout(freq,site)=vout(freq,site)+qromo(func_to_int_v,lambda,lambda+(stuv-2*lambda)/2._sp,midpnt&
                !     ,0,0,freq,stuv,site)           
                !vout(freq,site)=vout(freq,site)+qromo(func_to_int_v,lambda+(stuv-2*lambda)/2._sp,stuv-lambda,midpnt&
                !     ,0,0,freq,stuv,site)
             end if
          end if

       end do
    end do

    chiout=-chiout ! The fermion loop

    rhoout=0._sp
    !print*, maxval(abs(vrho_in2))
    do site=1,svar
       do freq=1,nf
          omega=omv(freq)
          do k=1,4
             rhotemp=0._sp  
             IF (k==1) THEN
                nu=lambda !-> |nu+omega|>lambda
             ELSEIF (k==2) THEN
                nu=-lambda
             ELSEIF (k==3) THEN
                nu=lambda-omega
             ELSEIF (k==4) THEN
                nu=-lambda-omega !-> |nu|>lambda
             END IF
             
             call calc_rho_i(nu,rhotemp,vrho_in2,freq,site)
             !print*, rhotemp
             IF (k==1.OR.k==2) THEN
                rhotemp=rhotemp*pfunc_i(nu,nu+omega,lambda,g_in)
             ELSEif (k==3.or.k==4) then
                rhotemp=rhotemp*pfunc_i(nu+omega,nu,lambda,g_in)
             END IF
             
             rhoout(freq,site)=rhoout(freq,site)+rhotemp

             if (.false.) then !NOT IMPLEMENTED
             uplim=max(lambda,omega-lambda)+20*omv(nf)        
             uplim2=lambda+omega+20*omv(nf)
             chiout(freq,site)=chiout(freq,site)+qromo(func_to_int_chi,max(lambda,omega-lambda),uplim,midinf&
                  ,0,0,freq,omega,site)
             chiout(freq,site)=chiout(freq,site)+qromo(func_to_int_chi,-uplim2,-lambda-omega,midinf&
                  ,0,0,freq,omega,site)
             if (omega>2*lambda) then
                chiout(freq,site)=chiout(freq,site)+qromo(func_to_int_chi,lambda-omega,-lambda,midpnt&
                     ,0,0,freq,omega,site)              
                chiout(freq,site)=chiout(freq,site)+qromo(func_to_int_chi,lambda,omega-lambda,midpnt&
                     ,0,0,freq,omega,site)
             end if
          end if
          end do
       end do
    end do

    !print*, maxval(abs(rhoout))

    vrhoout=0._sp
    
    DO site=1,svar ! The flow eqs. for the fermion-boson vertex functions
       do freq=1,vfreqs
          omegai=vind(freq,1)
          omega=vomv(omegai)
          stuv=vomv(vind(freq,1))

          do k=1,4
             vtemp=0._sp
             IF (k==1) THEN
                nu=lambda
             ELSEIF (k==2) THEN
                nu=-lambda 
             ELSEIF (k==3) THEN
                nu=lambda-omega 
             ELSEIF (k==4) THEN
                nu=-lambda-omega
             END IF

             CALL calc_vrho_i(nu,vtemp,s_in2,d_in2,vrho_in2,site,freq)

             IF (k==1.OR.k==2) THEN
                vtemp=vtemp*pfunc_i(nu,nu+omega,lambda,g_in)
             ELSE
                vtemp=vtemp*pfunc_i(nu+omega,nu,lambda,g_in)
             END IF
             vrhoout(freq,site)=vrhoout(freq,site)+vtemp
          end do


          if (.false.) then
             uplim=max(lambda,stuv-lambda)+20*omv(nf)        
             uplim2=lambda+stuv+20*omv(nf)
             vrhoout(freq,site)=vrhoout(freq,site)+qromo(func_to_int_vrho,max(lambda,stuv-lambda),uplim,midinf&
                  ,flag_stu,flag_sd,freq,stuv,site)
             vrhoout(freq,site)=vrhoout(freq,site)+qromo(func_to_int_vrho,-uplim2,-lambda-stuv,midinf&
                  ,flag_stu,flag_sd,freq,stuv,site)
             if (stuv>2*lambda) then
                vrhoout(freq,site)=vrhoout(freq,site)+qromo(func_to_int_vrho,lambda-stuv,-lambda,midpnt&
                     ,0,0,freq,stuv,site)              
                vrhoout(freq,site)=vrhoout(freq,site)+qromo(func_to_int_vrho,lambda,stuv-lambda,midpnt&
                     ,0,0,freq,stuv,site)
             end if
          end if

       end do
    end do

    rhoout=-rhoout ! fermion loop

    lask=0

    dydx=0._sp
 
    DO j=1,gnf
       lask=lask+1
       dydx(j)=gout(j)
    END DO
  
    DO j=1,svar
       DO i=1,freqs
          lask=lask+1
          dydx(lask)=sout(i,j)
       END DO
    END DO
    
    DO j=1,svar
       DO i=1,freqs
          lask=lask+1
          dydx(lask)=dout(i,j)
       END DO
    END DO

    DO j=1,svar
       do i=1,nf
          lask=lask+1
          dydx(lask)=chiout(i,j)
       end do
    end do

    DO j=1,svar
       do i=1,vfreqs
          lask=lask+1
          dydx(lask)=vout(i,j)
       end do
    end do

    DO j=1,svar
       do i=1,nf
          lask=lask+1
          dydx(lask)=rhoout(i,j)
       end do
    end do

    DO j=1,svar
       do i=1,vfreqs
          lask=lask+1
          dydx(lask)=vrhoout(i,j)
       end do
    end do



    ! DO j=1,svar
    !   do i=1,vfreqs
          lask=lask+1
          dydx(lask)=4*log(1+interpg(lambda,g_in)/lambda)
    !   end do
    !end do
    
    
    dydx=dydx/(2*PI) ! The overall prefactor
    dydx=dydx*(-lambda) ! Exponential flow


    if (abs(x)<tiny) then
       open(55,file='derivs.dat')
       do i=1,nvar
          write(55,*) dydx(i)
       end do
       write(*,*) 'wrote to file derivs.dat, max(y)=',maxval(abs(y))
       close(55)  
    end if

    
  END SUBROUTINE derivs

  SUBROUTINE CALC_CHI_I(nu,chitemp,v_in,freq,sud)
    use nrtype
    use globalvariables, only: vind2,sites_x,sites_y,nf,vomv,site_ind2,site_ind
    implicit none
    real(sp), intent(in) :: nu,  v_in(:,:)
    integer, intent(in) :: freq,sud
    real(sp), intent(out) :: chitemp
    real(sp) :: omega,r1,v1,v2,s2v
    integer :: s2i,s2in,i1,i2,m,omegai
     integer :: abstand1x,abstand1y,abstand2x,abstand2y,n,sudx,sudy

    omegai=freq
    omega=vomv(freq)

    s2v=2*nu+omega ! The s-argument for the vertex functions
    s2i=find_ind_down(s2v)
         
    s2in=min(s2i+1,nf) ! The number of frequencies in vomv is nf
    s2i=max(s2i,1)
    
    if (s2in>s2i) then
       r1=(abs(s2v)-vomv(s2i))/(vomv(s2in)-vomv(s2i))
    else
       r1=0._sp
    end if
    
    chitemp=0._sp

    sudy=site_ind(sud,1)
    sudx=site_ind(sud,2)
             
!    do m=1,sites
    do m=1,sites_y ! In the y-direction
       abstand1y=min(m-1,sites_y-(m-1))
       abstand2y=min(abs(sudy-m),sites_y-abs(sudy-m))
       do n=1,sites_x ! In the x-direction
          abstand1x=min(n-1,sites_x-(n-1))
          abstand2x=min(abs(sudx-n),sites_x-abs(sudx-n))
          i1=site_ind2(abstand1y+1,abstand1x+1)
          i2=site_ind2(abstand2y+1,abstand2x+1)
          !      i1=min(m-1,sites-(m-1))+1
          !      i2=min(abs(m-sud),sites-abs(m-sud))+1
          ! Argument of the first is -omega<0
          call blend_101(r1,v_in(vind2(omegai,s2i),i1),v_in(vind2(omegai,s2in),i1),v1)
          ! First argument omega>0
          call blend_101(r1,v_in(vind2(omegai,s2i),i2),v_in(vind2(omegai,s2in),i2),v2)
          chitemp=chitemp+2*v1*v2
       end do
    END do
    
    
  end SUBROUTINE CALC_CHI_I

  SUBROUTINE CALC_RHO_I(nu,rhotemp,vrho_in,freq,sud)
    use nrtype
    use globalvariables, only: omv,vind2,sites_x,sites_y,nf,vomv,site_ind2,site_ind
    implicit none
    real(sp), intent(in) :: nu,  vrho_in(:,:)
    integer, intent(in) :: freq,sud
    real(sp), intent(out) :: rhotemp
    real(sp) :: omega,r1,v1,v2,s2v
    integer :: s2i,s2in,i1,i2,m,omegai
    integer :: abstand1x,abstand1y,abstand2x,abstand2y,n,sudx,sudy

    omegai=freq
    omega=vomv(freq)

    s2v=2*nu+omega ! The s-argument for the vertex functions
    s2i=find_ind_down(s2v)
         
    s2in=min(s2i+1,nf) ! The number of frequencies in vomv is nf
    s2i=max(s2i,1)
    
    if (s2in>s2i) then
       r1=(abs(s2v)-vomv(s2i))/(vomv(s2in)-vomv(s2i))
    else
       r1=0._sp
    end if
    
    rhotemp=0._sp

    sudy=site_ind(sud,1)
    sudx=site_ind(sud,2)
             
    do m=1,sites_y ! In the y-direction
       abstand1y=min(m-1,sites_y-(m-1))
       abstand2y=min(abs(sudy-m),sites_y-abs(sudy-m))
       do n=1,sites_x ! In the x-direction
          abstand1x=min(n-1,sites_x-(n-1))
          abstand2x=min(abs(sudx-n),sites_x-abs(sudx-n))
          i1=site_ind2(abstand1y+1,abstand1x+1)
          i2=site_ind2(abstand2y+1,abstand2x+1)
          ! Argument of the first is -omega<0
          call blend_101(r1,vrho_in(vind2(omegai,s2i),i1),vrho_in(vind2(omegai,s2in),i1),v1)
          ! First argument omega>0
          call blend_101(r1,vrho_in(vind2(omegai,s2i),i2),vrho_in(vind2(omegai,s2in),i2),v2)
          rhotemp=rhotemp+2*v1*v2
       end do
    END do
    
  end SUBROUTINE CALC_RHO_I

    SUBROUTINE CALC_VCH_I(nu,vtemp,s_in,d_in,v_in,sud,freq)
      USE nrtype
      USE globalvariables, ONLY : omv, ind2, vind, vind2, sites_x,sites_y,&
           nf,vomv, site_ind2, site_ind
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: nu
      REAL(SP), INTENT(OUT) :: vtemp
      REAL(SP), INTENT(IN) :: s_in(:,:), d_in(:,:), v_in(:,:)
      INTEGER, INTENT(IN) :: sud, freq ! Site and frequency under study
      INTEGER :: s1i, u1i,nuin,s1in,u1in
      REAL(SP) :: omega, sw, sw2,s1v,u1v
      INTEGER :: omegai,sw2i,sw2in
      INTEGER :: m, i1, i2
      REAL(SP) :: v1,s1,d1,rn,rr,ss
      integer :: abstand1x,abstand1y,abstand2x,abstand2y,n,sudy,sudx

      vtemp=0._sp

      omega=vomv(vind(freq,1))
      sw=vomv(vind(freq,2))

      omegai=vind(freq,1)

      sw2=-2*nu-omega
      sw2i=find_ind_down(sw2)

      sw2in=min(sw2i+1,nf)
      sw2i=max(sw2i,1)

      if (sw2in>sw2i) then ! the scaled interpolation point
         rn=(abs(sw2)-vomv(sw2i))/(vomv(sw2in)-vomv(sw2i))
      else
         rn=0._sp
      end if

      if (abs(rn-0.5_sp)>1) then
         stop
      end if

      s1v=(-omega+sw)/2._sp-nu
      u1v=(omega+sw)/2._sp+nu
      
      s1i=find_ind_down(s1v)
      s1in=min(s1i+1,nf)
      s1i=max(s1i,1)

      u1i=find_ind_down(u1v)
      u1in=min(u1i+1,nf)
      u1i=max(u1i,1)

      if (s1in>s1i) then
         rr=(abs(s1v)-omv(s1i))/(omv(s1in)-omv(s1i))
      else
         rr=0._sp ! arbitrary
      end if

      if (u1in>u1i) then
         ss=(abs(u1v)-omv(u1i))/(omv(u1in)-omv(u1i))
      else
         ss=0._sp ! arbitrary
      end if

      sudy=site_ind(sud,1)
      sudx=site_ind(sud,2)
      do m=1,sites_y ! In the y-direction
         abstand1y=min(m-1,sites_y-(m-1))
         abstand2y=min(abs(sudy-m),sites_y-abs(sudy-m))
         do n=1,sites_x ! In the x-direction
            abstand1x=min(n-1,sites_x-(n-1))
            abstand2x=min(abs(sudx-n),sites_x-abs(sudx-n))
            ! First the y-direction (smaller), then the x-direction (bigger)
            i1=site_ind2(abstand1y+1,abstand1x+1)
            i2=site_ind2(abstand2y+1,abstand2x+1)
            v1=0._sp
            s1=0._sp
            call blend_101(rn,v_in(vind2(omegai,sw2i),i1),v_in(vind2(omegai,sw2in),i1),v1)
            call blend_102(rr,ss,s_in(ind2(s1i,omegai,u1i),i2),s_in(ind2(s1in,omegai,u1i),i2),&
                 s_in(ind2(s1i,omegai,u1in),i2),s_in(ind2(s1in,omegai,u1in),i2),s1)
            vtemp=vtemp+2*v1*s1
         end do
      END do
      v1=0._sp
      call blend_101(rn,v_in(vind2(omegai,sw2i),sud),v_in(vind2(omegai,sw2in),sud),v1)
      s1=0._sp
      d1=0._sp
      call blend_102(rr,ss,s_in(ind2(s1i,u1i,omegai),1),s_in(ind2(s1in,u1i,omegai),1),&
           s_in(ind2(s1i,u1in,omegai),1),s_in(ind2(s1in,u1in,omegai),1),s1)
      call blend_102(rr,ss,d_in(ind2(s1i,u1i,omegai),1),d_in(ind2(s1in,u1i,omegai),1),&
           d_in(ind2(s1i,u1in,omegai),1),d_in(ind2(s1in,u1in,omegai),1),d1)
      
      vtemp=vtemp+v1*(s1-d1)
      
    end SUBROUTINE CALC_VCH_I

    SUBROUTINE CALC_VRHO_I(nu,vtemp,s_in,d_in,vrho_in,sud,freq)
      USE nrtype
      USE globalvariables, ONLY : omv, ind2, vind, vind2, sites_x,sites_y,&
           nf,vomv, site_ind2, site_ind
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: nu
      REAL(SP), INTENT(OUT) :: vtemp
      REAL(SP), INTENT(IN) :: s_in(:,:), d_in(:,:), vrho_in(:,:)
      INTEGER, INTENT(IN) :: sud, freq ! Site and frequency under study
      INTEGER :: s1i, u1i,nuin,s1in,u1in
      REAL(SP) :: omega, sw, sw2,s1v,u1v
      INTEGER :: omegai,sw2i,sw2in
      INTEGER :: m, i1, i2
      REAL(SP) :: v1,s1,d1,rn,rr,ss
      integer :: abstand1x,abstand1y,abstand2x,abstand2y,n,sudy,sudx

      vtemp=0._sp

      omega=vomv(vind(freq,1))
      sw=vomv(vind(freq,2))

      omegai=vind(freq,1)

      sw2=-2*nu-omega
      sw2i=find_ind_down(sw2)

      sw2in=min(sw2i+1,nf)
      sw2i=max(sw2i,1)

      if (sw2in>sw2i) then ! the scaled interpolation point
         rn=(abs(sw2)-vomv(sw2i))/(vomv(sw2in)-vomv(sw2i))
      else
         rn=0._sp
      end if

      if (abs(rn-0.5_sp)>1) then
         stop
      end if

      s1v=(-omega+sw)/2._sp-nu
      u1v=(omega+sw)/2._sp+nu
      
      s1i=find_ind_down(s1v)
      s1in=min(s1i+1,nf)
      s1i=max(s1i,1)

      u1i=find_ind_down(u1v)
      u1in=min(u1i+1,nf)
      u1i=max(u1i,1)

      if (s1in>s1i) then
         rr=(abs(s1v)-omv(s1i))/(omv(s1in)-omv(s1i))
      else
         rr=0._sp ! arbitrary
      end if

      if (u1in>u1i) then
         ss=(abs(u1v)-omv(u1i))/(omv(u1in)-omv(u1i))
      else
         ss=0._sp ! arbitrary
      end if

      sudy=site_ind(sud,1)
      sudx=site_ind(sud,2)
      do m=1,sites_y ! In the y-direction
         abstand1y=min(m-1,sites_y-(m-1))
         abstand2y=min(abs(sudy-m),sites_y-abs(sudy-m))
         do n=1,sites_x ! In the x-direction
            abstand1x=min(n-1,sites_x-(n-1))
            abstand2x=min(abs(sudx-n),sites_x-abs(sudx-n))
            ! First the y-direction (smaller), then the x-direction (bigger)
            i1=site_ind2(abstand1y+1,abstand1x+1)
            i2=site_ind2(abstand2y+1,abstand2x+1)
            v1=0._sp
            d1=0._sp
            call blend_101(rn,vrho_in(vind2(omegai,sw2i),i1),vrho_in(vind2(omegai,sw2in),i1),v1)
            call blend_102(rr,ss,d_in(ind2(s1i,omegai,u1i),i2),d_in(ind2(s1in,omegai,u1i),i2),&
                 d_in(ind2(s1i,omegai,u1in),i2),d_in(ind2(s1in,omegai,u1in),i2),d1)
            vtemp=vtemp+2*v1*d1
           ! print*, 'hello'
         end do
      END do
      v1=0._sp
      call blend_101(rn,vrho_in(vind2(omegai,sw2i),sud),vrho_in(vind2(omegai,sw2in),sud),v1)
      s1=0._sp
      d1=0._sp
      call blend_102(rr,ss,s_in(ind2(s1i,u1i,omegai),1),s_in(ind2(s1in,u1i,omegai),1),&
           s_in(ind2(s1i,u1in,omegai),1),s_in(ind2(s1in,u1in,omegai),1),s1)
      call blend_102(rr,ss,d_in(ind2(s1i,u1i,omegai),1),d_in(ind2(s1in,u1i,omegai),1),&
           d_in(ind2(s1i,u1in,omegai),1),d_in(ind2(s1in,u1in,omegai),1),d1)

      vtemp=vtemp+v1*(-3*s1-d1)

    end SUBROUTINE CALC_VRHO_I

    SUBROUTINE CALC_SCH_I(oms,stemp,dtemp,s,d,freq)
      USE nrtype
      USE globalvariables, ONLY : nind, ind2, omv,nf
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: oms
      REAL(SP), INTENT(OUT) :: stemp, dtemp
      REAL(SP), INTENT(IN) :: s(:,:,:), d(:,:,:)
      INTEGER, INTENT(IN) :: freq ! Site and frequency under study
      REAL(SP) :: sv, tv, uv, om1s, om1, om2s, om2
      INTEGER :: si, ti, ui, i, j, k
      INTEGER :: t1i,t2i,u1i,u2i, t1in,t2in,u1in,u2in
      REAL(SP) :: t1v,t2v,u1v,u2v,rr1,ss1,rr2,ss2,s1,s2,d1,d2

      s1=0._sp;s2=0._sp;d1=0._sp;d2=0._sp
      
      stemp=0._sp
      dtemp=0._sp

      si=nind(freq,1)          
      sv=omv(si)
      ti=nind(freq,2)
      tv=omv(ti)
      ui=nind(freq,3)    
      uv=omv(ui)
          
      om1s=(sv+tv+uv)/2._sp  ! omega_1 strich etc.
      om1=(sv-tv+uv)/2._sp
      om2s=(sv-tv-uv)/2._sp
      om2=(sv+tv-uv)/2._sp
      
      t1v=-om2s-oms  
      t1i=FIND_IND_DOWN(abs(t1v))
      t1in=min(t1i+1,nf) ! The next index or the largest one
      
      u1v=om1s+oms
      u1i=FIND_IND_DOWN(abs(u1v))
      u1in=min(u1i+1,nf) 
      
      t2v=om2+oms
      t2i=FIND_IND_DOWN(abs(t2v))
      t2in=min(t2i+1,nf) 
      
      u2v=om1+oms
      u2i=FIND_IND_DOWN(abs(u2v))
      u2in=min(u2i+1,nf) 
      
      if ((t1in>t1i)) then ! The boundaries
         rr1=(abs(t1v)-omv(t1i))/(omv(t1in)-omv(t1i)) ! The scaled variable
      else
         rr1=0._sp ! Evaluate simply with x==0 (should be arbitrary)
      end if
      if ((u1in>u1i)) then
         ss1=(abs(u1v)-omv(u1i))/(omv(u1in)-omv(u1i)) ! The scaled variable
      else
         ss1=0._sp
      end if
      
      call blend_102(rr1,ss1,s(si,t1i,u1i),s(si,t1in,u1i),s(si,t1i,u1in),s(si,t1in,u1in),s1)
      call blend_102(rr1,ss1,d(si,t1i,u1i),d(si,t1in,u1i),d(si,t1i,u1in),d(si,t1in,u1in),d1)

      if ((t2in>t2i)) then ! The boundaries
         rr2=(abs(t2v)-omv(t2i))/(omv(t2in)-omv(t2i)) ! The scaled variable
      else
         rr2=0._sp ! Evaluate simply with x==0 (arbitrary)
      end if
      if ((u2in>u2i)) then
         ss2=(abs(u2v)-omv(u2i))/(omv(u2in)-omv(u2i)) ! The scaled variable
      else
         ss2=0._sp
      end if

      if (rr1<0) then
         write(*,*) t1i,t1in,abs(t1v),omv(t1i),omv(t1in),FIND_IND_DOWN(abs(t1v))
         stop
      elseif (rr2<0) then
         write(*,*) t2i,t2in,abs(t2v),omv(t2i),omv(t2in),FIND_IND_DOWN(abs(t2v))
         stop
      end if

      call blend_102(rr2,ss2,s(si,t2i,u2i),s(si,t2in,u2i),s(si,t2i,u2in),s(si,t2in,u2in),s2)
      call blend_102(rr2,ss2,d(si,t2i,u2i),d(si,t2in,u2i),d(si,t2i,u2in),d(si,t2in,u2in),d2)

      stemp=s1*(-2*s2+d2)
      stemp=stemp+d1*s2          
      
      dtemp=3*s1*s2
      dtemp=dtemp+d1*d2
      
    END SUBROUTINE CALC_SCH_I

    SUBROUTINE CALC_TCH_I(oms,stemp,dtemp,s_in,d_in,sud,freq)
      USE nrtype
      USE globalvariables, ONLY : nind, ind2, omv, sites_x, sites_y,nf,&
           site_ind2,site_ind
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: oms
      REAL(SP), INTENT(OUT) :: stemp, dtemp
      REAL(SP), INTENT(IN) :: s_in(:,:), d_in(:,:)
      INTEGER, INTENT(IN) :: sud, freq ! Site and frequency under study
      REAL(SP) :: sv, tv, uv, om1s, om1, om2s, om2
      INTEGER :: si, ti, ui
      INTEGER :: s1i,s2i,u1i,u2i,s1in,s2in,u1in,u2in
      REAL(SP) :: s1v,s2v,u1v,u2v,s1,s2,d1,d2
      INTEGER :: m, i1, i2
      REAL(SP) :: rr1,ss1,rr2,ss2
      integer :: abstand1x,abstand1y,abstand2x,abstand2y,n,sudy,sudx
      
      stemp=0._sp
      dtemp=0._sp
      
      si=nind(freq,1)          
      sv=omv(si)
      ti=nind(freq,2)
      tv=omv(ti)
      ui=nind(freq,3)    
      uv=omv(ui)
      
      om1s=(sv+tv+uv)/2._sp  ! omega_1 strich etc.
      om1=(sv-tv+uv)/2._sp
      om2s=(sv-tv-uv)/2._sp
      om2=(sv+tv-uv)/2._sp

      s1v=om1s+oms
      s1i=find_ind_down(s1v)
      s1in=min(s1i+1,nf)
      
      u1v=om1-oms
      u1i=find_ind_down(u1v)
      u1in=min(u1i+1,nf)
      
      s2v=om2+oms
      s2i=find_ind_down(s2v)
      s2in=min(s2i+1,nf)
            
      u2v=-om2s+oms
      u2i=find_ind_down(u2v)
      u2in=min(u2i+1,nf)
       
      if ((s1in>s1i)) then ! The boundaries
         rr1=(abs(s1v)-omv(s1i))/(omv(s1in)-omv(s1i)) ! The scaled variable
      else
         rr1=0._sp ! Evaluate simply with x==0 (should be arbitrary)
      end if
      if ((u1in>u1i)) then
         ss1=(abs(u1v)-omv(u1i))/(omv(u1in)-omv(u1i)) ! The scaled variable
      else
         ss1=0._sp
      end if
      

      if ((s2in>s2i)) then ! The boundaries
         rr2=(abs(s2v)-omv(s2i))/(omv(s2in)-omv(s2i)) ! The scaled variable
      else
         rr2=0._sp ! Evaluate simply with x==0 (should be arbitrary)
      end if
      if ((u2in>u2i)) then
         ss2=(abs(u2v)-omv(u2i))/(omv(u2in)-omv(u2i)) ! The scaled variable
      else
         ss2=0._sp
      end if
      
!      DO m=1,sites 
!         i1=min(m-1,sites-(m-1))+1
!         i2=min(abs(sud-m),sites-abs(sud-m))+1
      sudy=site_ind(sud,1)
      sudx=site_ind(sud,2)
      do m=1,sites_y ! In the y-direction
         abstand1y=min(m-1,sites_y-(m-1))
         abstand2y=min(abs(sudy-m),sites_y-abs(sudy-m))
         do n=1,sites_x ! In the x-direction
            abstand1x=min(n-1,sites_x-(n-1))
            abstand2x=min(abs(sudx-n),sites_x-abs(sudx-n))
            ! First the y-direction (smaller), then the x-direction (bigger)
            i1=site_ind2(abstand1y+1,abstand1x+1)
            i2=site_ind2(abstand2y+1,abstand2x+1)
            call blend_102(rr1,ss1,s_in(ind2(s1i,ti,u1i),i1),s_in(ind2(s1in,ti,u1i),i1),&
                 s_in(ind2(s1i,ti,u1in),i1),s_in(ind2(s1in,ti,u1in),i1),s1)
            call blend_102(rr2,ss2,s_in(ind2(s2i,ti,u2i),i2),s_in(ind2(s2in,ti,u2i),i2),&
                 s_in(ind2(s2i,ti,u2in),i2),s_in(ind2(s2in,ti,u2in),i2),s2)
            stemp=stemp+2*s1*s2
            call blend_102(rr1,ss1,d_in(ind2(s1i,ti,u1i),i1),d_in(ind2(s1in,ti,u1i),i1),&
                 d_in(ind2(s1i,ti,u1in),i1),d_in(ind2(s1in,ti,u1in),i1),d1)
            call blend_102(rr2,ss2,d_in(ind2(s2i,ti,u2i),i2),d_in(ind2(s2in,ti,u2i),i2),&
                 d_in(ind2(s2i,ti,u2in),i2),d_in(ind2(s2in,ti,u2in),i2),d2)
            dtemp=dtemp+2*d1*d2
         END DO
      END do

      call blend_102(rr1,ss1,s_in(ind2(s1i,ti,u1i),sud),s_in(ind2(s1in,ti,u1i),sud),&
           s_in(ind2(s1i,ti,u1in),sud),s_in(ind2(s1in,ti,u1in),sud),s1)
      call blend_102(rr2,ss2,s_in(ind2(s2i,u2i,ti),1),s_in(ind2(s2in,u2i,ti),1),&
           s_in(ind2(s2i,u2in,ti),1),s_in(ind2(s2in,u2in,ti),1),s2)
      call blend_102(rr2,ss2,d_in(ind2(s2i,u2i,ti),1),d_in(ind2(s2in,u2i,ti),1),&
           d_in(ind2(s2i,u2in,ti),1),d_in(ind2(s2in,u2in,ti),1),d2)
      stemp=stemp+s1*(s2-d2)
      call blend_102(rr1,ss1,d_in(ind2(s1i,ti,u1i),sud),d_in(ind2(s1in,ti,u1i),sud),&
           d_in(ind2(s1i,ti,u1in),sud),d_in(ind2(s1in,ti,u1in),sud),d1)
      dtemp=dtemp-d1*(3*s2+d2)
      call blend_102(rr1,ss1,s_in(ind2(s1i,u1i,ti),1),s_in(ind2(s1in,u1i,ti),1),&
           s_in(ind2(s1i,u1in,ti),1),s_in(ind2(s1in,u1in,ti),1),s1)
      call blend_102(rr1,ss1,d_in(ind2(s1i,u1i,ti),1),d_in(ind2(s1in,u1i,ti),1),&
           d_in(ind2(s1i,u1in,ti),1),d_in(ind2(s1in,u1in,ti),1),d1)
      call blend_102(rr2,ss2,s_in(ind2(s2i,ti,u2i),sud),s_in(ind2(s2in,ti,u2i),sud),&
           s_in(ind2(s2i,ti,u2in),sud),s_in(ind2(s2in,ti,u2in),sud),s2)
      stemp=stemp+(s1-d1)*s2
      call blend_102(rr2,ss2,d_in(ind2(s2i,ti,u2i),sud),d_in(ind2(s2in,ti,u2i),sud),&
           d_in(ind2(s2i,ti,u2in),sud),d_in(ind2(s2in,ti,u2in),sud),d2)
      dtemp=dtemp-(3*s1+d1)*d2
      
    END SUBROUTINE CALC_TCH_I 

    SUBROUTINE CALC_UCH_I(oms,stemp,dtemp,s,d,freq)
      USE nrtype
      USE globalvariables, ONLY : nind, ind2, omv,nf
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: oms
      REAL(SP), INTENT(OUT) :: stemp, dtemp
      REAL(SP), INTENT(IN) :: s(:,:,:), d(:,:,:)
      INTEGER, INTENT(IN) :: freq ! Site and frequency under study
      REAL(SP) :: sv, tv, uv, om1s, om1, om2s, om2
      INTEGER :: si, ti, ui
      INTEGER :: s1i,s2i,t1i,t2i,s1in,s2in,t1in,t2in
      REAL(SP) :: s1v,s2v,t1v,t2v, s1,s2,d1,d2
      REAL(SP) :: rr1, ss1,rr2,ss2

      s1=0._sp;s2=0._sp;d1=0._sp;d2=0._sp
      
      stemp=0._sp
      dtemp=0._sp
      
      si=nind(freq,1)          
      sv=omv(si)
      ti=nind(freq,2)
      tv=omv(ti)
      ui=nind(freq,3)    
      uv=omv(ui)
      
      om1s=(sv+tv+uv)/2._sp  ! omega_1 strich etc.
      om1=(sv-tv+uv)/2._sp
      om2s=(sv-tv-uv)/2._sp
      om2=(sv+tv-uv)/2._sp

      s1v=om2s-oms  
      s1i=FIND_IND_DOWN(abs(s1v))
      s1in=min(s1i+1,nf) ! The next index or the largest one
      
      t1v=-om1-oms
      t1i=FIND_IND_DOWN(abs(t1v))
      t1in=min(t1i+1,nf) 
      
      s2v=om2-oms
      s2i=FIND_IND_DOWN(abs(s2v))
      s2in=min(s2i+1,nf) 
      
      t2v=om1s+oms
      t2i=FIND_IND_DOWN(abs(t2v))
      t2in=min(t2i+1,nf) 
      
      if ((s1in>s1i)) then ! The boundaries
         rr1=(abs(s1v)-omv(s1i))/(omv(s1in)-omv(s1i)) ! The scaled variable
      else
         rr1=0._sp ! Evaluate simply with x==0 (should be arbitrary)
      end if
      if ((t1in>t1i)) then
         ss1=(abs(t1v)-omv(t1i))/(omv(t1in)-omv(t1i)) ! The scaled variable
      else
         ss1=0._sp
      end if
      
      call blend_102(rr1,ss1,s(s1i,t1i,ui),s(s1in,t1i,ui),s(s1i,t1in,ui),s(s1in,t1in,ui),s1)
      call blend_102(rr1,ss1,d(s1i,t1i,ui),d(s1in,t1i,ui),d(s1i,t1in,ui),d(s1in,t1in,ui),d1)
      
      if ((s2in>s2i)) then ! The boundaries
         rr2=(abs(s2v)-omv(s2i))/(omv(s2in)-omv(s2i)) ! The scaled variable
      else
         rr2=0._sp ! Evaluate simply with x==0 (should be arbitrary)
      end if
      if ((t2in>t2i)) then
         ss2=(abs(t2v)-omv(t2i))/(omv(t2in)-omv(t2i)) ! The scaled variable
      else
         ss2=0._sp
      end if

      call blend_102(rr2,ss2,s(s2i,t2i,ui),s(s2in,t2i,ui),s(s2i,t2in,ui),s(s2in,t2in,ui),s2)
      call blend_102(rr2,ss2,d(s2i,t2i,ui),d(s2in,t2i,ui),d(s2i,t2in,ui),d(s2in,t2in,ui),d2)
      
      stemp=-s1*(2*s2+d2)
      stemp=stemp-d1*s2
       
      dtemp=-3*s1*s2
      dtemp=dtemp-d1*d2
       
    END SUBROUTINE CALC_UCH_I

    SUBROUTINE CALC_G_I(gout,g_in,s_in,d_in,lambda)
      use nrtype
      use globalvariables, only: omv, nf, ind2, gnf, gomv, site_ind2, sites_y, sites_x
      implicit none
      REAL(SP), INTENT(IN) :: g_in(:), s_in(:,:), d_in(:,:), lambda
      REAL(SP), INTENT(OUT) :: gout(size(g_in))
      REAL(SP) :: gtemp, om1, om2, rr, ss, s1, s2, d1, d2
      INTEGER :: om1i,om1in,om2i,om2in, abstand ,j
      INTEGER :: m,n,abstandx,abstandy,i1

      s1=0._sp;s2=0._sp;d1=0._sp;d2=0._sp

      DO j=2,gnf ! Loop over the frequencies for gout (the first fixed at zero)
         gtemp=0.0_sp
         om1=gomv(j)+lambda
         om2=gomv(j)-lambda
         om1i=FIND_IND_down(om1)
         om1in=min(om1i+1,nf)
         om2i=FIND_IND_down(om2)
         om2in=min(om2i+1,nf)
         if (om1in>om1i) then
            rr=(abs(om1)-omv(om1i))/(omv(om1in)-omv(om1i))
         else
            rr=0._sp ! arbitrary
         end if
         if (om2in>om2i) then
            ss=(abs(om2)-omv(om2i))/(omv(om2in)-omv(om2i))
         else
            ss=0._sp ! arbitrary
         end if
         
!         DO abstand=0,sites-1
!            i=MIN(abstand,sites-abstand) ! The actual relevant distance for a 1D chain
         DO m=1,sites_y
            abstandy=min(m-1,sites_y-(m-1))
            DO n=1,sites_x
               abstandx=min(n-1,sites_x-(n-1))
               i1=site_ind2(abstandy+1,abstandx+1)
               call blend_102(rr,ss,d_in(ind2(om1i,1,om2i),i1),d_in(ind2(om1in,1,om2i),i1),&
                    d_in(ind2(om1i,1,om2in),i1),d_in(ind2(om1in,1,om2in),i1),d1)
               call blend_102(ss,rr,d_in(ind2(om2i,1,om1i),i1),d_in(ind2(om2in,1,om1i),i1),&
                    d_in(ind2(om2i,1,om1in),i1),d_in(ind2(om2in,1,om1in),i1),d2)
               gtemp=gtemp-2*(d1-d2)
            END DO
         END DO
         CALL blend_102(rr,ss,s_in(ind2(om1i,om2i,1),1),s_in(ind2(om1in,om2i,1),1),&
              s_in(ind2(om1i,om2in,1),1),s_in(ind2(om1in,om2in,1),1),s1)
         CALL blend_102(ss,rr,s_in(ind2(om2i,om1i,1),1),s_in(ind2(om2in,om1i,1),1),&
              s_in(ind2(om2i,om1in,1),1),s_in(ind2(om2in,om1in,1),1),s2)
         gtemp=gtemp+3*(s1-s2)
         CALL blend_102(rr,ss,d_in(ind2(om1i,om2i,1),1),d_in(ind2(om1in,om2i,1),1),&
              d_in(ind2(om1i,om2in,1),1),d_in(ind2(om1in,om2in,1),1),d1)
         CALL blend_102(ss,rr,d_in(ind2(om2i,om1i,1),1),d_in(ind2(om2in,om1i,1),1),&
              d_in(ind2(om2i,om1in,1),1),d_in(ind2(om2in,om1in,1),1),d2)
         gtemp=gtemp+(d1-d2)
         gout(j)=gtemp
      END DO
      
      gout=gout/(lambda+interpg2(lambda,g_in)) ! With full interpolation
      
    END SUBROUTINE CALC_G_I
    
    SUBROUTINE vec2mat(vec,mat,m,n)
      USE nrtype
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m,n
      REAL(SP),INTENT(IN) :: vec(m*n)
      REAL(SP),INTENT(OUT) :: mat(m,n)
      INTEGER :: lask, i,j
      
      lask=0
      DO j=1,n
         DO i=1,m
            lask=lask+1
            mat(i,j)=vec(lask)
         END DO
      END DO
    END SUBROUTINE VEC2MAT
    
    FUNCTION find_ind_DOWN(omega) result(index)
      USE nrtype
      USE globalvariables, ONLY: omv, nf, aa, om0, logaa
      use nrutil, only: assert
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: omega
      INTEGER :: index
      INTRINSIC LOG

      ! The automatic flooring of a _positive_ real number to integer
      index=floor(LOG(1+(aa-1)*abs(omega)/om0)/logaa)+1 ! Note the plus one!
      index=min(index,nf)
     
      if (index>1 .and. (omv(index)>abs(omega))) then
         index=index-1 ! This should never happen, would lead to a (small) negative rr=-tiny
      elseif (index<nf.and.(omv(index+1)<abs(omega))) then
         index=index+1 ! This neither, would lead to rr=1+tiny
      end if
      
    END FUNCTION FIND_IND_DOWN   

     FUNCTION find_ind_DOWN_g(omega) result(index)
      USE nrtype
      USE globalvariables, ONLY: gomv, gnf, gaa, gom0, log_gaa
      use nrutil, only: assert
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: omega
      INTEGER :: index
      INTRINSIC LOG

      ! The automatic flooring of a _positive_ real number to integer
      index=floor(LOG(1+(gaa-1)*abs(omega)/gom0)/log_gaa)+1 ! Note the plus one!
      index=min(index,gnf)
     
      if (index>1 .and. (gomv(index)>abs(omega))) then
         index=index-1 ! This should never happen, would lead to a (small) negative rr=-tiny
      elseif (index<gnf.and.(gomv(index+1)<abs(omega))) then
         index=index+1 ! This neither, would lead to rr=1+tiny
      end if
      
    END FUNCTION FIND_IND_DOWN_g 
   
    FUNCTION pfunc_i(om1,om2,lambda,g_in) result(value)
      USE nrtype 
      USE globalvariables, ONLY: omv, nf
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: om1,om2,lambda
      REAL(SP) :: value
      REAL(SP), INTENT(IN):: g_in(:) 
      REAL(SP) :: g1, g2, rr, ga, gb
      INTEGER :: om1i, om2i, om1in, om2in
      REAL(SP) :: tiny=1e-15_sp
      value=0._sp
      if (abs(abs(om2)-abs(om1))<tiny) then ! The arguments are equal
         g1=interpg(abs(om1),g_in)
         g2=interpg(abs(om2),g_in)
         value=1.0_sp/((abs(om1)+g1)*&
              (abs(om2)+g2))/2 ! The one-half from the combination of delta and theta
         value=value*abs(om1)/om1*abs(om2)/om2
      ELSEIF (abs(om2)>lambda) THEN
         g1=interpg(abs(om1),g_in)
         g2=interpg(abs(om2),g_in)
         value=1.0_sp/((ABS(om1)+g1)&
              *(abs(om2)+g2))
         value=value*abs(om1)/om1*abs(om2)/om2
      ELSE
         value=0.0_sp
      END IF
      
    END FUNCTION PFUNC_I

    function interpg(om,g_in) result(val)
      ! Returns the interpolated value for _positive_ omegas
      USE nrtype
      USE globalvariables, ONLY: gomv, gnf
      IMPLICIT NONE
      REAL(SP),INTENT(IN) :: om, g_in(:)
      REAL(SP) :: val, rr, ga, gb
      INTEGER :: omi, omin

      omi=FIND_IND_down_g(abs(om))
      omin=min(omi+1,gnf)
      if ((omi==1).and.(.true.)) then ! Do not interpolate from zero
         !ga=0._sp
         !gb=g_in(omin)
         !call blend_101(abs(om)/omv(1),ga,gb,value)
         val=g_in(2) ! No interpolation from zero => stronger damping for lambda->0
         return
      end if
      
      if (omin>omi) then
         rr=(abs(om)-gomv(omi))/(gomv(omin)-gomv(omi))
      else
         rr=0._sp ! arbitrary
      end if
      ga=g_in(omi)
      gb=g_in(omin)
      call blend_101(rr,ga,gb,val)
 
    end function interpg

    function interpg2(om,g_in) result(val) ! For the self-energy derivative, uses interpolation
      ! Returns the interpolated value for _positive_ omegas
      USE nrtype
      USE globalvariables, ONLY: gomv, gnf
      IMPLICIT NONE
      REAL(SP),INTENT(IN) :: om, g_in(:)
      REAL(SP) :: val, rr, ga, gb
      INTEGER :: omi, omin

      omi=FIND_IND_down_g(abs(om))
      omin=min(omi+1,gnf)
      if ((omi==1).and.(.false.)) then ! Interpolate from zero
         !ga=0._sp
         !gb=g_in(omin)
         !call blend_101(abs(om)/omv(1),ga,gb,value)
         val=g_in(2) ! No interpolation from zero => stronger damping for lambda->0
         return
      end if
      
      if (omin>omi) then
         rr=(abs(om)-gomv(omi))/(gomv(omin)-gomv(omi))
      else
         rr=0._sp ! arbitrary
      end if
      ga=g_in(omi)
      gb=g_in(omin)
      call blend_101(rr,ga,gb,val)
 
    end function interpg2
      

    function func_to_int_v(x,flag_stu,flag_sd,fval,freq,site) result(value)
      USE nrtype
      USE globalvariables, only: s_in=>s_in2, d_in=>d_in2, g_in=>g_in2, gout=>gout2,&
           omv, v_in=>v_in2
      IMPLICIT NONE
      REAL(SP),INTENT(IN) :: x(:), fval
      INTEGER, INTENT(IN) :: flag_stu, flag_sd, freq, site
      REAL(SP) :: value(size(x))
      REAL(SP) :: vtemp
      INTEGER :: k
      REAL(SP) :: a,b, g1, g2

      value=0._sp
      
      do k=1,size(x)
         vtemp=0._sp
         a=x(k) ! NU
         b=x(k)+fval ! nu+omega
         CALL CALC_VCH_I(a,vtemp,s_in,d_in,v_in,site,freq)
            
         g1=interpg(abs(a),g_in)
         g2=interpg(abs(b),g_in)
         vtemp=vtemp*(interpg(abs(a),gout)/(abs(a)+g1)+&
              interpg(abs(b),gout)/(abs(b)+g2))
         vtemp=vtemp/((abs(a)+g1)*(abs(b)+g2))
         vtemp=vtemp*abs(a)/a*abs(b)/b ! Odd in both variables
         value(k)=vtemp
      end do
    

    end function func_to_int_v

 function func_to_int_vrho(x,flag_stu,flag_sd,fval,freq,site) result(value)
      USE nrtype
      USE globalvariables, only: s_in=>s_in2, d_in=>d_in2, g_in=>g_in2, gout=>gout2,&
           omv, v_in=>v_in2
      IMPLICIT NONE
      REAL(SP),INTENT(IN) :: x(:), fval
      INTEGER, INTENT(IN) :: flag_stu, flag_sd, freq, site
      REAL(SP) :: value(size(x))
      REAL(SP) :: vtemp
      INTEGER :: k
      REAL(SP) :: a,b, g1, g2

      value=0._sp
      
      do k=1,size(x)
         vtemp=0._sp
         a=x(k) ! NU
         b=x(k)+fval ! nu+omega
         CALL CALC_VRHO_I(a,vtemp,s_in,d_in,v_in,site,freq)
            
         g1=interpg(abs(a),g_in)
         g2=interpg(abs(b),g_in)
         vtemp=vtemp*(interpg(abs(a),gout)/(abs(a)+g1)+&
              interpg(abs(b),gout)/(abs(b)+g2))
         vtemp=vtemp/((abs(a)+g1)*(abs(b)+g2))
         vtemp=vtemp*abs(a)/a*abs(b)/b ! Odd in both variables
         value(k)=vtemp
      end do
    

    end function func_to_int_vrho
    

    function func_to_int_chi(x,flag_stu,flag_sd,fval,freq,site) result(val)
      USE nrtype
      USE globalvariables, only: g_in=>g_in2, gout=>gout2,&
           omv, v_in=>v_in2
      IMPLICIT NONE
      REAL(SP),INTENT(IN) :: x(:), fval
      INTEGER, INTENT(IN) :: flag_stu, flag_sd, freq, site
      REAL(SP) :: val(size(x))
      REAL(SP) :: chitemp
      INTEGER :: k
      REAL(SP) :: a,b, g1, g2

      val=0._sp
      
      do k=1,size(x)
         chitemp=0._sp
         a=x(k) ! NU
         b=x(k)+fval ! nu+omega
         CALL CALC_CHI_I(a,chitemp,v_in,freq,site)
            
         g1=interpg(abs(a),g_in)
         g2=interpg(abs(b),g_in)
         chitemp=chitemp*(interpg(abs(a),gout)/(abs(a)+g1)+&
              interpg(abs(b),gout)/(abs(b)+g2))
         chitemp=chitemp/((abs(a)+g1)*(abs(b)+g2))
         chitemp=chitemp*abs(a)/a*abs(b)/b ! Odd in both variables
         val(k)=chitemp
      end do
    

    end function func_to_int_chi


 
    function func_to_int(x,flag_stu,flag_sd,fval,freq,site) result(value)
      USE nrtype
      USE globalvariables, only: s_in=>s_in2, d_in=>d_in2, g_in=>g_in2, gout=>gout2,&
           omv, smat=>smat2, dmat=>dmat2
      ! smat and dmat relevant for a constant site (no parallelization over sites!)
      IMPLICIT NONE
      REAL(SP),INTENT(IN) :: x(:), fval
      INTEGER, INTENT(IN) :: flag_stu, flag_sd, freq, site
      REAL(SP) :: value(size(x))
      REAL(SP) :: stemp, dtemp, vtemp
      INTEGER :: k
      REAL(SP) :: a,b, g1, g2

      value=0._sp
      
      do k=1,size(x)
         stemp=0._sp
         dtemp=0._sp
         vtemp=0._sp
         a=x(k) ! omega'
         b=x(k)+fval ! s+omega' or t+omega' or u+omega'
         if (flag_stu==1) then
            !CALL CALC_SCH(a,stemp,dtemp,s_in,d_in,site,freq)
            CALL CALC_SCH_I(a,stemp,dtemp,smat,dmat,freq)
         elseif (flag_stu==2) then
            CALL CALC_TCH_I(a,stemp,dtemp,s_in,d_in,site,freq)
         elseif (flag_stu==3) then
            !CALL CALC_UCH(a,stemp,dtemp,s_in,d_in,site,freq)
            CALL CALC_UCH_I(a,stemp,dtemp,smat,dmat,freq)
         else
             stop
         end if
                  
         if (flag_sd==1) then
            vtemp=stemp
         elseif (flag_sd==2) then
            vtemp=dtemp
         else
            stop
         end if
         vtemp=vtemp*(interpg(abs(a),gout)/(abs(a)+interpg(abs(a),g_in))+&
              interpg(abs(b),gout)/(abs(b)+interpg(abs(b),g_in)))
         vtemp=vtemp/((abs(a)+interpg(abs(a),g_in))*&
              (abs(b)+interpg(abs(b),g_in)))
         vtemp=vtemp*abs(a)/a*abs(b)/b ! Odd in both variables
         value(k)=vtemp

      end do
    

    end function func_to_int
    
    FUNCTION qromo(func,a,b,choose,flag_stu,flag_sd,freq,fval,site)
      USE nrtype; USE nrutil, ONLY : nrerror
      USE integrators, ONLY : polint
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: a,b
      REAL(SP) :: qromo
      INTEGER, INTENT(IN) :: flag_stu, flag_sd, freq, site
      REAL(SP), INTENT(IN) :: fval
      INTERFACE
         FUNCTION func(x,flag_stu,flag_sd,fval,freq,site)
           USE nrtype
           IMPLICIT NONE
           REAL(SP), DIMENSION(:), INTENT(IN) :: x
           REAL(SP), DIMENSION(size(x)) :: func
           INTEGER, INTENT(IN) :: flag_stu, flag_sd, freq, site
           REAL(SP), INTENT(IN) :: fval
         END FUNCTION func
         !BL
         SUBROUTINE choose(funk,aa,bb,s,n,flag_stu,flag_sd,freq,fval, site)
           USE nrtype
           IMPLICIT NONE
           REAL(SP), INTENT(IN) :: aa,bb
           REAL(SP), INTENT(INOUT) :: s
           INTEGER(I4B), INTENT(IN) :: n
           REAL(SP), INTENT(IN) :: fval
           INTEGER, INTENT(IN) :: flag_stu, flag_sd, freq, site
           INTERFACE
              FUNCTION funk(x,flag_stu,flag_sd,fval,freq,site)
                USE nrtype
                IMPLICIT NONE
                REAL(SP), DIMENSION(:), INTENT(IN) :: x
                REAL(SP), INTENT(IN) :: fval
                INTEGER, INTENT(IN) :: flag_stu, flag_sd, freq,site
                REAL(SP), DIMENSION(size(x)) :: funk
              END FUNCTION funk
           END INTERFACE
         END SUBROUTINE choose
      END INTERFACE
      INTEGER(I4B), PARAMETER :: JMAX=14,JMAXP=JMAX+1,K=3,KM=K-1
      REAL(SP), PARAMETER :: EPS=1.0e-1
      REAL(SP), DIMENSION(JMAXP) :: h,s
      REAL(SP) :: dqromo
      INTEGER(I4B) :: j
      h(1)=1.0
      do j=1,JMAX
         call choose(func,a,b,s(j),j,flag_stu,flag_sd,freq,fval,site)
         if (j >= K) then
            call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromo,dqromo)
            if (abs(dqromo) <= EPS*max(abs(qromo),1e-2_sp)) then
               !if (freq==1.and.(abs(min(abs(a),abs(b))-0.0369)<1e-2).and.site==2) then
               !   print*, min(abs(a),abs(b)), qromo,flag_sd,flag_stu
               !end if
               RETURN
            end if
          end if
         s(j+1)=s(j)

         h(j+1)=h(j)/9.0_sp
      end do
      write(*,*) 'Problems in the integration, error=', abs(dqromo)/abs(qromo), ', value=', qromo
       open(99,file='integrand_fail.dat')
       do j=1,100
         write(99,*)  a+(b-a)*j/100._sp, func( a+(b-a)* (/ (j/100._sp) /),flag_stu,flag_sd,fval,freq,site )
      end do
      write(*,*) 'Wrote to file integrand_fail.dat.'
      close(99)
    END FUNCTION qromo  

    function funktio(a) result(val)
      use nrtype
      use globalvariables, only: om0, om_max, nf
      implicit none
      REAL(SP), intent(in):: a
      REAL(SP):: val
      val=om0/om_max*(a**(nf-1)-1)/(a-1)-1
    end function funktio

    function gfunktio(ga) result(val)
      use nrtype
      use globalvariables, only: gom0, gom_max, gnf
      implicit none
      REAL(SP), intent(in):: ga
      REAL(SP):: val
      val=gom0/gom_max*(ga**(gnf-1)-1)/(ga-1)-1
    end function gfunktio
    

  END MODULE derivatives
