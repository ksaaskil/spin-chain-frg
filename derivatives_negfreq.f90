MODULE derivatives
PRIVATE
PUBLIC:: derivs
CONTAINS

  SUBROUTINE derivs2(x,y,dydx)
    USE nrtype
    IMPLICIT NONE
    REAL(SP),INTENT(IN) :: x
    REAL(SP),DIMENSION(:),INTENT(IN) :: y
    REAL(SP),DIMENSION(:),INTENT(OUT):: dydx

    dydx=0._sp
    dydx(1)=y(1)
    dydx(2)=3*y(2)

  END SUBROUTINE derivs2

  SUBROUTINE derivs(x,y,dydx)
    ! Derivatives for a linear frequency mesh symmetric around zero
    USE nrtype
    USE globalvariables, ONLY: nf,sites,freqs,FO,LO
    IMPLICIT NONE
    !INTEGER,INTENT(IN) :: nvar,n,sites,freqs
    REAL(SP),INTENT(IN) :: x
    REAL(SP),DIMENSION(:),INTENT(IN) :: y
    REAL(SP),DIMENSION(:),INTENT(OUT):: dydx
    INTEGER :: nvar
    INTEGER :: indg(nf),inds(sites*freqs),indd(sites*freqs)
    INTEGER :: i,j,k,lask,m
    INTEGER :: ind(freqs,3),ind2(nf,nf,nf)
    REAL(SP) :: gout(nf), sout(freqs,sites),dout(freqs,sites) ! The outgoing derivatives of selfenergies and vertex functions
    REAL(SP) :: stemp, dtemp
    REAL(SP) :: omv(nf),sv(nf),tv(nf),uv(nf) ! The frequency meshes
    REAL(SP) :: g_in(nf),s_in(freqs,sites),d_in(freqs,sites) ! vertex function matrices
    REAL(SP) :: oms, s1v, t1v, u1v, s2v, t2v, u2v, om1s, om1, om2s, om2, sval,tval,uval
    INTEGER :: s1i,s2i,t1i,t2i,u1i,u2i,om1i,om2i,lambda_ind,si,ti,ui
    INTRINSIC EXP, ATAN, ISNAN
    REAL(SP) :: lambda
    REAL(SP), PARAMETER :: TINY=1.0e-15_sp
    REAL(SP), PARAMETER :: hf=2*FO/(nf-1)
    REAL(SP) :: incr

    gout=0._sp
    sout=0._sp
    dout=0._sp

    g_in=0._sp
    s_in=0._sp
    d_in=0._sp

    
    nvar=size(y)
    dydx=0.0_sp
    DO j=1,nvar
       IF (ISNAN(y(j))) THEN
          write(*,*) 'error, encountered NAN!'
          write(*,*) x,j
          stop
       END IF
    END DO
    !write(*,*) maxval(abs(y))
          
    lambda=LO*exp(-x)
    !lambda=LO-x
    !write(*,*) lambda
    
    if (mod(nf,2)>0) THEN
       write(*,*) 'error, odd nf!'
       stop
    end if

    omv=0.0_sp
    DO i=1,nf
       indg(i)=i
       omv(i)=-FO+(i-1)*hf
    END DO
    g_in=y(indg)

    !omv=omv+1e-3_sp
    sv=omv
    tv=omv
    uv=omv

    incr=maxval(omv)

    j=0
    DO i=nf+1,nf+sites*freqs
       j=j+1
       inds(j)=i
    END DO

    ! Construct sin, the incoming matrix for spin interactions
    CALL VEC2MAT(y(inds),s_in,freqs,sites)
    !write(*,*) size(y(inds))

    j=0
    DO i=nf+sites*freqs+1,nvar
       j=j+1
       indd(j)=i
    END DO
    
    ! Construct din, the incoming matrix for density interactions
    CALL VEC2MAT(y(indd),d_in,freqs,sites)
    
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

    !i=NINT(0.50_sp)
    
    !---------------------------------------------------------------------

    DO j=1,nf ! Loop over the frequencies for gout
       stemp=0.0_sp
       om1i=MAX(MIN(NINT((omv(j)+lambda+incr)/hf)+1,nf),1)
       !index1=ind2(om1i,1,om2i)
       om2i=MAX(MIN(NINT((omv(j)-lambda+incr)/hf)+1,nf),1)
       !index2=ind2(omi,1,om1i)
       DO i=1,sites
          stemp=stemp-2*(d_in(ind2(om1i,1,om2i),i)-d_in(ind2(om2i,1,om1i),i))
       END DO
       stemp=stemp+3*(s_in(ind2(om1i,om2i,1),1)-s_in(ind2(om2i,om1i,1),1))
       stemp=stemp+d_in(ind2(om1i,om2i,1),1)-d_in(ind2(om2i,om1i,1),1)
       gout(j)=stemp
    END DO
    lambda_ind=MIN(NINT((lambda+incr)/hf)+1,nf)
    !gout=gout/(2*pi)
    gout=gout/(lambda+g_in(lambda_ind))
    
    

    DO j=1,sites
       DO i=1,freqs
          !write(*,*) 1
          si=ind(i,1)
          
          !write(*,*) si
          sval=sv(si)
          ti=ind(i,2)
          
          !write(*,*) ti
          tval=tv(ti)
          ui=ind(i,3)
          !write(*,*) ui
          uval=uv(ui)

          !write(*,*) uval
          
          om1s=(sval+tval+uval)/2._sp  ! omega_1 strich etc.
          om1=(sval-tval+uval)/2._sp
          om2s=(sval-tval-uval)/2._sp
          om2=(sval+tval-uval)/2._sp
          
          !---------------------------------------------------------------------
          !                      s-channel
          !---------------------------------------------------------------------
          !if (.false.) then
          DO k=1,4
             stemp=0.0_sp
             dtemp=0.0_sp
             IF (k==1) THEN
                oms=lambda
             ELSEIF (k==2) THEN
                oms=-lambda
             ELSEIF (k==3) THEN
                oms=lambda-sval
             ELSEIF (k==4) THEN
                oms=-lambda-sval
             END IF
             !write(*,*) 2
             t1v=-om2s-oms
             t1i=MAX(MIN(NINT((t1v+incr)/hf)+1,nf),1)
             !write(*,*) t1i
             u1v=om1s+oms
             u1i=MAX(MIN(NINT((u1v+incr)/hf)+1,nf),1)
             !write(*,*) u1i
             t2v=om2+oms
             t2i=MAX(MIN(NINT((t2v+incr)/hf)+1,nf),1)
             !write(*,*) t2i
             u2v=om1+oms
             u2i=MAX(MIN(NINT((u2v+incr)/hf)+1,nf),1)
             !write(*,*) u2i
             
             stemp=s_in(ind2(si,t1i,u1i),j)*(-2*s_in(ind2(si,t2i,u2i),j)+d_in(ind2(si,t2i,u2i),j))
             stemp=stemp+d_in(ind2(si,t1i,u1i),j)*s_in(ind2(si,t2i,u2i),j)
             
             IF (k==1.OR.k==2) THEN
                stemp=stemp*pfunc(oms,sval+oms,lambda,g_in,nf)
             ELSE
                stemp=stemp*pfunc(sval+oms,oms,lambda,g_in,nf)
             END IF
             
             sout(i,j)=sout(i,j)+stemp
             
             dtemp=3*s_in(ind2(si,t1i,u1i),j)*s_in(ind2(si,t2i,u2i),j)
             dtemp=dtemp+d_in(ind2(si,t1i,u1i),j)*d_in(ind2(si,t2i,u2i),j)
             
             IF (k==1.OR.k==2) THEN
                dtemp=dtemp*pfunc(oms,sval+oms,lambda,g_in,nf)
             ELSE
                dtemp=dtemp*pfunc(sval+oms,oms,lambda,g_in,nf)
             END IF
             
             dout(i,j)=dout(i,j)+dtemp
             
          END DO
       !end if
          
          !---------------------------------------------------------------------
          !                      t-channel
          !---------------------------------------------------------------------
          
          !if (.false.) then
          DO k=1,4
             stemp=0.0_sp
             dtemp=0.0_sp
             IF (k==1) THEN
                oms=lambda
             ELSEIF (k==2) THEN
                oms=-lambda
             ELSEIF (k==3) THEN
                oms=lambda-tval
             ELSEIF (k==4) THEN
                oms=-lambda-tval
             END IF
             !write(*,*) 2
             s1v=om1s+oms
             s1i=MAX(MIN(NINT((s1v+incr)/hf)+1,nf),1)
             !write(*,*) t1i
             u1v=om1-oms
             u1i=MAX(MIN(NINT((u1v+incr)/hf)+1,nf),1)
             !write(*,*) u1i
             s2v=om2+oms
             s2i=MAX(MIN(NINT((s2v+incr)/hf)+1,nf),1)
             !write(*,*) t2i
             u2v=-om2s+oms
             u2i=MAX(MIN(NINT((u2v+incr)/hf)+1,nf),1)
             !write(*,*) u2i
             
             DO m=1,sites
                stemp=stemp+2*s_in(ind2(s1i,ti,u1i),m)*s_in(ind2(s2i,ti,u2i),abs(j-m)+1)
                dtemp=dtemp+2*d_in(ind2(s1i,ti,u1i),m)*d_in(ind2(s2i,ti,u2i),abs(j-m)+1)
             END DO
             stemp=stemp+s_in(ind2(s1i,ti,u1i),j)*(s_in(ind2(s2i,u2i,ti),1)-d_in(ind2(s2i,u2i,ti),1))
             stemp=stemp+(s_in(ind2(s1i,u1i,ti),1)-d_in(ind2(s1i,u1i,ti),1))*s_in(ind2(s2i,ti,u2i),j)
             
             IF (k==1.OR.k==2) THEN
                stemp=stemp*pfunc(oms,tval+oms,lambda,g_in,nf)
             ELSE
                stemp=stemp*pfunc(tval+oms,oms,lambda,g_in,nf)
             END IF
             
             sout(i,j)=sout(i,j)+stemp
             
             dtemp=dtemp-d_in(ind2(s1i,ti,u1i),j)*(3*s_in(ind2(s2i,u2i,ti),1)+d_in(ind2(s2i,u2i,ti),1))
             dtemp=dtemp-(3*s_in(ind2(s1i,u1i,ti),1)+d_in(ind2(s1i,u1i,ti),1))*d_in(ind2(s2i,ti,u2i),j)
             
             IF (k==1.OR.k==2) THEN
                dtemp=dtemp*pfunc(oms,tval+oms,lambda,g_in,nf)
             ELSE
                dtemp=dtemp*pfunc(tval+oms,oms,lambda,g_in,nf)
             END IF
             
             dout(i,j)=dout(i,j)+dtemp
             
          END DO
          
          !end if
          
          
          !---------------------------------------------------------------------
          !                      u-channel
          !---------------------------------------------------------------------      
    
         ! if (.false.) then
          DO k=1,4
             stemp=0.0_sp
             dtemp=0.0_sp
             IF (k==1) THEN
                oms=lambda
             ELSEIF (k==2) THEN
                oms=-lambda
             ELSEIF (k==3) THEN
                oms=lambda-uval
             ELSEIF (k==4) THEN
                oms=-lambda-uval
             END IF
             !write(*,*) 2
             s1v=om2s-oms
             s1i=MAX(MIN(NINT((s1v+incr)/hf)+1,nf),1)
             !write(*,*) t1i
             t1v=-om1-oms
             t1i=MAX(MIN(NINT((t1v+incr)/hf)+1,nf),1)
             !write(*,*) u1i
             s2v=om2-oms
             s2i=MAX(MIN(NINT((s2v+incr)/hf)+1,nf),1)
             !write(*,*) t2i
             t2v=om1s+oms
             t2i=MAX(MIN(NINT((t2v+incr)/hf)+1,nf),1)
             !write(*,*) u2i
             
             stemp=-s_in(ind2(s1i,t1i,ui),j)*(2*s_in(ind2(s2i,t2i,ui),j)+d_in(ind2(s2i,t2i,ui),j))
             stemp=stemp-d_in(ind2(s1i,t1i,ui),j)*s_in(ind2(s2i,t2i,ui),j)

             
             
             IF (k==1.OR.k==2) THEN
                stemp=stemp*pfunc(oms,uval+oms,lambda,g_in,nf)
             ELSE
                stemp=stemp*pfunc(uval+oms,oms,lambda,g_in,nf)
             END IF
             
             sout(i,j)=sout(i,j)+stemp
             
             dtemp=-3*s_in(ind2(s1i,t1i,ui),j)*s_in(ind2(s2i,t2i,ui),j)
             dtemp=dtemp-d_in(ind2(s1i,t1i,ui),j)*d_in(ind2(s2i,t2i,ui),j)
             
             IF (k==1.OR.k==2) THEN
                dtemp=dtemp*pfunc(oms,uval+oms,lambda,g_in,nf)
             ELSE
                dtemp=dtemp*pfunc(uval+oms,oms,lambda,g_in,nf)
             END IF
             
             !if (isnan(dtemp)) THEN
             !   write(*,*) 'error, NaN encountered in dtemp!'
             !   stop
             !end if
             dout(i,j)=dout(i,j)+dtemp
             
             
          END DO
          !end if
          !if (i==1) THEN
          !write(*,*) 'here',i,j
          !end if
          !write(*,*) 'here2'
          !write(*,*) ind(912-freqs-9,1), ind(912-freqs-9,2),ind(912-freqs-9,3)
       END DO
    END DO
    

    
    
    
    !DO j=1,nvar
    !  IF (ISNAN(dydx(j))) THEN
    !     write(*,*) 'error at the end4!'
    !     write(*,*) x,j
    !     stop
    !  END IF
    !END DO
    
    lask=0
    DO j=1,nf
       lask=lask+1
       dydx(j)=gout(j)
    END DO
    
    !DO j=1,nvar
    !   IF (ISNAN(dydx(j))) THEN
    !      write(*,*) 'error at the end3!'
    !      write(*,*) x,j
    !      stop
    !   END IF
    !END DO
    
    

    DO j=1,sites
       DO i=1,freqs
          lask=lask+1
          dydx(lask)=sout(i,j)
       END DO
    END DO
    
    ! DO j=1,nvar
    !   IF (ISNAN(dydx(j))) THEN
    !      write(*,*) 'error at the end2!'
    !      write(*,*) x,j
    !      stop
    !   END IF
    !END DO
    
    DO j=1,sites
       DO i=1,freqs
          lask=lask+1
          dydx(lask)=dout(i,j)
       END DO
    END DO


    DO j=1,nvar
       IF (ISNAN(dydx(j))) THEN
          write(*,*) 'error, encountered Nan in dydx!'
          write(*,*) 'x=',x,'var=',j
          stop
       END IF
    END DO
    
    dydx=dydx/(2*PI)
    dydx=dydx*(-lambda)
    !dydx=dydx*(-1)

    !write(*,*) omv

   

    if (abs(x)<TINY) then
       do i=1,nvar
          write(55,*) dydx(i)
       end do
       write(*,*) 'wrote to file, max(y)=',maxval(abs(y))
    end if
    
    !write(*,*) 'stopping in derivs'
    !stop

    ! stop
    !write(*,*) nvar

    !write(77,*) dydx
    
    
    !DO j=1,nvar
    !   IF (ISNAN(dydx(j))) THEN
    !      write(*,*) 'error at the end!'
    !      write(*,*) x,j
    !      stop
    !   END IF
    !END DO


  END SUBROUTINE derivs


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

  SUBROUTINE mat2vec(vec,mat,m,n)
    USE nrtype
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m,n
    REAL(SP),INTENT(OUT) :: vec(m*n)
    REAL(SP),INTENT(IN) :: mat(m,n)
    INTEGER :: lask, i,j

    lask=0
    DO j=1,n
       DO i=1,m
          lask=lask+1
          vec(lask)=mat(i,j)
       END DO
    END DO
  END SUBROUTINE MAT2VEC

  FUNCTION pfunc(om1,om2,lambda,g_in,nf) result(value)
    USE nrtype 
    USE globalvariables, ONLY: hf
    INTEGER :: nf
    REAL(SP):: om1,om2,lambda,value
    REAL(SP):: g_in(nf)
    value=0._sp
    !write(*,*) om2
    !REAL(SP):: hf
    !hf=FO/(nf-1)
    !write(*,*) 'hf is in pfunc ',hf
    IF (abs(om2)>lambda) THEN
       value=1.0_sp/((om1+g_in(MAX(MIN(NINT((om1+incr)/hf)+1,nf),1)))*(om2+g_in(MAX(MIN(NINT((om2+incr)/hf)+1,nf),1))))
       !value=value/(om2+g_in(MIN(NINT(abs(om2)/hf)+1,nf)))
    ELSE
       value=0.0_sp
    END IF
    
  END FUNCTION PFUNC

END MODULE derivatives
