MODULE globalvariables
  USE nrtype
  IMPLICIT NONE
  !REAL(SP), PARAMETER :: aa=2.065_sp ! The multiplier for the increase of distance
  REAL(SP) :: aa, logaa,gaa,log_gaa
  REAL(SP), PARAMETER :: om_max=250._sp, gom_max=250._sp
  REAL(SP), PARAMETER :: om0=1e-4_sp, gom0=1e-4
  !REAL(SP), PARAMETER :: aa=1.35
  !REAL(SP), PARAMETER :: logaa=log(aa)
  INTEGER, PARAMETER :: nf=37, gnf=37 ! The number of frequencies
  !REAL(SP), PARAMETER :: om0=(aa-1)*om_max/(aa**(nf-1)-1)
 ! REAL(SP), PARAMETER :: gom0=(aa-1)*om_max/(aa**(gnf-1)-1)
 ! INTEGER, PARAMETER :: sites=16
  INTEGER, PARAMETER :: sites_x=1, sites_y=2
  INTEGER, PARAMETER :: svar=(sites_x/2+1)*(sites_y/2+1) 
  !INTEGER, PARAMETER :: svar=sites/2+1 ! The number of needed site parameters
  INTEGER, PARAMETER :: freqs=nf*nf*(nf+1)/2, vfreqs=(nf)**2
  INTEGER, PARAMETER :: nvar=gnf+2*svar*freqs+2*nf*svar+2*(svar*vfreqs)+1 ! Not used that much
  REAL(SP), PARAMETER :: LO=10000.0_sp ! Initial cut-off
  INTEGER :: ind(nf*nf*nf,3), ind2(nf,nf,nf)
  INTEGER :: nind(freqs,3), nind2(nf,nf,nf)
  INTEGER :: vind(vfreqs,2), vind2(nf,nf)
  REAL(SP) :: omv(nf), gomv(gnf)
  REAL(SP) :: vomv(size(omv))
  REAL(SP) :: s_in2(nf*nf*nf,svar), d_in2(nf*nf*nf,svar), g_in2(gnf), gout2(gnf)
  REAL(SP) :: smat2(nf,nf,nf), dmat2(nf,nf,nf)
  REAL(SP) :: v_in2(vfreqs,svar)
  REAL(SP):: vrho_in2(vfreqs,svar)
  INTEGER :: site2, freq2
  INTEGER :: site_ind(svar,2), site_ind2(sites_y/2+1,sites_x/2+1)
  INTEGER :: flag_stu, flag_sd ! Some flags to help the integrator
  LOGICAL:: symmetry_flag
  
END MODULE globalvariables
