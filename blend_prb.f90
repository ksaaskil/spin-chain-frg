program blend_prb
!
!*******************************************************************************
!
!! BLEND_PRB tests routines from BLEND.
!
  implicit none
!
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BLEND_PRB'
  write ( *, '(a)' ) '  Test routines in BLEND.'
 
  call test01
  call test02
  call test03
  call test04
  call test05
  call test06
  call test07
  call test08
  call test09
  call test10
  call test11
  call test12
  call test13

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BLEND_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
subroutine test01
!
!*******************************************************************************
!
!! TEST01 checks for a gross error in the blend coefficients.
!
  implicit none
!
  integer, parameter :: nmax = 3
!
  integer i
  integer n
  real r
  real s
  real t
  real x(nmax)
!
  external identity_r
  external identity_rs
  external identity_rst
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Simple identity test to detect gross errors.'
!
!  Set N to 1.
!
  n = 1
!
!  Test BLEND_R_0DN on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Identity test for BLEND_R_0DN:'
  write ( *, '(a)' ) ' '
  r = 0.0E+00
  call blend_r_0dn ( r, x, n, identity_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  r = 1.0E+00
  call blend_r_0dn ( r, x, n, identity_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  r = 0.5E+00
  call blend_r_0dn ( r, x, n, identity_r )
  write ( *, '(2f8.4)' ) r, x(1:n)
!
!  Set N to 2.
!
  n = 2
!
!  Test BLEND_RS_0DN on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Identity test for BLEND_RS_0DN:'
  write ( *, '(a)' ) ' '
  r = 0.0E+00
  s = 0.0E+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0E+00
  s = 0.0E+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.0E+00
  s = 1.0E+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0E+00
  s = 1.0E+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.5E+00
  s = 0.5E+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)
!
!  Test BLEND_RS_1DN on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Identity test for BLEND_RS_1DN:'
  write ( *, '(a)' ) ' '
  r = 0.0E+00
  s = 0.0E+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0E+00
  s = 0.0E+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.0E+00
  s = 1.0E+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0E+00
  s = 1.0E+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.5E+00
  s = 0.5E+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)
!
!  Set N to 3.
!
  n = 3
!
!  Test BLEND_RST_0DN on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Identity test for BLEND_RST_0DN:'
  write ( *, '(a)' ) ' '
  r = 0.0E+00
  s = 0.0E+00
  t = 0.0E+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0E+00
  s = 0.0E+00
  t = 0.0E+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0E+00
  s = 1.0E+00
  t = 0.0E+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0E+00
  s = 0.0E+00
  t = 1.0E+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0E+00
  s = 1.0E+00
  t = 1.0E+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5E+00
  s = 0.5E+00
  t = 0.5E+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)
!
!  Test BLEND_RST_1DN on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Identity test for BLEND_RST_1DN:'
  write ( *, '(a)' ) ' '
  r = 0.0E+00
  s = 0.0E+00
  t = 0.0E+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0E+00
  s = 0.0E+00
  t = 0.0E+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0E+00
  s = 1.0E+00
  t = 0.0E+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0E+00
  s = 0.0E+00
  t = 1.0E+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0E+00
  s = 1.0E+00
  t = 1.0E+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5E+00
  s = 0.5E+00
  t = 0.5E+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)
!
!  Test BLEND_RST_2DN on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Identity test for BLEND_RST_2DN:'
  write ( *, '(a)' ) ' '
  r = 0.0E+00
  s = 0.0E+00
  t = 0.0E+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0E+00
  s = 0.0E+00
  t = 0.0E+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0E+00
  s = 1.0E+00
  t = 0.0E+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0E+00
  s = 0.0E+00
  t = 1.0E+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0E+00
  s = 1.0E+00
  t = 1.0E+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5E+00
  s = 0.5E+00
  t = 0.5E+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  return
end
subroutine test02
!
!*******************************************************************************
!
!! TEST02 checks for simple errors in the blend coefficients.
!
  implicit none
!
  integer, parameter :: nmax = 3
!
  integer i
  integer n
  real r
  real s
  real t
  real x(nmax)
!
  external stretch_r
  external stretch_rs
  external stretch_rst
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Shift and stretch test to detect simple errors.'
!
!  Set N to 1.
!
  n = 1
!
!  Test BLEND_R_0DN on shift by 1, stretch by 2.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Shift and stretch test for BLEND_R_0DN:'
  write ( *, '(a)' ) ' '
  r = 0.0E+00
  call blend_r_0dn ( r, x, n, stretch_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  r = 1.0E+00
  call blend_r_0dn ( r, x, n, stretch_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  r = 0.5E+00
  call blend_r_0dn ( r, x, n, stretch_r )
  write ( *, '(2f8.4)' ) r, x(1:n)
!
!  Set N to 2.
!
  n = 2
!
!  Test BLEND_RS_0DN on shift by (1,2), stretch by (3,4).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Shift and stretch test for BLEND_RS_0DN:'
  write ( *, '(a)' ) ' '
  r = 0.0E+00
  s = 0.0E+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0E+00
  s = 0.0E+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.0E+00
  s = 1.0E+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0E+00
  s = 1.0E+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.5E+00
  s = 0.5E+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)
!
!  Test BLEND_RS_1D on shift by (1,2), stretch by (3,4).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Shift and stretch test for BLEND_RS_1DN:'
  write ( *, '(a)' ) ' '
  r = 0.0E+00
  s = 0.0E+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0E+00
  s = 0.0E+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.0E+00
  s = 1.0E+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0E+00
  s = 1.0E+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.5E+00
  s = 0.5E+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)
!
!  Set N to 3.
!
  n = 3
!
!  Test BLEND_RST_0DN on shift by (1,2,3), stretch by (4,5,6).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Shift and stretch test for BLEND_RST_0DN:'
  write ( *, '(a)' ) ' '
  r = 0.0E+00
  s = 0.0E+00
  t = 0.0E+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0E+00
  s = 0.0E+00
  t = 0.0E+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0E+00
  s = 1.0E+00
  t = 0.0E+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0E+00
  s = 0.0E+00
  t = 1.0E+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0E+00
  s = 1.0E+00
  t = 1.0E+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5E+00
  s = 0.5E+00
  t = 0.5E+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)
!
!  Test BLEND_RST_1DN on shift by (1,2,3), stretch by (4,5,6).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Shift and stretch test for BLEND_RST_1DN:'
  write ( *, '(a)' ) ' '
  r = 0.0E+00
  s = 0.0E+00
  t = 0.0E+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0E+00
  s = 0.0E+00
  t = 0.0E+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0E+00
  s = 1.0E+00
  t = 0.0E+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0E+00
  s = 0.0E+00
  t = 1.0E+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0E+00
  s = 1.0E+00
  t = 1.0E+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)
  r = 0.5E+00
  s = 0.5E+00
  t = 0.5E+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)
!
!  Test BLEND_RST_2DN on shift by (1,2,3), stretch by (4,5,6).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Shift and stretch test for BLEND_RST_2DN:'
  write ( *, '(a)' ) ' '
  r = 0.0E+00
  s = 0.0E+00
  t = 0.0E+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0E+00
  s = 0.0E+00
  t = 0.0E+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0E+00
  s = 1.0E+00
  t = 0.0E+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0E+00
  s = 0.0E+00
  t = 1.0E+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0E+00
  s = 1.0E+00
  t = 1.0E+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5E+00
  s = 0.5E+00
  t = 0.5E+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  return
end
subroutine test03
!
!*******************************************************************************
!
!! TEST03 checks out BLEND_I_0D1.
!
  implicit none
!
  integer, parameter :: m = 5
!
  integer i
  real x(m)
!
  x(1) = 100.0E+00
  x(m) = 100.0 + real ( ( m - 1 ) * 5 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  BLEND_I_0D1 interpolates data in a vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  X(1) = ', x(1)
  write ( *, '(a,g14.6)' ) '  X(', m, ')= ', x(m)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Interpolated values:'
  write ( *, '(a)' ) ' '

  call blend_i_0d1 ( x, m )

  do i = 1, m
    write ( *, '(i6,g14.6)' ) i, x(i)
  end do

  return
end
subroutine test04
!
!*******************************************************************************
!
!! TEST04 checks out BLEND_IJ_0D1 and BLEND_IJ_1D1.
!
  implicit none
!
  integer, parameter :: m1 = 5
  integer, parameter :: m2 = 4
!
  integer i
  integer j
  integer n
  real r
  real s
  real x(m1,m2)
!
  external cubic_rs
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  BLEND_IJ_0D1 interpolates data in a table,'
  write ( *, '(a)' ) '  from corner data.'
  write ( *, '(a)' ) '  BLEND_IJ_1D1 interpolates data in a table,'
  write ( *, '(a)' ) '  from edge data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a,i6,a)' ) '  The table is ', m1, ' rows by ', m2, ' columns.'
!
!  Load data in the corners only.
!
  i = 1
  j = 1
  r = real ( i - 1 ) / real ( m1 - 1 )
  s = real ( j - 1 ) / real ( m2 - 1 )
  call cubic_rs ( r, s, 1, x(i,j) )

  i = m1
  j = 1
  r = real ( i - 1 ) / real ( m1 - 1 )
  s = real ( j - 1 ) / real ( m2 - 1 )
  call cubic_rs ( r, s, 1, x(i,j) )

  i = 1
  j = m2
  r = real ( i - 1 ) / real ( m1 - 1 )
  s = real ( j - 1 ) / real ( m2 - 1 )
  call cubic_rs ( r, s, 1, x(i,j) )

  i = m1
  j = m2
  r = real ( i - 1 ) / real ( m1 - 1 )
  s = real ( j - 1 ) / real ( m2 - 1 )
  call cubic_rs ( r, s, 1, x(i,j) )

  call blend_ij_0d1 ( x, m1, m2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by BLEND_IJ_0D1:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do
!
!  Load data in the edges.
!
  j = 1
  do i = 1, m1
    r = real ( i - 1 ) / real ( m1 - 1 )
    s = real ( j - 1 ) / real ( m2 - 1 )
    call cubic_rs ( r, s, 1, x(i,j) )
  end do

  j = m2
  do i = 1, m1
    r = real ( i - 1 ) / real ( m1 - 1 )
    s = real ( j - 1 ) / real ( m2 - 1 )
    call cubic_rs ( r, s, 1, x(i,j) )
  end do

  i = 1
  do j = 2, m2 - 1
    r = real ( i - 1 ) / real ( m1 - 1 )
    s = real ( j - 1 ) / real ( m2 - 1 )
    call cubic_rs ( r, s, 1, x(i,j) )
  end do

  i = m1
  do j = 2, m2 - 1
    r = real ( i - 1 ) / real ( m1 - 1 )
    s = real ( j - 1 ) / real ( m2 - 1 )
    call cubic_rs ( r, s, 1, x(i,j) )
  end do

  call blend_ij_1d1 ( x, m1, m2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by BLEND_IJ_1D1:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do
!
!  Compare with BLEND_RS_1D1
!
  do i = 1, m1
    r = real ( i - 1 ) / real ( m1 - 1 )
    do j = 1, m2
      s = real ( j - 1 ) / real ( m2 - 1 )
      call blend_rs_1dn ( r, s, x(i,j), n, cubic_rs )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data blended by BLEND_RS_1DN:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do
!
!  Load all data.
!
  do i = 1, m1
    do j = 1, m2
      r = real ( i - 1 ) / real ( m1 - 1 )
      s = real ( j - 1 ) / real ( m2 - 1 )
      call cubic_rs ( r, s, 1, x(i,j) )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact data:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do

  return
end
subroutine test05
!
!*******************************************************************************
!
!! TEST05 checks out BLEND_IJK_0D1
!
  implicit none
!
  integer, parameter :: m1 = 4
  integer, parameter :: m2 = 3
  integer, parameter :: m3 = 3
!
  integer i
  integer j
  integer k
  integer num_extreme
  real r
  real s
  real t
  real x(m1,m2,m3)
!
  external quad_rst
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  BLEND_IJK_0D1 interpolates data in a table,'
  write ( *, '(a)' ) '  from corner data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a,i6,a,i6,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns by ', m3, ' layers.'
!
!  Load data on the faces.
!
  do i = 1, m1
    r = real ( i - 1 ) / real ( m1 - 1 )
    do j = 1, m2
      s = real ( j - 1 ) / real ( m2 - 1 )
      do k = 1, m3
        t = real ( k - 1 ) / real ( m3 - 1 )

        num_extreme = 0
        if ( i == 1 .or. i == m1 ) then
          num_extreme = num_extreme + 1
        end if
        if ( j == 1 .or. j == m2 ) then
          num_extreme = num_extreme + 1
        end if
        if ( k == 1 .or. k == m3 ) then
          num_extreme = num_extreme + 1
        end if

        if ( num_extreme == 3 ) then
          call quad_rst ( r, s, t, 1, x(i,j,k) )
        else
          x(i,j,k) = 0.0E+00
        end if

      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data given to BLEND_IJK_0D1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  call blend_ijk_0d1 ( x, m1, m2, m3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by BLEND_IJK_0D1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do
!
!  Load all data.
!
  do i = 1, m1
    r = real ( i - 1 ) / real ( m1 - 1 )
    do j = 1, m2
      s = real ( j - 1 ) / real ( m2 - 1 )
      do k = 1, m3
        t = real ( k - 1 ) / real ( m3 - 1 )
        call quad_rst ( r, s, t, 1, x(i,j,k) )
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact data:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  return
end
subroutine test06
!
!*******************************************************************************
!
!! TEST06 checks out BLEND_IJK_1D1
!
  implicit none
!
  integer, parameter :: m1 = 4
  integer, parameter :: m2 = 3
  integer, parameter :: m3 = 3
!
  integer i
  integer j
  integer k
  integer num_extreme
  real r
  real s
  real t
  real x(m1,m2,m3)
!
  external quad_rst
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  BLEND_IJK_1D1 interpolates data in a table,'
  write ( *, '(a)' ) '  from edge data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a,i6,a,i6,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns by ', m3, ' layers.'
!
!  Load data on the faces.
!
  do i = 1, m1
    r = real ( i - 1 ) / real ( m1 - 1 )
    do j = 1, m2
      s = real ( j - 1 ) / real ( m2 - 1 )
      do k = 1, m3
        t = real ( k - 1 ) / real ( m3 - 1 )

        num_extreme = 0
        if ( i == 1 .or. i == m1 ) then
          num_extreme = num_extreme + 1
        end if
        if ( j == 1 .or. j == m2 ) then
          num_extreme = num_extreme + 1
        end if
        if ( k == 1 .or. k == m3 ) then
          num_extreme = num_extreme + 1
        end if

        if ( num_extreme >= 2 ) then
          call quad_rst ( r, s, t, 1, x(i,j,k) )
        else
          x(i,j,k) = 0.0E+00
        end if

      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data given to BLEND_IJK_1D1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  call blend_ijk_1d1 ( x, m1, m2, m3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by BLEND_IJK_1D1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do
!
!  Load all data.
!
  do i = 1, m1
    r = real ( i - 1 ) / real ( m1 - 1 )
    do j = 1, m2
      s = real ( j - 1 ) / real ( m2 - 1 )
      do k = 1, m3
        t = real ( k - 1 ) / real ( m3 - 1 )
        call quad_rst ( r, s, t, 1, x(i,j,k) )
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact data:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  return
end
subroutine test07
!
!*******************************************************************************
!
!! TEST07 checks out BLEND_IJK_2D1
!
  implicit none
!
  integer, parameter :: m1 = 4
  integer, parameter :: m2 = 3
  integer, parameter :: m3 = 3
!
  integer i
  integer j
  integer k
  integer num_extreme
  real r
  real s
  real t
  real x(m1,m2,m3)
!
  external quad_rst
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  BLEND_IJK_2D1 interpolates data in a table,'
  write ( *, '(a)' ) '  from face data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a,i6,a,i6,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns by ', &
    m3, ' layers.'
!
!  Load data on the faces.
!
  do i = 1, m1
    r = real ( i - 1 ) / real ( m1 - 1 )
    do j = 1, m2
      s = real ( j - 1 ) / real ( m2 - 1 )
      do k = 1, m3
        t = real ( k - 1 ) / real ( m3 - 1 )

        num_extreme = 0
        if ( i == 1 .or. i == m1 ) then
          num_extreme = num_extreme + 1
        end if
        if ( j == 1 .or. j == m2 ) then
          num_extreme = num_extreme + 1
        end if
        if ( k == 1 .or. k == m3 ) then
          num_extreme = num_extreme + 1
        end if

        if ( num_extreme >= 1 ) then
          call quad_rst ( r, s, t, 1, x(i,j,k) )
        else
          x(i,j,k) = 0.0E+00
        end if

      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data given to BLEND_IJK_2D1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  call blend_ijk_2d1 ( x, m1, m2, m3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by BLEND_IJK_2D1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do
!
!  Load all data.
!
  do i = 1, m1
    r = real ( i - 1 ) / real ( m1 - 1 )
    do j = 1, m2
      s = real ( j - 1 ) / real ( m2 - 1 )
      do k = 1, m3
        t = real ( k - 1 ) / real ( m3 - 1 )
        call quad_rst ( r, s, t, 1, x(i,j,k) )
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact data:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  return
end
subroutine cubic_rs ( r, s, i, x )
!
!*******************************************************************************
!
!! CUBIC_RS evaluates a function of R and S used for some tests.
!
  implicit none
!
  integer i
  real r
  real s
  real x
!
  x = 20.0E+00 * ( r**2 * s**3 )

  return
end
subroutine quad_rst ( r, s, t, i, x )
!
!*******************************************************************************
!
!! QUAD_RS evaluates a function of R and S used for some tests.
!
  implicit none
!
  integer i
  real r
  real s
  real t
  real x
!
  x = 18.0E+00 * ( r**2 + s + t )

  return
end
subroutine identity_r ( r, i, x )
!
!*******************************************************************************
!
!! IDENTITY_R returns a data component given (R).
!
!
!  Discussion:
!
!    This is a dummy routine, which simply returns (R).
!
!  Modified:
!
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the coordinate of a point that lies on the
!    boundary of the cube.
!
!    Input, integer I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R), which is on an endpoint of the unit line segment.
!
  implicit none
!
  integer i
  real r
  real x
!
  if ( i == 1 ) then
    x = r
  else 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IDENTITY_R - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal component index I = ', i
    stop
  end if

  return
end
subroutine identity_rs ( r, s, i, x )
!
!*******************************************************************************
!
!! IDENTITY_RS returns a data component given (R,S).
!
!
!  Discussion:
!
!    This is a dummy routine, which simply returns (R,S).
!
!  Modified:
!
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, the coordinates of a point that lies on the
!    boundary of the square.
!
!    Input, integer I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, or edge, of the unit square.
!
  implicit none
!
  integer i
  real r
  real s
  real x
!
  if ( i == 1 ) then
    x = r
  else if ( i == 2 ) then
    x = s
  else 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IDENTITY_RS - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal component index I = ', i
    stop
  end if

  return
end
subroutine identity_rst ( r, s, t, i, x )
!
!*******************************************************************************
!
!! IDENTITY_RST returns a data component given (R,S,T).
!
!
!  Discussion:
!
!    This is a dummy routine, which simply returns (R,S,T).
!
!  Modified:
!
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, T, the coordinates of a point that lies on the
!    boundary of the cube.
!
!    Input, integer I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, edge or face of the unit cube.
!
  implicit none
!
  integer i
!
  real r
  real s
  real t
  real x
!
  if ( i == 1 ) then
    x = r
  else if ( i == 2 ) then
    x = s
  else if ( i == 3 ) then
    x = t
  else 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IDENTITY_RST - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal component index I = ', i
    stop
  end if

  return
end
subroutine stretch_r ( r, i, x )
!
!*******************************************************************************
!
!! STRETCH_R returns a data component given (R).
!
!
!  Discussion:
!
!    This routine shifts by 1 and stretches by 2.
!
!  Modified:
!
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the coordinate of a point that lies on the
!    boundary of the cube.
!
!    Input, integer I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R), which is on an endpoint of the unit line segment.
!
  implicit none
!
  integer i
  real r
  real x
!
  if ( i == 1 ) then
    x = 2.0E+00 * r + 1.0E+00
  else 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRETCH_R - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal component index I = ', i
    stop
  end if

  return
end
subroutine stretch_rs ( r, s, i, x )
!
!*******************************************************************************
!
!! STRETCH_RS returns a data component given (R,S).
!
!
!  Discussion:
!
!    This routine shifts by (1,2) and stretches by (3,4).
!
!  Modified:
!
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, the coordinates of a point that lies on the
!    boundary of the square.
!
!    Input, integer I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, or edge, of the unit square.
!
  implicit none
!
  integer i
  real r
  real s
  real x
!
  if ( i == 1 ) then
    x = 3.0E+00 * r + 1.0E+00
  else if ( i == 2 ) then
    x = 4.0E+00 * s + 2.0E+00
  else 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRETCH_RS - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal component index I = ', i
    stop
  end if

  return
end
subroutine stretch_rst ( r, s, t, i, x )
!
!*******************************************************************************
!
!! STRETCH_RST returns a data component given (R,S,T).
!
!
!  Discussion:
!
!    This routine shifts by (1,2,3) and stretches by (4,5,6)
!
!  Modified:
!
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, T, the coordinates of a point that lies on the
!    boundary of the cube.
!
!    Input, integer I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, edge or face of the unit cube.
!
  implicit none
!
  integer i
!
  real r
  real s
  real t
  real x
!
  if ( i == 1 ) then
    x = 4.0E+00 * r + 1.0E+00
  else if ( i == 2 ) then
    x = 5.0E+00 * s + 2.0E+00
  else if ( i == 3 ) then
    x = 6.0E+00 * t + 3.0E+00
  else 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRETCH_RST - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal component index I = ', i
    stop
  end if

  return
end
subroutine ellipse_rs ( r, s, i, x )
!
!*******************************************************************************
!
!! ELLIPSE_RS maps the boundary of the unit square to an ellipse.
!
!
!  Discussion:
!
!    The ellipse is ( 3 * cos ( THETA ), 2 * sin ( THETA ) ).
!
!  Modified:
!
!    15 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, the coordinates of a point that lies on the
!    boundary of the square.
!
!    Input, integer I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, or edge, of the unit square.
!
  implicit none
!
  integer i
  real, parameter :: pi = 3.14159265E+00
  real r
  real s
  real theta
  real x
!
  if ( r == 0.0E+00 ) then
    theta = 0.25E+00 * pi * ( 5.0E+00 * ( 1.0E+00 - s ) + 3.0E+00 * s  )
  else if ( r == 1.0E+00 ) then
    theta = 0.25E+00 * pi * ( - 1.0E+00 * ( 1.0E+00 - s ) + 1.0E+00 * s )
  else if ( s == 0.0E+00 ) then
    theta = 0.25E+00 * pi * ( 5.0E+00 * ( 1.0E+00 - r ) + 7.0E+00 * r )
  else if ( s == 1.0E+00 ) then
    theta = 0.25E+00 * pi * ( 3.0E+00 * ( 1.0E+00 - r ) + 1.0E+00 * r )
  end if 

  if ( i == 1 ) then

    x = 3.0E+00 * cos ( theta )

  else if ( i == 2 ) then

    x = 2.0E+00 * sin ( theta )

  else 

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ELLIPSE_RS - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal component index I = ', i
    stop

  end if

  return
end
subroutine sphere_rst ( r, s, t, i, x )
!
!*******************************************************************************
!
!! SPHERE_RST maps the boundary of the unit cube to a sphere.
!
!
!  Discussion:
!
!    The sphere is
!      x = cos ( theta ) * cos ( phi )
!      y = sin ( theta ) * cos ( phi )
!      z = sin ( phi )
!
!  Modified:
!
!    15 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, the coordinates of a point that lies on the
!    boundary of the square.
!
!    Input, integer I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, or edge, of the unit square.
!
  implicit none
!
  integer i
  real norm
  real r
  real s
  real t
  real x
!
!  Compute length of vector from ( 0.5, 0.5, 0.5 ) to ( r, s, t )
!
  norm = sqrt ( ( r - 0.5E+00 )**2 + ( s - 0.5E+00 )**2 + ( t - 0.5E+00 )**2 )
!
!  Compute ( x, y, z ) coordinates of a point on the sphere 
!  ( x - 0.5 )**2 + ( y - 0.5 )**2 + ( z - 0.5 )**2 = 0.25 that is
!  the projection of the point ( r, s, t ) on the unit cube.
!
  if ( i == 1 ) then

    x = 0.5E+00 + 0.5E+00 * ( r - 0.5E+00 ) / norm

  else if ( i == 2 ) then

    x = 0.5E+00 + 0.5E+00 * ( s - 0.5E+00 ) / norm

  else if ( i == 3 ) then

    x = 0.5E+00 + 0.5E+00 * ( t - 0.5E+00 ) / norm

  else 

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_RST - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal component index I = ', i
    stop

  end if

  return
end
subroutine test08
!
!*******************************************************************************
!
!! TEST08 tests BLEND_IJ_W_1D1.
!
  implicit none
!
  integer, parameter :: n1 = 5
  integer, parameter :: n2 = 5
!
  integer i
  integer j
  real, parameter :: pi = 3.14159265E+00
  real r(n1)
  real rad
  real rr
  real s(n2)
  real ss
  real x(n1,n2)
  real y(n1,n2)
!
  rad = 3.0E+00
 
  x(1:n1,1:n2) = 0.0E+00
  y(1:n1,1:n2) = 0.0E+00
!
!  Set the boundary values.
!
!  It turns out that our values correspond to the X and Y
!  coordinates of a quarter circle of radius 3, although
!  it is by no means necessary that a formula for the data
!  be known.
!
  do i = 1, n1
    rr = ( real ( i - 1 ) / real ( n1 - 1 ) )**2
    r(i) = rr
    x(i,1) = 0.0E+00
    y(i,1) = 0.0E+00
    x(i,n2) = rad * cos ( 0.5E+00 * pi * ( 1.0E+00 - rr ) )
    y(i,n2) = rad * sin ( 0.5E+00 * pi * ( 1.0E+00 - rr ) )
  end do
 
  do j = 1, n2
    ss = ( real ( j - 1 ) / real ( n2 - 1 ) )**2
    s(j) = ss
    x(1,j) = 0.0E+00
    y(1,j) = rad * ss
    x(n1,j) = rad * ss
    y(n1,j) = 0.0E+00
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  BLEND_IJ_W_1D1 uses blending to fill in the'
  write ( *, '(a)' ) '    interior of a table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     R           S           X           Y'
  write ( *, '(a)' ) ' '
  call blend_ij_w_1d1 ( x, r, s, n1, n2 )
  call blend_ij_w_1d1 ( y, r, s, n1, n2 )
 
  do i = 1, n1
    write ( *, '(a)' ) ' '
    do j = 1, n2
      write ( *, '(4g12.4)' ) r(i), s(j), x(i,j), y(i,j)
    end do
  end do
 
  return
end
subroutine test09
!
!*******************************************************************************
!
!! TEST09 tests BLEND_102.
!
  implicit none
!
  integer, parameter :: n1 = 5
  integer, parameter :: n2 = 5
!
  real d(n1,n2)
  integer i
  integer j
  real r
  real s
!
  d(1:n1,1:n2) = 0.0E+00

  do i = 1, n1
    do j = 1, n2
      if ( &
        ( i == 1 .or. i == n1 ) .and. &
        ( j == 1 .or. j == n2 ) ) then
        d(i,j) = real ( i + j )
      end if
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  BLEND_102 blends corner values into a table.'
  write ( *, '(a)' ) ' '

  call rmat_print ( n1, n1, n2, d, '  Initial data array' )

  do i = 1, n1

    r = real ( i - 1 ) / real ( n1 - 1 )

    do j = 1, n2

      s = real ( j - 1 ) / real ( n2 - 1 )

      if ( &
        ( i == 1 .or. i == n1 ) .and. &
        ( j == 1 .or. j == n2 ) ) then
        cycle
      end if

      call blend_102 ( r, s, d(1,1), d(1,n2), d(n1,1), d(n1,n2), d(i,j) )

    end do

  end do
 
  call rmat_print ( n1, n1, n2, d, '  Interpolated data array' )

  return
end
subroutine test10
!
!*******************************************************************************
!
!! TEST10 tests BLEND_112.
!
  implicit none
!
  integer, parameter :: n1 = 5
  integer, parameter :: n2 = 5
!
  real d(n1,n2)
  integer i
  integer j
  real r
  real s
!
  d(1:n1,1:n2) = 0.0E+00

  do i = 1, n1
    do j = 1, n2
      if ( i == 1 .or. i == n1 .or. j == 1 .or. j == n2 ) then
        d(i,j) = real ( i + j )
      end if
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  BLEND_112 blends side values into a table.'
  write ( *, '(a)' ) ' '

  call rmat_print ( n1, n1, n2, d, '  Initial data array' )

  do i = 1, n1

    r = real ( i - 1 ) / real ( n1 - 1 )

    do j = 1, n2

      s = real ( j - 1 ) / real ( n2 - 1 )

      if ( i == 1 .or. i == n1 .or. j == 1 .or. j == n2 ) then
        cycle
      end if

      call blend_112 ( r, s, d(1,1), d(1,n2), d(n1,1), d(n1,n2), &
        d(i,1), d(i,n2), d(1,j), d(n1,j), d(i,j) )

    end do

  end do
 
  call rmat_print ( n1, n1, n2, d, '  Interpolated data array' )

  return
end
subroutine test11
!
!*******************************************************************************
!
!! TEST11 tests BLEND_103.
!
  implicit none
!
  integer, parameter :: n1 = 3
  integer, parameter :: n2 = 5
  integer, parameter :: n3 = 4
!
  real d(n1,n2,n3)
  integer i
  integer j
  integer k
  real r
  real s
  real t
!
  d(1:n1,1:n2,1:n3) = 0.0E+00

  do i = 1, n1
    do j = 1, n2
      do k = 1, n3
        if ( &
          ( i == 1 .or. i == n1 ) .and. &
          ( j == 1 .or. j == n2 ) .and. &
          ( k == 1 .or. k == n3  ) ) then
          d(i,j,k) = real ( i + j + k )
        end if
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  BLEND_103 blends corner values into a table.'
  write ( *, '(a)' ) ' '

  call rblock_print ( n1, n2, n3, d, '  Initial data array' )

  do i = 1, n1

    r = real ( i - 1 ) / real ( n1 - 1 )

    do j = 1, n2

      s = real ( j - 1 ) / real ( n2 - 1 )

      do k = 1, n3

        t = real ( k - 1 ) / real ( n3 - 1 )

        if ( &
          ( i == 1 .or. i == n1 ) .and. &
          ( j == 1 .or. j == n2 ) .and. &
          ( k == 1 .or. k == n3  ) ) then
          cycle
        end if

        call blend_103 ( r, s, t, &
          d(1,1,1), d(1,1,n3), d(1,n2,1), d(1,n2,n3), &
          d(n1,1,1), d(n1,1,n3), d(n1,n2,1), d(n1,n2,n3), d(i,j,k) )

      end do

    end do

  end do

  call rblock_print ( n1, n2, n3, d, '  Interpolated data array' )

  return
end
subroutine test12
!
!*******************************************************************************
!
!! TEST12 tests BLEND_113.
!
  implicit none
!
  integer, parameter :: n1 = 3
  integer, parameter :: n2 = 5
  integer, parameter :: n3 = 4
!
  real d(n1,n2,n3)
  integer i
  integer j
  integer k
  real r
  real s
  real t
!
  d(1:n1,1:n2,1:n3) = 0.0E+00

  do i = 1, n1
    do j = 1, n2
      do k = 1, n3
        if ( &
          ( ( i == 1 .or. i == n1 ) .and. ( j == 1 .or. j == n2 ) ) .or. &
          ( ( i == 1 .or. i == n1 ) .and. ( k == 1 .or. k == n3 ) ) .or. &
          ( ( j == 1 .or. j == n2 ) .and. ( k == 1 .or. k == n3 ) ) ) then
          d(i,j,k) = real ( i + j + k )
        end if
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  BLEND_113 blends edge values into a table.'
  write ( *, '(a)' ) ' '

  call rblock_print ( n1, n2, n3, d, '  Initial data array' )

  do i = 1, n1

    r = real ( i - 1 ) / real ( n1 - 1 )

    do j = 1, n2

      s = real ( j - 1 ) / real ( n2 - 1 )

      do k = 1, n3

        t = real ( k - 1 ) / real ( n3 - 1 )

        if ( &
          ( ( i == 1 .or. i == n1 ) .and. ( j == 1 .or. j == n2 ) ) .or. &
          ( ( i == 1 .or. i == n1 ) .and. ( k == 1 .or. k == n3 ) ) .or. &
          ( ( j == 1 .or. j == n2 ) .and. ( k == 1 .or. k == n3 ) ) ) then
          cycle
        end if

        call blend_113 ( r, s, t, &
          d(1,1,1), d(1,1,n3), d(1,n2,1), d(1,n2,n3), &
          d(n1,1,1), d(n1,1,n3), d(n1,n2,1), d(n1,n2,n3), &
          d(i,1,1), d(i,1,n3), d(i,n2,1), d(i,n2,n3), &
          d(1,j,1), d(1,j,n3), d(n1,j,1), d(n1,j,n3), &
          d(1,1,k), d(1,n2,k), d(n1,1,k), d(n1,n2,k), &
          d(i,j,k) )

      end do

    end do

  end do

  call rblock_print ( n1, n2, n3, d, '  Interpolated data array' )

  return
end
subroutine test13
!
!*******************************************************************************
!
!! TEST11 tests BLEND_123.
!
  implicit none
!
  integer, parameter :: n1 = 3
  integer, parameter :: n2 = 5
  integer, parameter :: n3 = 4
!
  real d(n1,n2,n3)
  integer i
  integer j
  integer k
  real r
  real s
  real t
!
  d(1:n1,1:n2,1:n3) = 0.0E+00

  do i = 1, n1
    do j = 1, n2
      do k = 1, n3
        if ( &
          ( i == 1 .or. i == n1 ) .or. &
          ( j == 1 .or. j == n2 ) .or. &
          ( k == 1 .or. k == n3  ) ) then
          d(i,j,k) = real ( i + j + k )
        end if
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  BLEND_123 blends face values into a table.'
  write ( *, '(a)' ) ' '

  call rblock_print ( n1, n2, n3, d, '  Initial data array' )

  do i = 1, n1

    r = real ( i - 1 ) / real ( n1 - 1 )

    do j = 1, n2

      s = real ( j - 1 ) / real ( n2 - 1 )

      do k = 1, n3

        t = real ( k - 1 ) / real ( n3 - 1 )

        if ( &
          ( i == 1 .or. i == n1 ) .or. &
          ( j == 1 .or. j == n2 ) .or. &
          ( k == 1 .or. k == n3  ) ) then
          cycle
        end if

        call blend_123 ( r, s, t, &
          d(1,1,1), d(1,1,n3), d(1,n2,1), d(1,n2,n3), &
          d(n1,1,1), d(n1,1,n3), d(n1,n2,1), d(n1,n2,n3), &
          d(i,1,1), d(i,1,n3), d(i,n2,1), d(i,n2,n3), &
          d(1,j,1), d(1,j,n3), d(n1,j,1), d(n1,j,n3), &
          d(1,1,k), d(1,n2,k), d(n1,1,k), d(n1,n2,k), &
          d(1,j,k), d(n1,j,k), d(i,1,k), d(i,n2,k), d(i,j,1), d(i,j,n3), &
          d(i,j,k) )

      end do

    end do

  end do

  call rblock_print ( n1, n2, n3, d, '  Interpolated data array' )

  return
end
