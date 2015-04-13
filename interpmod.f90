module interpmod
use nrtype
use nrutil, only: assert
implicit none

contains

  subroutine blend_101 ( r, x0, x1, x )
    !
    !*******************************************************************************
    !
    !! BLEND_101 extends scalar endpoint data to a line.
    !
    !
    !  Diagram:
    !
    !    0-----r-----1
    !
    !  Reference:
    !
    !    William Gordon,
    !    Blending-Function Methods of Bivariate and Multivariate Interpolation
    !      and Approximation,
    !    SIAM Journal on Numerical Analysis,
    !    Volume 8, Number 1, March 1971, pages 158-177.
    !
    !    William Gordon and Charles Hall,
    !    Transfinite Element Methods: Blending-Function Interpolation over
    !      Arbitrary Curved Element Domains,
    !    Numerische Mathematik,
    !    Volume 21, Number 1, 1973, pages 109-129.
    !
    !    William Gordon and Charles Hall,
    !    Construction of Curvilinear Coordinate Systems and Application to
    !      Mesh Generation,
    !    International Journal of Numerical Methods in Engineering,
    !    Volume 7, 1973, pages 461-477.
    !
    !    Joe Thompson, Bharat Soni, Nigel Weatherill,
    !    Handbook of Grid Generation,
    !    CRC Press, 1999.
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
    !    Input, real R, the coordinate where an interpolated value is desired.  
    !
    !    Input, real X0, X1, the data values at the ends of the line.
    !
    !    Output, real X, the interpolated data value at (R).
    !
    implicit none
    !
    real(SP), INTENT(IN):: r,x0,x1
    real(SP), INTENT(OUT):: x
    
    call assert(floor(r*.99999)==0,'out of range in blend_101')

    
    !
    x = ( 1.0E+00 - r ) * x0 + r * x1

  end subroutine blend_101
  subroutine blend_102 ( r, s, x00, x10, x01, x11, x )
    !
    !*******************************************************************************
    !
    !! BLEND_102 extends scalar point data into a square.
    !
    !
    !  Diagram:
    !
    !    01------------11
    !     |      .      |
    !     |      .      |
    !     |.....rs......|
    !     |      .      |
    !     |      .      |
    !    00------------10
    !
    !  Formula:
    !
    !    Written in terms of R and S, the map has the form:
    !
    !      X(R,S) =
    !               1     * ( + x00 )
    !             + r     * ( - x00 + x10 )
    !             + s     * ( - x00       + x01 )
    !             + r * s * ( + x00 - x10 - x01 + x11 )
    !
    !    Written in terms of the coefficients, the map has the form:
    !
    !      X(R,S) = x00 * ( 1 - r - s + r * s )
    !             + x01 * (         s - r * s )
    !             + x10 * (     r     - r * s )
    !             + x11 * (             r * s )
    !
    !             = x00 * ( 1 - r ) * ( 1 - s )
    !             + x01 * ( 1 - r ) *       s
    !             + x10 *       r   * ( 1 - s )
    !             + x11 *       r           s
    !
    !    The nonlinear term ( r * s ) has an important role:
    !
    !      If ( x01 + x10 - x00 - x11 ) is zero, then the input data lies in
    !      a plane, and the mapping is affine.  All the interpolated data 
    !      will lie on the plane defined by the four corner values.  In 
    !      particular, on any line through the square, data values at 
    !      intermediate points will lie between the values at the endpoints.  
    !
    !      If ( x01 + x10 - x00 - x11 ) is not zero, then the input data does
    !      not lie in a plane, and the interpolation map is nonlinear.  On
    !      any line through the square, data values at intermediate points
    !      may lie above or below the data values at the endpoints.  The
    !      size of the coefficient of r * s will determine how severe this
    !      effect is.
    !
    !  Reference:
    !
    !    William Gordon,
    !    Blending-Function Methods of Bivariate and Multivariate Interpolation
    !      and Approximation,
    !    SIAM Journal on Numerical Analysis,
    !    Volume 8, Number 1, March 1971, pages 158-177.
    !
    !    William Gordon and Charles Hall,
    !    Transfinite Element Methods: Blending-Function Interpolation over
    !      Arbitrary Curved Element Domains,
    !    Numerische Mathematik,
    !    Volume 21, Number 1, 1973, pages 109-129.
    !
    !    William Gordon and Charles Hall,
    !    Construction of Curvilinear Coordinate Systems and Application to
    !      Mesh Generation,
    !    International Journal of Numerical Methods in Engineering,
    !    Volume 7, 1973, pages 461-477.
    !
    !    Joe Thompson, Bharat Soni, Nigel Weatherill,
    !    Handbook of Grid Generation,
    !    CRC Press, 1999.
    !
    !  Modified:
    !
    !    11 October 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real R, S, the coordinates where an interpolated value is 
    !    desired.  
    !
    !    Input, real X00, X01, X10, X11, the data values at the corners.
    !
    !    Output, real X, the interpolated data value at (R,S).
    !
    implicit none
    !
    real(sp):: r
    real(sp):: s
    real(sp):: x
    real(sp):: x00
    real(sp):: x01
    real(sp):: x10
    real(sp):: x11

    call assert(floor(r*.99999)==0,'r out of range in blend_102')
    !if (floor(r*0.999)>0..or.floor(r*.0999)<0) then
    !   print*, r
    !   stop
    !end if
    call assert(floor(s*.99999)==0,'s out of range in blend_102')
    !
    x =             + x00 &
         + r *     ( - x00 + x10 ) & 
         + s *     ( - x00       + x01 ) &
         + r * s * ( + x00 - x10 - x01 + x11 )
    
    
  end subroutine blend_102

  subroutine blend_101_bf ( r, x0, x1, x )
    use nrtype
    implicit none
    real(SP), INTENT(IN):: r,x0,x1
    REAL(SP), INTENT(OUT) :: x

    call assert(floor(r)==0,'out of range in blend_101_bf')


    x=(1.0_sp-nint(r))*x0+nint(r)*x1
    
    
  end subroutine blend_101_bf
  
  subroutine blend_102_bf ( r, s, x00, x10, x01, x11, x )
    use nrtype
    implicit none
    real(sp),intent(IN):: r
    real(sp),intent(IN):: s
    real(sp),intent(out):: x
    real(sp):: x00,x01,x10,x11

    !if (floor(s)>0.or.floor(s)<0) then
    !   print*,s
    !end if

    call assert(floor(r*.99999)==0,'out of range in blend_102_bf')
    call assert(floor(s*.99999)==0,'out of range in blend_102_bf')
    
    x =             + x00 &
         + nint(r) *     ( - x00 + x10 ) & 
         + nint(s) *     ( - x00       + x01 ) &
         + nint(r) * nint(s) * ( + x00 - x10 - x01 + x11 )
    
  end subroutine blend_102_bf
  
end module interpmod
