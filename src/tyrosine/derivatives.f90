module derivatives
  !
  ! A module containing finite-difference approximations for derivatives
  !
  !    - First derivatives  -->  [centred; 2nd and 4th order, backward; 1st, 2nd, & 3rd order]
  !    - Second derivatives -->  [partial; 2nd and 4th order, mixed; 2nd and 4th order]
  !
  use options, only:c_double
  use periodic, only:apply_periodic,apply_periodic_fourth
  implicit none

contains

  ! ------------------------------------------
  ! ------------------------------------------
  !
  !       First derivatives: df/dx
  !
  ! ------------------------------------------
  ! ------------------------------------------

  ! =========== Centred difference ============

  !
  ! function to return the 2nd order approx of first derivative
  !
  real(c_double) function deriv1(fp1,fm1,h)
    real(c_double) :: fp1, fm1, h    ! values of the function at i+1,i-1 (or j, k) and the grid spacing

    deriv1 = (fp1 - fm1) / (2.d0 * h)

  end function deriv1

  !
  ! function to return 4th order approx of first derivative
  !
  real(c_double) function deriv1fourth(fp1,fp2,fm1,fm2,h)
    real(c_double) :: fp1,fm1,fp2,fm2,h

    deriv1fourth = (-fp2 + 8.d0 * fp1 - 8.d0 * fm1 + fm2) / (12.d0 * h)

  end function deriv1fourth


  ! =========== Backward difference ============

  !
  ! a function to return the first order approx of first deriv using backward difference
  !       --> can be forward diff with h --> -h
  !
  real(c_double) function deriv1_bckwrd(f,fm1,h)
    real(c_double) :: f,fm1,h

    deriv1_bckwrd = (f - fm1) / h

  end function deriv1_bckwrd

  !
  ! function to return the SECOND order approx of first deriv BACKWARD difference (good for time derivs)
  !       --> can be forward diff with h --> -h
  !
  real(c_double) function deriv1_bckwrd2nd(f,fm1,fm2,h)
    real(c_double) :: f,fm1,fm2,h

    deriv1_bckwrd2nd = (3.d0 * f - 4.d0 * fm1 + fm2) / (2.d0 * h)

  end function deriv1_bckwrd2nd

  !
  ! function to return the THIRD order approx of first deriv BACKWARD difference
  !       --> can be forward diff with h --> -h
  !
  real(c_double) function deriv1_bckwrd3rd(f,fm1,fm2,fm3,h)
    real(c_double) :: f,fm1,fm2,fm3,h
    real(c_double) :: c1,c2,c3 ! to simplify

    c1 = 11.d0 / 6.d0
    c2 = 3.d0 / 2.d0
    c3 = 1.d0 / 3.d0

    deriv1_bckwrd3rd = (c1 * f - 3.d0 * fm1 + c2 * fm2 - c3 * fm3) /  h

  end function deriv1_bckwrd3rd

  !
  ! function to return the FOURTH order approx of first deriv BACKWARD difference
  !       --> can be forward diff with h --> -h
  !
  real(c_double) function deriv1_bckwrd4th(f,fm1,fm2,fm3,fm4,h)
    real(c_double) :: f,fm1,fm2,fm3,fm4,h
    real(c_double) :: c,c1,c2,c3,c4

    c  = 25.d0 / 12.d0
    c1 = - 4.d0
    c2 = 3.d0
    c3 = - 4.d0 / 3.d0
    c4 = 1.d0 / 4.d0

    deriv1_bckwrd4th = (c * f + c1 * fm1 + c2 * fm2 + c3 * fm3 + c4 * fm4) / h

  end function deriv1_bckwrd4th


  ! =========== Forward difference ============

  !
  ! function to return the SECOND order approx of first deriv FORWARD difference (good for time derivs)
  !       --> can be backward diff with h --> -h
  !
  real(c_double) function deriv1_forward2nd(f,fp1,fp2,h)
    real(c_double) :: f,fp1,fp2,h

    deriv1_forward2nd = (-3.d0 * f + 4.d0 * fp1 - fp2) / (2.d0 * h)

  end function deriv1_forward2nd


  !
  ! function to return the THIRD order approx of first deriv FORWARD difference
  !       --> can be backward diff with h --> -h
  !
  real(c_double) function deriv1_forward3rd(f,fp1,fp2,fp3,h)
    real(c_double) :: f,fp1,fp2,fp3,h
    real(c_double) :: c1,c2,c3 ! to simplify

    c1 = - 11.d0 / 6.d0
    c2 = - 3.d0 / 2.d0
    c3 = 1.d0 / 3.d0

    deriv1_forward3rd = (c1 * f + 3.d0 * fp1 + c2 * fp2 + c3 * fp3) /  h

  end function deriv1_forward3rd

  !
  ! function to return the FOURTH order approx of first deriv FORWARD difference
  !       --> can be backward diff with h --> -h
  !
  real(c_double) function deriv1_forward4th(f,fp1,fp2,fp3,fp4,h)
    real(c_double) :: f,fp1,fp2,fp3,fp4,h
    real(c_double) :: c,c1,c2,c3,c4

    c  = - 25.d0 / 12.d0
    c1 = 4.d0
    c2 = - 3.d0
    c3 = 4.d0 / 3.d0
    c4 = - 1.d0 / 4.d0

    deriv1_forward4th = (c * f + c1 * fp1 + c2 * fp2 + c3 * fp3 + c4 * fp4) / h

  end function deriv1_forward4th


  ! ------------------------------------------
  ! ------------------------------------------
  !
  !       Second derivatives: d^2f/dx^2
  !
  ! ------------------------------------------
  ! ------------------------------------------



  ! =========== Centred difference ============

  !
  ! function to return the 2nd order approx of partial 2nd deriv one variable e.g. d2/dx2(f)
  !
  real(c_double) function deriv2(fp1,f,fm1,h)
    real(c_double) :: fp1, f, fm1  !! values of the function at i+1,i,i-1 (or j, k)
    real(c_double) :: h !! the spacing in whatever dimension we're in

    deriv2 = (fp1 - 2.d0 * f + fm1) / h**2

  end function deriv2

  !
  ! function to return the 4th order approx of partial 2nd deriv of one variable
  !
  real(c_double) function deriv2fourth(fp1,fp2,f,fm1,fm2,h)
    real(c_double) :: fp1,fp2,f,fm1,fm2,h

    deriv2fourth = (-fp2 + 16.d0 * fp1 - 30.d0 * f + 16.d0 * fm1 - fm2) / (12.d0 * h**2 )

  end function deriv2fourth


  ! =========== Backward difference ============

  !
  ! function to return the 2nd order accurate 2ND DERIVATIVE using a
  !      backward stencil -- deriv w.r.t ONE variable
  !
  !     --> EQUIV to forward diff
  !
  real(c_double) function deriv2_bckwrd2nd(f,fm1,fm2,fm3,h)
    real(c_double), intent(in) :: f,fm1,fm2,fm3,h

    deriv2_bckwrd2nd = (2.d0*f - 5.d0*fm1 + 4.d0*fm2 - fm3) / h**2

  end function deriv2_bckwrd2nd

  !
  ! function to return the 4th order accurate 2ND DERIVATIVE using a
  !     backward stencil -- deriv w.r.t ONE variable
  !
  !    --> EQUIV to forward diff (fm1-->fp1, etc)
  !
  real(c_double) function deriv2_bckwrd4th(f,fm1,fm2,fm3,fm4,fm5,h)
    real(c_double), intent(in) :: f,fm1,fm2,fm3,fm4,fm5,h
    real(c_double) :: n,nm1,nm2,nm3,nm4,nm5,denom
    !
    ! coefficients
    n   = 45.d0
    nm1 = - 154.d0
    nm2 = 214.d0
    nm3 = - 156.d0
    nm4 = 61.d0
    nm5 = - 10.d0
    denom = 12.d0 * h**2

    deriv2_bckwrd4th = (n*f + nm1*fm1 + nm2*fm2 + nm3*fm3 + nm4*fm4 + nm5*fm5) / denom

  end function deriv2_bckwrd4th


  ! ------------------------------------------
  ! ------------------------------------------
  !
  !     Mixed second derivatives: d^2f/dxdy
  !
  ! ------------------------------------------
  ! ------------------------------------------


  !
  ! a function to return the secon *mixed* derivative, e.g. d2/dxdy(f) (general to any dim mix)
  ! * TESTED for linear solution phi = sin(x*y*z) *
  !
  real(c_double) function deriv2_mix_wholestencil(f,f_xp1yp1,f_xm1ym1,f_xp1,f_xm1,f_yp1,f_ym1,dx,dy)
    real(c_double), intent(in) :: dx, dy ! grid spacing in the two dimensions
    real(c_double), intent(in) :: f, f_xp1yp1, f_xm1ym1 ! f(x,y), f(x+dx,y+dy), f(x-dx,y-dy)
    real(c_double), intent(in) :: f_xp1, f_xm1, f_yp1, f_ym1 ! f(x+dx,y), f(x-dx,y), f(x,y+dy), f(x,y-dy)

    real(c_double) :: num, denom, d2fdx2, d2fdy2

    denom = 2.d0 * dx * dy ! denominator

    d2fdx2 = deriv2(f_xp1,f,f_xm1,dx) ! second deriv w.r.t x
    d2fdy2 = deriv2(f_yp1,f,f_ym1,dy) ! second deriv w.r.t y

    num = f_xp1yp1 + f_xm1ym1 - 2.d0 * f - dx*dx*d2fdx2 - dy*dy*d2fdy2

    deriv2_mix_wholestencil = num / denom

  end function deriv2_mix_wholestencil


  !
  ! a function for the mixed second derivative
  !      different to above, found it here: https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781119083405.app1
  !         (but also in other places too)
  !
  real(c_double) function deriv2_mix_2ndorder(f_xp1yp1,f_xp1ym1,f_xm1yp1,f_xm1ym1,dx,dy)
      real(c_double), intent(in) :: f_xp1yp1,f_xp1ym1,f_xm1yp1,f_xm1ym1,dx,dy
      real(c_double) :: num

      num = f_xp1yp1 - f_xm1yp1 - f_xp1ym1 + f_xm1ym1

      deriv2_mix_2ndorder = num / (4.d0 * dx * dy)


  end function deriv2_mix_2ndorder

  !
  ! a function to return the mixed 2nd derivative (same as above)
  !   -- takes function we want deriv. of in the full 2D plane of derivative
  !         (can be simpler to call in a loop)
  !
  real(c_double) function deriv2_mix(i,j,nx,func,dx)
    integer, intent(in) :: i,j,nx    ! = ny
    real(c_double), intent(in) :: dx ! = dy
    real(c_double), intent(in) :: func(nx,nx) ! function in the 2D plane derivative is taken in

    real(c_double) :: f, f_xp1yp1, f_xm1ym1      ! f(x,y), f(x+dx,y+dy), f(x-dx,y-dy)
    real(c_double) :: f_xp1, f_xm1, f_yp1, f_ym1 ! f(x+dx,y), f(x-dx,y), f(x,y+dy), f(x,y-dy)
    real(c_double) :: num, denom, d2fdx2, d2fdy2
    integer :: ip1,jp1,im1,jm1   ! for periodic boundaries

    ! get periodic stencil
    call apply_periodic(i,ip1,im1,nx)
    call apply_periodic(j,jp1,jm1,nx)

    f        = func(i,j)
    f_xp1yp1 = func(ip1,jp1)
    f_xm1ym1 = func(im1,jm1)
    f_xp1    = func(ip1,j)
    f_xm1    = func(im1,j)
    f_yp1    = func(i,jp1)
    f_ym1    = func(i,jm1)

    denom = 2.d0 * dx * dx ! denominator

    d2fdx2 = deriv2(f_xp1,f,f_xm1,dx) ! second deriv w.r.t x
    d2fdy2 = deriv2(f_yp1,f,f_ym1,dx) ! second deriv w.r.t y

    num = f_xp1yp1 + f_xm1ym1 - 2.d0 * f - dx*dx*d2fdx2 - dx*dx*d2fdy2

    deriv2_mix = num / denom

  end function deriv2_mix




  !
  ! a function to return the second *mixed* derivative to FOURTH order d2(f)/dxdy for any x,y dims
  ! this was calculated (by me) by taking the fourth order approx. for 1st deriv, and taking the fourth order
  ! approx of the y-deriv of that x-deriv approx...
  ! * TESTED for linear solution phi = sin(x*y*z) *
  !
  real(c_double) function deriv2_mixfourth_wholestencil(fxp2yp2,fxp2yp1,fxp2ym1,fxp2ym2,fxp1yp2,fxp1yp1,&
       fxp1ym1,fxp1ym2,fxm1yp2,fxm1yp1,fxm1ym1,fxm1ym2,fxm2yp2,fxm2yp1,fxm2ym1,fxm2ym2,dx,dy)
    real(c_double), intent(in) :: dx,dy
    real(c_double), intent(in) :: fxp2yp2,fxp2yp1,fxp2ym1,fxp2ym2
    real(c_double), intent(in) :: fxp1yp2,fxp1yp1,fxp1ym1,fxp1ym2
    real(c_double), intent(in) :: fxm1yp2,fxm1yp1,fxm1ym1,fxm1ym2
    real(c_double), intent(in) :: fxm2yp2,fxm2yp1,fxm2ym1,fxm2ym2
    real(c_double) :: t1,t2,t3,t4,denom

    denom = 12.d0 * dx * 12.d0 * dy ! denominator

    t1 = fxp2yp2 - 8.d0 * fxp2yp1 + 8.d0 * fxp2ym1 - fxp2ym2
    t2 = 8.d0 * (- fxp1yp2 + 8.d0 * fxp1yp1 - 8.d0 * fxp1ym1 + fxp1ym2)
    t3 = -8.d0 * (- fxm1yp2 + 8.d0 * fxm1yp1 - 8.d0 * fxm1ym1 + fxm1ym2)
    t4 = - fxm2yp2 + 8.d0 * fxm2yp1 - 8.d0 * fxm2ym1 + fxm2ym2

    deriv2_mixfourth_wholestencil = (t1 + t2 + t3 + t4) / denom

  end function deriv2_mixfourth_wholestencil



  !
  ! a function to return the second *mixed* derivative to FOURTH order d2(f)/dxdy for any x,y dims
  !
  ! --> copied from above but takes a (4,4) array with stencil instead of each component seperately
  !    (written for dL calc curlyh derivs)
  !
  real(c_double) function deriv2_mixfourth_stencilarray(f_xystencil,dx,dy)
    real(c_double), intent(in) :: dx,dy,f_xystencil(4,4)

    real(c_double) :: fxp2yp2,fxp2yp1,fxp2ym1,fxp2ym2
    real(c_double) :: fxp1yp2,fxp1yp1,fxp1ym1,fxp1ym2
    real(c_double) :: fxm1yp2,fxm1yp1,fxm1ym1,fxm1ym2
    real(c_double) :: fxm2yp2,fxm2yp1,fxm2ym1,fxm2ym2
    real(c_double) :: t1,t2,t3,t4,denom

    fxp2yp2 = f_xystencil(1,1)
    fxp2yp1 = f_xystencil(1,2)
    fxp2ym1 = f_xystencil(1,3)
    fxp2ym2 = f_xystencil(1,4)
    fxp1yp2 = f_xystencil(2,1)
    fxp1yp1 = f_xystencil(2,2)
    fxp1ym1 = f_xystencil(2,3)
    fxp1ym2 = f_xystencil(2,4)
    fxm1yp2 = f_xystencil(3,1)
    fxm1yp1 = f_xystencil(3,2)
    fxm1ym1 = f_xystencil(3,3)
    fxm1ym2 = f_xystencil(3,4)
    fxm2yp2 = f_xystencil(4,1)
    fxm2yp1 = f_xystencil(4,2)
    fxm2ym1 = f_xystencil(4,3)
    fxm2ym2 = f_xystencil(4,4)

    denom = 12.d0 * dx * 12.d0 * dy ! denominator

    t1 = fxp2yp2 - 8.d0 * fxp2yp1 + 8.d0 * fxp2ym1 - fxp2ym2
    t2 = 8.d0 * (- fxp1yp2 + 8.d0 * fxp1yp1 - 8.d0 * fxp1ym1 + fxp1ym2)
    t3 = -8.d0 * (- fxm1yp2 + 8.d0 * fxm1yp1 - 8.d0 * fxm1ym1 + fxm1ym2)
    t4 = - fxm2yp2 + 8.d0 * fxm2yp1 - 8.d0 * fxm2ym1 + fxm2ym2

    deriv2_mixfourth_stencilarray = (t1 + t2 + t3 + t4) / denom

  end function deriv2_mixfourth_stencilarray




  !
  ! a function to return the secon MIXED derivative (copied from above)
  !   but removes the ugly call to the function by putting the indexing of the function in here
  !
  real(c_double) function deriv2_mixfourth(xi,yi,nx,func,dx)
    integer, intent(in) :: xi,yi,nx ! = ny
    real(c_double), intent(in) :: dx! = dy
    real(c_double), intent(in) :: func(nx,nx) ! function in the 2D plane derivative is taken in

    real(c_double) :: fxp2yp2,fxp2yp1,fxp2ym1,fxp2ym2
    real(c_double) :: fxp1yp2,fxp1yp1,fxp1ym1,fxp1ym2
    real(c_double) :: fxm1yp2,fxm1yp1,fxm1ym1,fxm1ym2
    real(c_double) :: fxm2yp2,fxm2yp1,fxm2ym1,fxm2ym2
    real(c_double) :: t1,t2,t3,t4,denom
    integer :: ip1,jp1,im1,jm1,ip2,im2,jp2,jm2 ! for periodic boundaries

    ! get periodic stencil
    call apply_periodic(xi,ip1,im1,nx)
    call apply_periodic(yi,jp1,jm1,nx)
    call apply_periodic_fourth(xi,ip2,im2,nx)
    call apply_periodic_fourth(yi,jp2,jm2,nx)

    ! define stencil from function
    fxp2yp2 = func(ip2,jp2)
    fxp2yp1 = func(ip2,jp1)
    fxp2ym1 = func(ip2,jm1)
    fxp2ym2 = func(ip2,jm2)
    fxp1yp2 = func(ip1,jp2)
    fxp1yp1 = func(ip1,jp1)
    fxp1ym1 = func(ip1,jm1)
    fxp1ym2 = func(ip1,jm2)
    fxm1yp2 = func(im1,jp2)
    fxm1yp1 = func(im1,jp1)
    fxm1ym1 = func(im1,jm1)
    fxm1ym2 = func(im1,jm2)
    fxm2yp2 = func(im2,jp2)
    fxm2yp1 = func(im2,jp1)
    fxm2ym1 = func(im2,jm1)
    fxm2ym2 = func(im2,jm2)

    denom = 12.d0 * dx * 12.d0 * dx ! denominator

    t1 = fxp2yp2 - 8.d0 * fxp2yp1 + 8.d0 * fxp2ym1 - fxp2ym2
    t2 = 8.d0 * (- fxp1yp2 + 8.d0 * fxp1yp1 - 8.d0 * fxp1ym1 + fxp1ym2)
    t3 = -8.d0 * (- fxm1yp2 + 8.d0 * fxm1yp1 - 8.d0 * fxm1ym1 + fxm1ym2)
    t4 = - fxm2yp2 + 8.d0 * fxm2yp1 - 8.d0 * fxm2ym1 + fxm2ym2

    deriv2_mixfourth = (t1 + t2 + t3 + t4) / denom

  end function deriv2_mixfourth






end module derivatives
