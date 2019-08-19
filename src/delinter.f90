  !!<summary>Contains scripts used to correct or 'delint' a POSCAR.</summary>

Module delinter

  use num_types
  use numerical_utilities
  use vector_matrix_utilities, only: matrix_inverse, norm
  
  implicit none
  private
  public delint

CONTAINS

  !!<summary>Converts cell parameters into a matrix of basis
  !!vectors.</summary>
  !!<parameter name="params" regular="true">The cell parameters, a, b,
  !!c, alpha, beta, gamma.</parameter>
  !!<parameter name="basis" regular="true">The basis vectors of the
  !!cell.</parameter>
  SUBROUTINE cell2basis(params, basis)
    real(dp), intent(in) :: params(6)
    real(dp), intent(out) :: basis(3,3)

    real(dp) :: a, b, c, alpha, beta, gamma, temp

    a = params(1)
    b = params(2)
    c = params(3)
    alpha = params(4)
    beta = params(5)
    gamma = params(6)

    temp = (cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)

    basis(1,1) = a
    basis(1,2) = b*cos(gamma)
    basis(1,3) = c*cos(beta)
    basis(2,1) = 0.0_dp
    basis(2,2) = b*sin(gamma)
    basis(2,3) = c*temp
    basis(3,1) = 0.0_dp
    basis(3,2) = 0.0_dp
    basis(3,3) = c*sqrt(1-cos(beta)**2 -temp**2)
  end SUBROUTINE cell2basis

  !!<summary>Returns the cosine of the angle between two
  !!vectors.</summary>
  !!<parameter name="vec1" regular="true">Vector 1.</parameter>
  !!<parameter name="vec2" regular="true">Vecotr 2.</parameter>
  !!<parameter name="ang" regular="true">The cosine of the angle
  !!between vecs 1 and 2.</parameter>
  SUBROUTINE cosangs(vec1, vec2, ang)
    real(dp), intent(in) :: vec1(3), vec2(3)
    real(dp), intent(out) :: ang

    ang = dot_product(vec1, vec2)/(sqrt(dot_product(vec1,vec1))*sqrt(dot_product(vec2,vec2)))
  end SUBROUTINE cosangs
  
  !!<summary>Converts crystall basis to cell parameters.</summary>
  !!<parameter name="vecs" regular="true">The input basis
  !!vectors.</parameter>
  !!<parameter name="params" regular="true">The output cell
  !!parameters.</parameter>
  SUBROUTINE basis2cell(vecs, params)
    real(dp), intent(in) :: vecs(3,3)
    real(dp), intent(out) :: params(6)

    integer :: i
    real(dp) :: alpha, beta, gamma

    do i=1,3
       params(i) = sqrt(dot_product(vecs(:,i),vecs(:,i)))
    end do

    call cosangs(vecs(:,2), vecs(:,3), alpha)
    call cosangs(vecs(:,1), vecs(:,3), beta)
    call cosangs(vecs(:,1), vecs(:,2), gamma)
    params(4) = acos(alpha)
    params(5) = acos(beta)
    params(6) = acos(gamma)    
  end SUBROUTINE basis2cell

  !!<summary>Enforces the symmetries of the niggli case onto the
  !!cell parameters.</summary>
  !!<parameter name="in_params" regular="true">The input cell
  !!parameters.</parameter>
  !!<parameter name="case" regular="true">The niggli case
  !!number.</parameter>
  !!<parameter name="out_params" regular="true">The output symmetrized
  !!parameters.</parameter>
  !!<parameter name="eps" regular="true">Floating point tolerance.</parameter>
  SUBROUTINE symmetrize(in_params, case, out_params, eps)
    real(dp), intent(in) :: in_params(6)
    integer, intent(in) :: case
    real(dp), intent(out) :: out_params(6)
    real(dp), intent(in) :: eps

    real(dp) :: a0, b0, c0, alpha0, beta0, gamma0
    real(dp) :: a, b, c, alpha, beta, gamma
    real(dp) :: temp, pi

    pi = 3.1415926535897931_dp

    a0 = in_params(1)
    b0 = in_params(2)
    c0 = in_params(3)
    alpha0 = in_params(4)
    beta0 = in_params(5)
    gamma0 = in_params(6)
    
    if (case==1) then
       temp = (a0+b0+c0)/3.0_dp
       a = temp
       b = temp
       c = temp
       alpha = pi/3.0_dp
       beta= pi/3.0_dp
       gamma = pi/3.0_dp
    else if ((case==2) .or. (case==4)) then
       temp = (a0+b0+c0)/3.0_dp
       a = temp
       b = temp
       c = temp
       temp = (alpha0+beta0+gamma0)/3.0_dp
       alpha = temp
       beta = temp
       gamma = temp
    else if (case==3) then
       temp = (a0+b0+c0)/3.0_dp
       a = temp
       b = temp
       c = temp
       alpha = pi/2.0_dp
       beta = pi/2.0_dp
       gamma= pi/2.0_dp
    else if (case==5) then
       temp = (a0+b0+c0)/3.0_dp
       a = temp
       b = temp
       c = temp
       alpha = 1.910633236249019_dp
       beta = 1.910633236249019_dp
       gamma = 1.910633236249019_dp
    else if (case==6) then
       temp = (a0+b0+c0)/3.0_dp
       a = temp
       b = temp
       c = temp
       temp = acos(((((a**2+b**2)/2.0_dp)-cos(gamma0)*a*c)/2.0_dp)/(b+c))
       alpha = temp
       beta = temp
       gamma = gamma0
    else if (case==7) then
       temp = (a0+b0+c0)/3.0_dp
       a = temp
       b = temp
       c = temp
       alpha = acos(((((a**2+b**2)/2.0_dp)-2*cos(gamma0)*a*c))/(b+c))
       temp = (beta0+gamma0)/2.0_dp
       beta = temp
       gamma = temp
    else if (case==8) then
       temp = (a0+b0+c0)/3.0_dp
       a = temp
       b = temp
       c = temp
       alpha = acos((((a**2+b**2)/2.0_dp)-cos(gamma0)*a*c-cos(beta0)*a*b)/(b*c))
       beta = beta0
       gamma = gamma0
    else if (case==9) then
       temp = (a0+b0)/2.0_dp
       a = temp
       b = temp
       c = c0
       alpha = acos((a*a)/(2*b*c))
       beta = acos((a*a)/(2*a*c))
       gamma = pi/3.0_dp
    else if (case==10) then
       temp = (a0+b0)/2.0_dp
       a = temp
       b = temp
       c = c0
       temp = (alpha0+beta0)/2.0_dp
       alpha = temp
       beta = temp
       gamma = gamma0
    else if (case==11) then
       temp = (a0+b0)/2.0_dp
       a = temp
       b = temp
       c = c0
       alpha = pi/2.0_dp
       beta = pi/2.0_dp
       gamma = pi/2.0_dp
    else if (case==12) then
       temp = (a0+b0)/2.0_dp
       a = temp
       b = temp
       c = c0
       alpha = pi/2.0_dp
       beta = pi/2.0_dp
       gamma = pi/3.0_dp
    else if (case==13) then
       temp = (a0+b0)/2.0_dp
       a = temp
       b = temp
       c = c0
       alpha = pi/2.0_dp
       beta = pi/2.0_dp
       gamma = gamma0
    else if (case==15) then
       temp = (a0+b0)/2.0_dp
       a = temp
       b = temp
       c = c0
       alpha = acos(-(a*a)/(2*b*c))
       beta = acos(-(a*a)/(2*a*c))
       gamma = pi/2.0_dp
    else if (case==16) then
       temp = (a0+b0)/2.0_dp
       a = temp
       b = temp
       c = c0
       gamma = pi/2.0_dp
       alpha = acos((((a*a+b*b)/2.0_dp-cos(gamma)*a*c)/2.0_dp)/(b*c))
       beta = acos((((a*a+b*b)/2.0_dp-cos(gamma)*a*c)/2.0_dp)/(a*c))
    else if (case==14) then
       temp = (a0+b0)/2.0_dp
       a = temp
       b = temp
       c = c0
       gamma = gamma0
       temp = (alpha0+beta0)/2.0_dp
       alpha = temp
       beta = temp
    else if (case==17) then
       temp = (a0+b0)/2.0_dp
       a = temp
       b = temp
       c = c0
       gamma = gamma0
       beta= beta0
       alpha = acos(((a*a+b*b)/2.0_dp-cos(gamma)*a*c-cos(beta)*a*b)/(b*c))
    else if (case==18) then
       temp = (c0+b0)/2.0_dp
       a = a0
       b = temp
       c = temp
       alpha = acos((a*a)/(4*b*c))
       beta = acos((a*a)/(2*a*c))
       gamma = acos((a*a)/(2*a*b))
    else if (case==19) then
       temp = (c0+b0)/2.0_dp
       a = a0
       b = temp
       c = temp
       alpha = alpha0
       beta = acos((a*a)/(2*a*c))
       gamma = acos((a*a)/(2*a*b))
    else if ((case==20) .or. (case==25)) then
       temp = (c0+b0)/2.0_dp
       a = a0
       b = temp
       c = temp
       alpha = alpha0
       temp = (beta0+gamma0)/2.0_dp
       beta = temp
       gamma = temp
    else if (case==21) then
       temp = (c0+b0)/2.0_dp
       a = a0
       b = temp
       c = temp
       alpha = pi/2.0_dp
       beta = pi/2.0_dp
       gamma = pi/2.0_dp
    else if (case==22) then
       temp = (c0+b0)/2.0_dp
       a = a0
       b = temp
       c = temp
       alpha = acos(-(b*b)/(2*b*c))
       beta = pi/2.0_dp
       gamma = pi/2.0_dp
    else if (case==23) then
       temp = (c0+b0)/2.0_dp
       a = a0
       b = temp
       c = temp
       alpha = alpha0
       beta = pi/2.0_dp
       gamma = pi/2.0_dp
    else if (case==24) then
       temp = (c0+b0)/2.0_dp
       a = a0
       b = temp
       c = temp
       beta = acos(-(a*a)/(3*a*c))
       gamma = acos(-(a*a)/(3*a*b))
       alpha = acos(((a*a+b*b)*2.0_dp -cos(gamma)*a*c-cos(beta)*a*b)/(b*c))
    else if (case==26) then
       a = a0
       b = b0
       c = c0
       alpha = acos((a*a)/(4*b*c))
       beta = acos(-(a*a)/(2*a*c))
       gamma = acos(-(a*a)/(2*a*b))
    else if (case==27) then
       a = a0
       b = b0
       c = c0
       alpha = alpha0
       beta = acos(-(a*a)/(2*a*c))
       gamma = acos(-(a*a)/(2*a*b))
    else if (case==28) then
       a = a0
       b = b0
       c = c0
       alpha = alpha0
       beta = acos(-(a*a)/(2*a*c))
       gamma = acos((2*cos(alpha)*b*c)/(a*b))
    else if (case==29) then
       a = a0
       b = b0
       c = c0
       alpha = alpha0
       beta = acos((2*cos(alpha)*b*c)/(a*c))
       gamma = acos((a*a)/(2*a*b))
    else if (case==30) then
       a = a0
       b = b0
       c = c0
       alpha = acos((b*b)/(2*b*c))
       beta = beta0
       gamma = acos(2*cos(beta)*a*c/(a*b))
    else if ((case==31) .or. (case==44)) then
       a = a0
       b = b0
       c = c0
       alpha = alpha0
       beta = beta0
       gamma = gamma0
    else if (case==32) then
       a = a0
       b = b0
       c = c0
       alpha = pi/2.0_dp
       beta = pi/2.0_dp
       gamma = pi/2.0_dp
    else if (case==40) then
       a = a0
       b = b0
       c = c0
       alpha = acos(-(b*b)/(2*b*c))
       beta = pi/2.0_dp
       gamma = pi/2.0_dp
    else if (case==35) then
       a = a0
       b = b0
       c = c0
       alpha = alpha0
       beta = pi/2.0_dp
       gamma = pi/2.0_dp
    else if (case==36) then
       a = a0
       b = b0
       c = c0
       alpha = pi/2.0_dp
       beta = acos(-(a*a)/(2*a*c))
       gamma = pi/2.0_dp
    else if (case==33) then
       a = a0
       b = b0
       c = c0
       alpha = pi/2.0_dp
       beta = beta0
       gamma = pi/2.0_dp
    else if (case==38) then
       a = a0
       b = b0
       c = c0
       alpha = pi/2.0_dp
       beta = pi/2.0_dp
       gamma = acos(-(a*a)/(2*a*b))
    else if (case==34) then
       a = a0
       b = b0
       c = c0
       alpha = pi/2.0_dp
       beta = pi/2.0_dp
       gamma = gamma0
    else if (case==42) then
       a = a0
       b = b0
       c = c0
       alpha = acos(-(b*b)/(2*b*c))
       beta = acos(-(a*a)/(2*a*c))
       gamma = pi/2.0_dp
              
    else if (case==41) then
       a = a0
       b = b0
       c = c0
       alpha = acos(-(b*b)/(2*b*c))
       beta = beta0
       gamma = pi/2.0_dp
    else if (case==37) then
       a = a0
       b = b0
       c = c0
       alpha = alpha0
       beta = acos(-(a*a)/(2*a*c))
       gamma = pi/2.0_dp
     else if (case==39) then
       a = a0
       b = b0
       c = c0
       alpha = alpha0
       beta = pi/2.0_dp
       gamma = acos(-(a*a)/(2*a*b))
     else if (case==43) then
       a = a0
       b = b0
       c = c0
       gamma = gamma0
       beta = beta0
       temp = ((a*a+b*b)/2.0_dp)-cos(gamma)*a*c-cos(beta)*a*b
       if (equal(2*temp*cos(gamma)*a*b,b*b,eps)) then
          alpha = acos(temp/(b*c))
       else
          temp = (-(a*a+b*b)/2.0_dp)-cos(gamma)*a*c-cos(beta)*a*b
          alpha = acos(temp/(b*c))
       end if
    end if

    out_params(1) = a
    out_params(2) = b
    out_params(3) = c
    out_params(4) = alpha
    out_params(5) = beta
    out_params(6) = gamma
  end SUBROUTINE symmetrize
  
  !!<summary>Fixes a POSCAR so that it perfectly matches the expected
  !!symmetries.</summary>
  !!<parameter name="in_vecs" regular="true">The original lattice
  !!vectors.</parameter>
  !!<parameter name="case" regular="true">The niggli id for this
  !!lattice.</parameter>
  !!<paramter name="out_vecs" regular="true">The corrected lattice
  !!vectors.</parameter>
  !!<parameter name="eps_" regular="true">The floating point
  !!tollerance.</parameter>
  SUBROUTINE delint(in_vecs, case, out_vecs, eps_)
    real(dp), intent(in) :: in_vecs(3,3)
    integer, intent(in) :: case
    real(dp), intent(out) :: out_vecs(3,3)
    real(dp), optional, intent(in) :: eps_

    real(dp) :: eps
    real(dp) :: in_params(6), canonical_basis(3,3), rot_mat(3,3)
    real(dp) :: rot_test(3,3), cb_inv(3,3), ident(3,3), c_params(6)
    real(dp) :: fixed_basis(3,3)
    integer :: i

    ident = reshape((/1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp/),(/3,3/))

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-3
    end if

    call basis2cell(in_vecs, in_params)
    call cell2basis(in_params, canonical_basis)
    call matrix_inverse(canonical_basis, cb_inv)
    
    rot_mat = matmul(in_vecs, cb_inv)
    rot_test = matmul(rot_mat, transpose(rot_mat))
    if (.not. equal(rot_test, ident, eps)) then
       write(*,*) "Error generating rotation matrix in delinter."
       stop
    end if

    call symmetrize(in_params, case, c_params, eps)
    call cell2basis(c_params, fixed_basis)
    out_vecs = matmul(rot_mat, fixed_basis)
    
  end SUBROUTINE delint
  
  
end Module delinter  
