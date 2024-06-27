!!<summary>Contains the subroutines needed to find a niggli reduced
!!cell and to identify it.</summary>

Module niggli
  use vector_matrix_utilities, only: determinant
  use numerical_utilities
  use num_types
  implicit none
  private
  public reduce_cell, id_cell


CONTAINS

  !!<summary>Finds the niggli reduced cell from the users given
  !!cell.</summary>
  !!<parameter name="IN" regular="true">A 3 by 3 matrix with the basis
  !!vectors as columns.</parameter>
  !!<parameter name="eps_" regular="true">Optional floating point
  !!tolerance.</parameter>
  !!<parameter name="path_" regular="true">An optional bool, true if
  !!the reduction path is to be printed at the end.</parameter>
  !!<parameter name="n_cell" regular="true">The niggli reduced
  !!cell.</parameter>
  !!<parameter name="trans" regular="true">The transformation
  !!matrix.</parameter>
  SUBROUTINE reduce_cell(IN,n_cell,trans,eps_,path_)
    real(dp), intent(in) :: IN(3,3)
    real(dp), intent(in), optional :: eps_
    logical, intent(in), optional :: path_
    real(dp), intent(out) :: n_cell(3,3)
    integer, intent(out) :: trans(3,3)

    !!<local name="vol">The vol of the input matrix.</local>
    !!<local name="eps">The floating point tolerance.</local>
    !!<local name="path">The path taken by the reduction.</local>
    !!<local name="A">The square of the first lattice vector.</local>
    !!<local name="B">The square of the second lattice vector.</local>
    !!<local name="C">The square of the third lattice vector.</local>
    !!<local name="xi">Twice the angle between the 2nd and 3rd lattice
    !!vectors.</local>
    !!<local name="eta">Twice the angle between the 3nd and 1st lattice
    !!vectors.</local>
    !!<local name="zeta">Twice the angle between the 1st and 2nd lattice
    !!vectors.</local>
    !!<local name="l">The sign of xi.</local>
    !!<local name="m">The sign of eta.</local>
    !!<local name="n">The sign of zeta.</local>
    !!<local name="reduced">True if the cell is reduced.</local>
    !!<local name="count">The number of passes through the steps that
    !!have been taken.</local>
    !!<local name="temp_lat">A temporary lattice.</local>
    !!<local name="temp_M">Temporary transoframtion matrix.</local>
    real(dp) :: vol, eps, temp_lat(3,3)
    character(len=2000):: path
    real(dp) :: A, B, C, xi, eta, zeta
    integer :: l, m, n, count, temp_M(3,3)
    logical :: reduced

    vol = determinant(IN)

    if (equal(vol,0._dp,0.000001_dp)) stop "Input matrix is linearly dependent."

    if (present(eps_)) then
       eps = eps_*abs(vol)**(1.0_dp/3.0_dp)
    else
       eps = (10.0_dp**(-5.0_dp))*abs(vol)**(1.0_dp/3.0_dp)
    end if

    path = ''
    count = 0
    reduced = .False.
    trans = transpose(reshape((/ 1, 0, 0, 0, 1, 0, 0, 0, 1 /),(/3,3/)))

    call get_params(IN,eps,A,B,C,xi,eta,zeta,l,m,n)

    do while ((.not. reduced) .and. (count <= 1000))
       reduced = .True.
       count = count + 1
       ! step #1
       if ((A>(B+eps)) .or. ((.not. (ABS(A-B)>eps)) .and. (ABS(xi)>(ABS(eta)+eps))))then
          path = trim(path) // "1"
          reduced = .False.
          trans = matmul(trans,transpose(reshape((/0, -1, 0, -1, 0, 0, 0, 0, -1/),(/3,3/))))
          temp_lat = matmul(IN,trans)
          call get_params(temp_lat,eps,A,B,C,xi,eta,zeta,l,m,n)
       end if

       ! step #2
       if ((B>(C+eps)) .or. ((.not. (ABS(C-B)>eps) .and. (ABS(eta)>ABS(zeta)+eps)))) then
          path = trim(path) // "2"
          reduced = .False.
          trans = matmul(trans, transpose(reshape((/-1, 0, 0, 0, 0, -1, 0, -1, 0/),(/3,3/))))
          temp_lat = matmul(IN,trans)
          call get_params(temp_lat,eps,A,B,C,xi,eta,zeta,l,m,n)
          cycle
       end if

       ! step #3
       if (l*m*n==1) then
          path = trim(path) // "3"
          call find_C3(l,m,n,temp_M)
          if (.not. (all(abs(temp_M-transpose(reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/),(/3,3/))))==0))) then
             reduced = .False.
          end if
          trans = matmul(trans,temp_M)
          temp_lat = matmul(IN,trans)
          call get_params(temp_lat,eps,A,B,C,xi,eta,zeta,l,m,n)
       end if

       ! step #4
       if ((l*m*n==0) .or. (l*m*n==(-1))) then
          path = trim(path) // "4"
          call find_C4(l,m,n,temp_M)
          if (.not. (all(abs(temp_M-transpose(reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/),(/3,3/))))==0))) then
             reduced = .False.
          end if
          trans = matmul(trans,temp_M)
          temp_lat = matmul(IN,trans)
          call get_params(temp_lat,eps,A,B,C,xi,eta,zeta,l,m,n)
       end if

       ! step #5
       if ((ABS(xi)>(B+eps)) .or. ((.not. (ABS(B-xi)>eps)) .and. (2.0_dp*eta<(zeta-eps))) .or. ((.not. (ABS(B+xi)>eps)) .and. (zeta<(-eps)))) then
          path = trim(path) // "5"
          reduced = .False.
          trans = matmul(trans,transpose(reshape((/1, 0, 0, 0, 1, -get_sign(xi), 0, 0, 1/),(/3,3/))))
          temp_lat = matmul(IN,trans)
          call get_params(temp_lat,eps,A,B,C,xi,eta,zeta,l,m,n)
          cycle
       end if

       ! step #6
       if ((ABS(eta)>(A+eps)) .or. ((.not. (ABS(A-eta)>eps)) .and. (2.0_dp*xi<(zeta-eps))) .or. ((.not. (ABS(A+eta)>eps)) .and. (zeta<(-eps)))) then
          path = trim(path) // "6"
          reduced = .False.
          trans = matmul(trans,transpose(reshape((/1, 0, -get_sign(eta), 0, 1, 0, 0, 0, 1/),(/3,3/))))
          temp_lat = matmul(IN,trans)
          call get_params(temp_lat,eps,A,B,C,xi,eta,zeta,l,m,n)
          cycle
       end if

       ! step #7
       if ((ABS(zeta)>(A+eps)) .or. ((.not. (ABS(A-zeta)>eps)) .and. (2.0_dp*xi<(eta-eps))) .or. ((.not. ABS(A+zeta)>eps) .and. (eta<(-eps)))) then
          path = trim(path) // "7"
          reduced = .False.
          trans = matmul(trans,transpose(reshape((/1, -get_sign(zeta), 0, 0, 1, 0, 0, 0, 1/),(/3,3/))))
          temp_lat = matmul(IN,trans)
          call get_params(temp_lat,eps,A,B,C,xi,eta,zeta,l,m,n)
          cycle
       end if

       ! step #8
       if (((xi+eta+zeta+A+B)<(-eps)) .or. ((.not. (ABS(xi+eta+zeta+A+B)>eps)) .and. ((2.0_dp*(A+eta)+zeta)>eps))) then
          path = trim(path) // "8"
          reduced = .False.
          trans = matmul(trans,transpose(reshape((/1, 0, 1, 0, 1, 1, 0, 0, 1/),(/3,3/))))
          temp_lat = matmul(IN,trans)
          call get_params(temp_lat,eps,A,B,C,xi,eta,zeta,l,m,n)
          cycle
       end if
    end do

    if (present(path_)) then
       if (path_) then
          print *, trim(path)
       end if
    end if

    if (count >= 1000) stop "Could not reduce the cell in 1000 iterations. This could be because of floating point error, try providing a smaller eps."

    if (.not. (condition_check(A,B,C,xi,eta,zeta,eps))) stop "Cell reduction failed to satisfy niggli conditions."

    n_cell = matmul(IN,trans)

  end SUBROUTINE reduce_cell

  !!<summary>Returns .True. if the niggli cell is completely reduced,
  !!.False. if otherwise.</summary>
  !!<parameter name="A" regular="true">The value of
  !!IN(:,1).IN(:,1)</parameter>
  !!<parameter name="B" regular="true">The value of
  !!IN(:,2).IN(:,2)</parameter>
  !!<parameter name="C" regular="true">The value of
  !!IN(:,3).IN(:,3)</parameter>
  !!<parameter name="xi" regular="true">The value of
  !!2*IN(:,2).IN(:,3)</parameter>
  !!<parameter name="eta" regular="true">The value of
  !!2*IN(:,3).IN(:,1)</parameter>
  !!<parameter name="zeta" regular="true">The value of
  !!2*IN(:,1).IN(:,2)</parameter>
  !!<parameter name="eps" regular="true">The floating point
  !!tolerance.</parameter>
  function condition_check(A,B,C,xi,eta,zeta,eps)
    logical :: condition_check
    real(dp), intent(in) :: A, B, C, xi, eta, zeta, eps

    condition_check = .True.

    if (.not. (((A-eps)>0.0_dp) .and. ((A<(B-eps)) .or. (ABS(A-B)<eps)) .and. ((B<(C-eps)) .or. (ABS(B-C)<eps)))) then
       condition_check = .False.
    end if

    if ((ABS(A-B)<eps) .and. (.not. ((ABS(xi)<(ABS(eta)-eps)) .or. (ABS(ABS(xi)-ABS(eta))<eps)))) then
       condition_check = .False.
    end if

    if ((ABS(B-C)<eps) .and. (.not. ((ABS(eta)<(ABS(zeta)-eps)) .or. (ABS(ABS(eta)-ABS(zeta))<eps)))) then
       condition_check = .False.
    end if

    if (.not. ((((xi-eps)>0.0_dp) .and. ((eta-eps)>0.0_dp) .and. ((zeta-eps)>0.0_dp)) .or. (((xi<0.0_dp) .or. (ABS(xi)<eps)) .and. ((eta<0.0_dp) .or. (ABS(eta)<eps)) .and. ((zeta<0.0_dp) .or. (ABS(zeta)<eps))))) then
       condition_check = .False.
    end if

    if (.not. ((ABS(xi)<(B-eps)) .or. (ABS(ABS(xi)-B)<eps))) then
       condition_check = .False.
    end if

    if (.not. (((ABS(eta)<(A-eps)) .or. (ABS(ABS(eta)-A)<eps)) .and. ((ABS(zeta)<(A-eps)) .or. (ABS(ABS(zeta)-A)<eps)))) then
       condition_check = .False.
    end if

    if (.not. ((C<(A+B+C+xi+eta+zeta-eps)) .or. (ABS(C-(A+B+C+xi+eta+zeta))<eps))) then
       condition_check = .False.
    end if

    if ((ABS(xi-B)<eps) .and. (.not. ((zeta<((2.0_dp*eta)-eps)) .or. (ABS(zeta-(2.0_dp*eta))<eps)))) then
       condition_check = .False.
    end if

    if ((ABS(eta-A)<eps) .and. (.not. ((zeta<((2.0_dp*xi)-eps)) .or. (ABS(zeta-(2.0_dp*xi))<eps)))) then
       condition_check = .False.
    end if

    if ((ABS(zeta-A)<eps) .and. (.not. ((eta<((2.0_dp*xi)-eps)) .or. (ABS(eta-(2.0_dp*xi))<eps)))) then
       condition_check = .False.
    end if

    if ((ABS(xi+B)<eps) .and. (.not. (ABS(zeta)<eps))) then
       condition_check = .False.
    end if

    if ((ABS(eta+A)<eps) .and. (.not. (ABS(zeta)<eps))) then
       condition_check = .False.
    end if

    if ((ABS(zeta+A)<eps) .and. (.not. (ABS(eta)<eps))) then
       condition_check = .False.
    end if

    if ((ABS(C-(A+B+C+xi+eta+zeta))<eps) .and. (.not. (((2.0_dp*A+2.0_dp*eta+zeta)<(-eps)) .or. (ABS(2.0_dp*A+2.0_dp*eta+zeta)<eps)))) then
       condition_check = .False.
    end if

  end function condition_check

  !!<summary>Finds the sign (-1,0,0) of a real number.</summary>
  !!<parameter name="a" regular="true">A real number.</parameter>
  function get_sign(a)
    integer :: get_sign
    real(dp), intent(in) :: a

    if (a>0) then
       get_sign = 1
    else if (a<0) then
       get_sign = -1
    else
       get_sign = 0
    end if
  end function get_sign

  !!<summary>Finds the transformation matrix for the operations of
  !!step 4.</summary>
  !!<parameter name="l" regular="true">The sign of xi.</parameter>
  !!<parameter name="m" regular="true">The sign of eta.</parameter>
  !!<parameter name="n" regular="true">The sign of zeta.</parameter>
  !!<parameter name="trans" regular="true">The transformation
  !!matrix.</parameter>
  SUBROUTINE find_C4(l,m,n,trans)
    integer, intent(in) :: l,m,n
    integer, intent(out) :: trans(3,3)

    integer :: i, j, k, r

    if ((l==(-1)) .and. (m==(-1)) .and. (n==(-1))) then
       i = 1
       j = 1
       k = 1
    else

       i = 1
       j = 1
       k = 1
       r = -1

       if (l==1) then
          i = -1
       else if (l==0) then
          r = 0
       end if

       if (m==1) then
          j = -1
       else if (m==0) then
          r = 1
       end if

       if (n==1) then
          k = -1
       else if (n==0) then
          r = 2
       end if

       if (i*j*k == (-1)) then
          if (r==0) then
             i = -1
          else if (r==1) then
             j = -1
          else if (r==2) then
             k = -1
          end if
       end if
    end if

    trans = transpose(reshape((/i, 0, 0, 0, j, 0, 0, 0, k/),(/3,3/)))

  end SUBROUTINE find_C4

  !!<summary>Finds the transformation matrix for the operation of step
  !!3.</summary>
  !!<parameter name="l" regular="true">The sign of xi.</parameter>
  !!<parameter name="m" regular="true">The sign of eta.</parameter>
  !!<parameter name="n" regular="true">The sign of zeta.</parameter>
  !!<parameter name="trans" regular="true">The transformation
  !!matrix.</parameter>
  SUBROUTINE find_C3(l,m,n,trans)
    integer, intent(in) :: l,m,n
    integer, intent(out) :: trans(3,3)

    integer :: i, j, k

    if (l==(-1)) then
       i = -1
    else
       i = 1
    end if

    if (m==(-1)) then
       j = -1
    else
       j = 1
    end if

    if (n==(-1)) then
       k = -1
    else
       k = 1
    end if

    trans = transpose(reshape((/i, 0, 0, 0, j, 0, 0, 0, k/),(/3,3/)))

  end SUBROUTINE find_C3

  !!<summary>Gets the G vector (A,B,C,xi,eta,zeta) for a given input
  !!matrix.</summary>
  !!<parameter name="IN" regular="true">The matrix being
  !!reduced.</parameter>
  !!<parameter name="eps" regular="true">The floating point
  !!tolerance.</parameter>
  !!<parameter name="A" regular="true">The value of
  !!IN(:,1).IN(:,1)</parameter>
  !!<parameter name="B" regular="true">The value of
  !!IN(:,2).IN(:,2)</parameter>
  !!<parameter name="C" regular="true">The value of
  !!IN(:,3).IN(:,3)</parameter>
  !!<parameter name="xi" regular="true">The value of
  !!2*IN(:,2).IN(:,3)</parameter>
  !!<parameter name="eta" regular="true">The value of
  !!2*IN(:,3).IN(:,1)</parameter>
  !!<parameter name="zeta" regular="true">The value of
  !!2*IN(:,1).IN(:,2)</parameter>
  !!<parameter name="l" regular="true">1 if xi is positive, -1 if it
  !!is negatev, 0 otherwise.</parameter>
  !!<parameter name="m" regular="true">1 if eta is positive, -1 if it
  !!is negatev, 0 otherwise.</parameter>
  !!<parameter name="n" regular="true">1 if zeta is positive, -1 if it
  !!is negatev, 0 otherwise.</parameter>
  SUBROUTINE get_params(IN,eps,A,B,C,xi,eta,zeta,l,m,n)
    real(dp), intent(in) :: IN(3,3), eps
    real(dp), intent(out) :: A, B, C, xi, eta, zeta
    integer, intent(out) :: l, m, n

    A  = dot_product(IN(:,1),IN(:,1))
    B  = dot_product(IN(:,2),IN(:,2))
    C  = dot_product(IN(:,3),IN(:,3))
    xi  = 2.0_dp*dot_product(IN(:,2),IN(:,3))
    eta  = 2.0_dp*dot_product(IN(:,3),IN(:,1))
    zeta  = 2.0_dp*dot_product(IN(:,1),IN(:,2))

    if (xi<-eps) then
       l = -1
    else if (xi>eps) then
       l = 1
    else
       l = 0
    end if

    if (eta<-eps) then
       m = -1
    else if (eta>eps) then
       m = 1
    else
       m = 0
    end if

    if (zeta<-eps) then
       n = -1
    else if (zeta>eps) then
       n = 1
    else
       n = 0
    end if

  end SUBROUTINE get_params

  !!<summary>Identifies which of the niggli cells is currently being
  !!used.</summary>
  !!<parameter name="U" regular="true">The users basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli cell
  !!basis.</parameter>
  !!<parameter name="No" regular="true">The niggli basis for our
  !!cell.</parameter>
  !!<parameter name="O" regular="true">Our basis for the
  !!cell.</parameter>
  !!<parameter name="Cu" regular="true">The transformation for the
  !!users basis.</parameter>
  !!<parameter name="Co" regular="true">The transformaiton for our
  !!basis.</parameter>
  !!<parameter name="id" regular="true">The niggli cell
  !!number.</parameter>
  !!<parameter name="s_range" regular="true">The range of determinant
  !!sizes to search over.</parameter>
  !!<parameter name="eps_" regular="true">Optional floating point
  !!tolerance.</parameter>
  SUBROUTINE id_cell(U, Nu, Cu, O, No, Co, id, s_range, eps_)
    real(dp), intent(in) :: U(3,3)
    real(dp), intent(in), optional :: eps_
    real(dp), intent(out) :: O(3,3), Nu(3,3), No(3,3)
    integer, intent(out) :: Cu(3,3), Co(3,3)
    integer, intent(out) :: id, s_range

    real(dp) :: eps
    real(dp) :: temp_a(3), temp_b(3), temp_c(3), A, B, C, D, E, F
    logical :: positive

    if (present(eps_)) then
       eps = eps_
    else
       eps = (10.0_dp**(-5.0_dp))
    end if

    call reduce_cell(U,Nu,Cu,eps_=eps)
    !! Integer transformation check TODO

    temp_a = Nu(:,1)
    temp_b = Nu(:,2)
    temp_c = Nu(:,3)

    A = dot_product(temp_a,temp_a)
    B = dot_product(temp_b,temp_b)
    C = dot_product(temp_c,temp_c)
    D = dot_product(temp_b,temp_c)
    E = dot_product(temp_a,temp_c)
    F = dot_product(temp_a,temp_b)

    positive = .False.
    if (((D-eps)>0) .and. ((E-eps)>0) .and. ((F-eps)>0)) then
       positive = .True.
    end if

    id = (-1)
    if (isclose(A,B,atol_=eps) .and. isclose(B,C,atol_=eps)) then !!lengths the same
       if (positive .eqv. .True.) then !! all angles acute?
          if (isclose(D,E,atol_=eps) .and. isclose(D,F,atol_=eps)) then !!all angles the same?
             if (isclose((A/2.0_dp),D,atol_=eps)) then
                id = 1 !! fcc
                s_range = 1
                O = reshape((/0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
                     1.0_dp, 0.0_dp/),(/3,3/))
             else
                id = 2 !! Acute Rhombohedral, angles not 90 or 60 but all the same
                s_range = 25
                O = reshape((/-1.11652_dp, -0.610985_dp, 0.616515_dp, 0.0_dp, &
                     -1.32288_dp, -0.5_dp, 1.0_dp, 1.32288_dp, 1.5_dp/),(/3,3/))
             end if
          end if
       else !! not acute angles.
          if (isclose(D,E,atol_=eps) .and. isclose(D,F,atol_=eps)) then !! all angles the same.
             if (isclose(0.0_dp,D,atol_=eps)) then !! all angles close to 0.
                id = 3 !! simple cubic.
                s_range = 1
                O = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                     0.0_dp, 1.0_dp/),(/3,3/))
             else if (isclose((-A/3.0_dp),D,atol_=eps)) then
                id = 5 !! bcc
                s_range = 1
                O = reshape((/-1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, &
                     1.0_dp, 1.0_dp, -1.0_dp/),(/3,3/))
             else
                id = 4 !! Obtuse Rhombohedral, angles not 90 or 60 but all the same
                s_range = 25
                O = reshape((/-0.548584_dp, 0.774292_dp, 1.04858_dp, 0.0_dp, &
                     -1.32288_dp, 0.5_dp, 1.0_dp, 1.32288_dp, 0.5_dp/),(/3,3/))
             end if
          else if (isclose((2.0_dp*ABS(D+E+F)),(A+B),atol_=eps)) then !! not all angles the same.
             if (isclose(D,E,atol_=eps)) then
                id = 6 !! body-centered Tetragonal
                s_range = 25
                O = reshape((/1.80278_dp, -1.47253_dp, 0.762655_dp, 2.80278_dp, &
                     0.13535_dp, -0.791285_dp, 0.80278_dp, -0.47253_dp, 2.762655_dp/),(/3,3/))
             else if (isclose(E,F,atol_=eps)) then
                id = 7 !! body-centered Tetragonal
                s_range = 25
                O = reshape((/1.95095_dp, 1.19163_dp, 0.879663_dp, 0.0_dp, &
                     2.60788_dp, 0.44606_dp, 0.95095_dp, -0.41625_dp, 2.433603_dp/),(/3,3/))
             else
                id = 8 !! body-centered Orthorhombic
                s_range = 15
                O = reshape((/1.41144_dp, 0.0885622_dp, -2.0_dp, -0.99868_dp, 2.21232_dp,&
                     1.268178_dp, 3.41012_dp, -1.1237578_dp, -1.268178_dp/),(/3,3/))
             end if
          end if !! 2|D+E+F|==(A+B)
       end if !! are the angles acute?
    end if !! are all the lattice vector lengths the same?

    if (isclose(A,B,atol_=eps) .and. (id==(-1))) then !! A==B!=C
       if (positive .eqv. .True.) then !! all angles acute.
          if (isclose(D,E,atol_=eps) .and. isclose(D,F,atol_=eps) .and. &
               isclose((A/2.0_dp),D,atol_=eps)) then !! all angles the same and 60 degrees.
             id = 9 !! Hexagonal Rhombohedral
             s_range = 25
             O = reshape((/1.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, &
                  3.0_dp, 3.0_dp/),(/3,3/))
          else if (isclose(D,E,atol_=eps)) then
             id = 10 !! base-centered Monoclinic
             s_range = 10
             O = reshape((/1.0_dp, -1.0_dp, 1.0_dp, -1.46391_dp, 0.0_dp, 1.96391_dp, &
                  0.0_dp, 2.0_dp, 0.0_dp/),(/3,3/))
          end if
       else !! some angles obtuse.
          if (isclose(D,E,atol_=eps) .and. isclose(D,F,atol_=eps) .and. &
               isclose(0.0_dp,D,atol_=eps)) then !! all angles close to 0.
             id = 11 !! Primitive Tetragonal
             s_range = 25
             O = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                  0.0_dp, 2.0_dp/),(/3,3/))
          else if (isclose(D,E,atol_=eps)) then !! two angles the same
             if (isclose(0.0_dp,D,atol_=eps) .and. isclose((-A/2.0_dp),F,atol_=eps)) then !! D==E==0 and F==120
                id = 12 !! Primitive Hexagonal 120 degrees.
                s_range = 25
                O = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, -0.8660254037844386_dp, &
                     0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp/),(/3,3/))
             else if (isclose((-A/2.0_dp),D,atol_=eps) .and. isclose(0.0_dp,F,atol_= eps)) then !! D==E==120 and F==0
                id = 15 !! body-centered Tetragonal
                s_range = 25
                O = reshape((/-1.0_dp, -1.0_dp, 2.0_dp, 0.0_dp, -2.0_dp, 0.0_dp, -2.0_dp, &
                     0.0_dp, 0.0_dp/),(/3,3/))
             else if (isclose(0.0_dp,D,atol_=eps)) then !! D==E==0 and no constraint on F.
                id = 13 !! base-centered Orthorhombic
                s_range = 15
                O = reshape((/1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, 0.0_dp, &
                     -1.73205_dp, 1.73205_dp/),(/3,3/))
             else if (isclose((2.0_dp*ABS(D+E+F)),(A+B),atol_=eps)) then !! 2|D+E+F|==(A+B)
                id = 16 !! face-centered Orthorhombic
                s_range = 20
                O = reshape((/1.04442_dp, 1.43973_dp, 1.68415_dp, 0.779796_dp, -1.1789_dp, &
                     1.0_dp, 1.779796_dp, -0.1789_dp, 0.0_dp/),(/3,3/))
             else
                id = 14 !! base-centered Monoclinic
                s_range = 10
                O = reshape((/-1.0_dp, 1.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 2.0_dp, 0.0_dp, &
                     -2.0_dp, 0.0_dp/),(/3,3/))
             end if
          else if (isclose((2.0_dp*ABS(D+E+F)),(A+B),atol_=eps)) then
             id = 17 !! base-centered Monoclinic
             s_range = 10
             O = reshape((/-1.05387_dp, -1.61088_dp, 1.51474_dp, -0.244302_dp, -2.77045_dp, &
                  0.51474_dp, 1.809568_dp, -0.15957_dp, 0.0_dp/),(/3,3/))
          end if
       end if !! acute or obtuse angles.
    end if !! A==B!=C

    if (isclose(B,C,atol_=eps) .and. (id==(-1))) then !! B==C!=A
       if (positive .eqv. .True.) then !! all angles acute
          if (isclose(E,F,atol_=eps)) then !! E==F
             if (isclose((A/4.0_dp),D,atol_=eps) .and. isclose((A/2.0_dp),E,atol_=eps)) then
                id = 18
                s_range = 25
                O = reshape((/-2.0_dp, -1.0_dp, 1.0_dp, -3.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, &
                     -3.0_dp, 0.0_dp/),(/3,3/))
             else if (isclose((A/2.0_dp),E,atol_=eps)) then
                id = 19
                s_range = 15
                O = reshape((/0.5_dp, 1.0_dp, 1.5_dp, 0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, &
                     0.0_dp, 3.0_dp/),(/3,3/))
             else
                id = 20
                s_range = 10
                O = reshape((/1.0_dp, 1.0_dp, 1.0_dp, 1.70119_dp, -1.45119_dp, 1.0_dp, &
                     0.69779_dp, -1.4322505_dp, 3.23446_dp/),(/3,3/))
             end if
          end if
       else
          if (isclose(E,F,atol_= eps)) then
             if (isclose(0.0_dp,D,atol_=eps) .and. isclose(0.0_dp,E,atol_=eps)) then
                id = 21
                s_range = 25
                O = reshape((/0.0_dp, 0.0_dp, 0.5_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                     1.0_dp, 0.0_dp/),(/3,3/))
             else if (isclose((-B/2.0_dp),D,atol_=eps) .and. isclose(0.0_dp,E,atol_=eps)) then
                id = 22
                s_range = 25
                O = reshape((/0.0_dp, 0.0_dp, -0.5_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                     -0.5_dp, 0.8660254037844386_dp, 0.0_dp/),(/3,3/))
             else if (isclose(0.0_dp,E,atol_=eps)) then
                id = 23
                s_range = 15
                O = reshape((/-0.3333333_dp, -1.54116_dp, 1.87449_dp, 1.0_dp, 1.0_dp, &
                     1.0_dp, 2.0_dp, -1.0_dp, -1.0_dp/),(/3,3/))
             else if (isclose((2.0_dp*ABS(D+E+F)),(A+B),atol_=eps) .and. &
                  isclose((-A/3.0_dp),E,atol_=eps)) then
                id = 24
                s_range = 25
                O = reshape((/-0.255922_dp, -1.44338_dp, 0.92259_dp, 1.51184_dp, &
                     0.0_dp, -0.845178_dp, 1.255922_dp, 1.44338_dp, 0.07741_dp/),(/3,3/))
             else
                id = 25
                s_range = 10
                O = reshape((/1.0_dp, 1.0_dp, 1.0_dp, 1.45119_dp, -1.70119_dp, &
                     -1.0_dp, 0.28878_dp, -3.26895_dp, 0.48018_dp/),(/3,3/))
             end if
          end if
       end if
    end if

    if (id==(-1)) then
       if (positive .eqv. .True.) then
          if (isclose(E,F,atol_=eps)) then
             !! Why did we not use isclose()? TODO
             if ((ABS((A/4.0_dp)-D)<= ABS(eps*D)) .and. (ABS((A/2.0_dp)-E)<= ABS(eps*E))) then
                id = 26
                s_range = 20
                O = reshape((/0.0_dp, 1.0_dp, 1.5_dp, 0.5_dp, 0.0_dp, 1.5_dp, 0.0_dp, &
                     0.0_dp, 3.0_dp/),(/3,3/))
             else if (isclose((A/2.0_dp),E,atol_=eps)) then
                id = 27
                s_range = 10
                O = reshape((/-1.464824_dp,0.464824_dp,1.907413_dp,-0.153209_dp,0.153209_dp, &
                     -2.907413_dp,1.0_dp,1.0_dp,0.0_dp/),(/3,3/))
             else ! Default to Triclinic (same behavior as VASP)
                id = 31
                s_range = 1
                O = U
             end if
          else !! E==F
             if (isclose((A/2.0_dp),E,atol_=eps) .and. isclose((2.0_dp*D),F,atol_=eps)) then
                id = 28
                s_range = 10
                O = reshape((/-1.44896_dp, 0.948958_dp, -1.0_dp, -1.0_dp, -1.0_dp, &
                     0.0_dp, 0.342424_dp, -1.342424_dp, -2.02006_dp/),(/3,3/))
             else if (isclose((A/2.0_dp),F,atol_=eps) .and. &
                  isclose((2.0_dp*D),E,atol_=eps)) then
                id = 29
                s_range = 10
                O = reshape((/-0.666125_dp, 1.16613_dp, 2.04852_dp, 1.0_dp, 1.0_dp, &
                     0.0_dp, 1.61803_dp, -0.618034_dp, 1.0_dp/),(/3,3/))
             else if (isclose((B/2.0_dp),D,atol_=eps) .and. &
                  isclose((2.0_dp*E),F,atol_=eps)) then
                id = 30
                s_range = 10
                O = reshape((/1.0_dp, 1.0_dp, 0.0_dp, 1.61803_dp, -0.618034_dp, 1.0_dp, &
                     -0.0361373_dp, 0.536137_dp, 2.38982_dp/),(/3,3/))
             else
                id = 31
                s_range = 1
                O = U
             end if
          end if
       else
          if (isclose(E,F,atol_=eps) .and. isclose(0.0_dp,E,atol_=eps)) then
             if (isclose(0.0_dp,D,atol_=eps)) then
                id = 32
                s_range = 15
                O = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, &
                     0.0_dp, 3.0_dp/),(/3,3/))
             else if (isclose(D,(-B/2.0_dp),atol_=eps)) then
                id = 40
                s_range = 15
                O = reshape((/1.0_dp, 1.0_dp, 1.0_dp, 1.61803_dp, -0.618034_dp, &
                     -1.0_dp, -1.05557_dp, 1.99895_dp, -0.943376_dp/),(/3,3/))
             else
                id = 35
                s_range = 10
                O = reshape((/-0.668912_dp, 1.96676_dp, -1.29785_dp, 1.61803_dp, &
                     -0.618034_dp, -1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp/),(/3,3/))
             end if
          else if (isclose(D,F,atol_=eps) .and. isclose(0.0_dp,D,atol_=eps)) then
             if (isclose((-A/2.0_dp),E,atol_=eps)) then
                id = 36
                s_range = 15
                O = reshape((/1.0_dp, 1.0_dp, 1.0_dp, 1.41421_dp, -1.41421_dp, &
                     0.0_dp, -1.43541_dp, -1.43541_dp, 1.37083_dp/),(/3,3/))
             else
                id = 33
                s_range = 10
                O = reshape((/2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 0.5_dp, &
                     0.0_dp, 2.0_dp/),(/3,3/))
             end if
          else if (isclose(D,E,atol_=eps) .and. isclose(0.0_dp,D,atol_=eps)) then
             if (ABS((-A/2.0_dp)-F)<= ABS(eps*F)) then
                id = 38
                s_range = 15
                O = reshape((/0.5_dp, 1.0_dp, 0.0_dp, 0.5_dp, -1.0_dp, 0.0_dp, 0.0_dp, &
                     0.0_dp, 3.0_dp/),(/3,3/))
             else
                id = 34
                s_range = 10
                O = reshape((/1.0_dp, 1.0_dp, 1.0_dp, 1.22474487_dp, -1.22474487_dp, &
                     -1.0_dp, -0.16598509_dp, -1.64308297_dp, 1.80906806_dp/),(/3,3/))
             end if
          else
             if (isclose((-B/2.0_dp),D,atol_=eps) .and. isclose((-A/2.0_dp),E,atol_=eps) &
                  .and. isclose(0.0_dp,F,atol_=eps)) then
                id = 42
                s_range = 20
                O = reshape((/-1.53633_dp, 1.36706_dp, -1.33073_dp, 1.0_dp, 1.0_dp, &
                     1.0_dp, 1.61803_dp, -0.61803_dp, -1.0_dp/),(/3,3/))
             else if (isclose((-B/2.0_dp),D,atol_=eps) .and. isclose(0.0_dp,F,atol_=eps)) then
                id = 41
                s_range = 10
                O = reshape((/-1.85397_dp, -0.854143_dp, 1.35397_dp, 1.0_dp, 0.0_dp, &
                     1.0_dp, 1.0_dp, -1.41421_dp, -1.0_dp/),(/3,3/))
             else if (isclose((-A/2.0_dp),E,atol_=eps) .and. isclose(0.0_dp,F,atol_=eps)) then
                id = 37
                s_range = 10
                O = reshape((/-1.79092_dp, -1.47209_dp, 0.790922_dp, 1.0_dp, -1.41421_dp, &
                     -1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp/),(/3,3/))
             else if (isclose(0.0_dp,E,atol_=eps) .and. isclose((-A/2.0_dp),F,atol_=eps)) then
                id = 39
                s_range = 10
                O = reshape((/0.0_dp, -1.73205_dp, -1.0_dp, -1.66542_dp, -0.672857_dp, &
                     1.66542_dp, 1.0_dp, 0.0_dp, 1.0_dp/),(/3,3/))
             else if (isclose((2.0_dp*ABS(D+E+F)),(A+B),atol_=eps) .and. &
                  isclose(ABS(2.0_dp*D+F),B,atol_=eps)) then
                id = 43
                s_range = 10
                O = reshape((/-0.39716_dp, -0.34718_dp, 2.49434_dp, 2.64194_dp, &
                     -0.14194_dp, 0.0_dp, -1.39716_dp, -1.34718_dp, 1.49434_dp/),(/3,3/))
             else
                id = 44
                s_range = 1
                O = U
             end if
          end if
       end if
    end if

    call reduce_cell(O,No,Co,eps_=eps)

  end SUBROUTINE id_cell

  !!<summary>Matches the numpy allclose function for single floats.</summary>
  !!<parameter name="A" regular="true">The first value to compare.</parameter>
  !!<parameter name="B" regular="true">The second value to compare.</parameter>
  !!<parameter name="rtol_" regular="true">Relative tolerance.</parameter>
  !!<parameter name="atol_" regular="true">Absolute tolerance.</parameter>
  Function isclose(A,B,rtol_,atol_)
    logical :: isclose
    real(dp), intent(in) :: A, B
    real(dp), optional, intent(in) :: rtol_, atol_

    real(dp) :: rtol, atol

    if (present(rtol_)) then
       rtol = rtol_
    else
       rtol = 10.0_dp**(-5.0_dp)
    end if

    if (present(atol_)) then
       atol = atol_
    else
       atol = 10.0_dp**(-8.0_dp)
    end if

    if (ABS(A-B) <= (atol+rtol*ABS(B))) then
       isclose = .True.
    else
       isclose = .False.
    end if

  end Function isclose

end Module niggli
