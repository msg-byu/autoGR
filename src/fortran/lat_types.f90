!!<summary>Code used to identify the lattice type from the basis vectors.</summary>

Module lat_types
  use symmetry, only: get_lattice_pointGroup
  use vector_matrix_utilities, only: minkowski_reduce_basis, norm
  use num_types
  implicit none
  private
  public identify_lattice, canonical_basis

CONTAINS

  !!<summary>Finds the lattice type for the vectors
  !!provided.</summary>
  !!<parameter name="lat_vecs" regular="true">The lattice vectors in
  !!the users basis.</parameter>
  !!<parameter name="lat_id" regular="true">The integer that
  !!corresponds to the lattice type ('simple cubic' =&gt 1, 'body
  !!centered cubic' =&gt 2, 'face centered cubic' =&gt 3, 'hexagonal'
  !!=&gt 4, 'trigonal' =&gt 5, 'simple tetragonal' =&gt 6, 'body
  !!centered tetragonal' =&gt 7, 'simple orthorhombic' =&gt 8, 'base
  !!centered orthorhombic' =&gt 9, 'body centered orthorhombic'
  !!=&gt10, 'face centered orthorhombic' =&gt 11, 'simplo monoclinic'
  !!=&gt 12, 'base centered monoclinic' =&gt 13, 'triclinic' =&gt
  !!14).</parameter>
  !!<parameter name="eps" regular="true">Finite precision
  !!tolerance.</parameter>
  SUBROUTINE identify_lattice(lat_vecs, lat_id, eps)
    real(dp), intent(in) :: lat_vecs(3,3), eps
    integer, intent(out) :: lat_id

    real(dp) :: r_vecs(3,3)
    real(dp) :: a1(3), a2(3), a3(3)
    real(dp) :: a1n, a2n, a3n, a1ta2, a1ta3, a2ta3
    real(dp) :: pa2ta1n, pa2ta3n, pa1ta2n, pa1ta3n, pa3ta1n, pa3ta2n
    integer :: pcount, fam
    real(dp) :: tiny
    real(dp), pointer :: lattpg(:,:,:)
    
    call minkowski_reduce_basis(lat_vecs,r_vecs,eps)
    
    a1 = r_vecs(:,1)
    a2 = r_vecs(:,2)
    a3 = r_vecs(:,3)

    a1n = norm(a1)
    a2n = norm(a2)
    a3n = norm(a3)
    a1ta2 = abs(dot_product(a1,a2)/abs(a1n*a2n))
    a1ta3 = abs(dot_product(a1,a3)/abs(a1n*a3n))
    a2ta3 = abs(dot_product(a2,a3)/abs(a2n*a3n))
    
    tiny = 1E-6
    pcount = 0
    if (abs(a1ta2)<tiny) then
       pcount = pcount + 1
    end if
    if (abs(a1ta3)<tiny) then
       pcount = pcount + 1
    end if
    if (abs(a2ta3)<tiny) then
       pcount = pcount + 1
    end if

    call get_lattice_pointGroup(r_vecs,lattpg,eps_=eps)
    fam = size(lattpg,3)

    if (fam==48) then
       if (pcount==3) then
          lat_id = 1
       else if ((abs(a1ta2 -(1.0_dp/3.0_dp)) < tiny) .or. (abs(a1ta3 -(1.0_dp/3.0_dp)) <tiny) .or. (abs(a2ta3 -(1.0_dp/3.0_dp)) < tiny)) then
          lat_id = 2
       else if ((abs(a1ta2 -0.5_dp) <tiny) .or. (abs(a1ta3-0.5_dp)<tiny) .or. (abs(a2ta3-0.5_dp)<tiny)) then
          lat_id = 3
       else
          stop "Could not identify lattice type (cubic) in lat_id."
       end if
    else if (fam==24) then
       lat_id = 4
    else if (fam==12) then
       lat_id = 5
    else if (fam==16) then
       if (pcount==3) then
          lat_id = 6
       else
          lat_id = 7
       end if
    else if (fam==8) then
       pa2ta1n = abs(dot_product(a2,a1)/dot_product(a1,a1))
       pa2ta3n = abs(dot_product(a2,a3)/dot_product(a3,a3))
       pa1ta2n = abs(dot_product(a1,a2)/dot_product(a2,a2))
       pa1ta3n = abs(dot_product(a1,a3)/dot_product(a3,a3))
       pa3ta1n = abs(dot_product(a3,a1)/dot_product(a1,a1))
       pa3ta2n = abs(dot_product(a3,a2)/dot_product(a2,a2))
       if (pcount==3) then
          lat_id = 8
       else if (pcount ==2) then
          if (((abs(abs(a1ta2)-0.5_dp)<tiny) .and. (abs(a1n-a2n)<tiny)) .or. &
               ((abs(abs(a1ta3)-0.5_dp)<tiny) .and. (abs(a1n-a3n)<tiny)) .or. &
               ((abs(abs(a2ta3)-0.5_dp)<tiny) .and. (abs(a2n-a3n)<tiny))) then
             lat_id = 4
          else
             lat_id = 9
          end if
       else if (pcount==1) then
          if (abs(a1ta2)<tiny) then
             if (((pa3ta1n-(a1n/2.0_dp))<tiny) .and. ((pa3ta2n-(a2n/2.0_dp))<tiny)) then
                lat_id = 10
             else
                lat_id = 11
             end if
          else if (abs(a1ta3)<tiny) then
             if (((pa2ta1n-(a1n/2.0_dp))<tiny) .and. ((pa2ta3n-(a3n/2.0_dp))<tiny)) then
                lat_id = 10
             else
                lat_id = 11
             end if
          else
             if (((pa1ta2n-(a2n/2.0_dp))<tiny) .and. ((pa1ta3n-(a3n/2.0_dp))<tiny)) then
                lat_id = 10
             else
                lat_id = 11
             end if
          end if
       else if (((abs(a1n - a2n)<tiny) .and. (abs(a2n - a3n)<tiny)) .or. &
            ((abs(pa3ta1n - pa2ta1n)<tiny) .and. (abs(a3n - a2n)<tiny)) .or. &
            ((abs(pa3ta2n - pa1ta2n)<tiny) .and. (abs(a3n - a1n)<tiny)) .or. &
            ((abs(pa1ta3n - pa2ta3n)<tiny) .and. (abs(a1n - a2n)<tiny))) then
          lat_id =10
       else if ((abs(a1n - a2n)>tiny) .and. (abs(a1n - a3n)>tiny) .and. (abs(a2n - a3n)>tiny)) then
          lat_id = 11
       else 
          stop "Could not identify lattice (orthorhombic) in lat_id."
       end if
    else if (fam == 4) then
       if (pcount==2) then
          lat_id = 12
       else
          lat_id = 13
       end if
    else if (fam == 2) then
       lat_id = 14
    else
       stop "Could not identify lattice in lat_id."
    end if
  end SUBROUTINE identify_lattice

  !!<summary>Get's the canonical basis for the lattice type
  !!provided.</summary>
  !!<parameter name="lat_id" regular="true">The lattice type id from
  !!identify_lattice ('simple cubic' =&gt 1, 'body
  !!centered cubic' =&gt 2, 'face centered cubic' =&gt 3, 'hexagonal'
  !!=&gt 4, 'trigonal' =&gt 5, 'simple tetragonal' =&gt 6, 'body
  !!centered tetragonal' =&gt 7, 'simple orthorhombic' =&gt 8, 'base
  !!centered orthorhombic' =&gt 9, 'body centered orthorhombic'
  !!=&gt10, 'face centered orthorhombic' =&gt 11, 'simplo monoclinic'
  !!=&gt 12, 'base centered monoclinic' =&gt 13, 'triclinic' =&gt
  !!14).</parameter>
  !!<parameter name="lat_vecs" regula="true">The users basis vectors
  !!for the lattice.</parameter>
  !!<parameter name="c_basis" regula="true">The canonical basis for
  !!the lattice.</parameter>
  SUBROUTINE canonical_basis(lat_id,lat_vecs,c_basis)
    integer, intent(in) :: lat_id
    real(dp), intent(in) :: lat_vecs(3,3)
    real(dp), intent(out) :: c_basis(3,3)

    if (lat_id==1) then
       c_basis(:,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       c_basis(:,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       c_basis(:,3) = (/0.0_dp,0.0_dp,1.0_dp/)
    else if (lat_id==2) then
       c_basis(:,1) = (/-1.0_dp,1.0_dp,1.0_dp/)
       c_basis(:,2) = (/1.0_dp,-1.0_dp,1.0_dp/)
       c_basis(:,3) = (/1.0_dp,1.0_dp,-1.0_dp/)
    else if (lat_id==3) then
       c_basis(:,1) = (/0.0_dp,1.0_dp,1.0_dp/)
       c_basis(:,2) = (/1.0_dp,0.0_dp,1.0_dp/)
       c_basis(:,3) = (/1.0_dp,1.0_dp,0.0_dp/)
    else if (lat_id==4) then
       c_basis(:,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       c_basis(:,2) = (/0.5_dp,-0.86602540378443860_dp,0.0_dp/)
       c_basis(:,3) = (/0.0_dp,0.0_dp,2.0_dp/)
    else if (lat_id==5) then
       c_basis(:,1) = (/1.0_dp,2.0_dp,2.0_dp/)
       c_basis(:,2) = (/2.0_dp,1.0_dp,2.0_dp/)
       c_basis(:,3) = (/4.0_dp,3.0_dp,3.0_dp/)
    else if (lat_id==6) then
       c_basis(:,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       c_basis(:,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       c_basis(:,3) = (/0.0_dp,0.0_dp,2.0_dp/)
    else if (lat_id==7) then
       c_basis(:,1) = (/-1.0_dp,1.0_dp,2.0_dp/)
       c_basis(:,2) = (/1.0_dp,-1.0_dp,2.0_dp/)
       c_basis(:,3) = (/1.0_dp,1.0_dp,-2.0_dp/)
    else if (lat_id==8) then
       c_basis(:,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       c_basis(:,2) = (/0.0_dp,2.0_dp,0.0_dp/)
       c_basis(:,3) = (/0.0_dp,0.0_dp,3.0_dp/)
    else if (lat_id==9) then
       c_basis(:,1) = (/0.5_dp,1.0_dp,0.0_dp/)
       c_basis(:,2) = (/0.5_dp,-1.0_dp,0.0_dp/)
       c_basis(:,3) = (/0.0_dp,0.0_dp,3.0_dp/)
    else if (lat_id==10) then
       c_basis(:,1) = (/0.5_dp,1.0_dp,1.5_dp/)
       c_basis(:,2) = (/0.0_dp,2.0_dp,0.0_dp/)
       c_basis(:,3) = (/0.0_dp,0.0_dp,3.0_dp/)
    else if (lat_id==11) then
       c_basis(:,1) = (/0.0_dp,1.0_dp,1.5_dp/)
       c_basis(:,2) = (/0.5_dp,0.0_dp,1.5_dp/)
       c_basis(:,3) = (/0.0_dp,0.0_dp,3.0_dp/)
    else if (lat_id==12) then
       c_basis(:,1) = (/2.0_dp,0.0_dp,0.0_dp/)
       c_basis(:,2) = (/0.0_dp,2.0_dp,0.0_dp/)
       c_basis(:,3) = (/0.5_dp,0.0_dp,2.0_dp/)
    else if (lat_id==13) then
       c_basis(:,1) = (/1.0_dp,1.0_dp,0.0_dp/)
       c_basis(:,2) = (/0.0_dp,2.0_dp,0.0_dp/)
       c_basis(:,3) = (/0.5_dp,0.0_dp,2.0_dp/)
    else if (lat_id==14) then
       c_basis = lat_vecs
    else
       stop "Only the integers 1-14 can be mapped to a canonical basis."
    end if
  end SUBROUTINE canonical_basis  
end Module lat_types
