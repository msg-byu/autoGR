!!<summary>Subroutines used to find the kpoint grids.</summary>

Module find_kgrids
  use num_types
  use niggli
  use sp_hnfs
  use vector_matrix_utilities, only: matrix_inverse, minkowski_reduce_basis, norm, cross_product
  use numerical_utilities, only: equal
  use kpointgeneration, only: generateIrredKpointList
  use symmetry, only : get_lattice_pointGroup

  implicit none
  private
  public find_grids, grid_selection

CONTAINS

  !!<summary>Finds the symmetry preserving supercell in the users
  !!basis.</summary>
  !!<parameter name="spHNF" regular="true">One of the symmetry
  !!preserving HNFs in our basis</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="UB" regular="true">The supercell in the users
  !!basis.</parameter>
  SUBROUTINE transform_supercell(spHNF,No,Nu,Co,Cu,O,UB)
    integer, intent(in) :: spHNF(3,3)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3)
    real(dp), intent(out) :: UB(3,3)

    integer :: F(3,3)
    real(dp) :: Oinv(3,3), Noinv(3,3), Cuinv(3,3), L(3,3)

    call matrix_inverse(O,Oinv)
    call matrix_inverse(No,Noinv)
    call matrix_inverse(real(Cu,dp),Cuinv)
    L = matmul(matmul(O,spHNF),Oinv)
    F = nint(matmul(Noinv,matmul(matmul(L,O),Co)))
    UB = matmul(matmul(Nu,F),Cuinv)

  end SUBROUTINE transform_supercell

  !!<summary>Determines the symmetry preserving kpoint grids, with
  !!target denstiy, for the provided lattice.</summary>
  !!<parameter name="lat_vecs" regular="true">The parent cell lattice
  !!vectors.</parameter>
  !!<parameter name="kpd" regular="true">The target kpoint
  !!density.</parameter>
  !!<parameter name="grids" regular="true">The kpoint grids
  !!found.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  SUBROUTINE find_grids(lat_vecs,kpd,grids,eps_)
    real(dp), intent(in) :: lat_vecs(3,3)
    integer, intent(in) :: kpd
    real(dp), optional, intent(in) :: eps_
    real(dp), allocatable, intent(out) :: grids(:,:,:)

    integer :: lat_id, a_kpd, c_kpd(3), i, status, count, old, news
    integer, allocatable :: sp_hnfs(:,:,:), temp_hnfs(:,:,:), temp_hnfs2(:,:,:), n_irr_kp(:)
    real(dp) :: O(3,3), Nu(3,3), No(3,3), UB(3,3)
    integer :: Cu(3,3), Co(3,3), mult
    real(dp) :: eps

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-5
    end if

    call id_cell(lat_vecs,Nu,Cu,O,No,Co,lat_id,eps_=eps)

    if ((lat_id==3) .or. (lat_id==5) .or. (lat_id==1)) then
       call get_kpd_cubic(lat_id,kpd,c_kpd)
       allocate(sp_hnfs(3,3,3), n_irr_kp(3))
       do i=1,3
          a_kpd = c_kpd(i)
          if (lat_id==3) then
             call sc_3(a_kpd,temp_hnfs)
          else if (lat_id==5) then
             call bcc_5(a_kpd,temp_hnfs)
          else if (lat_id==1) then
             call fcc_1(a_kpd,temp_hnfs)
          end if
          sp_hnfs(:,:,i) = temp_hnfs(:,:,1)
          deallocate(temp_hnfs)
       end do
    else if ((lat_id==44) .or. (lat_id==31)) then
       call get_kpd_tric(kpd,a_kpd,mult)
       call tric_31_44(a_kpd,sp_hnfs)
       sp_hnfs = sp_hnfs*mult
    else
       count = 0
       a_kpd = kpd
       do while ((count <5) .and. (a_kpd-kpd<=10))
          if ((lat_id==2) .or. (lat_id==4)) then
             call rhom_4_2(a_kpd,temp_hnfs)
          else if (lat_id==6 .or.lat_id==7 .or. lat_id==15 .or.lat_id==18) then
             call bct_6_7_15_18(a_kpd,temp_hnfs)
          else if (lat_id==8) then
             call bco_8(a_kpd,temp_hnfs)
          else if (lat_id==9 .or. lat_id==24) then
             call rhom_9_24(a_kpd,temp_hnfs)
          else if ((lat_id==10) .or. (lat_id==14) .or. (lat_id==17) .or. (lat_id==27) &
                  .or. (lat_id==37) .or. (lat_id==39) .or. (lat_id == 41)) then
             call basecm_10_14_17_27_37_39_41(a_kpd,temp_hnfs)
          else if (lat_id==11) then
             call st_11(a_kpd,temp_hnfs)
          else if (lat_id==12) then
             call hex_12(a_kpd,temp_hnfs)
          else if ((lat_id==13) .or. (lat_id==38)) then
             call baseco_38_13(a_kpd,temp_hnfs)
          else if (lat_id==16) then
             call fco_16(a_kpd,temp_hnfs)
          else if (lat_id==19) then
             call bco_19(a_kpd,temp_hnfs)
          else if ((lat_id==20) .or. (lat_id==25)) then
             call basecm_20_25(a_kpd,temp_hnfs)
          else if (lat_id==21) then
             call st_21(a_kpd,temp_hnfs)
          else if (lat_id==22) then
             call hex_22(a_kpd,temp_hnfs)
          else if (lat_id==23) then
             call baseco_23(a_kpd,temp_hnfs)
          else if (lat_id==26) then
             call fco_26(a_kpd,temp_hnfs)
          else if (lat_id==28) then
             call basecm_28(a_kpd,temp_hnfs)
          else if ((lat_id==29) .or. (lat_id==30)) then
             call basecm_29_30(a_kpd,temp_hnfs)
          else if (lat_id==32) then
             call so_32(a_kpd,temp_hnfs)
          else if (lat_id==33) then
             call sm_33(a_kpd,temp_hnfs)
          else if (lat_id==34 .or. lat_id==35) then
             call sm_34_35(a_kpd,temp_hnfs)
          else if (lat_id==36) then
             call baseco_36(a_kpd,temp_hnfs)
          else if (lat_id==40) then
             call baseco_40(a_kpd,temp_hnfs)
          else if (lat_id==42) then
             call bco_42(a_kpd,temp_hnfs)
          else if (lat_id==43) then
             call basecm_43(a_kpd,temp_hnfs)
          end if
          if (size(temp_hnfs,3)>1) then
             if (count==0) then
                allocate(sp_hnfs(3,3,size(temp_hnfs,3)))
                sp_hnfs(:,:,1:size(temp_hnfs,3)) = temp_hnfs(:,:,1:size(temp_hnfs,3))
             else
                old = size(sp_hnfs,3)
                news = old + 1
                allocate(temp_hnfs2(3,3,old))
                temp_hnfs2(:,:,1:old) = sp_hnfs(:,:,1:old)
                deallocate(sp_hnfs)
                allocate(sp_hnfs(3,3,old+size(temp_hnfs,3)))
                sp_hnfs(:,:,1:old) = temp_hnfs2(:,:,1:old)
                sp_hnfs(:,:,news:size(temp_hnfs,3)+old) = temp_hnfs(:,:,1:size(temp_hnfs,3))
                deallocate(temp_hnfs2)
             end if
             count = count + 1
             a_kpd = a_kpd + 1
          else
             a_kpd = a_kpd + 1
          end if
          deallocate(temp_hnfs)
       end do
    end if

    allocate(grids(3,3,size(sp_hnfs,3)),STAT=status)
    if (status /=0) stop "Failed to allocate memory in find_kgrids."

    do i=1,size(sp_hnfs,3)
       call transform_supercell(sp_HNFs(:,:,i),No,Nu,Co,Cu,O,UB)
       grids(:,:,i) = transpose(UB)
       call matrix_inverse(grids(:,:,i),grids(:,:,i))
    end do
    deallocate(sp_hnfs)
  end SUBROUTINE find_grids

  !!<summary>Gets the trial range of kpoint densities for cubic
  !!lattices.</summary>
  !!<parameter name="lat_id" regular="true">The integer lattice
  !!identifier.</parameter>
  !!<parameter name="kpd" regular="true">The target kpoint
  !!density.</parameter>
  !!<parameter name="densities" regular="true">The corrected
  !!kpoint density.</parameter>
  SUBROUTINE get_kpd_cubic(lat_id,kpd,densities)
    integer, intent(in) :: lat_id, kpd
    integer, intent(out) :: densities(3)

    integer :: j, a, b, c , nb, nmax, nc, temp
    integer :: mults(3)

    nb = int(real(kpd,dp)**(1.0_dp/3.0_dp))-1
    nmax  = kpd+1000
    a = 0
    b = 0
    c = 0

    if (lat_id==1) then
       mults = (/1,4,16/)
    else
       mults = (/1,2,4/)
    end if

    do while (nb<nmax)
       nc = nb**3
       do j=1,3
          temp = mults(j)*nc
          if (temp >= kpd) then
             if (((temp-kpd) < abs(a-kpd)) .or. (a==0)) then
                c = b
                b = a
                a = temp
             else if (((temp-kpd) < abs(b-kpd)) .or. (b==0)) then
                c = b
                b = temp
             else if (((temp-kpd) < abs(c-kpd)) .or. (c==0)) then
                c = temp
             end if
          end if
       end do
       nb = nb + 1
    end do

    densities = (/a,b,c/)

  end SUBROUTINE get_kpd_cubic

  !!<summary>Finds the k-point density for triclinic crystals.</summary>
  !!<parameter name="kpd" regular="true">The target k-point density.</parameter>
  !!<parameter name="r_kpd" regular="true">The returned k-point
  !!density.</parameter>
  !!<parameter name="mult" regular="true">The multiplication factor
  !!between kpd and r_kpd.</parameter>
  SUBROUTINE get_kpd_tric(kpd,r_kpd,mult)
    integer, intent(in) :: kpd
    integer, intent(out) :: r_kpd, mult

    integer :: tmult, test

    mult = 0

    tmult = 1

    do while (mult==0)
       test = kpd/(tmult**3)
       if (test <= 100) then
          mult = tmult
          r_kpd = test
       else
          tmult = tmult + 1
       end if
    end do

  end SUBROUTINE get_kpd_tric

  !!<summary>Selects the best grid from the list of grids.</summary>
  !!<parameter name="lat_vecs" regular="true">The lattice vectors for
  !!the crystal.</parameter>
  !!<parameter name="grids" regular="true">The list of generating
  !!vectors for the candidate grids.</parameter>
  !!<parameter name="best_grid" regular="true">The best grid given the
  !!criteria.</parameter>
  !!<parameter name="shift" regular="true">The shift off gamma to be
  !!used.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  SUBROUTINE grid_selection(lat_vecs,B_vecs,at, grids, shift, best_grid,eps_)
    real(dp), intent(in) :: lat_vecs(3,3), shift(3)
    real(dp), allocatable, intent(in) :: grids(:,:,:)
    real(dp), optional, intent(in) :: eps_
    real(dp), intent(out) :: best_grid(3,3)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)

    real(dp) :: reduced_grid(3,3), norms(3)
    integer :: i, n_irreducible
    real(dp) :: r_min, r_min_best, eps

    real(dp)              :: R(3,3), invLat(3,3)
    real(dp), pointer     :: rdKlist(:,:)
    integer, pointer      :: weights(:)


    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call matrix_inverse(lat_vecs, invLat)
    R = transpose(invLat)

    r_min_best = 0
    n_irreducible = 0
    do i=1,size(grids,3)
       call minkowski_reduce_basis(grids(:,:,i),reduced_grid,eps)
       norms(1) = sqrt(dot_product(reduced_grid(:,1),reduced_grid(:,1)))
       norms(2) = sqrt(dot_product(reduced_grid(:,2),reduced_grid(:,2)))
       norms(3) = sqrt(dot_product(reduced_grid(:,3),reduced_grid(:,3)))
       r_min = min(norms(1),norms(2),norms(3))
       if ((r_min_best==0) .or. (r_min>r_min_best)) then
          r_min_best = r_min
          best_grid = grids(:,:,i)
          call generateIrredKpointList(lat_vecs,B_vecs,at,best_grid, R, shift, rdKlist, weights, eps_=eps)
          n_irreducible = size(rdKlist,1)
       else if (abs(r_min-r_min_best)<eps) then
          call generateIrredKpointList(lat_vecs,B_vecs,at,best_grid, R, shift, rdKlist, weights, eps_=eps)
          if (size(rdKlist,1) < n_irreducible) then
             r_min_best = r_min
             best_grid = grids(:,:,i)
             n_irreducible = size(rdKlist,1)
          end if
       end if
    end do

  end SUBROUTINE grid_selection
end Module find_kgrids
