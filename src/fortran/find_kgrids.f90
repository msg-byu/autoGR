!!<summary>Subroutines used to find the kpoint grids.</summary>

Module find_kgrids
  use num_types
  use niggli
  use sp_hnfs
  use kpointgeneration, only: generateIrredKpointList

  use symmetry
  use vector_matrix_utilities
  use grid_utils, only: transform_supercell

  implicit none
  private
  public find_grid

CONTAINS
  
  !!<summary>Determines the best k-point grid to use near the target
  !!density.</summary>
  !!<parameter name="lat_vecs" regular="true">The parent cell lattice
  !!vectors.</parameter>
  !!<parameter name="kpd" regular="true">The target kpoint
  !!density.</parameter>
  !!<parameter name="best_grid" regular="true">The best k-point grid
  !!found.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="True">The atom types in the
  !!basis.</parameter>
  SUBROUTINE find_grid(lat_vecs, kpd, B_vecs, at, best_grid, eps_)
    real(dp), intent(in) :: lat_vecs(3,3)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(in) :: kpd
    integer, intent(inout) :: at(:)
    real(dp), optional, intent(in) :: eps_
    real(dp), intent(out) :: best_grid(3,3)

    integer :: lat_id, a_kpd, c_kpd(3), i, count, mult
    integer, allocatable :: sp_hnfs(:,:,:), temp_hnfs(:,:,:)
    integer, allocatable :: n_irr_kp(:), nhnfs(:), nt_kpts(:), ratio(:)
    real(dp), allocatable :: grids(:,:,:), rmin(:)
    real(dp) :: O(3,3), Nu(3,3), No(3,3), temp_grid(3,3), temp_rmin
    integer :: Cu(3,3), Co(3,3), min_kpn_loc(1), temp_nirr, temp_nhnfs
    real(dp) :: eps

    real(dp), pointer :: O_pg_ops(:,:,:), U_pg_ops(:,:,:)
    integer :: j, k, z
    real(dp) :: N_u(3,3), N_o(3,3), Hinv(3,3), Uinv(3,3), Oinv(3,3), X(3,3), S(3,3), H(3,3), St(3,3)
    

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call id_cell(lat_vecs,Nu,Cu,O,No,Co,lat_id,eps_=eps)
    count = 0

    call get_lattice_pointGroup(lat_vecs, U_pg_ops)
    call get_lattice_pointGroup(O, O_pg_ops, eps_=real(1E-3,dp))
    ! print *, "Niggli N", lat_id
    ! print *, "User group order", size(U_pg_ops,3)
    ! print *, "Our group order", size(O_pg_ops,3)
    
    if ((lat_id==3) .or. (lat_id==5) .or. (lat_id==1)) then
       call get_kpd_cubic(lat_id,kpd,c_kpd)
       allocate(sp_hnfs(3,3,3), n_irr_kp(3), rmin(3), nhnfs(3), grids(3,3,3), nt_kpts(3))
       sp_hnfs = 0
       n_irr_kp = 0
       rmin = 0
       nhnfs = 0
       grids = 0
       do i=1,3
          a_kpd = c_kpd(i)
          if (lat_id==3) then
             call sc_3(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==5) then
             call bcc_5(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==1) then
             call fcc_1(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          end if
          sp_hnfs(:,:,i) = temp_hnfs(:,:,1)
          grids(:,:,i) = temp_grid
          rmin(i) = temp_rmin
          n_irr_kp(i) = temp_nirr
          nhnfs(i) = temp_nhnfs
          nt_kpts(i) = a_kpd
          count = count + 1
          deallocate(temp_hnfs)
       end do
    else if ((lat_id==44) .or. (lat_id==31)) then
       call get_kpd_tric(kpd, a_kpd, mult)
       allocate(sp_hnfs(3,3,1), n_irr_kp(1), rmin(1), nhnfs(1), grids(3,3,1), nt_kpts(1))
       call tric_31_44(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, mult, temp_hnfs, &
                  grids(:,:,1), rmin(1), n_irr_kp(1), nhnfs(1), eps_=eps)
       sp_hnfs = temp_hnfs*mult
       nt_kpts(1) = a_kpd*mult
       count = 1
    else
       allocate(sp_hnfs(3,3,5), n_irr_kp(5), rmin(5), nhnfs(5), grids(3,3,5), nt_kpts(5))
       a_kpd = kpd
       sp_hnfs = 0
       n_irr_kp = 0
       rmin = 0
       nhnfs = 0
       grids = 0
       do while ((count <5) .and. (a_kpd-kpd<=10))
          if ((lat_id==2) .or. (lat_id==4)) then
             call rhom_4_2(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==6 .or.lat_id==7 .or. lat_id==15 .or.lat_id==18) then
             call bct_6_7_15_18(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==8) then
             call bco_8(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==9 .or. lat_id==24) then
             call rhom_9_24(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if ((lat_id==10) .or. (lat_id==14) .or. (lat_id==17) .or. (lat_id==27) &
               .or. (lat_id==37) .or. (lat_id==39) .or. (lat_id == 41)) then
             call basecm_10_14_17_27_37_39_41(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==11) then
             call st_11(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==12) then
             call hex_12(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if ((lat_id==13) .or. (lat_id==38)) then
             call baseco_38_13(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==16) then
             call fco_16(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==19) then
             call bco_19(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if ((lat_id==20) .or. (lat_id==25)) then
             call basecm_20_25(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==21) then
             call st_21(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==22) then
             call hex_22(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==23) then
             call baseco_23(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==26) then
             call fco_26(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==28) then
             call basecm_28(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if ((lat_id==29) .or. (lat_id==30)) then
             call basecm_29_30(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==32) then
             call so_32(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==33) then
             call sm_33(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==34 .or. lat_id==35) then
             call sm_34_35(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==36) then
             call baseco_36(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==40) then
             call baseco_40(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==42) then
             call bco_42(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==43) then
             call basecm_43(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, temp_hnfs, &
                  temp_grid, temp_rmin, temp_nirr, temp_nhnfs, eps_=eps)
          end if
          if (temp_nhnfs>1) then
             sp_hnfs(:,:,count+1) = temp_hnfs(:,:,1)
             grids(:,:,count+1) = temp_grid
             rmin(count+1) = temp_rmin
             n_irr_kp(count+1) = temp_nirr
             nhnfs(count+1) = temp_nhnfs
             nt_kpts(count+1) = a_kpd
             count = count + 1
             a_kpd = a_kpd + 1
             ! print *, "temp_nirr", temp_nirr
             ! print *, "n_irr_kp", n_irr_kp
             
          else
             a_kpd = a_kpd + 1
          end if
          deallocate(temp_hnfs)
       end do
    end if

    allocate(ratio(count))
    do i=1,size(ratio,1)
       ratio(i) = real(n_irr_kp(i),dp)/real(nt_kpts(i),dp)
    end do
    
    ! call matrix_inverse(lat_vecs, Uinv)
    ! call matrix_inverse(O, Oinv)
    
    ! print *, "Checking the HNFs generated to ensure they preserve the symmetries of our basis."
    ! do i = 1, count
    !    call matrix_inverse(real(sp_hnfs(:,:,i),dp), Hinv)
    !    do k=1, size(O_pg_ops,3)
    !       X = matmul(Oinv, matmul(O_pg_ops(:,:,k),O))
    !       N_u = matmul(Hinv, matmul(X, sp_hnfs(:,:,i)))
    !       if (any(abs(N_u-nint(N_u))>1E-3) .or. any(abs(X-nint(X))>1E-3)) then
    !          print *, "Failed to preserve our symmetry."
    !          print *, "i", i, "k", k
    !          print *, "HNF"
    !          do z=1,3
    !             print *, sp_hnfs(z,:,i)
    !          end do
    !          print *, "N"
    !          do z=1,3
    !             print *, N_u(:,z)
    !          end do
    !          print *, "X"
    !          do z=1,3
    !             print *, X(:,z)
    !          end do
    !       end if
    !    end do
    ! end do
 
    
    ! print *, "Checking the grids generated to ensure they preserve the symmetries of the users basis."
    ! do i = 1, count
    !    call transform_supercell(sp_hnfs(:,:,i), No, Nu, Co, Cu, O, S)
    !    H = matmul(Uinv,S)
    !    call matrix_inverse(H, Hinv)
    !    do k=1, size(U_pg_ops, 3)
    !       X = matmul(Uinv, matmul(U_pg_ops(:,:,k), lat_vecs))
    !       N_u = matmul(Hinv, matmul(X, H))
    !       if (any(abs(N_u-nint(N_u))>1E-3) .or. any(abs(X-nint(X))>1E-3)) then
    !          print *, "Failed to user symmetry."
    !          print *, "i", i, "k", k
    !          print *, "grid"
    !          do z=1,3
    !             print *, grids(z,:,i)
    !          end do
    !          ! print *, "S"
    !          ! do z=1,3
    !          !    print *, S(z,:)
    !          ! end do
    !          print *, "Uinv"
    !          do z=1,3
    !             print *, Uinv(z,:)
    !          end do
    !          print *, "HNF"
    !          do z=1,3
    !             print *, H(z,:)
    !          end do
    !          print *, "N"
    !          do z=1,3
    !             print *, N_u(z,:)
    !          end do
    !          print *, "X"
    !          do z=1,3
    !             print *, X(z,:)
    !          end do
    !       end if
    !    end do
    ! end do

    ! print *, "check if grid and supercell match."
    ! do j=1, count
    !    call transform_supercell(sp_hnfs(:,:,j), No, Nu, Co, Cu, O, S)
    !    call matrix_inverse(transpose(S),St)
    !    call get_lattice_pointGroup(St,U_pg_ops)
    !    if (any(abs(grids(:,:,j)-St)>1E-3)) then
    !       print *, "grids and supercell don't match."
    !       print *, "j",j, size(U_pg_ops,3)
    !       print *, "grids"
    !       do i=1,3
    !          print *, grids(i,:,j)
    !       end do
    !       print *, "S"
    !       do i=1,3
    !          print *, S(i,:)
    !       end do
    !       print *, "St"
    !       do i=1,3
    !          print *, St(i,:)
    !       end do
    !    end if
    ! end do
    ! do i=1,count
    !    call get_lattice_pointGroup(grids(:,:,i),U_pg_ops,eps_=real(1E-3,dp))
    !    print *, "i", i, size(U_pg_ops,3)
    ! end do
    ! ! Select the grid that has the fewest irreducible k-points as the
    ! ! best grid for this system.
    ! print *, n_irr_kp
    min_kpn_loc = MINLOC(ratio)
    best_grid = grids(:,:, min_kpn_loc(1))
    ! print *, "selected", min_kpn_loc
    ! print *, "count", count
  end SUBROUTINE find_grid

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

    nb = ((int(real(kpd,dp)**(1.0_dp/3.0_dp)))/16)-1
    nmax  = nb+500
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

end Module find_kgrids
