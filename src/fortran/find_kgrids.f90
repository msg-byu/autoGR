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

  !!<summary>Returns the possible symmetry preserving offsets given
  !!the niggli number.</summary>
  !!<parameter name="lat_id" regular="true">The niggli lattice
  !!id.</parameter>
  !!<parameter name="offsets" regular="true">The symmetry preserving
  !!offsets for this case.</parameter>
  SUBROUTINE get_offsets(lat_id, offsets)
    integer, intent(in) :: lat_id
    real(dp), allocatable, intent(out) :: offsets(:,:)

    if ((lat_id==1) .or. (lat_id==3) .or. (lat_id==5) .or. (lat_id==16) .or. &
         (lat_id==26) .or. (lat_id==11) .or. (lat_id==21) .or. (lat_id==32) .or. &
         (lat_id==33) .or. (lat_id==34) .or. (lat_id==35) .or. (lat_id==10) .or. &
         (lat_id==14) .or. (lat_id==17) .or. (lat_id==27) .or. (lat_id==37) .or. &
         (lat_id==39) .or. (lat_id==41) .or. (lat_id==43) .or. (lat_id==28) .or. &
         (lat_id==29) .or. (lat_id==30) .or. (lat_id==20) .or. (lat_id==25)) then
       allocate(offsets(8,3))
       offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
       offsets(2,:) = (/0.5_dp, 0.0_dp, 0.0_dp/)
       offsets(3,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
       offsets(4,:) = (/0.5_dp, 0.5_dp, 0.0_dp/)
       offsets(5,:) = (/0.0_dp, 0.0_dp, 0.5_dp/)
       offsets(6,:) = (/0.0_dp, 0.5_dp, 0.5_dp/)
       offsets(7,:) = (/0.5_dp, 0.0_dp, 0.5_dp/)
       offsets(8,:) = (/0.5_dp, 0.5_dp, 0.5_dp/)
    elseif ((lat_id==12) .or. (lat_id==9) .or. (lat_id==15) .or. (lat_id==22) .or. &
         (lat_id==24) .or. (lat_id==18) .or. (lat_id==4) .or. (lat_id==2) .or. &
         (lat_id==6) .or. (lat_id==7) .or. (lat_id==19) .or. (lat_id==8) .or. &
         (lat_id==42)) then
       allocate(offsets(4,3))
       offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
       offsets(2,:) = (/0.0_dp, 0.0_dp, 0.5_dp/)
       offsets(3,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
       offsets(4,:) = (/0.5_dp, 0.0_dp, 0.0_dp/)
    elseif ((lat_id==40) .or. (lat_id==38) .or. (lat_id==36) .or. (lat_id==23) .or. &
         (lat_id==13)) then
       allocate(offsets(7,3))
       offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
       offsets(2,:) = (/0.5_dp, 0.0_dp, 0.0_dp/)
       offsets(3,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
       offsets(4,:) = (/0.5_dp, 0.5_dp, 0.0_dp/)
       offsets(5,:) = (/0.0_dp, 0.0_dp, 0.5_dp/)
       offsets(6,:) = (/0.0_dp, 0.5_dp, 0.5_dp/)
       offsets(7,:) = (/0.5_dp, 0.0_dp, 0.5_dp/)
    ! elseif ((lat_id==10) .or. (lat_id==14) .or. (lat_id==17) .or. (lat_id==27) .or. &
    !      (lat_id==37) .or. (lat_id==39) .or. (lat_id==41) .or. (lat_id==43) .or. &
    !      (lat_id==28) .or. (lat_id==29) .or. (lat_id==30)) then
    !    allocate(offsets(8,3))
    !    offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
    !    offsets(2,:) = (/0.0_dp, 0.0_dp, 0.5_dp/)
    !    offsets(3,:) = (/-0.25_dp, 0.25_dp, 0.0_dp/)
    !    offsets(4,:) = (/-0.25_dp, 0.25_dp, 0.5_dp/)
    !    offsets(5,:) = (/0.25_dp, 0.25_dp, 0.0_dp/)
    !    offsets(6,:) = (/0.25_dp, 0.25_dp, 0.5_dp/)
    !    offsets(7,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
    !    offsets(8,:) = (/0.0_dp, 0.5_dp, 0.5_dp/)
    ! elseif ((lat_id==20) .or. (lat_id==25)) then
    !    allocate(offsets(8,3))
    !    offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
    !    offsets(2,:) = (/0.5_dp, 0.0_dp, 0.0_dp/)
    !    offsets(3,:) = (/0.0_dp, 0.25_dp, -0.25_dp/)
    !    offsets(4,:) = (/0.5_dp, 0.25_dp, -0.25_dp/)
    !    offsets(5,:) = (/0.0_dp, 0.25_dp, 0.25_dp/)
    !    offsets(6,:) = (/0.5_dp, 0.25_dp, 0.25_dp/)
    !    offsets(7,:) = (/0.0_dp, 0.5_dp, 0.0_dp/)
    !    offsets(8,:) = (/0.5_dp, 0.5_dp, 0.0_dp/)
    elseif ((lat_id==31) .or. (lat_id==44)) then
       allocate(offsets(8,3))
       offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
       offsets(2,:) = (/0.5_dp, 0.5_dp, 0.5_dp/)       
    else
       write(*,*) "Failed to find offsets for this case. Please report this occurance. Attempting to use defaults."
       allocate(offsets(2,3))
       offsets(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
       offsets(2,:) = (/0.5_dp, 0.5_dp, 0.5_dp/)       
    end if

  end SUBROUTINE get_offsets
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
  !!<parameter name="at" regular="true">The atom types in the
  !!basis.</parameter>
  !!<parameter name="offset" regular="true">The offset for the
  !!k-points grid.</parameter>
  !!<parameter name="find_offset" regular="true">'True' if the offset
  !!needs to be determined by the algorithm.</parameter>
  !!<parameter name="best_offset" regular="true">The best offset from
  !!those checked for the final grid.</parameter>
  SUBROUTINE find_grid(lat_vecs, kpd, B_vecs, at, offset, find_offset, best_grid, &
       best_offset, eps_)
    real(dp), intent(in) :: lat_vecs(3,3), offset(3)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(in) :: kpd
    integer, intent(inout) :: at(:)
    real(dp), optional, intent(in) :: eps_
    real(dp), intent(out) :: best_grid(3,3), best_offset(3)
    logical, intent(in) :: find_offset

    integer :: lat_id, a_kpd, c_kpd(3), i, count, mult
    integer, allocatable :: sp_hnfs(:,:,:), temp_hnfs(:,:,:)
    integer, allocatable :: n_irr_kp(:), nhnfs(:), nt_kpts(:)
    real(dp), allocatable :: grids(:,:,:), ratio(:), offsets(:,:), grid_offsets(:,:)
    real(dp) :: O(3,3), Nu(3,3), No(3,3), temp_grid(3,3)
    integer :: Cu(3,3), Co(3,3), min_kpn_loc(1), temp_nirr, temp_nhnfs, s_range
    real(dp) :: eps
    
    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-3
    end if

    call id_cell(lat_vecs, Nu, Cu, O, No, Co, lat_id, s_range,eps_=eps)
    count = 0

    if (find_offset .eqv. .False.) then
       allocate(offsets(1,3))
       offsets(1,:) = offset
    else
       call get_offsets(lat_id, offsets)
    end if
    if ((lat_id==3) .or. (lat_id==5) .or. (lat_id==1)) then
       call get_kpd_cubic(lat_id,kpd,c_kpd)
       allocate(sp_hnfs(3,3,3), n_irr_kp(3), nhnfs(3), grids(3,3,3), nt_kpts(3))
       allocate(grid_offsets(3,3))
       sp_hnfs = 0
       n_irr_kp = 0
       nhnfs = 0
       grids = 0
       do i=1,3
          a_kpd = c_kpd(i)
          if (lat_id==3) then
             call sc_3(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==5) then
             call bcc_5(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==1) then
             call fcc_1(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          end if
          sp_hnfs(:,:,i) = temp_hnfs(:,:,1)
          grids(:,:,i) = temp_grid
          grid_offsets(i,:) = best_offset
          n_irr_kp(i) = temp_nirr
          nhnfs(i) = temp_nhnfs
          nt_kpts(i) = a_kpd
          count = count + 1
          deallocate(temp_hnfs)
       end do
    else if ((lat_id==44) .or. (lat_id==31)) then
       call get_kpd_tric(kpd, a_kpd, mult)
       allocate(sp_hnfs(3,3,1), n_irr_kp(1), nhnfs(1), grids(3,3,1), nt_kpts(1))
       allocate(grid_offsets(1,3))
       call tric_31_44(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, mult, offsets, &
            best_offset, temp_hnfs, grids(:,:,1), n_irr_kp(1), nhnfs(1), eps_=eps)
       sp_hnfs = temp_hnfs*mult
       nt_kpts(1) = a_kpd*mult
       grid_offsets(1,:) = best_offset
       count = 1
    else
       allocate(sp_hnfs(3,3,s_range), n_irr_kp(s_range), nhnfs(s_range))
       allocate(grids(3,3,s_range), nt_kpts(s_range))
       allocate(grid_offsets(5,3))
       a_kpd = kpd
       sp_hnfs = 0
       n_irr_kp = 0
       nhnfs = 0
       grids = 0
       do count = 1, s_range
       ! do while ((count < s_range) .and. (a_kpd-kpd<=s_range))
          if ((lat_id==2) .or. (lat_id==4)) then
             call rhom_4_2(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==6 .or.lat_id==7 .or. lat_id==15 .or.lat_id==18) then
             call bct_6_7_15_18(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, &
                  offsets, best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, &
                  eps_=eps)
          else if (lat_id==8) then
             call bco_8(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==9 .or. lat_id==24) then
             call rhom_9_24(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if ((lat_id==10) .or. (lat_id==14) .or. (lat_id==17) .or. (lat_id==27) &
               .or. (lat_id==37) .or. (lat_id==39) .or. (lat_id == 41)) then
             call basecm_10_14_17_27_37_39_41(a_kpd, No, Nu, Co, Cu, O, lat_vecs, &
                  B_vecs, at, offsets, best_offset, temp_hnfs, temp_grid, temp_nirr, &
                  temp_nhnfs, eps_=eps)
          else if (lat_id==11) then
             call st_11(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==12) then
             call hex_12(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if ((lat_id==13) .or. (lat_id==38)) then
             call baseco_38_13(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, &
                  offsets, best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, &
                  eps_=eps)
          else if (lat_id==16) then
             call fco_16(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==19) then
             call bco_19(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if ((lat_id==20) .or. (lat_id==25)) then
             call basecm_20_25(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, &
                  offsets, best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, &
                  eps_=eps)
          else if (lat_id==21) then
             call st_21(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==22) then
             call hex_22(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==23) then
             call baseco_23(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==26) then
             call fco_26(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==28) then
             call basecm_28(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if ((lat_id==29) .or. (lat_id==30)) then
             call basecm_29_30(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, &
                  offsets, best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, &
                  eps_=eps)
          else if (lat_id==32) then
             call so_32(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==33) then
             call sm_33(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==34 .or. lat_id==35) then
             call sm_34_35(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==36) then
             call baseco_36(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==40) then
             call baseco_40(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==42) then
             call bco_42(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          else if (lat_id==43) then
             call basecm_43(a_kpd, No, Nu, Co, Cu, O, lat_vecs, B_vecs, at, offsets, &
                  best_offset, temp_hnfs, temp_grid, temp_nirr, temp_nhnfs, eps_=eps)
          end if
          sp_hnfs(:,:,count) = temp_hnfs(:,:,1)
          grids(:,:,count) = temp_grid
          grid_offsets(count,:) = best_offset
          n_irr_kp(count) = temp_nirr
          nhnfs(count) = temp_nhnfs
          nt_kpts(count) = a_kpd
          ! count = count + 1
          a_kpd = a_kpd + 1
          deallocate(temp_hnfs)
       end do
    end if

    allocate(ratio(count))
    do i=1, count
       ratio(i) = real(n_irr_kp(i),dp)/real(nt_kpts(i),dp)
    end do

    min_kpn_loc = MINLOC(ratio)
    best_grid = grids(:,:, min_kpn_loc(1))
    best_offset = grid_offsets(min_kpn_loc(1),:)
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
