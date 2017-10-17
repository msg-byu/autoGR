!!<summary>Subroutines used to find the kpoint grids.</summary>

Module find_kgrids
  use num_types
  use niggli
  use sp_hnfs
  use vector_matrix_utilities, only: matrix_inverse, minkowski_reduce_basis, norm, cross_product
  use numerical_utilities, only: equal

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
    integer, intent(in) :: No(3,3), Nu(3,3), Co(3,3), Cu(3,3), O(3,3)
    integer, intent(out) :: UB(3,3)

    integer :: i
    integer :: L(3,3), F(3,3), B(3,3)
    real(dp) :: Oinv(3,3), Noinv(3,3), Cuinv(3,3)

    call matrix_inverse(O,Oinv)
    call matrix_inverse(No,Noinv)
    call matrix_inverse(Cu,Cuinv)
    L = matmul(matmul(O,spHNF),Oinv)
    F = matmul(Noinv,matmul(matmul(L,O),Co))
    UB = int(matmul(matmul(Nu,F),Cuinv))

  end SUBROUTINE transform_supercells
  
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

    integer :: lat_id, a_kpd, c_kpd(3), i, status, larger, smaller, old, news, j, k
    integer, allocatable :: sp_hnfs(:,:,:), temp_hnfs(:,:,:), temp_hnfs2(:,:,:)
    real(dp) :: O(3,3), Nu(3,3), No(3,3)
    integer :: Cu(3,3), Co(3,3), lat_id
    real(dp) :: eps
    
    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call id_cell(lat_vecsU,Nu,Cu,O,No,Co,lat_id,eps_=eps)

    if ((lat_id==3) .or. (lat_id==5) .or. (lat_id==1)) then
       call get_kpd(lat_id,kpd,c_kpd)
       allocate(sp_hnfs(3,3,3))
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
    else
       larger = 0
       smaller = 0
       a_kpd = kpd
       do while (larger <2)
          if ((lat_id==2) .or. (lat_id==4)) then
             call rhom_4_2(a_kpd,temp_hnfs)
          else if (lat_id==6) then
             call bct_6(a_kpd,temp_hnfs)
          else if (lat_id==7) then
             call bct_7(a_kpd,temp_hnfs)
          else if (lat_id==8) then
             call bco_8(a_kpd,temp_hnfs)
          else if (lat_id==9) then
             call rhom_9(a_kpd,temp_hnfs)
          else if ((lat_id==10) .or. (lat_id==17)) then
             call basecm_10_17(a_kpd,temp_hnfs)
          else if (lat_id==11) then
             call st_11(a_kpd,temp_hnfs)
          else if (lat_id==12) then
             call hex_12(a_kpd,temp_hnfs)
          else if ((lat_id==13) .or. (lat_id==38)) then
             call baseco_38_13(a_kpd,temp_hnfs)
          else if (lat_id==14) then
             call basecm_14(a_kpd,temp_hnfs)
          else if (lat_id==15) then
             call bct_15(a_kpd,temp_hnfs)
          else if (lat_id==16) then
             call fco_16(a_kpd,temp_hnfs)
          else if (lat_id==18) then
             call bct_18(a_kpd,temp_hnfs)
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
          else if (lat_id==24) then
             call rhom_24(a_kpd,temp_hnfs)
          else if (lat_id==26) then
             call fco_26(a_kpd,temp_hnfs)
          else if (lat_id==27) then
             call basecm_27(a_kpd,temp_hnfs)
          else if (lat_id==28) then
             call basecm_28(a_kpd,temp_hnfs)
          else if ((lat_id==21) .or. (lat_id==30)) then
             call basecm_29_30(a_kpd,temp_hnfs)
          else if ((lat_id==31) .or. (lat_id==44)) then
             call tric_31_44(a_kpd,temp_hnfs)
          else if (lat_id==32) then
             call so_32(a_kpd,temp_hnfs)
          else if (lat_id==33) then
             call sm_33(a_kpd,temp_hnfs)
          else if (lat_id==34) then
             call sm_34(a_kpd,temp_hnfs)
          else if (lat_id==35) then
             call sm_35(a_kpd,temp_hnfs)
          else if (lat_id==36) then
             call baseco_36(a_kpd,temp_hnfs)
          else if ((lat_id==37) .or. (lat_id==39)) then
             call basecm_37_39(a_kpd,temp_hnfs)
          else if (lat_id==40) then
             call baseco_40(a_kpd,temp_hnfs)
          else if (lat_id==41) then
             call basecm_41(a_kpd,temp_hnfs)
          else if (lat_id==42) then
             call bco_42(a_kpd,temp_hnfs)
          else if (lat_id==43) then
             call basecm_43(a_kpd,temp_hnfs)
          end if
          if (size(temp_hnfs,3)>1) then
             if (larger==0) then
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
             larger = larger + 1
          else
             a_kpd = a_kpd + 1
          end if
          deallocate(temp_hnfs)
       end do
       
       a_kpd = kpd -1
       
       do while (smaller < 1)
          if ((lat_id==2) .or. (lat_id==4)) then
             call rhom_4_2(a_kpd,temp_hnfs)
          else if (lat_id==6) then
             call bct_6(a_kpd,temp_hnfs)
          else if (lat_id==7) then
             call bct_7(a_kpd,temp_hnfs)
          else if (lat_id==8) then
             call bco_8(a_kpd,temp_hnfs)
          else if (lat_id==9) then
             call rhom_9(a_kpd,temp_hnfs)
          else if ((lat_id==10) .or. (lat_id==17)) then
             call basecm_10_17(a_kpd,temp_hnfs)
          else if (lat_id==11) then
             call st_11(a_kpd,temp_hnfs)
          else if (lat_id==12) then
             call hex_12(a_kpd,temp_hnfs)
          else if ((lat_id==13) .or. (lat_id==38)) then
             call baseco_38_13(a_kpd,temp_hnfs)
          else if (lat_id==14) then
             call basecm_14(a_kpd,temp_hnfs)
          else if (lat_id==15) then
             call bct_15(a_kpd,temp_hnfs)
          else if (lat_id==16) then
             call fco_16(a_kpd,temp_hnfs)
          else if (lat_id==18) then
             call bct_18(a_kpd,temp_hnfs)
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
          else if (lat_id==24) then
             call rhom_24(a_kpd,temp_hnfs)
          else if (lat_id==26) then
             call fco_26(a_kpd,temp_hnfs)
          else if (lat_id==27) then
             call basecm_27(a_kpd,temp_hnfs)
          else if (lat_id==28) then
             call basecm_28(a_kpd,temp_hnfs)
          else if ((lat_id==21) .or. (lat_id==30)) then
             call basecm_29_30(a_kpd,temp_hnfs)
          else if ((lat_id==31) .or. (lat_id==44)) then
             call tric_31_44(a_kpd,temp_hnfs)
          else if (lat_id==32) then
             call so_32(a_kpd,temp_hnfs)
          else if (lat_id==33) then
             call sm_33(a_kpd,temp_hnfs)
          else if (lat_id==34) then
             call sm_34(a_kpd,temp_hnfs)
          else if (lat_id==35) then
             call sm_35(a_kpd,temp_hnfs)
          else if (lat_id==36) then
             call baseco_36(a_kpd,temp_hnfs)
          else if ((lat_id==37) .or. (lat_id==39)) then
             call basecm_37_39(a_kpd,temp_hnfs)
          else if (lat_id==40) then
             call baseco_40(a_kpd,temp_hnfs)
          else if (lat_id==41) then
             call basecm_41(a_kpd,temp_hnfs)
          else if (lat_id==42) then
             call bco_42(a_kpd,temp_hnfs)
          else if (lat_id==43) then
             call basecm_43(a_kpd,temp_hnfs)
          end if
          if (size(temp_hnfs,3)>1) then
             if (larger==0) then
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
             smaller = smaller + 1
          else
             a_kpd = a_kpd - 1
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
  SUBROUTINE get_kpd(lat_id,kpd,densities)
    integer, intent(in) :: lat_id, kpd
    integer, intent(out) :: densities(3)

    integer :: i, j, a, b, c, temp1, temp2
    integer :: mults(3)

    temp1= int(real(kpd,dp)**(1.0_dp/3.0_dp))
    a = 0
    b = temp1**3
    c = 0
    
    if (lat_id==3) then
       mults = (/1,4,16/)
    else
       mults = (/1,2,4/)
    end if

    do i=1,temp1
       do j=1,3
          temp2 = mults(j)*(i**3)
          if (temp2 /= b) then
             if (temp2<kpd) then
                if ((a==0) .or. ((kpd-a)>(kpd-temp2))) then
                   a = int(temp2)
                end if
             else
                if ((c==0) .or. ((c-kpd)>(temp2-kpd))) then
                   c = int(temp2)
                end if
             end if
          end if
       end do
    end do

    densities = (/a,b,c/)

  end SUBROUTINE get_kpd

  !!<summary>Selects the best grid from the list of grids.</summary>
  !!<parameter name="grids" regular="true">The list of generating
  !!vectors for the candidate grids.</parameter>
  !!<parameter name="best_grid" regular="true">The best grid given the
  !!criteria.</parameter>
  SUBROUTINE grid_selection(grids, best_grid)
    real(dp), allocatable, intent(in) :: grids(:,:,:)
    real(dp), intent(out) :: best_grid(3,3)

    real(dp) :: reduced_grid(3,3), norms(3)
    integer :: i, j
    real(dp) :: r_min, pf, r_min_best, pf_best, pi, eps
    

    pi = 3.14159265358979323846
    eps = 1E-6
    
    r_min_best = 0
    pf_best = 0
    
    do i=1,size(grids,3)
       call minkowski_reduce_basis(grids(:,:,i),reduced_grid,eps)
       r_min = max(norms(1),norms(2),norms(3))
       pf = (4.0_dp/3.0_dp)*pi*(min(norms(1),norms(2),norms(3))**2)/(dot_product(reduced_grid(:,1),cross_product(reduced_grid(:,2),reduced_grid(:,3))))
       if (r_min_best==0) then
          r_min_best = r_min
          pf_best = pf
          best_grid = grids(:,:,i)
       else
          if (r_min<r_min_best) then
             r_min_best = r_min
             pf_best = pf
             best_grid = grids(:,:,i)
          elseif (abs(r_min-r_min_best)<1E-3) then
             if (pf <pf_best) then
                r_min_best = r_min
                pf_best = pf
                best_grid = grids(:,:,i)
             end if
          end if
       end if
    end do

  end SUBROUTINE grid_selection
end Module find_kgrids
