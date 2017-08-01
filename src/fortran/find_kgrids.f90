!!<summary>Subroutines used to find the kpoint grids.</summary>

Module find_kgrids
  use num_types
  use lat_types
  use sp_hnfs
  use vector_matrix_utilities, only: matrix_inverse, minkowski_reduce_basis, norm, cross_product
  use numerical_utilities, only: equal
  use kpointgeneration, only: generateIrredKpointList

  implicit none
  private
  public find_grids, grid_selection

CONTAINS

  !!<summary>Determines the symmetry preserving kpoint grids, with
  !!target denstiy, for the provided lattice.</summary>
  !!<parameter name="lat_vecs" regular="true">The parent cell lattice
  !!vectors.</parameter>
  !!<parameter name="kpd" regular="true">The target kpoint
  !!density.</parameter>
  !!<parameter name="grids" regular="true">The kpoint grids
  !!found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets allowed for
  !!this system.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  SUBROUTINE find_grids(lat_vecs,kpd,grids,offsets,eps_)
    real(dp), intent(in) :: lat_vecs(3,3)
    integer, intent(in) :: kpd
    real(dp), optional, intent(in) :: eps_
    real(dp), allocatable, intent(out) :: grids(:,:,:), offsets(:,:)

    integer :: lat_id, a_kpd, c_kpd(3), i, status, larger, smaller, old, news, j, k, mult
    integer, allocatable :: sp_hnfs(:,:,:), temp_hnfs(:,:,:), temp_hnfs2(:,:,:)
    real(dp) :: c_basis(3,3), M(3,3), pLVinv(3,3), Minv(3,3), temp(3,3)
    real(dp) :: eps
    
    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call identify_lattice(lat_vecs,lat_id,eps)
    call canonical_basis(lat_id,lat_vecs,c_basis)

    call matrix_inverse(lat_vecs,pLVinv)
    M = matmul(c_basis,pLVinv)
    call matrix_inverse(M,Minv)

    if ((lat_id==1) .or. (lat_id==2) .or. (lat_id==3)) then
       call get_kpd(lat_id,kpd,c_kpd)
       allocate(sp_hnfs(3,3,3))
       do i=1,3
          a_kpd = c_kpd(i)
          if (lat_id==1) then
             call sc(a_kpd,temp_hnfs)
          else if (lat_id==2) then
             call bcc(a_kpd,temp_hnfs)
          else if (lat_id==3) then
             call fcc(a_kpd,temp_hnfs)
          end if
          sp_hnfs(:,:,i) = temp_hnfs(:,:,1)
          deallocate(temp_hnfs)
       end do

       if (lat_id ==2) then
          allocate(offsets(1,3))
          offsets(1,:) = (/0.0_dp,0.0_dp,0.0_dp/)
       else 
          allocate(offsets(2,3))
          offsets(1,:) = (/0.5_dp,0.5_dp,0.5_dp/)
          offsets(2,:) = (/0.0_dp,0.0_dp,0.0_dp/)
       end if
    else if (lat_id /=14) then
       larger = 0
       smaller = 0
       a_kpd = kpd
       do while (larger <2)
          if (lat_id==4) then
             call hex(a_kpd,temp_hnfs)
          else if (lat_id==5) then
             call trig(a_kpd,temp_hnfs)
          else if (lat_id==6) then
             call st(a_kpd,temp_hnfs)
          else if (lat_id==7) then
             call bct(a_kpd,temp_hnfs)
          else if (lat_id==8) then
             call so(a_kpd,temp_hnfs)
          else if (lat_id==9) then
             call baseco(a_kpd,temp_hnfs)
          else if (lat_id==10) then
             call bco(a_kpd,temp_hnfs)
          else if (lat_id==11) then
             call fco(a_kpd,temp_hnfs)
          else if (lat_id==12) then
             call sm(a_kpd,temp_hnfs)
          else if (lat_id==13) then
             call basecm(a_kpd,temp_hnfs)
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
          if (lat_id==4) then
             call hex(a_kpd,temp_hnfs)
          else if (lat_id==5) then
             call trig(a_kpd,temp_hnfs)
          else if (lat_id==6) then
             call st(a_kpd,temp_hnfs)
          else if (lat_id==7) then
             call bct(a_kpd,temp_hnfs)
          else if (lat_id==8) then
             call so(a_kpd,temp_hnfs)
          else if (lat_id==9) then
             call baseco(a_kpd,temp_hnfs)
          else if (lat_id==10) then
             call bco(a_kpd,temp_hnfs)
          else if (lat_id==11) then
             call fco(a_kpd,temp_hnfs)
          else if (lat_id==12) then
             call sm(a_kpd,temp_hnfs)
          else if (lat_id==13) then
             call basecm(a_kpd,temp_hnfs)
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

       if ((lat_id==4) .or. (lat_id==5) .or. (lat_id==7)) then
          allocate(offsets(2,3))
          offsets(1,:) = (/0.0_dp,0.0_dp,0.5_dp/)
          offsets(2,:) = (/0.0_dp,0.0_dp,0.0_dp/)
       else if (lat_id==6) then
          allocate(offsets(4,3))
          offsets(1,:) = (/0.5_dp,0.5_dp,0.5_dp/)
          offsets(2,:) = (/0.0_dp,0.0_dp,0.0_dp/)
          offsets(3,:) = (/0.0_dp,0.0_dp,0.5_dp/)
          offsets(4,:) = (/0.5_dp,0.5_dp,0.0_dp/)
       else if (lat_id==11) then
          allocate(offsets(2,3))
          offsets(1,:) = (/0.5_dp,0.5_dp,0.5_dp/)
          offsets(2,:) = (/0.0_dp,0.0_dp,0.0_dp/)
       else if (lat_id==10) then
          allocate(offsets(4,3))
          offsets(1,:) = (/0.5_dp,0.0_dp,0.0_dp/)
          offsets(2,:) = (/0.0_dp,0.0_dp,0.0_dp/)
          offsets(3,:) = (/0.0_dp,0.5_dp,0.0_dp/)
          offsets(4,:) = (/0.0_dp,0.0_dp,0.5_dp/)
       else if (lat_id==9) then
          allocate(offsets(4,3))
          offsets(1,:) = (/0.0_dp,0.5_dp,0.5_dp/)
          offsets(2,:) = (/0.0_dp,0.0_dp,0.0_dp/)
          offsets(3,:) = (/0.0_dp,0.5_dp,0.0_dp/)
          offsets(4,:) = (/0.0_dp,0.0_dp,0.5_dp/)
       else if (lat_id==13) then
          allocate(offsets(8,3))
          offsets(1,:) = (/-0.25_dp,0.25_dp,0.5_dp/)
          offsets(2,:) = (/0.0_dp,0.0_dp,0.0_dp/)
          offsets(3,:) = (/0.0_dp,0.0_dp,0.5_dp/)
          offsets(4,:) = (/-0.25_dp,0.25_dp,0.0_dp/)
          offsets(5,:) = (/0.0_dp,0.5_dp,0.0_dp/)
          offsets(6,:) = (/0.25_dp,0.25_dp,0.0_dp/)
          offsets(7,:) = (/0.25_dp,0.25_dp,0.5_dp/)
          offsets(8,:) = (/0.0_dp,0.5_dp,0.5_dp/)
       else
          allocate(offsets(8,3))
          offsets(1,:) = (/0.5_dp,0.5_dp,0.5_dp/)
          offsets(2,:) = (/0.0_dp,0.0_dp,0.0_dp/)
          offsets(3,:) = (/0.0_dp,0.5_dp,0.0_dp/)
          offsets(4,:) = (/0.0_dp,0.0_dp,0.5_dp/)
          offsets(5,:) = (/0.5_dp,0.0_dp,0.0_dp/)
          offsets(6,:) = (/0.5_dp,0.5_dp,0.0_dp/)
          offsets(7,:) = (/0.5_dp,0.0_dp,0.5_dp/)
          offsets(8,:) = (/0.0_dp,0.5_dp,0.5_dp/)
       end if
    else

       mult = 1
       do while ((kpd/(mult**3)) > 500)
          mult = mult + 1
       end do

       a_kpd = kpd/(mult**3)
       
       call tric(a_kpd,temp_hnfs)
       
       allocate(sp_hnfs(3,3,size(temp_hnfs,3)),STAT= status)
       if (status /=0) stop "Failed to allocate memory in find_kgrids."

       sp_hnfs = temp_hnfs*mult
       deallocate(temp_hnfs)
       
       allocate(offsets(8,3))
       offsets(1,:) = (/0.5_dp,0.5_dp,0.5_dp/)
       offsets(2,:) = (/0.0_dp,0.0_dp,0.0_dp/)
       offsets(3,:) = (/0.0_dp,0.5_dp,0.0_dp/)
       offsets(4,:) = (/0.0_dp,0.0_dp,0.5_dp/)
       offsets(5,:) = (/0.5_dp,0.0_dp,0.0_dp/)
       offsets(6,:) = (/0.5_dp,0.5_dp,0.0_dp/)
       offsets(7,:) = (/0.5_dp,0.0_dp,0.5_dp/)
       offsets(8,:) = (/0.0_dp,0.5_dp,0.5_dp/)
       
    end if

    allocate(grids(3,3,size(sp_hnfs,3)),STAT=status)
    if (status /=0) stop "Failed to allocate memory in find_kgrids."

    do i=1,size(sp_hnfs,3)
       temp = matmul(c_basis,sp_hnfs(:,:,i))
       grids(:,:,i) = matmul(Minv,temp)
       grids(:,:,i) = transpose(grids(:,:,i))
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

    if (abs((temp1**3)-kpd) > abs(((temp1+1)**3)-kpd)) then
       temp1 = temp1 + 1
    end if
    
    a = 0
    b = temp1**3
    c = 0
    
    if (lat_id==3) then
       mults = (/1,4,16/)
    else
       mults = (/1,2,4/)
    end if

    do i=1,temp1+1
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
  !!<parameter name="best_grid" regular="true">The offset for the grid given.</parameter>
  !!<parameter name="r_lattice" regular="true">The reciprocal cell
  !!lattice vectors.</parameter>
  !!<parameter name="offsets" regular="true">The possible offsets for the
  !!grid.</parameter>
  !!<parameter name="eps_" regular="true">Floating point tolerance.</parameter>
  SUBROUTINE grid_selection(grids, best_grid, best_offset, r_lattice, offsets,eps_)
    real(dp), allocatable, intent(in) :: grids(:,:,:), offsets(:,:)
    real(dp), intent(out) :: best_grid(3,3), best_offset(3)
    real(dp), intent(in) :: r_lattice(3,3)
    real(dp), intent(in), optional :: eps_

    real(dp) :: reduced_grid(3,3), norms(3)
    integer :: i, j
    real(dp) :: r_min, pf, r_min_best, pf_best, pi, eps, nIrKpts, best_nirkpts
    

    real(dp), pointer :: IRKps(:,:)
    integer, pointer :: weights(:)

    if (.not. present(eps_)) then
       eps = 1E-6
    else
       eps = eps_
    end if
    pi = 3.14159265358979323846
    
    r_min_best = 0
    pf_best = 0
    best_nirkpts = 0
    
    do i=1,size(grids,3)
       ! print *, "checking grid ",i," of ",size(grids,3)
       call minkowski_reduce_basis(grids(:,:,i),reduced_grid,eps)
       norms(1) = norm(reduced_grid(:,1))
       norms(2) = norm(reduced_grid(:,2))
       norms(3) = norm(reduced_grid(:,3))
       r_min = max(norms(1),norms(2),norms(3))
       pf = (4.0_dp/3.0_dp)*pi*(min(norms(1),norms(2),norms(3))**2)/(dot_product(reduced_grid(:,1),cross_product(reduced_grid(:,2),reduced_grid(:,3))))
       ! call generateIrredKpointList(grids(:,:,i),r_lattice,(/0.0_dp,0.0_dp,0.0_dp/),IRKps,weights,eps_=eps)
       ! nIrKpts = size(IRKps,1)
       ! if (best_nirkpts==0) then
       !    best_nirkpts = nIrKpts
       !    best_grid = grids(:,:,i)
       ! elseif (nIrKpts < best_nirkpts) then
       !    best_nirkpts = nIrKpts
       !    best_grid = grids(:,:,i)
       ! ! elseif (nIrKpts == best_nirkpts) then
       ! !    stop "Grids with equal number of irreducibel k-poins found. Need heuristic to determine which to use."
       ! end if
       if (r_min_best==0) then
          r_min_best = r_min
          ! pf_best = pf
          best_grid = grids(:,:,i)
          call generateIrredKpointList(grids(:,:,i),r_lattice,(/0.0_dp,0.0_dp,0.0_dp/),IRKps,weights,eps_=eps)
          best_nirkpts = size(IRKps,1)
       else
          if (r_min<r_min_best) then
             r_min_best = r_min
             ! pf_best = pf
             best_grid = grids(:,:,i)
             call generateIrredKpointList(grids(:,:,i),r_lattice,(/0.0_dp,0.0_dp,0.0_dp/),IRKps,weights,eps_=eps)
             best_nirkpts = size(IRKps,1)
          elseif (abs(r_min-r_min_best)<eps) then
             call generateIrredKpointList(grids(:,:,i),r_lattice,(/0.0_dp,0.0_dp,0.0_dp/),IRKps,weights,eps_=eps)
             nIrKpts = size(IRKps,1)
             
             if (nIrKpts < best_nirkpts) then
                r_min_best = r_min
                ! pf_best = pf
                best_grid = grids(:,:,i)
                best_nirkpts = nIrKpts
             end if
          end if
       end if
    end do

    best_nirkpts = 0
    do i=1,size(offsets,1)
       call generateIrredKpointList(best_grid,r_lattice,offsets(i,:),IRKps,weights,eps_=eps)
       nIrKpts = size(IRKps,1)
       if (best_nirkpts==0) then
          best_nirkpts = nIrKpts
          best_offset = offsets(i,:)
       else if (nIrKpts < best_nirkpts) then
          best_nirkpts = nIrKpts
          best_offset = offsets(i,:)
       end if
    end do

  end SUBROUTINE grid_selection
end Module find_kgrids
