PROGRAM lat_id_driver
  use find_kgrids, only: find_grid
  use kpointgeneration, only: generateIrredKpointList, mapKptsIntoBZ
  use vector_matrix_utilities, only: matrix_inverse, volume
  use num_types
  ! use fortpy, only: pysave, fpy_read_f, fpy_read
  implicit none

  real(dp) :: lat_vecs(3,3), grid(3,3), offset(3), r_vecs(3,3), point(3), Rinv(3,3)
  ! real(dp) :: kpd
  integer :: i, nkpts, count, j, z, kpd
  real(dp), pointer :: IRKps(:,:)
  integer, pointer :: weights(:)

  character(300) :: line
  integer, allocatable :: at(:), concs(:)
  real(dp), pointer :: B_vecs(:,:)

  open(1,file="POSCAR",status="old")

  read(1,'(a300)') line
  read(1,'(a300)') line
  do i=1,3
     read(1,*) lat_vecs(:,i)
  end do
  read(1,'(a300)') line
  call parse_line_for_numbers(line,count,concs)
  allocate(at(sum(concs)),B_vecs(3,sum(concs)))
  read(1,'(a300)') line
  do i=1,sum(concs)
     read(1,*) B_vecs(:,i)     
  end do
  close(1)
  
  z = 1
  do i=1,count
     if (concs(i) > 0) then
        do j=1,concs(i)
           at(z) = i-1
           z = z+1        
        end do
     end if
  end do

  open(2,file="kpgen",status="old")
  read(2,*) kpd
  read(2,*) offset
  close(2)

  call matrix_inverse(transpose(lat_vecs),r_vecs)

  nkpts = kpd!int(abs(volume(r_vecs(:,1),r_vecs(:,2),r_vecs(:,3)))*kpd)

  call find_grid(lat_vecs, nkpts, B_vecs, at, grid)

  call generateIrredKpointList(lat_vecs,B_vecs,at,grid,r_vecs,offset,IRKps,weights,eps_=1E-6_dp)
  call mapKptsIntoBZ(r_vecs,IRKps)

  print *, "OUR GRID"
  do i=1,3
     print *, grid(:,i)
  end do
  
  open(4,file="KPOINTS")
  ! write(4,'("Our new kpoint method. ",i6)') sum(weights)
  write(4,'(A22, I12, A1, I6)')"Our new kpoint method.", sum(weights), '/', size(IRKps,1)
  write(4,*) size(IRKps,1)
  write(4,*) "Fractional"
  ! write(4,'(A4)')"cart"
  call matrix_inverse(r_vecs,Rinv)
  do i = 1,size(IRKps,1)
     point = matmul(Rinv,IRKps(i,:))
     write(4,'(F16.14,A1,F16.14,A1,F16.14,A1,I5.1)') point(1), " ",point(2), " ",point(3), " ", weights(i)
  end do

  ! do i=1,3
  !    write(4,'(F16.14,A1,F16.14,A1,F16.14)') grid(1,i), " ", grid(2,i), " ", grid(3,i)
  ! end do
  !    write(4,'(F16.14,A1,F16.14,A1,F16.14)') offset(1), " ", offset(2), " ", offset(3)
  
     
  
  ! print *, "choosen grid"
  ! do i=1,3
  !    print *, grid(:,i)
  ! end do
 
contains 
  subroutine split_str(instr,out1)
    character(300), intent(inout) :: instr
    character(300), intent(out) :: out1
    integer :: index

    instr = adjustl(TRIM(instr))

    index = SCAN(instr," ")
    out1 = adjustl(TRIM(instr(1:index-1)))
    instr = adjustl(TRIM(instr(index+1:)))
  end subroutine split_str
  
  subroutine parse_line_for_numbers(line,readEntries,values)
    character(300), intent(inout):: line
    integer, intent(out)          :: readEntries
    integer, allocatable, intent(out) :: values(:)

    character(300) :: temp
    integer iE, ierr
    logical :: done

    done = .False.
    if (allocated(values)) deallocate(values)
    allocate(values(len(line)))
    values = 0

    readEntries=0
    iE = 0

    line = trim(adjustl(line)) ! Remove preceding blanks                                         
    if (trim(line)=="") then ! make sure we didn't get a blank line
       readEntries = 0
       return
    endif
    
    do while (.not. done)
       iE = iE+1
       call split_str(line,temp)
       if (temp=="") then
          done = .True.
          exit
       end if
       read(temp,*,iostat=ierr) values(iE)
       if (ierr /=0) then
          readEntries=iE-1
          return 
          !stop "ERROR: line parsing for number failed. Maybe format mismatch?"
       endif
       if (line=="") then
          done = .True.
       end if
    enddo
    
    readEntries=iE
  end subroutine parse_line_for_numbers

end PROGRAM lat_id_driver
