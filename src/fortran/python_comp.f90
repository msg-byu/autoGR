PROGRAM py_comp
  use sp_hnfs
  
  implicit none

  integer :: max_size, i, count
  integer, allocatable :: sphnfs(:,:,:)
  
  max_size = 500

  count = 0
  do i=1,max_size
     call sc_3(i,sphnfs)
     if (allocated(sphnfs)) then
        count = count + size(sphnfs,3)
        deallocate(sphnfs)
     end if
  end do

  print *, count

end PROGRAM py_comp
