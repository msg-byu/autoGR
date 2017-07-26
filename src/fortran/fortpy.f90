!!<fortpy codeversion="1.7.6" />
!!<summary>Provides an interface for saving the values of multiple variable
!!types using a single call. Used as part of the FORTPY unit testing framework.</summary>
module fortpy
  implicit none
  private
  public pysave, fdp, fsp, fsi, fli, fpy_linevalue_count, fpy_newunit, fpy_read, &
       fpy_value_count, fpy_period_join_indices, fpy_linevalue_count_all, fpy_read_p, &
       fpy_read_f, fpy_vararray, autoclass_analyze, fpy_set_verbosity, fpy_verbose

  !!<member name="fileunit">I/O unit for the file to write the output to.</member>
  integer :: fileunit
  !!<member>Output/debugging verbosity for auxiliary reads/writes.</member>
  integer :: fpy_verbose
  
  !!<member name="seeded">Specifies whether the random number generator has
  !!already been seeded for this program execution.</member>
  logical :: seeded

  integer, parameter:: fdp=selected_real_kind(15,307)
  integer, parameter:: fsp=selected_real_kind(6,37)
  integer, parameter:: fsi=selected_int_kind(1) ! very short integer -10..10 range
  integer, parameter:: fli=selected_int_kind(18) ! Big integer -10^18..10^18 range

  !!<summary>Provides a single call interface for integer and real data types
  !!for single values, 1D and 2D arrays. Also handles outputting true/false for
  !!logical type variables.</summary>
  interface pysave
     module procedure pysave_realsp_0d, pysave_realsp_1d, pysave_realsp_2d, &
                       pysave_realsp_3d, pysave_realsp_4d, pysave_realsp_5d, &
                       pysave_realsp_6d, pysave_realsp_7d, pysave_realdp_0d, &
                       pysave_realdp_1d, pysave_realdp_2d, pysave_realdp_3d, &
                       pysave_realdp_4d, pysave_realdp_5d, pysave_realdp_6d, &
                       pysave_realdp_7d, pysave_integer_0d, pysave_integer_1d, &
                       pysave_integer_2d, pysave_integer_3d, pysave_integer_4d, &
                       pysave_integer_5d, pysave_integer_6d, pysave_integer_7d, &
                       pysave_integersp_0d, pysave_integersp_1d, pysave_integersp_2d, &
                       pysave_integersp_3d, pysave_integersp_4d, pysave_integersp_5d, &
                       pysave_integersp_6d, pysave_integersp_7d, pysave_integerdp_0d, &
                       pysave_integerdp_1d, pysave_integerdp_2d, pysave_integerdp_3d, &
                       pysave_integerdp_4d, pysave_integerdp_5d, pysave_integerdp_6d, &
                       pysave_integerdp_7d, pysave_complexsp_0d, pysave_complexsp_1d, &
                       pysave_complexsp_2d, pysave_complexsp_3d, pysave_complexsp_4d, &
                       pysave_complexsp_5d, pysave_complexsp_6d, pysave_complexsp_7d, &
                       pysave_complexdp_0d, pysave_complexdp_1d, pysave_complexdp_2d, &
                       pysave_complexdp_3d, pysave_complexdp_4d, pysave_complexdp_5d, &
                       pysave_complexdp_6d, pysave_complexdp_7d, pysave_character_0d, &
                       pysave_character_1d, pysave_character_2d, pysave_character_3d, &
                       pysave_character_4d, pysave_character_5d, pysave_character_6d, &
                       pysave_character_7d, pysave_logical_0d, pysave_logical_1d, &
                       pysave_logical_2d, pysave_logical_3d, pysave_logical_4d, &
                       pysave_logical_5d, pysave_logical_6d, pysave_logical_7d
  end interface pysave

  !!<summary>Reads values from a data file into a variable.</summary>
  interface fpy_read
     module procedure fpy_read_realsp_0d, fpy_read_realsp_1d, fpy_read_realsp_2d, &
                       fpy_read_realsp_3d, fpy_read_realsp_4d, fpy_read_realsp_5d, &
                       fpy_read_realsp_6d, fpy_read_realsp_7d, fpy_read_realdp_0d, &
                       fpy_read_realdp_1d, fpy_read_realdp_2d, fpy_read_realdp_3d, &
                       fpy_read_realdp_4d, fpy_read_realdp_5d, fpy_read_realdp_6d, &
                       fpy_read_realdp_7d, fpy_read_integer_0d, fpy_read_integer_1d, &
                       fpy_read_integer_2d, fpy_read_integer_3d, fpy_read_integer_4d, &
                       fpy_read_integer_5d, fpy_read_integer_6d, fpy_read_integer_7d, &
                       fpy_read_integersp_0d, fpy_read_integersp_1d, fpy_read_integersp_2d, &
                       fpy_read_integersp_3d, fpy_read_integersp_4d, fpy_read_integersp_5d, &
                       fpy_read_integersp_6d, fpy_read_integersp_7d, fpy_read_integerdp_0d, &
                       fpy_read_integerdp_1d, fpy_read_integerdp_2d, fpy_read_integerdp_3d, &
                       fpy_read_integerdp_4d, fpy_read_integerdp_5d, fpy_read_integerdp_6d, &
                       fpy_read_integerdp_7d, fpy_read_complexsp_0d, fpy_read_complexsp_1d, &
                       fpy_read_complexsp_2d, fpy_read_complexsp_3d, fpy_read_complexsp_4d, &
                       fpy_read_complexsp_5d, fpy_read_complexsp_6d, fpy_read_complexsp_7d, &
                       fpy_read_complexdp_0d, fpy_read_complexdp_1d, fpy_read_complexdp_2d, &
                       fpy_read_complexdp_3d, fpy_read_complexdp_4d, fpy_read_complexdp_5d, &
                       fpy_read_complexdp_6d, fpy_read_complexdp_7d, fpy_read_character_0d, &
                       fpy_read_character_1d, fpy_read_character_2d, fpy_read_character_3d, &
                       fpy_read_character_4d, fpy_read_character_5d, fpy_read_character_6d, &
                       fpy_read_character_7d, fpy_read_logical_0d, fpy_read_logical_1d, &
                       fpy_read_logical_2d, fpy_read_logical_3d, fpy_read_logical_4d, &
                       fpy_read_logical_5d, fpy_read_logical_6d, fpy_read_logical_7d
  end interface fpy_read

  !!<summary>Interface for fpy_read for variables with fixed dimensions (i.e. they
  !!don't have the 'pointer' or 'allocatable' attributes.</summary>
  interface fpy_read_f
     module procedure fpy_read_realsp_1df, fpy_read_realsp_2df, fpy_read_realsp_3df, &
                       fpy_read_realsp_4df, fpy_read_realsp_5df, fpy_read_realsp_6df, &
                       fpy_read_realsp_7df, fpy_read_realdp_1df, fpy_read_realdp_2df, &
                       fpy_read_realdp_3df, fpy_read_realdp_4df, fpy_read_realdp_5df, &
                       fpy_read_realdp_6df, fpy_read_realdp_7df, fpy_read_integer_1df, &
                       fpy_read_integer_2df, fpy_read_integer_3df, fpy_read_integer_4df, &
                       fpy_read_integer_5df, fpy_read_integer_6df, fpy_read_integer_7df, &
                       fpy_read_integersp_1df, fpy_read_integersp_2df, &
                       fpy_read_integersp_3df, fpy_read_integersp_4df, &
                       fpy_read_integersp_5df, fpy_read_integersp_6df, &
                       fpy_read_integersp_7df, fpy_read_integerdp_1df, &
                       fpy_read_integerdp_2df, fpy_read_integerdp_3df, &
                       fpy_read_integerdp_4df, fpy_read_integerdp_5df, &
                       fpy_read_integerdp_6df, fpy_read_integerdp_7df, &
                       fpy_read_complexsp_1df, fpy_read_complexsp_2df, &
                       fpy_read_complexsp_3df, fpy_read_complexsp_4df, &
                       fpy_read_complexsp_5df, fpy_read_complexsp_6df, &
                       fpy_read_complexsp_7df, fpy_read_complexdp_1df, &
                       fpy_read_complexdp_2df, fpy_read_complexdp_3df, &
                       fpy_read_complexdp_4df, fpy_read_complexdp_5df, &
                       fpy_read_complexdp_6df, fpy_read_complexdp_7df, &
                       fpy_read_character_1df, fpy_read_character_2df, &
                       fpy_read_character_3df, fpy_read_character_4df, &
                       fpy_read_character_5df, fpy_read_character_6df, &
                       fpy_read_character_7df, fpy_read_logical_1df, fpy_read_logical_2df, &
                       fpy_read_logical_3df, fpy_read_logical_4df, fpy_read_logical_5df, &
                       fpy_read_logical_6df, fpy_read_logical_7df
  end interface fpy_read_f
  
  !!<summary>Provides an interface for fpy_read for pointer-valued variables that
  !!need to be allocated before reading data.</summary>
  interface fpy_read_p
     module procedure fpy_read_realsp_1dp, fpy_read_realsp_2dp, fpy_read_realsp_3dp, &
                       fpy_read_realsp_4dp, fpy_read_realsp_5dp, fpy_read_realsp_6dp, &
                       fpy_read_realsp_7dp, fpy_read_realdp_1dp, fpy_read_realdp_2dp, &
                       fpy_read_realdp_3dp, fpy_read_realdp_4dp, fpy_read_realdp_5dp, &
                       fpy_read_realdp_6dp, fpy_read_realdp_7dp, fpy_read_integer_1dp, &
                       fpy_read_integer_2dp, fpy_read_integer_3dp, fpy_read_integer_4dp, &
                       fpy_read_integer_5dp, fpy_read_integer_6dp, fpy_read_integer_7dp, &
                       fpy_read_integersp_1dp, fpy_read_integersp_2dp, &
                       fpy_read_integersp_3dp, fpy_read_integersp_4dp, &
                       fpy_read_integersp_5dp, fpy_read_integersp_6dp, &
                       fpy_read_integersp_7dp, fpy_read_integerdp_1dp, &
                       fpy_read_integerdp_2dp, fpy_read_integerdp_3dp, &
                       fpy_read_integerdp_4dp, fpy_read_integerdp_5dp, &
                       fpy_read_integerdp_6dp, fpy_read_integerdp_7dp, &
                       fpy_read_complexsp_1dp, fpy_read_complexsp_2dp, &
                       fpy_read_complexsp_3dp, fpy_read_complexsp_4dp, &
                       fpy_read_complexsp_5dp, fpy_read_complexsp_6dp, &
                       fpy_read_complexsp_7dp, fpy_read_complexdp_1dp, &
                       fpy_read_complexdp_2dp, fpy_read_complexdp_3dp, &
                       fpy_read_complexdp_4dp, fpy_read_complexdp_5dp, &
                       fpy_read_complexdp_6dp, fpy_read_complexdp_7dp, &
                       fpy_read_character_1dp, fpy_read_character_2dp, &
                       fpy_read_character_3dp, fpy_read_character_4dp, &
                       fpy_read_character_5dp, fpy_read_character_6dp, &
                       fpy_read_character_7dp, fpy_read_logical_1dp, fpy_read_logical_2dp, &
                       fpy_read_logical_3dp, fpy_read_logical_4dp, fpy_read_logical_5dp, &
                       fpy_read_logical_6dp, fpy_read_logical_7dp
  end interface fpy_read_p

  !!<summary>Provides a structure for variable length arrays that need to have their
  !!cartesian product taken.</summary>
  !!<usage>
  !!type(fpy_vararray) single
  !!allocate(single)
  !!single%init(array)
  !!</usage>
  type fpy_vararray
     !!<member name="items">The array of items to take cartesian product over.</member>
     !!<member name="length">The number of items in @CREF[this.items].</member>
     integer, pointer :: items(:)
     integer :: length
  contains
    procedure, public :: init => vararray_init
  end type fpy_vararray
contains
  !!<summary>Sets the global verbosity for the auxiliary reads/writes.</summary>
  subroutine fpy_set_verbosity(v)
    integer, intent(in) :: v
    fpy_verbose = v
    if (v > 0) write (*, *) "Fortpy F90 verbosity set to", v, "."
  end subroutine fpy_set_verbosity

  !!<summary>Analyzes the specified .fortpy.analysis file to determine the
  !!actual dimensionality of the data being read in via auto-class.</summary>
  subroutine autoclass_analyze(filename, analysis)
    character(len=*), intent(in) :: filename
    class(fpy_vararray), allocatable, intent(out) :: analysis(:)

    !!<local name="line">The integer dimensionality specified by a single line of the file.</local>
    !!<local name="ragvals">The number of values on each line of the file.</local>
    integer, allocatable :: line(:), ragvals(:)
    integer :: nlines, nvalues, funit, i

    call fpy_linevalue_count_all(filename, '#', nlines, ragvals)
    allocate(analysis(nlines))

    open(fpy_newunit(funit), file=filename)
    do i=1, nlines
      allocate(line(ragvals(i)))
      read(funit, *) line
      call analysis(i)%init(line, alloc=.true.)
      deallocate(line)
    end do
    close(funit)
  end subroutine autoclass_analyze
  
  !!<summary>Initializes the array items and length property.</summary>
  subroutine vararray_init(self, array, length, alloc)
    class(fpy_vararray) :: self
    integer, target, optional, intent(in) :: array(:)
    integer, optional, intent(in) :: length
    logical, optional, intent(in) :: alloc

    logical :: nalloc
    !We need to see if we are *copying* the array, or just referencing it.
    if (present(alloc)) then
       nalloc = alloc
    else
       nalloc = .false.
    end if

    if (present(array)) then
       if (nalloc) then
          allocate(self%items(size(array, 1)))
          self%items = array
       else
          self%items => array
       end if
       self%length = size(self%items, 1)
    else
       allocate(self%items(length))
       self%length = length
    end if
  end subroutine vararray_init

    subroutine fpy_read_realsp_0d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), intent(inout) :: variable

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nvalues .gt. 1) .or. (nlines /= 1)) then
      write(*,*) "Cannot read a single value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_0d
    
  subroutine fpy_read_realsp_1d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), allocatable, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_1d
    
  subroutine fpy_read_realsp_2d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), allocatable, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_2d
    
  subroutine fpy_read_realsp_3d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), allocatable, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_3d
    
  subroutine fpy_read_realsp_4d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), allocatable, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_4d
    
  subroutine fpy_read_realsp_5d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), allocatable, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_5d
    
  subroutine fpy_read_realsp_6d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), allocatable, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_6d
    
  subroutine fpy_read_realsp_7d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), allocatable, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_7d
      subroutine fpy_read_realdp_0d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), intent(inout) :: variable

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nvalues .gt. 1) .or. (nlines /= 1)) then
      write(*,*) "Cannot read a single value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_0d
    
  subroutine fpy_read_realdp_1d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), allocatable, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_1d
    
  subroutine fpy_read_realdp_2d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), allocatable, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_2d
    
  subroutine fpy_read_realdp_3d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), allocatable, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_3d
    
  subroutine fpy_read_realdp_4d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), allocatable, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_4d
    
  subroutine fpy_read_realdp_5d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), allocatable, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_5d
    
  subroutine fpy_read_realdp_6d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), allocatable, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_6d
    
  subroutine fpy_read_realdp_7d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), allocatable, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_7d
    
  subroutine fpy_read_integer_0d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, intent(inout) :: variable

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nvalues .gt. 1) .or. (nlines /= 1)) then
      write(*,*) "Cannot read a single value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_0d
    
  subroutine fpy_read_integer_1d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, allocatable, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_1d
    
  subroutine fpy_read_integer_2d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, allocatable, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_2d
    
  subroutine fpy_read_integer_3d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, allocatable, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_3d
    
  subroutine fpy_read_integer_4d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, allocatable, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_4d
    
  subroutine fpy_read_integer_5d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, allocatable, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_5d
    
  subroutine fpy_read_integer_6d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, allocatable, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_6d
    
  subroutine fpy_read_integer_7d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, allocatable, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_7d
      subroutine fpy_read_integersp_0d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), intent(inout) :: variable

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nvalues .gt. 1) .or. (nlines /= 1)) then
      write(*,*) "Cannot read a single value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_0d
    
  subroutine fpy_read_integersp_1d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), allocatable, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_1d
    
  subroutine fpy_read_integersp_2d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), allocatable, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_2d
    
  subroutine fpy_read_integersp_3d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), allocatable, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_3d
    
  subroutine fpy_read_integersp_4d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), allocatable, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_4d
    
  subroutine fpy_read_integersp_5d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), allocatable, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_5d
    
  subroutine fpy_read_integersp_6d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), allocatable, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_6d
    
  subroutine fpy_read_integersp_7d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), allocatable, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_7d
      subroutine fpy_read_integerdp_0d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), intent(inout) :: variable

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nvalues .gt. 1) .or. (nlines /= 1)) then
      write(*,*) "Cannot read a single value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_0d
    
  subroutine fpy_read_integerdp_1d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), allocatable, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_1d
    
  subroutine fpy_read_integerdp_2d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), allocatable, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_2d
    
  subroutine fpy_read_integerdp_3d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), allocatable, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_3d
    
  subroutine fpy_read_integerdp_4d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), allocatable, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_4d
    
  subroutine fpy_read_integerdp_5d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), allocatable, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_5d
    
  subroutine fpy_read_integerdp_6d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), allocatable, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_6d
    
  subroutine fpy_read_integerdp_7d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), allocatable, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_7d
    
  subroutine fpy_read_complexsp_0d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), intent(inout) :: variable

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nvalues .gt. 1) .or. (nlines /= 1)) then
      write(*,*) "Cannot read a single value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_0d
    
  subroutine fpy_read_complexsp_1d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), allocatable, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_1d
    
  subroutine fpy_read_complexsp_2d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), allocatable, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_2d
    
  subroutine fpy_read_complexsp_3d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), allocatable, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_3d
    
  subroutine fpy_read_complexsp_4d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), allocatable, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_4d
    
  subroutine fpy_read_complexsp_5d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), allocatable, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_5d
    
  subroutine fpy_read_complexsp_6d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), allocatable, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_6d
    
  subroutine fpy_read_complexsp_7d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), allocatable, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_7d
      subroutine fpy_read_complexdp_0d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), intent(inout) :: variable

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nvalues .gt. 1) .or. (nlines /= 1)) then
      write(*,*) "Cannot read a single value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_0d
    
  subroutine fpy_read_complexdp_1d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), allocatable, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_1d
    
  subroutine fpy_read_complexdp_2d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), allocatable, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_2d
    
  subroutine fpy_read_complexdp_3d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), allocatable, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_3d
    
  subroutine fpy_read_complexdp_4d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), allocatable, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_4d
    
  subroutine fpy_read_complexdp_5d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), allocatable, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_5d
    
  subroutine fpy_read_complexdp_6d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), allocatable, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_6d
    
  subroutine fpy_read_complexdp_7d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), allocatable, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_7d
    
  subroutine fpy_read_character_0d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), intent(inout) :: variable

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues, .true.)
    if ((nvalues .gt. 1) .or. (nlines /= 1)) then
      write(*,*) "Cannot read a single value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    variable = ''


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_0d
    
  subroutine fpy_read_character_1d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), allocatable, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues, .true.)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = ''


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_1d
    
  subroutine fpy_read_character_2d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), allocatable, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues, .true.)
    allocate(variable(nlines, nvalues))
    variable = ''


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_2d
    
  subroutine fpy_read_character_3d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), allocatable, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_3d
    
  subroutine fpy_read_character_4d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), allocatable, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_4d
    
  subroutine fpy_read_character_5d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), allocatable, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_5d
    
  subroutine fpy_read_character_6d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), allocatable, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_6d
    
  subroutine fpy_read_character_7d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), allocatable, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_7d
    
  subroutine fpy_read_logical_0d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, intent(inout) :: variable

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nvalues .gt. 1) .or. (nlines /= 1)) then
      write(*,*) "Cannot read a single value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    variable = .false.


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_0d
    
  subroutine fpy_read_logical_1d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, allocatable, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = .false.


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_1d
    
  subroutine fpy_read_logical_2d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, allocatable, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = .false.


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_2d
    
  subroutine fpy_read_logical_3d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, allocatable, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_3d
    
  subroutine fpy_read_logical_4d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, allocatable, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_4d
    
  subroutine fpy_read_logical_5d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, allocatable, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_5d
    
  subroutine fpy_read_logical_6d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, allocatable, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_6d
    
  subroutine fpy_read_logical_7d(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, allocatable, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      return
    end if

    if (allocated(variable)) deallocate(variable)


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_7d
    

    subroutine fpy_read_realsp_1df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if

    
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 1))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_1df
    
  subroutine fpy_read_realsp_2df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 2))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_2df
    
  subroutine fpy_read_realsp_3df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_3df
    
  subroutine fpy_read_realsp_4df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_4df
    
  subroutine fpy_read_realsp_5df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_5df
    
  subroutine fpy_read_realsp_6df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_6df
    
  subroutine fpy_read_realsp_7df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_7df
      subroutine fpy_read_realdp_1df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if

    
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 1))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_1df
    
  subroutine fpy_read_realdp_2df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 2))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_2df
    
  subroutine fpy_read_realdp_3df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_3df
    
  subroutine fpy_read_realdp_4df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_4df
    
  subroutine fpy_read_realdp_5df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_5df
    
  subroutine fpy_read_realdp_6df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_6df
    
  subroutine fpy_read_realdp_7df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_7df
    
  subroutine fpy_read_integer_1df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if

    
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 1))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_1df
    
  subroutine fpy_read_integer_2df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 2))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_2df
    
  subroutine fpy_read_integer_3df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_3df
    
  subroutine fpy_read_integer_4df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_4df
    
  subroutine fpy_read_integer_5df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_5df
    
  subroutine fpy_read_integer_6df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_6df
    
  subroutine fpy_read_integer_7df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_7df
      subroutine fpy_read_integersp_1df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if

    
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 1))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_1df
    
  subroutine fpy_read_integersp_2df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 2))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_2df
    
  subroutine fpy_read_integersp_3df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_3df
    
  subroutine fpy_read_integersp_4df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_4df
    
  subroutine fpy_read_integersp_5df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_5df
    
  subroutine fpy_read_integersp_6df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_6df
    
  subroutine fpy_read_integersp_7df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_7df
      subroutine fpy_read_integerdp_1df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if

    
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 1))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_1df
    
  subroutine fpy_read_integerdp_2df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 2))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_2df
    
  subroutine fpy_read_integerdp_3df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_3df
    
  subroutine fpy_read_integerdp_4df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_4df
    
  subroutine fpy_read_integerdp_5df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_5df
    
  subroutine fpy_read_integerdp_6df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_6df
    
  subroutine fpy_read_integerdp_7df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_7df
    
  subroutine fpy_read_complexsp_1df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if

    
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 1))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_1df
    
  subroutine fpy_read_complexsp_2df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 2))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_2df
    
  subroutine fpy_read_complexsp_3df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_3df
    
  subroutine fpy_read_complexsp_4df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_4df
    
  subroutine fpy_read_complexsp_5df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_5df
    
  subroutine fpy_read_complexsp_6df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_6df
    
  subroutine fpy_read_complexsp_7df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_7df
      subroutine fpy_read_complexdp_1df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if

    
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 1))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_1df
    
  subroutine fpy_read_complexdp_2df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 2))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_2df
    
  subroutine fpy_read_complexdp_3df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_3df
    
  subroutine fpy_read_complexdp_4df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_4df
    
  subroutine fpy_read_complexdp_5df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_5df
    
  subroutine fpy_read_complexdp_6df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_6df
    
  subroutine fpy_read_complexdp_7df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = 0
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_7df
    
  subroutine fpy_read_character_1df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = ''
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues, .true.)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if

    
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 1))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_1df
    
  subroutine fpy_read_character_2df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = ''
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues, .true.)
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 2))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_2df
    
  subroutine fpy_read_character_3df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = ''
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_3df
    
  subroutine fpy_read_character_4df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = ''
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_4df
    
  subroutine fpy_read_character_5df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = ''
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_5df
    
  subroutine fpy_read_character_6df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = ''
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_6df
    
  subroutine fpy_read_character_7df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = ''
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_7df
    
  subroutine fpy_read_logical_1df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = .false.
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if

    
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 1))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_1df
    
  subroutine fpy_read_logical_2df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = .false.
      return
    end if

    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .ne. size(variable, 1)) .and. (nvalues .ne. size(variable, 2))) then
      write(*,*) "File data dimensions don't match fixed variable shape ", shape(variable)
      write(*,*) "Fortpy sees data dimensions in '", filename, "' as ", nlines, nvalues
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_2df
    
  subroutine fpy_read_logical_3df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = .false.
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_3df
    
  subroutine fpy_read_logical_4df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = .false.
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_4df
    
  subroutine fpy_read_logical_5df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = .false.
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_5df
    
  subroutine fpy_read_logical_6df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = .false.
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_6df
    
  subroutine fpy_read_logical_7df(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable = .false.
      return
    end if



    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_7df
    

    subroutine fpy_read_realsp_1dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), pointer, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_1dp
    
  subroutine fpy_read_realsp_2dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), pointer, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_2dp
    
  subroutine fpy_read_realsp_3dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), pointer, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_3dp
    
  subroutine fpy_read_realsp_4dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), pointer, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_4dp
    
  subroutine fpy_read_realsp_5dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), pointer, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_5dp
    
  subroutine fpy_read_realsp_6dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), pointer, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_6dp
    
  subroutine fpy_read_realsp_7dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fsp), pointer, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realsp_7dp
      subroutine fpy_read_realdp_1dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), pointer, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_1dp
    
  subroutine fpy_read_realdp_2dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), pointer, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_2dp
    
  subroutine fpy_read_realdp_3dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), pointer, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_3dp
    
  subroutine fpy_read_realdp_4dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), pointer, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_4dp
    
  subroutine fpy_read_realdp_5dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), pointer, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_5dp
    
  subroutine fpy_read_realdp_6dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), pointer, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_6dp
    
  subroutine fpy_read_realdp_7dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    real(fdp), pointer, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_realdp_7dp
    
  subroutine fpy_read_integer_1dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, pointer, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_1dp
    
  subroutine fpy_read_integer_2dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, pointer, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_2dp
    
  subroutine fpy_read_integer_3dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, pointer, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_3dp
    
  subroutine fpy_read_integer_4dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, pointer, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_4dp
    
  subroutine fpy_read_integer_5dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, pointer, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_5dp
    
  subroutine fpy_read_integer_6dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, pointer, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_6dp
    
  subroutine fpy_read_integer_7dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer, pointer, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integer_7dp
      subroutine fpy_read_integersp_1dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), pointer, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_1dp
    
  subroutine fpy_read_integersp_2dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), pointer, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_2dp
    
  subroutine fpy_read_integersp_3dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), pointer, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_3dp
    
  subroutine fpy_read_integersp_4dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), pointer, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_4dp
    
  subroutine fpy_read_integersp_5dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), pointer, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_5dp
    
  subroutine fpy_read_integersp_6dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), pointer, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_6dp
    
  subroutine fpy_read_integersp_7dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fsi), pointer, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integersp_7dp
      subroutine fpy_read_integerdp_1dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), pointer, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_1dp
    
  subroutine fpy_read_integerdp_2dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), pointer, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_2dp
    
  subroutine fpy_read_integerdp_3dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), pointer, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_3dp
    
  subroutine fpy_read_integerdp_4dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), pointer, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_4dp
    
  subroutine fpy_read_integerdp_5dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), pointer, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_5dp
    
  subroutine fpy_read_integerdp_6dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), pointer, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_6dp
    
  subroutine fpy_read_integerdp_7dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    integer(fli), pointer, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_integerdp_7dp
    
  subroutine fpy_read_complexsp_1dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), pointer, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_1dp
    
  subroutine fpy_read_complexsp_2dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), pointer, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_2dp
    
  subroutine fpy_read_complexsp_3dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), pointer, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_3dp
    
  subroutine fpy_read_complexsp_4dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), pointer, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_4dp
    
  subroutine fpy_read_complexsp_5dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), pointer, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_5dp
    
  subroutine fpy_read_complexsp_6dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), pointer, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_6dp
    
  subroutine fpy_read_complexsp_7dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fsp), pointer, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexsp_7dp
      subroutine fpy_read_complexdp_1dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), pointer, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_1dp
    
  subroutine fpy_read_complexdp_2dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), pointer, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = 0


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_2dp
    
  subroutine fpy_read_complexdp_3dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), pointer, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_3dp
    
  subroutine fpy_read_complexdp_4dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), pointer, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_4dp
    
  subroutine fpy_read_complexdp_5dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), pointer, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_5dp
    
  subroutine fpy_read_complexdp_6dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), pointer, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_6dp
    
  subroutine fpy_read_complexdp_7dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    complex(fdp), pointer, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = 0
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_complexdp_7dp
    
  subroutine fpy_read_character_1dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), pointer, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues, .true.)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = ''


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_1dp
    
  subroutine fpy_read_character_2dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), pointer, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues, .true.)
    allocate(variable(nlines, nvalues))
    variable = ''


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_2dp
    
  subroutine fpy_read_character_3dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), pointer, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_3dp
    
  subroutine fpy_read_character_4dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), pointer, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_4dp
    
  subroutine fpy_read_character_5dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), pointer, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_5dp
    
  subroutine fpy_read_character_6dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), pointer, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_6dp
    
  subroutine fpy_read_character_7dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    character(len=*), pointer, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = ''
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_character_7dp
    
  subroutine fpy_read_logical_1dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, pointer, intent(inout) :: variable(:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    if ((nlines .gt. 1 .and. nvalues .gt. 1) .or. (nlines .eq. 0 .or. nvalues .eq. 0)) then
      write(*,*) "Cannot read a vector value from ", filename
      write(*,*) "Found ", nlines, " lines and ", nvalues, " values"
      if (present(success_)) success_ = .false.
      if (strict) stop
    end if
    if (nlines .gt. 1) then
      allocate(variable(nlines))
    else
      allocate(variable(nvalues))
    end if
    variable = .false.


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              if (nlines .gt. 1) then
                read(line, *) variable(i)
                i = i+1
              else
                read(line, *) variable
              end if
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_1dp
    
  subroutine fpy_read_logical_2dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, pointer, intent(inout) :: variable(:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: nlines, nvalues, i
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()
    call fpy_linevalue_count(filename, commentchar, nlines, nvalues)
    allocate(variable(nlines, nvalues))
    variable = .false.


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      i=1
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i,:)
              i = i+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_2dp
    
  subroutine fpy_read_logical_3dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, pointer, intent(inout) :: variable(:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(3), i1, i2, i3, indices(3)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 3D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3)))
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,:)
              i2 = i2+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_3dp
    
  subroutine fpy_read_logical_4dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, pointer, intent(inout) :: variable(:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(4), i1, i2, i3, i4, indices(4)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 4D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4)))
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,:)
              i3 = i3+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_4dp
    
  subroutine fpy_read_logical_5dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, pointer, intent(inout) :: variable(:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(5), i1, i2, i3, i4, i5, indices(5)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 5D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5)))
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,:)
              i4 = i4+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_5dp
    
  subroutine fpy_read_logical_6dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, pointer, intent(inout) :: variable(:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(6), i1, i2, i3, i4, i5, i6, indices(6)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 6D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)))
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,:)
              i5 = i5+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_6dp
    
  subroutine fpy_read_logical_7dp(filename, commentchar, variable, success_, strict_)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    logical, optional, intent(out) :: success_
    logical, optional, intent(in) :: strict_
    logical, pointer, intent(inout) :: variable(:,:,:,:,:,:,:)

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    logical :: exists, strict
    integer :: dims(7), i1, i2, i3, i4, i5, i6, i7, indices(7)
    character(250000) :: line

    if (present(strict_)) then
      strict = strict_
    else
      strict = .false.
    end if

    inquire(file=filename, exist=exists)
    if (present(success_)) success_ = exists .or. .false.
    if (.not. exists) then
      if (fpy_verbose > 0) write (*,*) "Target file '", filename, "' does not exist."
      variable => null()
      return
    end if

    if (associated(variable)) variable => null()


    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 1) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:) ,*) dims
              exit
            end if
          end if
        else
          write(*,*) "No shape information found for 7D data file '", filename, "'"
          stop
        end if
      end do
      
      allocate(variable(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)))
      variable = .false.
      
      do
        read(funit, "(A)", iostat=ioerr) line
        if (ioerr == 0) then
          cleaned = trim(adjustl(line))
          if (len(cleaned) .gt. 0) then
            if (cleaned(1:2) .eq. '##') then
              read(cleaned(3:), *) indices
              i1 = indices(1)
              i2 = indices(2)
              i3 = indices(3)
              i4 = indices(4)
              i5 = indices(5)
              i6 = 1
            elseif (cleaned(1:1) /= commentchar) then
              read(line, *) variable(i1,i2,i3,i4,i5,i6,:)
              i6 = i6+1
            end if
          end if
        else
          exit
        end if
      end do
    end if
    close(funit)    
    if (present(success_)) success_ = .true.
  end subroutine fpy_read_logical_7dp
    

  
  !!<summary>Joins an array of integer indices as a period-separated string for
  !!concatenating to a file name.</summary>
  !!<parameter name="pslist">The period-separated result.</parameter>
  !!<parameter name="indices">An array of integer indices to join together.</parameter>
  !!<parameter name="n">The number of entries in the 'indices' array.</parameter>
  subroutine fpy_period_join_indices(pslist, indices, n)
    character(100), intent(out) :: pslist
    integer :: n
    integer, intent(in) :: indices(n)

    character(n*10) :: tempstr
    character(40) :: buffer
    integer :: i

    tempstr = ""
    pslist = ""
    do i=1, n
       write(buffer, '(I10)') indices(i)
       buffer = adjustl(buffer)
       if (i .eq. n) then
          tempstr = trim(pslist) // trim(buffer)
          pslist = tempstr
       else
          tempstr = trim(pslist) // trim(buffer) // "."
          pslist = tempstr
       end if
    end do
    pslist = trim(tempstr)
  end subroutine fpy_period_join_indices

  subroutine random_real(variable, min, max)
    real(fdp), intent(out) :: variable
    real(fdp), intent(in) :: min, max
    real(fdp) :: x
    !Since we only get real numbers in [0, 1], we have to set the value to be 
    !in the range manually.
    call random_number(x)    
    variable = min + x * (max - min)
  end subroutine random_real

  subroutine random_integer(variable, min, max)
    integer, intent(out) :: variable
    integer, intent(in) :: min, max
    real(fdp) :: x
    !Since we only get real numbers, we have to set the value to be in the range
    !and then turn it to an integer.
    call random_number(x)    
    variable = ceiling(min + x * (max - min))
  end subroutine random_integer

  !!<summary>Initializes the seed of the random number generator.</summary>
  subroutine random_init(seed)
    integer, intent(in), optional :: seed(:)
    if (.not. seeded) then
       if (present(seed)) then
          call random_seed(PUT=seed)
       else
          call random_seed()
       end if
       seeded = .true.
    end if
  end subroutine random_init

  !!<summary>Escapes all the apostrophes in each word and surrounds it
  !!with a set of apostrophes.</summary>
  subroutine char_escape_word(word, escaped)
    character(len=*), intent(in) :: word
    character(1000), intent(inout) :: escaped
    integer :: ichar, ibuff

    escaped = "'"
    ibuff = 2
    do ichar=1, len(trim(adjustl(word)))
       if (word(ichar:ichar) .eq. "'") then
          escaped(ibuff:ibuff+1) = "''"
          ibuff = ibuff + 2
       else
          escaped(ibuff:ibuff) = word(ichar:ichar)
          ibuff = ibuff + 1
       end if
    end do
    escaped(ibuff:ibuff) = "'"
  end subroutine char_escape_word

  !!<summary>Writes a vector array of *escaped* strings as a single line
  !!to the currently open file unit.</summary>
  subroutine char_write_trimmed(variable)
    character(len=*), intent(in) :: variable(:)
    !!<local name="fmt">The compiled format string that includes an 'A#' directive
    !!for each string in the variable array where # is the length of the string.</local>
    !!<local name="word">Holds the length of the current array element as a string.</local>
    character(len=:), allocatable :: fmt
    character(10) :: word
    integer :: i
    character(1000) :: escaped(size(variable, 1))

    !First we need to escape all the quotes in each word in the array
    !so that they can be read in correctly by Fortran later.
    do i=1, size(variable, 1)
       call char_escape_word(variable(i), escaped(i))
    end do

    fmt = "("
    do i=1, size(escaped, 1)
       write (word, "(i10)") (len(trim(adjustl(escaped(i)))) + 2)
       fmt = fmt // "A" // trim(adjustl(word))
       if (i .lt. size(escaped, 1)) then
          fmt = fmt // ","
       end if
    end do
    fmt = fmt // ")"
    write(fileunit, fmt) escaped
  end subroutine char_write_trimmed
  
    subroutine pysave_realsp_0d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fsp), intent(in) :: variable
    character(20) :: FMT
    
    write(FMT, *) 1

    call file_open(filename, len(filename), 'real')
    write(fileunit, *) variable
    call file_close()
  end subroutine pysave_realsp_0d
    
  subroutine pysave_realsp_1d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fsp), intent(in) :: variable(:)
    character(20) :: FMT
    integer :: dims(1)

    dims = shape(variable)
    write(FMT, *) dims(1)
    if (dims(1) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    write(fileunit, *) variable(:)
    call file_close()
  end subroutine pysave_realsp_1d
    
  subroutine pysave_realsp_2d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fsp), intent(in) :: variable(:,:)
    character(20) :: FMT
    integer :: dims(2), i1

    dims = shape(variable)
    write(FMT, *) dims(2)
    if (dims(2) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    do i1=1, dims(1)
      write(fileunit, *) variable(i1, :)
    end do
    call file_close()
  end subroutine pysave_realsp_2d
    
  subroutine pysave_realsp_3d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fsp), intent(in) :: variable(:,:,:)
    character(20) :: FMT
    integer :: dims(3), i1, i2

    dims = shape(variable)
    write(FMT, *) dims(3)
    if (dims(3) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    write(fileunit, "(A, 3i15)") "## ", dims
    do i1=1, dims(1)
      write(fileunit, "(A, 3i15)") "## ", i1, 0, 0
      do i2=1, dims(2)
        write(fileunit, *) variable(i1, i2, :)
      end do
    end do
    call file_close()
  end subroutine pysave_realsp_3d
    
  subroutine pysave_realsp_4d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fsp), intent(in) :: variable(:,:,:,:)
    character(20) :: FMT
    integer :: dims(4), i1, i2, i3

    dims = shape(variable)
    write(FMT, *) dims(4)
    if (dims(4) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    write(fileunit, "(A, 4i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        write(fileunit, "(A, 4i15)") "## ", i1, i2, 0, 0
        do i3=1, dims(3)
          write(fileunit, *) variable(i1, i2, i3, :)
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_realsp_4d
    
  subroutine pysave_realsp_5d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fsp), intent(in) :: variable(:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(5), i1, i2, i3, i4

    dims = shape(variable)
    write(FMT, *) dims(5)
    if (dims(5) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    write(fileunit, "(A, 5i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          write(fileunit, "(A, 5i15)") "## ", i1, i2, i3, 0, 0
          do i4=1, dims(4)
            write(fileunit, *) variable(i1, i2, i3, i4, :)
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_realsp_5d
    
  subroutine pysave_realsp_6d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fsp), intent(in) :: variable(:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(6), i1, i2, i3, i4, i5

    dims = shape(variable)
    write(FMT, *) dims(6)
    if (dims(6) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    write(fileunit, "(A, 6i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            write(fileunit, "(A, 6i15)") "## ", i1, i2, i3, i4, 0, 0
            do i5=1, dims(5)
              write(fileunit, *) variable(i1, i2, i3, i4, i5, :)
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_realsp_6d
    
  subroutine pysave_realsp_7d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fsp), intent(in) :: variable(:,:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(7), i1, i2, i3, i4, i5, i6

    dims = shape(variable)
    write(FMT, *) dims(7)
    if (dims(7) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    write(fileunit, "(A, 7i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            do i5=1, dims(5)
              write(fileunit, "(A, 7i15)") "## ", i1, i2, i3, i4, i5, 0, 0
              do i6=1, dims(6)
                write(fileunit, *) variable(i1, i2, i3, i4, i5, i6, :)
              end do
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_realsp_7d
      subroutine pysave_realdp_0d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fdp), intent(in) :: variable
    character(20) :: FMT
    
    write(FMT, *) 1

    call file_open(filename, len(filename), 'real')
    write(fileunit, *) variable
    call file_close()
  end subroutine pysave_realdp_0d
    
  subroutine pysave_realdp_1d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fdp), intent(in) :: variable(:)
    character(20) :: FMT
    integer :: dims(1)

    dims = shape(variable)
    write(FMT, *) dims(1)
    if (dims(1) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    write(fileunit, *) variable(:)
    call file_close()
  end subroutine pysave_realdp_1d
    
  subroutine pysave_realdp_2d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fdp), intent(in) :: variable(:,:)
    character(20) :: FMT
    integer :: dims(2), i1

    dims = shape(variable)
    write(FMT, *) dims(2)
    if (dims(2) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    do i1=1, dims(1)
      write(fileunit, *) variable(i1, :)
    end do
    call file_close()
  end subroutine pysave_realdp_2d
    
  subroutine pysave_realdp_3d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fdp), intent(in) :: variable(:,:,:)
    character(20) :: FMT
    integer :: dims(3), i1, i2

    dims = shape(variable)
    write(FMT, *) dims(3)
    if (dims(3) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    write(fileunit, "(A, 3i15)") "## ", dims
    do i1=1, dims(1)
      write(fileunit, "(A, 3i15)") "## ", i1, 0, 0
      do i2=1, dims(2)
        write(fileunit, *) variable(i1, i2, :)
      end do
    end do
    call file_close()
  end subroutine pysave_realdp_3d
    
  subroutine pysave_realdp_4d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fdp), intent(in) :: variable(:,:,:,:)
    character(20) :: FMT
    integer :: dims(4), i1, i2, i3

    dims = shape(variable)
    write(FMT, *) dims(4)
    if (dims(4) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    write(fileunit, "(A, 4i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        write(fileunit, "(A, 4i15)") "## ", i1, i2, 0, 0
        do i3=1, dims(3)
          write(fileunit, *) variable(i1, i2, i3, :)
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_realdp_4d
    
  subroutine pysave_realdp_5d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fdp), intent(in) :: variable(:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(5), i1, i2, i3, i4

    dims = shape(variable)
    write(FMT, *) dims(5)
    if (dims(5) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    write(fileunit, "(A, 5i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          write(fileunit, "(A, 5i15)") "## ", i1, i2, i3, 0, 0
          do i4=1, dims(4)
            write(fileunit, *) variable(i1, i2, i3, i4, :)
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_realdp_5d
    
  subroutine pysave_realdp_6d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fdp), intent(in) :: variable(:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(6), i1, i2, i3, i4, i5

    dims = shape(variable)
    write(FMT, *) dims(6)
    if (dims(6) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    write(fileunit, "(A, 6i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            write(fileunit, "(A, 6i15)") "## ", i1, i2, i3, i4, 0, 0
            do i5=1, dims(5)
              write(fileunit, *) variable(i1, i2, i3, i4, i5, :)
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_realdp_6d
    
  subroutine pysave_realdp_7d(variable, filename)
    character(len=*), intent(in) :: filename
    real(fdp), intent(in) :: variable(:,:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(7), i1, i2, i3, i4, i5, i6

    dims = shape(variable)
    write(FMT, *) dims(7)
    if (dims(7) .eq. 0) return

    call file_open(filename, len(filename), 'real')
    write(fileunit, "(A, 7i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            do i5=1, dims(5)
              write(fileunit, "(A, 7i15)") "## ", i1, i2, i3, i4, i5, 0, 0
              do i6=1, dims(6)
                write(fileunit, *) variable(i1, i2, i3, i4, i5, i6, :)
              end do
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_realdp_7d
    
  subroutine pysave_integer_0d(variable, filename)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: variable
    character(20) :: FMT
    
    write(FMT, *) 1

    call file_open(filename, len(filename), 'integer')
    write(fileunit, '('// adjustl(FMT) // 'i25)') variable
    call file_close()
  end subroutine pysave_integer_0d
    
  subroutine pysave_integer_1d(variable, filename)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: variable(:)
    character(20) :: FMT
    integer :: dims(1)

    dims = shape(variable)
    write(FMT, *) dims(1)
    if (dims(1) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, '('// adjustl(FMT) // 'i25)') variable(:)
    call file_close()
  end subroutine pysave_integer_1d
    
  subroutine pysave_integer_2d(variable, filename)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: variable(:,:)
    character(20) :: FMT
    integer :: dims(2), i1

    dims = shape(variable)
    write(FMT, *) dims(2)
    if (dims(2) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    do i1=1, dims(1)
      write(fileunit, '('// adjustl(FMT) // 'i25)') variable(i1, :)
    end do
    call file_close()
  end subroutine pysave_integer_2d
    
  subroutine pysave_integer_3d(variable, filename)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: variable(:,:,:)
    character(20) :: FMT
    integer :: dims(3), i1, i2

    dims = shape(variable)
    write(FMT, *) dims(3)
    if (dims(3) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 3i15)") "## ", dims
    do i1=1, dims(1)
      write(fileunit, "(A, 3i15)") "## ", i1, 0, 0
      do i2=1, dims(2)
        write(fileunit, '('// adjustl(FMT) // 'i25)') variable(i1, i2, :)
      end do
    end do
    call file_close()
  end subroutine pysave_integer_3d
    
  subroutine pysave_integer_4d(variable, filename)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: variable(:,:,:,:)
    character(20) :: FMT
    integer :: dims(4), i1, i2, i3

    dims = shape(variable)
    write(FMT, *) dims(4)
    if (dims(4) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 4i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        write(fileunit, "(A, 4i15)") "## ", i1, i2, 0, 0
        do i3=1, dims(3)
          write(fileunit, '('// adjustl(FMT) // 'i25)') variable(i1, i2, i3, :)
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_integer_4d
    
  subroutine pysave_integer_5d(variable, filename)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: variable(:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(5), i1, i2, i3, i4

    dims = shape(variable)
    write(FMT, *) dims(5)
    if (dims(5) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 5i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          write(fileunit, "(A, 5i15)") "## ", i1, i2, i3, 0, 0
          do i4=1, dims(4)
            write(fileunit, '('// adjustl(FMT) // 'i25)') variable(i1, i2, i3, i4, :)
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_integer_5d
    
  subroutine pysave_integer_6d(variable, filename)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: variable(:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(6), i1, i2, i3, i4, i5

    dims = shape(variable)
    write(FMT, *) dims(6)
    if (dims(6) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 6i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            write(fileunit, "(A, 6i15)") "## ", i1, i2, i3, i4, 0, 0
            do i5=1, dims(5)
              write(fileunit, '('// adjustl(FMT) // 'i25)') variable(i1, i2, i3, i4, i5, :)
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_integer_6d
    
  subroutine pysave_integer_7d(variable, filename)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: variable(:,:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(7), i1, i2, i3, i4, i5, i6

    dims = shape(variable)
    write(FMT, *) dims(7)
    if (dims(7) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 7i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            do i5=1, dims(5)
              write(fileunit, "(A, 7i15)") "## ", i1, i2, i3, i4, i5, 0, 0
              do i6=1, dims(6)
                write(fileunit, '('// adjustl(FMT) // 'i25)') variable(i1, i2, i3, i4, i5, i6, :)
              end do
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_integer_7d
      subroutine pysave_integersp_0d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fsi), intent(in) :: variable
    character(20) :: FMT
    
    write(FMT, *) 1

    call file_open(filename, len(filename), 'integer')
    write(fileunit, '('// adjustl(FMT) // 'i5)') variable
    call file_close()
  end subroutine pysave_integersp_0d
    
  subroutine pysave_integersp_1d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fsi), intent(in) :: variable(:)
    character(20) :: FMT
    integer :: dims(1)

    dims = shape(variable)
    write(FMT, *) dims(1)
    if (dims(1) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, '('// adjustl(FMT) // 'i5)') variable(:)
    call file_close()
  end subroutine pysave_integersp_1d
    
  subroutine pysave_integersp_2d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fsi), intent(in) :: variable(:,:)
    character(20) :: FMT
    integer :: dims(2), i1

    dims = shape(variable)
    write(FMT, *) dims(2)
    if (dims(2) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    do i1=1, dims(1)
      write(fileunit, '('// adjustl(FMT) // 'i5)') variable(i1, :)
    end do
    call file_close()
  end subroutine pysave_integersp_2d
    
  subroutine pysave_integersp_3d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fsi), intent(in) :: variable(:,:,:)
    character(20) :: FMT
    integer :: dims(3), i1, i2

    dims = shape(variable)
    write(FMT, *) dims(3)
    if (dims(3) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 3i15)") "## ", dims
    do i1=1, dims(1)
      write(fileunit, "(A, 3i15)") "## ", i1, 0, 0
      do i2=1, dims(2)
        write(fileunit, '('// adjustl(FMT) // 'i5)') variable(i1, i2, :)
      end do
    end do
    call file_close()
  end subroutine pysave_integersp_3d
    
  subroutine pysave_integersp_4d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fsi), intent(in) :: variable(:,:,:,:)
    character(20) :: FMT
    integer :: dims(4), i1, i2, i3

    dims = shape(variable)
    write(FMT, *) dims(4)
    if (dims(4) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 4i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        write(fileunit, "(A, 4i15)") "## ", i1, i2, 0, 0
        do i3=1, dims(3)
          write(fileunit, '('// adjustl(FMT) // 'i5)') variable(i1, i2, i3, :)
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_integersp_4d
    
  subroutine pysave_integersp_5d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fsi), intent(in) :: variable(:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(5), i1, i2, i3, i4

    dims = shape(variable)
    write(FMT, *) dims(5)
    if (dims(5) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 5i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          write(fileunit, "(A, 5i15)") "## ", i1, i2, i3, 0, 0
          do i4=1, dims(4)
            write(fileunit, '('// adjustl(FMT) // 'i5)') variable(i1, i2, i3, i4, :)
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_integersp_5d
    
  subroutine pysave_integersp_6d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fsi), intent(in) :: variable(:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(6), i1, i2, i3, i4, i5

    dims = shape(variable)
    write(FMT, *) dims(6)
    if (dims(6) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 6i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            write(fileunit, "(A, 6i15)") "## ", i1, i2, i3, i4, 0, 0
            do i5=1, dims(5)
              write(fileunit, '('// adjustl(FMT) // 'i5)') variable(i1, i2, i3, i4, i5, :)
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_integersp_6d
    
  subroutine pysave_integersp_7d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fsi), intent(in) :: variable(:,:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(7), i1, i2, i3, i4, i5, i6

    dims = shape(variable)
    write(FMT, *) dims(7)
    if (dims(7) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 7i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            do i5=1, dims(5)
              write(fileunit, "(A, 7i15)") "## ", i1, i2, i3, i4, i5, 0, 0
              do i6=1, dims(6)
                write(fileunit, '('// adjustl(FMT) // 'i5)') variable(i1, i2, i3, i4, i5, i6, :)
              end do
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_integersp_7d
      subroutine pysave_integerdp_0d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fli), intent(in) :: variable
    character(20) :: FMT
    
    write(FMT, *) 1

    call file_open(filename, len(filename), 'integer')
    write(fileunit, '('// adjustl(FMT) // 'i50)') variable
    call file_close()
  end subroutine pysave_integerdp_0d
    
  subroutine pysave_integerdp_1d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fli), intent(in) :: variable(:)
    character(20) :: FMT
    integer :: dims(1)

    dims = shape(variable)
    write(FMT, *) dims(1)
    if (dims(1) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, '('// adjustl(FMT) // 'i50)') variable(:)
    call file_close()
  end subroutine pysave_integerdp_1d
    
  subroutine pysave_integerdp_2d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fli), intent(in) :: variable(:,:)
    character(20) :: FMT
    integer :: dims(2), i1

    dims = shape(variable)
    write(FMT, *) dims(2)
    if (dims(2) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    do i1=1, dims(1)
      write(fileunit, '('// adjustl(FMT) // 'i50)') variable(i1, :)
    end do
    call file_close()
  end subroutine pysave_integerdp_2d
    
  subroutine pysave_integerdp_3d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fli), intent(in) :: variable(:,:,:)
    character(20) :: FMT
    integer :: dims(3), i1, i2

    dims = shape(variable)
    write(FMT, *) dims(3)
    if (dims(3) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 3i15)") "## ", dims
    do i1=1, dims(1)
      write(fileunit, "(A, 3i15)") "## ", i1, 0, 0
      do i2=1, dims(2)
        write(fileunit, '('// adjustl(FMT) // 'i50)') variable(i1, i2, :)
      end do
    end do
    call file_close()
  end subroutine pysave_integerdp_3d
    
  subroutine pysave_integerdp_4d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fli), intent(in) :: variable(:,:,:,:)
    character(20) :: FMT
    integer :: dims(4), i1, i2, i3

    dims = shape(variable)
    write(FMT, *) dims(4)
    if (dims(4) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 4i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        write(fileunit, "(A, 4i15)") "## ", i1, i2, 0, 0
        do i3=1, dims(3)
          write(fileunit, '('// adjustl(FMT) // 'i50)') variable(i1, i2, i3, :)
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_integerdp_4d
    
  subroutine pysave_integerdp_5d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fli), intent(in) :: variable(:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(5), i1, i2, i3, i4

    dims = shape(variable)
    write(FMT, *) dims(5)
    if (dims(5) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 5i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          write(fileunit, "(A, 5i15)") "## ", i1, i2, i3, 0, 0
          do i4=1, dims(4)
            write(fileunit, '('// adjustl(FMT) // 'i50)') variable(i1, i2, i3, i4, :)
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_integerdp_5d
    
  subroutine pysave_integerdp_6d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fli), intent(in) :: variable(:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(6), i1, i2, i3, i4, i5

    dims = shape(variable)
    write(FMT, *) dims(6)
    if (dims(6) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 6i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            write(fileunit, "(A, 6i15)") "## ", i1, i2, i3, i4, 0, 0
            do i5=1, dims(5)
              write(fileunit, '('// adjustl(FMT) // 'i50)') variable(i1, i2, i3, i4, i5, :)
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_integerdp_6d
    
  subroutine pysave_integerdp_7d(variable, filename)
    character(len=*), intent(in) :: filename
    integer(fli), intent(in) :: variable(:,:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(7), i1, i2, i3, i4, i5, i6

    dims = shape(variable)
    write(FMT, *) dims(7)
    if (dims(7) .eq. 0) return

    call file_open(filename, len(filename), 'integer')
    write(fileunit, "(A, 7i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            do i5=1, dims(5)
              write(fileunit, "(A, 7i15)") "## ", i1, i2, i3, i4, i5, 0, 0
              do i6=1, dims(6)
                write(fileunit, '('// adjustl(FMT) // 'i50)') variable(i1, i2, i3, i4, i5, i6, :)
              end do
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_integerdp_7d
    
  subroutine pysave_complexsp_0d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fsp), intent(in) :: variable
    character(20) :: FMT
    
    write(FMT, *) 1

    call file_open(filename, len(filename), 'complex')
    write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable
    call file_close()
  end subroutine pysave_complexsp_0d
    
  subroutine pysave_complexsp_1d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fsp), intent(in) :: variable(:)
    character(20) :: FMT
    integer :: dims(1)

    dims = shape(variable)
    write(FMT, *) dims(1)
    if (dims(1) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(:)
    call file_close()
  end subroutine pysave_complexsp_1d
    
  subroutine pysave_complexsp_2d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fsp), intent(in) :: variable(:,:)
    character(20) :: FMT
    integer :: dims(2), i1

    dims = shape(variable)
    write(FMT, *) dims(2)
    if (dims(2) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    do i1=1, dims(1)
      write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(i1, :)
    end do
    call file_close()
  end subroutine pysave_complexsp_2d
    
  subroutine pysave_complexsp_3d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fsp), intent(in) :: variable(:,:,:)
    character(20) :: FMT
    integer :: dims(3), i1, i2

    dims = shape(variable)
    write(FMT, *) dims(3)
    if (dims(3) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    write(fileunit, "(A, 3i15)") "## ", dims
    do i1=1, dims(1)
      write(fileunit, "(A, 3i15)") "## ", i1, 0, 0
      do i2=1, dims(2)
        write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(i1, i2, :)
      end do
    end do
    call file_close()
  end subroutine pysave_complexsp_3d
    
  subroutine pysave_complexsp_4d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fsp), intent(in) :: variable(:,:,:,:)
    character(20) :: FMT
    integer :: dims(4), i1, i2, i3

    dims = shape(variable)
    write(FMT, *) dims(4)
    if (dims(4) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    write(fileunit, "(A, 4i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        write(fileunit, "(A, 4i15)") "## ", i1, i2, 0, 0
        do i3=1, dims(3)
          write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(i1, i2, i3, :)
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_complexsp_4d
    
  subroutine pysave_complexsp_5d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fsp), intent(in) :: variable(:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(5), i1, i2, i3, i4

    dims = shape(variable)
    write(FMT, *) dims(5)
    if (dims(5) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    write(fileunit, "(A, 5i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          write(fileunit, "(A, 5i15)") "## ", i1, i2, i3, 0, 0
          do i4=1, dims(4)
            write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(i1, i2, i3, i4, :)
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_complexsp_5d
    
  subroutine pysave_complexsp_6d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fsp), intent(in) :: variable(:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(6), i1, i2, i3, i4, i5

    dims = shape(variable)
    write(FMT, *) dims(6)
    if (dims(6) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    write(fileunit, "(A, 6i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            write(fileunit, "(A, 6i15)") "## ", i1, i2, i3, i4, 0, 0
            do i5=1, dims(5)
              write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(i1, i2, i3, i4, i5, :)
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_complexsp_6d
    
  subroutine pysave_complexsp_7d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fsp), intent(in) :: variable(:,:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(7), i1, i2, i3, i4, i5, i6

    dims = shape(variable)
    write(FMT, *) dims(7)
    if (dims(7) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    write(fileunit, "(A, 7i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            do i5=1, dims(5)
              write(fileunit, "(A, 7i15)") "## ", i1, i2, i3, i4, i5, 0, 0
              do i6=1, dims(6)
                write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(i1, i2, i3, i4, i5, i6, :)
              end do
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_complexsp_7d
      subroutine pysave_complexdp_0d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fdp), intent(in) :: variable
    character(20) :: FMT
    
    write(FMT, *) 1

    call file_open(filename, len(filename), 'complex')
    write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable
    call file_close()
  end subroutine pysave_complexdp_0d
    
  subroutine pysave_complexdp_1d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fdp), intent(in) :: variable(:)
    character(20) :: FMT
    integer :: dims(1)

    dims = shape(variable)
    write(FMT, *) dims(1)
    if (dims(1) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(:)
    call file_close()
  end subroutine pysave_complexdp_1d
    
  subroutine pysave_complexdp_2d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fdp), intent(in) :: variable(:,:)
    character(20) :: FMT
    integer :: dims(2), i1

    dims = shape(variable)
    write(FMT, *) dims(2)
    if (dims(2) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    do i1=1, dims(1)
      write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(i1, :)
    end do
    call file_close()
  end subroutine pysave_complexdp_2d
    
  subroutine pysave_complexdp_3d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fdp), intent(in) :: variable(:,:,:)
    character(20) :: FMT
    integer :: dims(3), i1, i2

    dims = shape(variable)
    write(FMT, *) dims(3)
    if (dims(3) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    write(fileunit, "(A, 3i15)") "## ", dims
    do i1=1, dims(1)
      write(fileunit, "(A, 3i15)") "## ", i1, 0, 0
      do i2=1, dims(2)
        write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(i1, i2, :)
      end do
    end do
    call file_close()
  end subroutine pysave_complexdp_3d
    
  subroutine pysave_complexdp_4d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fdp), intent(in) :: variable(:,:,:,:)
    character(20) :: FMT
    integer :: dims(4), i1, i2, i3

    dims = shape(variable)
    write(FMT, *) dims(4)
    if (dims(4) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    write(fileunit, "(A, 4i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        write(fileunit, "(A, 4i15)") "## ", i1, i2, 0, 0
        do i3=1, dims(3)
          write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(i1, i2, i3, :)
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_complexdp_4d
    
  subroutine pysave_complexdp_5d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fdp), intent(in) :: variable(:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(5), i1, i2, i3, i4

    dims = shape(variable)
    write(FMT, *) dims(5)
    if (dims(5) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    write(fileunit, "(A, 5i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          write(fileunit, "(A, 5i15)") "## ", i1, i2, i3, 0, 0
          do i4=1, dims(4)
            write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(i1, i2, i3, i4, :)
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_complexdp_5d
    
  subroutine pysave_complexdp_6d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fdp), intent(in) :: variable(:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(6), i1, i2, i3, i4, i5

    dims = shape(variable)
    write(FMT, *) dims(6)
    if (dims(6) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    write(fileunit, "(A, 6i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            write(fileunit, "(A, 6i15)") "## ", i1, i2, i3, i4, 0, 0
            do i5=1, dims(5)
              write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(i1, i2, i3, i4, i5, :)
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_complexdp_6d
    
  subroutine pysave_complexdp_7d(variable, filename)
    character(len=*), intent(in) :: filename
    complex(fdp), intent(in) :: variable(:,:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(7), i1, i2, i3, i4, i5, i6

    dims = shape(variable)
    write(FMT, *) dims(7)
    if (dims(7) .eq. 0) return

    call file_open(filename, len(filename), 'complex')
    write(fileunit, "(A, 7i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            do i5=1, dims(5)
              write(fileunit, "(A, 7i15)") "## ", i1, i2, i3, i4, i5, 0, 0
              do i6=1, dims(6)
                write(fileunit, '('// adjustl(FMT) // 'e22.12)') variable(i1, i2, i3, i4, i5, i6, :)
              end do
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_complexdp_7d
    
  subroutine pysave_character_0d(variable, filename)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: variable
    character(20) :: FMT
    
    write(FMT, *) 1

    call file_open(filename, len(filename), 'character')
    write(fileunit, *) trim(adjustl(variable))
    call file_close()
  end subroutine pysave_character_0d
    
  subroutine pysave_character_1d(variable, filename)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: variable(:)
    character(20) :: FMT
    integer :: dims(1)

    dims = shape(variable)
    write(FMT, *) dims(1)
    if (dims(1) .eq. 0) return

    call file_open(filename, len(filename), 'character')
    call char_write_trimmed(variable(:))
    call file_close()
  end subroutine pysave_character_1d
    
  subroutine pysave_character_2d(variable, filename)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: variable(:,:)
    character(20) :: FMT
    integer :: dims(2), i1

    dims = shape(variable)
    write(FMT, *) dims(2)
    if (dims(2) .eq. 0) return

    call file_open(filename, len(filename), 'character')
    do i1=1, dims(1)
      call char_write_trimmed(variable(i1, :))
    end do
    call file_close()
  end subroutine pysave_character_2d
    
  subroutine pysave_character_3d(variable, filename)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: variable(:,:,:)
    character(20) :: FMT
    integer :: dims(3), i1, i2

    dims = shape(variable)
    write(FMT, *) dims(3)
    if (dims(3) .eq. 0) return

    call file_open(filename, len(filename), 'character')
    write(fileunit, "(A, 3i15)") "## ", dims
    do i1=1, dims(1)
      write(fileunit, "(A, 3i15)") "## ", i1, 0, 0
      do i2=1, dims(2)
        call char_write_trimmed(variable(i1, i2, :))
      end do
    end do
    call file_close()
  end subroutine pysave_character_3d
    
  subroutine pysave_character_4d(variable, filename)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: variable(:,:,:,:)
    character(20) :: FMT
    integer :: dims(4), i1, i2, i3

    dims = shape(variable)
    write(FMT, *) dims(4)
    if (dims(4) .eq. 0) return

    call file_open(filename, len(filename), 'character')
    write(fileunit, "(A, 4i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        write(fileunit, "(A, 4i15)") "## ", i1, i2, 0, 0
        do i3=1, dims(3)
          call char_write_trimmed(variable(i1, i2, i3, :))
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_character_4d
    
  subroutine pysave_character_5d(variable, filename)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: variable(:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(5), i1, i2, i3, i4

    dims = shape(variable)
    write(FMT, *) dims(5)
    if (dims(5) .eq. 0) return

    call file_open(filename, len(filename), 'character')
    write(fileunit, "(A, 5i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          write(fileunit, "(A, 5i15)") "## ", i1, i2, i3, 0, 0
          do i4=1, dims(4)
            call char_write_trimmed(variable(i1, i2, i3, i4, :))
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_character_5d
    
  subroutine pysave_character_6d(variable, filename)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: variable(:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(6), i1, i2, i3, i4, i5

    dims = shape(variable)
    write(FMT, *) dims(6)
    if (dims(6) .eq. 0) return

    call file_open(filename, len(filename), 'character')
    write(fileunit, "(A, 6i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            write(fileunit, "(A, 6i15)") "## ", i1, i2, i3, i4, 0, 0
            do i5=1, dims(5)
              call char_write_trimmed(variable(i1, i2, i3, i4, i5, :))
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_character_6d
    
  subroutine pysave_character_7d(variable, filename)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: variable(:,:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(7), i1, i2, i3, i4, i5, i6

    dims = shape(variable)
    write(FMT, *) dims(7)
    if (dims(7) .eq. 0) return

    call file_open(filename, len(filename), 'character')
    write(fileunit, "(A, 7i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            do i5=1, dims(5)
              write(fileunit, "(A, 7i15)") "## ", i1, i2, i3, i4, i5, 0, 0
              do i6=1, dims(6)
                call char_write_trimmed(variable(i1, i2, i3, i4, i5, i6, :))
              end do
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_character_7d
    
  subroutine pysave_logical_0d(variable, filename)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: variable
    character(20) :: FMT
    
    write(FMT, *) 1

    call file_open(filename, len(filename), 'logical')
    write(fileunit, *) variable
    call file_close()
  end subroutine pysave_logical_0d
    
  subroutine pysave_logical_1d(variable, filename)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: variable(:)
    character(20) :: FMT
    integer :: dims(1)

    dims = shape(variable)
    write(FMT, *) dims(1)
    if (dims(1) .eq. 0) return

    call file_open(filename, len(filename), 'logical')
    write(fileunit, *) variable(:)
    call file_close()
  end subroutine pysave_logical_1d
    
  subroutine pysave_logical_2d(variable, filename)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: variable(:,:)
    character(20) :: FMT
    integer :: dims(2), i1

    dims = shape(variable)
    write(FMT, *) dims(2)
    if (dims(2) .eq. 0) return

    call file_open(filename, len(filename), 'logical')
    do i1=1, dims(1)
      write(fileunit, *) variable(i1, :)
    end do
    call file_close()
  end subroutine pysave_logical_2d
    
  subroutine pysave_logical_3d(variable, filename)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: variable(:,:,:)
    character(20) :: FMT
    integer :: dims(3), i1, i2

    dims = shape(variable)
    write(FMT, *) dims(3)
    if (dims(3) .eq. 0) return

    call file_open(filename, len(filename), 'logical')
    write(fileunit, "(A, 3i15)") "## ", dims
    do i1=1, dims(1)
      write(fileunit, "(A, 3i15)") "## ", i1, 0, 0
      do i2=1, dims(2)
        write(fileunit, *) variable(i1, i2, :)
      end do
    end do
    call file_close()
  end subroutine pysave_logical_3d
    
  subroutine pysave_logical_4d(variable, filename)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: variable(:,:,:,:)
    character(20) :: FMT
    integer :: dims(4), i1, i2, i3

    dims = shape(variable)
    write(FMT, *) dims(4)
    if (dims(4) .eq. 0) return

    call file_open(filename, len(filename), 'logical')
    write(fileunit, "(A, 4i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        write(fileunit, "(A, 4i15)") "## ", i1, i2, 0, 0
        do i3=1, dims(3)
          write(fileunit, *) variable(i1, i2, i3, :)
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_logical_4d
    
  subroutine pysave_logical_5d(variable, filename)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: variable(:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(5), i1, i2, i3, i4

    dims = shape(variable)
    write(FMT, *) dims(5)
    if (dims(5) .eq. 0) return

    call file_open(filename, len(filename), 'logical')
    write(fileunit, "(A, 5i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          write(fileunit, "(A, 5i15)") "## ", i1, i2, i3, 0, 0
          do i4=1, dims(4)
            write(fileunit, *) variable(i1, i2, i3, i4, :)
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_logical_5d
    
  subroutine pysave_logical_6d(variable, filename)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: variable(:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(6), i1, i2, i3, i4, i5

    dims = shape(variable)
    write(FMT, *) dims(6)
    if (dims(6) .eq. 0) return

    call file_open(filename, len(filename), 'logical')
    write(fileunit, "(A, 6i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            write(fileunit, "(A, 6i15)") "## ", i1, i2, i3, i4, 0, 0
            do i5=1, dims(5)
              write(fileunit, *) variable(i1, i2, i3, i4, i5, :)
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_logical_6d
    
  subroutine pysave_logical_7d(variable, filename)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: variable(:,:,:,:,:,:,:)
    character(20) :: FMT
    integer :: dims(7), i1, i2, i3, i4, i5, i6

    dims = shape(variable)
    write(FMT, *) dims(7)
    if (dims(7) .eq. 0) return

    call file_open(filename, len(filename), 'logical')
    write(fileunit, "(A, 7i15)") "## ", dims
    do i1=1, dims(1)
      do i2=1, dims(2)
        do i3=1, dims(3)
          do i4=1, dims(4)
            do i5=1, dims(5)
              write(fileunit, "(A, 7i15)") "## ", i1, i2, i3, i4, i5, 0, 0
              do i6=1, dims(6)
                write(fileunit, *) variable(i1, i2, i3, i4, i5, i6, :)
              end do
            end do
          end do
        end do
      end do
    end do
    call file_close()
  end subroutine pysave_logical_7d
    
  
  subroutine file_open(filename, n, template_name)
    integer, intent(in) :: n
    character(n), intent(in) :: filename    
    character(*), intent(in) :: template_name
    character(len=:), allocatable :: delim
    integer :: ioerr
    logical :: exist

    if (template_name .eq. "character") then
       delim = "APOSTROPHE"
    else
       delim = "NONE"
    end if
    
    inquire(file=filename, exist=exist)
    if (exist) then
       open(fpy_newunit(fileunit), file=filename, status="old", position="append",  &
            action="write", iostat=ioerr)
    else
       open(fpy_newunit(fileunit), file=filename, status="new", action="write", iostat=ioerr)
       if (ioerr == 0) write(fileunit, "(A)") '# <fortpy version="1" template="' // template_name // '"></fortpy>'
    end if

    if (ioerr /= 0) then
       print *, "ERROR opening file ", filename, " for pysave in fortpy", ioerr
    end if
  end subroutine file_open

  subroutine file_close()
    !This is just a one line routine, but if we decide to add cleanup later
    !it makes it easier to do it.
    close(fileunit)
  end subroutine file_close

  !!<summary>Returns the number of values in the specified line assuming
  !!that they are separated by spaces or tabs.</summary>
  !!<parameter name="length">The number of characters in line.</parameter>
  !!<parameter name="line">The string for the line to count values in.</parameter>
  !!<parameter name="ischar">When true, the type of data being read in is character
  !!so that whitespace is not the separator, but rather isolated apostrophes or quotes.</parameter>
  integer function fpy_value_count(line, length, ischar)
    integer, intent(in) :: length
    character(length), intent(in) :: line
    logical, optional, intent(in) :: ischar
    
    character(2) :: whitespace
    integer :: i, ichar, indx, prev, cindx(2)
    logical :: ischar_, isquote(2), fquote, fquote_set, qopen
    character :: cchar(2), quote='"', apost="'"

    if (present(ischar)) then
       ischar_ = ischar
    else
       ischar_ = .false.
    end if

    !Initialize the whitespace array. We will cycle through all the characters
    !in the specified line looking for whitespace. Each time we find it, if the
    !character immediately preceding it was not whitespace, we have a value.
    whitespace = ' ' // char(9)
    fquote_set = .false.
    qopen = .false.

    fpy_value_count = 0
    prev = -1
    ichar = 1

    do i = 1, length
       !indx will be zero if the current character is not a whitespace character.
       cchar(ichar) = line(i:i)
       if (ischar_) then
          !We need to identify whether the current character is a quote or apostrophe.
          !We also need to remember which one it is to handle the escaping of '' or "".
          if (cchar(ichar) .eq. apost) then
             isquote(ichar) = .false.
             cindx(ichar) = 1
             if (.not. fquote_set) then
                fquote = .false.
                fquote_set = .true.
             end if
          else
             if (cchar(ichar) .eq. quote) then
                isquote(ichar) = .true.
                cindx(ichar) = 1
                if (.not. fquote_set) then
                   fquote = .true.
                   fquote_set = .true.
                end if
             end if
          end if
          
          !Next, we can analyze whether this character and the last are both the same
          !character. If they are, it is as good as a non-match. If they aren't we
          !check to see if either is a quote and whether they match the first type of
          !quote in the line.
          if (cindx(1) .gt. 0) then
             !First, make sure that the previous character isn't the same since that would
             !be an escaped quote.
             if (ichar .gt. 1) then
                if (cchar(ichar) .eq. cchar(ichar-1)) then
                   !Reset the sequencer so that triple quotes count as an escaped quote
                   !followed by a non-escaped one.
                   indx = 0
                else
                   !If the previous one was a quote and this one is not, then we may have
                   !a match; in that case, make sure the quote matches the first one in
                   !the line.
                   if (cindx(1) .eq. 1 .and. cindx(2) .eq. 0) then
                      if (isquote(1) .eqv. fquote) then
                         indx = 1
                         qopen = (.not. qopen)
                      else
                         indx = 0
                      end if
                   end if
                end if
                ichar = 1
                cindx = 0
                cchar = ""
             else
                !we can't make decisions using an isolated character because of the
                !double character escaping standard in Fortran.
                ichar = ichar + 1
             end if
          else
             !Since neither was a quote, it is the equivalent of matching a non-whitespace
             !character type for the numeric counters.
             indx = 0
             cindx(1) = 0
             cchar(1) = ""
          end if
       else
          indx = index(whitespace, cchar(ichar))
       end if
       
       if (indx > 0 .and. prev == 0) then
          !We found the first whitespace/quote after the end of a value we want.
          if (ischar_) then
             if (qopen) fpy_value_count = fpy_value_count + 1
          else
             fpy_value_count = fpy_value_count + 1
          end if
       end if

       prev = indx
    end do

    if (ischar_) then
       !Check if the last character was a quote and of the same variety as the first
       !quote character on the line.
       if (ichar .eq. 1 .and. cindx(1) .eq. 1 .and. isquote(1) .eqv. fquote .and. qopen) &
            fpy_value_count = fpy_value_count + 1
    else
       !If the last value on the line ends right before \n, then we wouldn't have
       !picked it up; add an extra one.
       if (indx == 0) fpy_value_count = fpy_value_count + 1
    end if
  end function fpy_value_count

  !!<summary>Gets the lengths of ragged-array structured lines for each line
  !!in the specified data file.</summary>
  subroutine fpy_linevalue_count_all(filename, commentchar, nlines, nvalues, ischar)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    integer, intent(out) :: nlines
    integer, allocatable, intent(out) :: nvalues(:)
    logical, optional :: ischar

    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit, i, firstnval
    character(250000) :: line
    logical :: ischar_

    if (present(ischar)) then
       ischar_ = ischar
    else
       ischar_ = .false.
    end if

    !Initialize the value for the result; if we get an error during the read, we
    !end the loop. It can be caused by badly formatted data or the EOF marker.
    call fpy_linevalue_count(filename, commentchar, nlines, firstnval)
    allocate(nvalues(nlines))
    nvalues = 0
    i = 0

    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
       do
          read(funit, "(A)", iostat=ioerr) line
          if (ioerr == 0) then
             cleaned = trim(adjustl(line))
             if (len(cleaned) .gt. 0) then
                if (cleaned(1:1) /= commentchar) then
                   i = i + 1
                   nvalues(i) = fpy_value_count(cleaned, len(cleaned), ischar_)
                end if
             end if
          else
             exit
          end if
       end do
    end if
    close(funit)
  end subroutine fpy_linevalue_count_all

  !!<summary>Returns the number of lines in the file that aren't comments and
  !!the number of whitespace-separated values on the first non-comment line.</summary>
  !!<parameter name="filename">The name of the file to pass to open.</parameter>
  !!<parameter name="n">The number of characters in 'filename'.</parameter>
  !!<parameter name="commentchar">A single character which, when present at the start
  !!of a line designates it as a comment.</parameter>
  subroutine fpy_linevalue_count(filename, commentchar, nlines, nvalues, ischar)
    character(len=*), intent(in) :: filename
    character(1), intent(in) :: commentchar
    integer, intent(out) :: nlines, nvalues
    logical, optional :: ischar
    
    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    character(250000) :: line
    logical :: ischar_, exists

    if (present(ischar)) then
       ischar_ = ischar
    else
       ischar_ = .false.
    end if

    !Initialize the value for the result; if we get an error during the read, we
    !end the loop. It can be caused by badly formatted data or the EOF marker.
    nlines = 0
    nvalues = 0
    inquire(file=filename, exist=exists)
    if (.not. exists) then
       write(*,*) "The file ", filename, " does not exist."
    end if
    
    open(fpy_newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
       do
          read(funit, "(A)", iostat=ioerr) line
          if (ioerr == 0) then
             cleaned = trim(adjustl(line))
             if (abs(len(cleaned) - len(line)) < 10) write (*,*) "Number of characters in line ", len(cleaned), &
                  & " likely exceeds the hard-coded limit of ", len(line), " in file '", filename, "'."

             if (len(cleaned) .gt. 0) then
                if (cleaned(1:1) /= commentchar) then
                   nlines = nlines + 1
                   !We only need to get the number of values present in a line once.
                   !We restrict the file structure to have rectangular arrays.
                   if (nvalues == 0) then
                      nvalues = fpy_value_count(cleaned, len(cleaned), ischar_)
                   end if
                end if
             end if
          else
             if (ioerr .ne. -1) then
                write(*,*) "IO error (", ioerr, ") counting file lines and values. Found ", &
                     & nlines, " and ", nvalues, " in '", filename, "'."
             end if
             exit
          end if
       end do
    end if
    if (fpy_verbose > 0) write (*,*) "Found ", nlines, " lines and ", nvalues, " values in '", filename, "'."
    close(funit)
  end subroutine fpy_linevalue_count

  !!<summary>Returns lowest i/o unit number not in use.</summary>
  !!<parameter name="unit">Out parameter that will contain the lowest i/o number.</parameter>
  integer function fpy_newunit(unit)
    integer, intent(out), optional :: unit
    logical inuse
    integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
    integer, parameter :: nmax=999  ! may be system-dependent
    integer :: n

    do n = nmin, nmax
       inquire(unit=n, opened=inuse)
       if (.not. inuse) then
          if (present(unit)) unit=n
          fpy_newunit = n
          return
       end if
    end do
    stop "newunit ERROR: available unit not found."
  end function fpy_newunit
end module fortpy
