!>  \brief  Useful functions which are not methods and don't belong anywhere else
module UsefulFunctions
    use Globals
    implicit none

    interface FileSize
        module procedure FileSizeFromFilename
        module procedure FileSizeFromUnitNumber
    end interface

    interface RadiansToDegrees
        module procedure DegreesSingle
        module procedure DegreesDouble
    end interface

    interface DegreesToRadians
        module procedure RadiansSingle
        module procedure RadiansDouble
    end interface

    interface QuickSort
        module procedure QuickSortDouble
        module procedure QuickSortSingle
    end interface

    interface FileCopy
        module procedure FileCopyRaw
        !module procedure FileCopySystem
    end interface

    contains

        !>  Convert degrees to radians
        elemental pure real function RadiansSingle(degrees)
            real, intent(in) :: degrees
            RadiansSingle = degrees * pi / 180.
        end function RadiansSingle

        !>  convert radians to degrees
        elemental pure real function DegreesSingle(radians)
            real, intent(in) :: radians
            DegreesSingle = radians / pi * 180.
        end function DegreesSingle

        !>  convert degrees to radians (double precision)
        pure elemental function RadiansDouble(degrees)
            real(kind=8), intent(in) :: degrees
            real(kind=8)    ::  RadiansDouble
            RadiansDouble = degrees * dpi / 180.
        end function RadiansDouble

        !>  convert radians to degrees (double precision)
        pure elemental function DegreesDouble(radians)
            real(kind=8), intent(in) :: radians
            real(kind=8)    ::  DegreesDouble
            DegreesDouble = radians / dpi * 180.
        end function DegreesDouble

        pure elemental function DegreesBetween0And360(degrees)
            real(kind=8),   intent(in)  ::  degrees
            real(kind=8)                ::  DegreesBetween0And360
            DegreesBetween0And360 = mod(degrees,360.0d0)
            if (DegreesBetween0And360 .lt. 0.0d0) DegreesBetween0And360 = DegreesBetween0And360 + 360.0d0
        end function DegreesBetween0And360

        !>  \brief  Normalised cross-correlation between two arrays of reals
        !!  \todo   Optimize  - possibly replace with library function from say BLAS or MKL
        pure function NormalizedCrossCorrelation(array_1,array_2) result(cc)
            real(kind=8),   intent(in)  ::  array_1(:)
            real(kind=8),   intent(in)  ::  array_2(:)
            ! result
            real(kind=8)                ::  cc
            ! private variables
            real(kind=8)                ::  array_1_mean, array_2_mean
            real(kind=8)                ::  array_1_sigma, array_2_sigma
            ! start work
            array_1_mean = sum(array_1) / real(size(array_1))
            array_2_mean = sum(array_2) / real(size(array_2))
            array_1_sigma = sum(array_1**2) / real(size(array_1)) - (sum(array_1) / real(size(array_1)))**2
            array_2_sigma = sum(array_2**2) / real(size(array_2)) - (sum(array_2) / real(size(array_2)))**2
            if (array_1_sigma .gt. 0.0d0 .and. array_2_sigma .gt. 0.0d0) then
                cc = sum((array_1-array_1_mean)*(array_2-array_2_mean)) &
                     / sqrt(array_1_sigma*array_2_sigma) / real(size(array_1))
            else
                cc = 0.0d0
            endif
        end function NormalizedCrossCorrelation

        !>  \brief  Rank cross-correlation
        pure function RankCrossCorrelation(array_1,array_2) result(cc)
            use napack_sort2, only : sort2
            ! arguments
            real(kind=8),   intent(in)  ::  array_1(:)
            real(kind=8),   intent(in)  ::  array_2(:)
            ! result
            real(kind=8)                ::  cc
            ! private variables
            real(kind=8)                ::  work_array(size(array_1))
            real(kind=8)                ::  order_1(size(array_1)), order_2(size(array_2))
            real(kind=8)                ::  sorted_array_1(size(array_1)), sorted_array_2(size(array_2))
            ! start work
            sorted_array_1 = array_1
            sorted_array_2 = array_2
            call sort2(sorted_array_1,order_1,work_array,size(array_1))
            call sort2(sorted_array_2,order_2,work_array,size(array_2))

            cc = NormalizedCrossCorrelation(order_1,order_2)
        end function RankCrossCorrelation



        !>  \brief  Check a file exists on disk
        logical function FileExists(filename)
            character(len=*), intent(in)    ::  filename
            inquire(file=trim(adjustl(filename)), exist=FileExists)
        end function FileExists

        !>  \brief  Copy a file by using a system call via a shell. Only support *nix.
!        subroutine FileCopySystem(filename_src,filename_dest)
!            character(len=*),   intent(in)  ::  filename_src,filename_dest
!            ! Private variables
!            integer ::  exit_status
!            integer ::  command_status
!            character(len=:),   allocatable ::  command_return_message
!            character(len=:),   allocatable ::  command_line
!            ! Start work
!
!            command_line = 'cp '//trim(adjustl(filename_src))//' '//trim(adjustl(filename_dest))
!
!            ! First try the *nix command
!            call execute_command_line(command_line, &
!                                      wait=.true., &
!                                      exitstat=exit_status, &
!                                      cmdstat=command_status, &
!                                      cmdmsg=command_return_message &
!                                      )
!            select case (command_status)
!                case (0)
!                    ! Copy was successful
!                case (-1)
!                    write(*,'(a)') 'Copy failed because the processor does not support command line execution'
!                case (-2)
!                    write(*,'(a)') 'Copy failed because the processor does not support asynchronous execution'
!                case default
!                    write(*,'(a,i0)') 'Copy failed with exit status: ', exit_status
!                    write(*,'(a,a)') 'Command line was: ', command_line
!                    write(*,'(a)') 'Error message was: ', command_return_message
!                    call this_program%TerminateWithFatalError('FileCopyShell','Failed file copy')
!            end select
!
!        end subroutine FileCopySystem

        !>  \brief  Copy a file by reading in and writing out its content. Do not use for large files.
        subroutine FileCopyRaw(filename_src,filename_dest,handle_src,handle_dest)
            character(len=*),   intent(in)  ::  filename_src,filename_dest
            integer,  optional, intent(in)  ::  handle_src,handle_dest          !<  If the files are already open with stream access, supply the unit numbers here
            ! Private variables
            integer                         ::  lun_src, lun_dest ! logical unit numbers
            integer                         ::  number_of_bytes
            integer(kind=1), allocatable    ::  bytes(:)
            integer                         ::  io_status
            character(len=512)              ::  io_msg
            ! Start work
            if (.not. FileExists(filename_src)) then
                call this_program%TerminateWithFatalError('FileCopyRaw','Source file does not exist: '//trim(adjustl(filename_src)))
            endif
            ! Open the files
            if (present(handle_src)) then
                lun_src = handle_src
            else
                lun_src = this_program%GetAvailableUnit()
                open(unit=lun_src, file=filename_src, action='read', status='old', access='stream', iostat=io_status, iomsg=io_msg)

                if (io_status .ne. 0) then
                    close(lun_src)
                    call this_program%ReleaseUnit(lun_src)
                    call this_program%TerminateWithFatalError('UsefulFunctions::FileCopy', &
                              'File exists but cannot be opened: '//trim(adjustl(filename_src))//' ; '//trim(io_msg))
                endif
            endif

            if (present(handle_dest)) then
                lun_dest = handle_dest
            else
                lun_dest = this_program%GetAvailableUnit()
                open(unit=lun_dest,file=filename_dest,action='write',status='replace',access='stream', iostat=io_status)

                if (io_status .ne. 0) then
                    close(lun_dest)
                    call this_program%ReleaseUnit(lun_dest)
                endif
            endif

            ! Find out size of file to copy
            number_of_bytes = FileSize(lun_src)
            if (number_of_bytes .gt. 0) then
                ! Read file data in
                allocate(bytes(number_of_bytes))
                read(unit=lun_src,pos=1,iostat=io_status,iomsg=io_msg) bytes
                if (io_status .ne. 0) then
                    write(*,'(a,i0,4a)') '**error(file_copy): io error ', io_status, ' when reading from: ', filename_src, &
                                                                                                    ' ; ', trim(io_msg)
                    call this_program%TerminateWithFatalError('file_copy','Read error')
                endif
                ! Write file data to output
                write(unit=lun_dest,pos=1,iostat=io_status) bytes
                if (io_status .ne. 0) then
                    write(*,'(a,i0,2a)') '**error(file_copy): io error ', io_status, ' when writing to: ', filename_dest
                    call this_program%TerminateWithFatalError('file_copy','Write error')
                endif
            endif

            ! close files..

            if (UnitIsOpen(lun_src)) then
                close(lun_src)
                call this_program%ReleaseUnit(lun_src)
            endif

            if (UnitIsOpen(lun_dest)) then
                close(lun_dest)
                call this_program%ReleaseUnit(lun_dest)
            endif



            ! Deallocate memory
            if (allocated(bytes)) deallocate(bytes)
        end subroutine FileCopyRaw

        !>  \brief  Delete a file which not currently open
        subroutine FileDelete(filename)
            character(len=*),   intent(in)  ::  filename
            ! Private variables
            integer ::  tmp_lun, io_status
            character(len=512)  ::  io_msg
            ! Start work
            if (FileExists(filename)) then
                tmp_lun = this_program%GetAvailableUnit()
                open(unit=tmp_lun,file=filename,action='write',iostat=io_status,iomsg=io_msg)
                if (io_status .ne. 0) then
                    call this_program%TerminateWithFatalError('FileDelete','Could not open file '//trim(filename)//&
                                                                                ': '//trim(io_msg))
                endif
                close(tmp_lun, status='delete',iostat=io_status)
                if (io_status .ne. 0) then
                    write(*,'(a,i0,2a)') '**warning(FileDelete): error ', io_status, ' when trying to delete ', trim(filename)
                endif
                call this_program%ReleaseUnit(tmp_lun)
            else
                write(*,'(2a)') '**warning(FileDelete): attempt to delete file which does not exist: ', trim(filename)
            endif
            if (FileExists(filename)) then
                write(*,'(2a)') '**warning(FileDelete): could not delete ', trim(filename)
            endif
        end subroutine FileDelete

        !> \brief   Find file size in bytes
        function FileSizeFromFilename(filename) result(file_size)
            character(len=*), intent(in)    ::  filename
            integer(kind=8)                 ::  file_size
            inquire(file=trim(adjustl(filename)),size=file_size)
        end function FileSizeFromFilename

        !> \brief   Find file size in bytes
        function FileSizeFromUnitNumber(lun) result(file_size)
            integer, intent(in)    ::  lun
            integer(kind=8)                 ::  file_size
            inquire(unit=lun,size=file_size)
        end function FileSizeFromUnitNumber

        !>  \brief  Check whether a IO unit is current open
        logical function UnitIsOpen(unit_number)
            integer, intent(in)      ::  unit_number
            integer :: io_status
            character(len=100) :: io_message
            io_status = 0
            inquire(unit=unit_number, opened=UnitIsOpen,iostat=io_status,iomsg=io_message)
            if (io_status .ne. 0) then
                print *, 'UnitIsOpen: IO error ', io_status, ': ', trim(adjustl(io_message))
                call this_program%TerminateWithFatalError('UnitIsOpen','IO error: '//trim(adjustl(io_message)))
            endif
        end function UnitIsOpen

        !>  \brief  Check whether a IO unit is current open
        logical function FileIsOpen(filename)
            character(len=*), intent(in)      ::  filename
            integer :: io_status
            character(len=100) :: io_message
            io_status = 0
            inquire(file=filename, opened=FileIsOpen,iostat=io_status,iomsg=io_message)
            if (io_status .ne. 0) then
                print *, 'FileIsOpen: IO error ', io_status, ': ', trim(adjustl(io_message))
                call this_program%TerminateWithFatalError('FileIsOpen','IO error: '//trim(adjustl(io_message)))
            endif
        end function FileIsOpen

        !>    \brief    returns true if the argument is even
        pure elemental logical function IsEven(int)
            integer,    intent(in)    ::    int
            ! test bit 0 of number. if it is 0, the number is even
            IsEven = .not. btest(int,0)
        end function IsEven

        !>    \brief    returns true if the argument is odd
        pure elemental logical function IsOdd(int)
            integer,    intent(in)    ::    int
            ! start work
            IsOdd = btest(int,0)
        end function IsOdd

        !>  \brief  Compute the phase to be applied to an element distance_from_origin away from the origin so that the real-space shift
        !!          real_space_shift is applied
        pure function ReturnPhaseFromShift(real_space_shift, distance_from_origin, dimension_size)

            ! Arguments

            real,           intent(in)      ::  real_space_shift          !<  Real space shift (in pixels)
            integer,        intent(in)      ::  distance_from_origin      !<  Distance from origin in Fourier transform (in pixels) along the axis of the shift
            integer,        intent(in)      ::  dimension_size            !<  The total extent of the image in the dimension relevant to the shift

            ! Result
            real                            ::  ReturnPhaseFromShift
            ! Private variables

            ! Start work
            ReturnPhaseFromShift = real_space_shift * real(distance_from_origin) * 2.0e0 * pi / real(dimension_size)

        end function ReturnPhaseFromShift

        !>  \brief  Compute the phase shift, given the phase in each dimension
        pure complex function Return3DPhaseFromIndividualDimensions(phase_x, phase_y, phase_z)
            !use LookUpTables
            ! Arguments
            real,   intent(in)        ::  phase_x
            real,   intent(in)        ::  phase_y
            real,   intent(in)        ::  phase_z

            ! Private Variables
            real                      ::  temp_phase
            logical,       parameter  ::  use_look_up_tables = .false.

            ! Start work
            temp_phase  = - phase_x - phase_y - phase_z
            !if (use_look_up_tables) then
            !    Return3DPhaseFromIndividualDimensions = cmplx(cos_lookup(temp_phase), sin_lookup(temp_phase))
            !else
                Return3DPhaseFromIndividualDimensions = cmplx(cos(temp_phase), sin(temp_phase))
            !endif

        end function Return3DPhaseFromIndividualDimensions

        !>  \brief  Sort the elements of an array using the quick sort algorithm
        recursive subroutine QuickSortDouble(array)
            real(kind=8),   intent(inout)   ::  array(:)
            ! private variables
            integer ::  pivot_index
            ! start work
            if (size(array) .gt. 1) then
                call PartitionDouble(array,pivot_index)
                call QuickSortDouble(array(:pivot_index-1))
                call QuickSortDouble(array(pivot_index:))
            endif
        end subroutine QuickSortDouble

        recursive subroutine QuickSortSingle(array)
            real(kind=4),   intent(inout)   ::  array(:)
            ! private variables
            integer ::  pivot_index
            ! start work
            if (size(array) .gt. 1) then
                call PartitionSingle(array,pivot_index)
                call QuickSortSingle(array(:pivot_index-1))
                call QuickSortSingle(array(pivot_index:))
            endif
        end subroutine QuickSortSingle

        !>  \brief  Partition an array so that all elements to "the left" and to "the right" are below and above the pivot point
        pure subroutine PartitionDouble(array,pivot_index)
            real(kind=8),   intent(inout)   ::  array(:)
            integer,        intent(out)     ::  pivot_index
            ! private variables
            integer         ::  i,j
            real(kind=8)    ::  temp
            real(kind=8)    ::  pivot_value
            ! start work
            pivot_value = array(1)
            i = 0
            j = size(array) + 1

            do
                j=j-1
                do
                    if(array(j) .le. pivot_value) exit
                    j=j-1
                enddo
                i=i+1
                do
                    if (array(i) .ge. pivot_value) exit
                    i=i+1
                enddo
                if (i .lt. j) then
                    ! exchange
                    temp = array(i)
                    array(i) = array(j)
                    array(j) = temp
                else if (i .eq. j) then
                    pivot_index = i+1
                    return
                else
                    pivot_index = i
                    return
                endif
            enddo
        end subroutine PartitionDouble

        pure subroutine PartitionSingle(array,pivot_index)
            real(kind=4),   intent(inout)   ::  array(:)
            integer,        intent(out)     ::  pivot_index
            ! private variables
            integer         ::  i,j
            real(kind=8)    ::  temp
            real(kind=8)    ::  pivot_value
            ! start work
            pivot_value = array(1)
            i = 0
            j = size(array) + 1

            do
                j=j-1
                do
                    if(array(j) .le. pivot_value) exit
                    j=j-1
                enddo
                i=i+1
                do
                    if (array(i) .ge. pivot_value) exit
                    i=i+1
                enddo
                if (i .lt. j) then
                    ! exchange
                    temp = array(i)
                    array(i) = array(j)
                    array(j) = temp
                else if (i .eq. j) then
                    pivot_index = i+1
                    return
                else
                    pivot_index = i
                    return
                endif
            enddo
        end subroutine PartitionSingle

        !>  \brief  Equivalent to the unix grep command, except it returns only the first line that contains pattern
        function grep(filename,pattern,shift) result(line)
            ! Arguments
            character(len=*),       intent(in)  ::  filename
            character(len=*),       intent(in)  ::  pattern
            integer,    optional,   intent(in)  ::  shift   !<  If shift is n, the line returned will the line n lines after that with the pattern
            ! Result
            character(len=512)  ::  line
            ! Private variables
            character(len=512)  ::  buffer
            integer ::  ios
            logical ::  pattern_found
            integer ::  sshift
            integer ::  shifted
            integer ::  tmp_file_handle
            ! Start work
            ! Shift?
            if (present(shift)) then
                if (shift .lt. 0) then
                    write(*,'(a,i3)') '**ERROR(grep): shift parameter has to be positive ', shift
                    call this_program%TerminateWithFatalError('grep','Fatal error in GREP')
                endif
                sshift = shift
            else
                sshift = 0
            endif
            pattern_found = .false.
            ! Open the flat-text file
            tmp_file_handle = this_program%GetAvailableUnit()
            open(tmp_file_handle,file=FILENAME,status='OLD',iostat=ios)
            if (ios .gt. 0) then
                write(*,'(a,a,a)')  '**ERROR(grep): Error opening file ', trim(filename), '.'
            else if (ios .lt. 0) then
                write(*,'(a,a,a)')  '**ERROR(grep): Error opening file ', trim(filename), '. eof or EOR.'
            endif
            ! Make sure we are at beginning of the file
            rewind (tmp_file_handle,iostat=ios)
            if (ios .ne. 0) then
                write(*,*) '**ERROR(grep): Error when REWINDing file ', trim(filename)
                call this_program%TerminateWithFatalError('grep','Fatal error in GREP')
            endif
            ! Read through the file, line by line
            ios = 0
            shifted = 0
            do while (ios == 0)
                read(tmp_file_handle, '(a)', iostat=ios) buffer
                if (ios .eq. 0) then
                    if (index(buffer,pattern) .ne. 0 .or. pattern_found) then
                        ! Found pattern in buffer
                        line = trim(buffer)
                        pattern_found = .true.
                        if (shifted .eq. sshift) then
                            exit
                        endif
                        shifted = shifted + 1
                    endif
                endif
            enddo
            ! If we get here, it means the pattern was not found, so we return a blank
            if (pattern_found) then

            else
                line = ' '
            endif
            call this_program%ReleaseUnit(tmp_file_handle)
            close(tmp_file_handle)
        end function grep

        pure elemental function IntegerFromUnsignedToSigned(self_unsigned) result(self_signed)
            integer(kind=4),    intent(in)  ::  self_unsigned
            integer(kind=4)                 ::  self_signed
            self_signed = iand(self_unsigned,huge(int(1,kind=4)))
        end function IntegerFromUnsignedToSigned




end module UsefulFunctions
