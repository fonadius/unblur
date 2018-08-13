!>  \brief  Class to deal with text files of numbers
module NumericTextFiles

    use Globals

    implicit none

    integer,            parameter ::  OPEN_TO_READ = 1
    integer,            parameter ::  OPEN_TO_WRITE = 2

    type,         public  ::  NumericTextFile

            integer                        :: file_unit
            character(len=:), allocatable  :: filename
            integer                        :: records_per_line = 0
            integer                        :: number_of_data_lines = 0
            integer                        :: access_type


            contains

            procedure                      :: Init
            procedure                      :: ReadNextDataLine
            procedure                      :: WriteDataLineReal
            procedure                      :: WriteDataLineInteger
            generic, public                :: WriteDataLine => WriteDataLineReal, WriteDataLineInteger
            procedure                      :: WriteCommentLine
            procedure, public               ::  GetNumberOfRecordsPerLine
            procedure, public               ::  GetNumberOfDataLines
            !procedure                      :: Rewind
            !procedure                      :: Flush
            final                          :: Destroy

    end type

    contains

    pure integer function GetNumberOfRecordsPerLine(self)
        class(NumericTextFile), intent(in)  ::  self
        GetNumberOfRecordsPerLine = self%records_per_line
    end function GetNumberOfRecordsPerLine

    pure integer function GetNumberOfDataLines(self)
        class(NumericTextFile), intent(in)  ::  self
        GetNumberOfDataLines = self%number_of_data_lines
    end function GetNumberOfDataLines

    subroutine Init(self, filename, access_type, wanted_records_per_line)

        use StringManipulations
        ! arguments
        class(NumericTextFile), intent(inout)  ::  self
        character(len=*), intent(in)           ::  filename
        integer, intent(in)                    ::  access_type                  !<  Either OPEN_TO_READ or OPEN_TO_WRITE
        integer, optional, intent(in)          ::  wanted_records_per_line

        ! variables

        character(len=line_max_len)            ::  buffer      !< will hold a line from the file
        integer                                ::  records_on_current_line
        integer                                ::  total_number_of_records
        integer                                ::  ios         !< ios is negative if an end of record condition is encountered or if
                                                               !! an endfile condition was detected.  it is positive if an error was
                                                               !! detected.  ios is zero otherwise.

        character(len=512)                      ::  io_msg
        !start

        total_number_of_records = 0

        self%filename = trim(adjustl(filename))
        self%access_type = access_type

        !open the file

        ! if we are opening to read, work out the file details..

        if (self%access_type .eq. OPEN_TO_READ) then

            self%file_unit = this_program%GetAvailableUnit()
            open(unit=self%file_unit, file=self%filename, iostat=ios, status='old')

            if (ios .ne. 0) then
                close(self%file_unit)
                call this_program%ReleaseUnit(self%file_unit)

                call this_program%TerminateWithFatalError('NumericTextFiles::Init', 'Error when opening file for reading: '// &
                                                            trim(adjustl(self%filename)))
            endif

            do
                ! work out records per line, and number_of_lines

                read(self%file_unit, '(a)', iostat=ios) buffer

                if ( ios == 0 ) then
                    if (StringIsAComment(buffer) .or. StringIsBlank(buffer)) then
                        ! don't do anything
                    else
                        self%number_of_data_lines = self%number_of_data_lines + 1

                        ! work out how many records on this line..

                        records_on_current_line = CountRecordsPerLine(buffer)
                        total_number_of_records = total_number_of_records + records_on_current_line
                    endif
                else
                    exit !we found the end..
                endif
            enddo

            ! if this file is correct, then the value left in records_on_current_line
            ! should be the same as the total_number_of_records / number_of_lines.
            ! This is possibly a bit of a nasty way to check things, and may need to
            ! be changed.

            if (records_on_current_line * self%number_of_data_lines .eq. total_number_of_records) then
                !All seems to be ok..

                self%records_per_line = records_on_current_line
                rewind(self%file_unit)
            else
                !Something went wrong..
                call this_program%TerminateWithFatalError('NumericTextFile::Init', &
                                                          'Not all lines contain the same number of records?')
            endif

        else if (self%access_type .eq. OPEN_TO_WRITE) then

            self%file_unit = this_program%GetAvailableUnit()
            open(unit=self%file_unit, file=self%filename, iostat=ios, status='replace',iomsg=io_msg)

            if (ios .ne. 0) then
                close(self%file_unit)
                call this_program%ReleaseUnit(self%file_unit)

                call this_program%TerminateWithFatalError('NumericTextFiles::Init', 'Error when opening file for writing. '// &
                                                            trim(self%filename)//' ; '//trim(io_msg))
            endif

            if (present(wanted_records_per_line)) then
                self%records_per_line = wanted_records_per_line
            else
                self%records_per_line = 1
            endif

            self%number_of_data_lines = 0

        endif

end subroutine Init


subroutine Destroy(self)

    use UsefulFunctions

    type(NumericTextFile), intent(inout)   ::  self

    if (allocated(self%filename)) deallocate(self%filename)

    if (UnitIsOpen(self%file_unit)) then
        close(self%file_unit)
        call this_program%ReleaseUnit(self%file_unit)
    endif

end subroutine Destroy


subroutine ReadNextDataLine(self, read_data)

        use StringManipulations
        ! arguments
        class(NumericTextFile), intent(inout)  ::  self
        real, intent(inout)                    ::  read_data(:)

        !variables

        character(len=line_max_len)            ::  buffer
        !integer                                ::  record_counter
        integer                                ::  ios
        character(len=256)                      ::  io_message


        ! Start

        ! Check we are open to read

        if (self%access_type .ne. OPEN_TO_READ) then
            call this_program%TerminateWithFatalError('NumericTextFiles::ReadNextDataLine','File is not OPEN_TO_READ')
        endif

        ! Check the passed array is big enough

        if (size(read_data) .lt. self%records_per_line) then
            call this_program%TerminateWithFatalError('NumericTextFiles::ReadNextDataLine', &
                                                      'Supplied array is smaller than records per line')
        endif

        ! read the next line data line
        do
            read(self%file_unit, '(a)', iostat=ios, iomsg=io_message) buffer

            if ( ios .ne. 0 ) then
                call this_program%TerminateWithFatalError('NumericTextFiles::ReadNextDataLine', &
                                                          'Encountered iostat error: '//trim(io_message))
            endif

            if (.not. StringIsAComment(buffer) .and. .not. StringIsBlank(buffer)) then
                buffer = trim(adjustl(buffer))

                read(buffer, *) read_data(1:self%records_per_line)
                exit

            endif

        enddo


end subroutine ReadNextDataLine

subroutine WriteDataLineReal(self, data_to_write)

        use StringManipulations
        ! arguments
        class(NumericTextFile), intent(inout)  ::  self
        real, intent(in)                    ::  data_to_write(:)

        !variables

        integer                                ::  record_counter

        ! Start

        ! Check we are open to write

        if (self%access_type .ne. OPEN_TO_WRITE) then
            call this_program%TerminateWithFatalError('NumericTextFiles::WriteDataLineReal', &
                                                      'File is not OPEN_TO_WRITE')
        endif

        ! Check the passed array is big enough

        if (size(data_to_write) .lt. self%records_per_line) then
            call this_program%TerminateWithFatalError('NumericTextFiles::WriteDataLineReal', &
                                                      'Supplied array is smaller than records per line')
        endif

        ! write out the line..

        do record_counter = 1, self%records_per_line
            write(self%file_unit, '(g14.7,a)', advance='no') data_to_write(record_counter), ' '
        enddo

        ! finish the line

        write(self%file_unit,*)

        ! increase the number of data lines

        self%number_of_data_lines = self%number_of_data_lines + 1

end subroutine WriteDataLineReal

subroutine WriteDataLineInteger(self, data_to_write)

        use StringManipulations
        ! arguments
        class(NumericTextFile), intent(inout)  ::  self
        integer, intent(inout)                    ::  data_to_write(:)

        !variables

        integer                                ::  record_counter

        ! Start

        ! Check we are open to write

        if (self%access_type .ne. OPEN_TO_WRITE) then
            call this_program%TerminateWithFatalError('NumericTextFiles::WriteDataLineInteger','File is not OPEN_TO_WRITE')
        endif

        ! Check the passed array is big enough

        if (size(data_to_write) .lt. self%records_per_line) then
            call this_program%TerminateWithFatalError('NumericTextFiles::WriteDataLineInteger', &
                                                      'Supplied array is smaller than records per line')
        endif

        ! write out the line..

        do record_counter = 1, self%records_per_line
            write(self%file_unit, '(g14.7,a)', advance='no') real(data_to_write(record_counter)), ' '
        enddo

        ! finish the line

        write(self%file_unit,*)

        ! increase the number of data lines

        self%number_of_data_lines = self%number_of_data_lines + 1

end subroutine WriteDataLineInteger


subroutine WriteCommentLine(self, comment_to_write)

        use StringManipulations
        ! arguments
        class(NumericTextFile), intent(inout)  ::  self
        character(len=*), intent(in)           ::  comment_to_write

        ! Start

        ! Check we are open to write

        if (self%access_type .ne. OPEN_TO_WRITE) then
            call this_program%TerminateWithFatalError('NumericTextFiles::WriteCommentLine','File is not OPEN_TO_WRITE')
        endif

        ! write out the line..

        write(self%file_unit, '(2a)') '# ', trim(adjustl(comment_to_write))

end subroutine WriteCommentLine







end module


