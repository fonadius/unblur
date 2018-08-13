!>  \brief  Set of routines to manipulate strings of characters
module StringManipulations
    use iso_c_binding
    use Globals

    implicit none
    public

    ! Useful constants
    character(len=*),               private,    parameter :: LOWER_CASE_LETTERS = 'abcdefghijklmnopqrstuvwxyz'
    character(len=*),               private,    parameter :: UPPER_CASE_LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len=*),               private,    parameter :: INTEGERS = '1234567890'
    character(kind=c_char,len=*),   private,    parameter :: NEW_LINES_C = C_FORM_FEED // C_NEW_LINE // C_CARRIAGE_RETURN &
                                                                                   // C_VERTICAL_TAB
    character(kind=c_char,len=*),   private,    parameter :: BLANK_C_CHARACTERS = C_NULL_CHAR // C_HORIZONTAL_TAB
    character(len=*),               private,    parameter :: BLANK_CHARACTERS =  ' '//BLANK_C_CHARACTERS
    character(len=*),               private,    parameter :: BLANKS_AND_NEW_LINES = BLANK_CHARACTERS // NEW_LINES_C

    character(len=*),               private,    parameter :: INTEGER_DIGITS = '10'   ! Maximum number of digits expected when reading an integer


    contains
        !>  \brief  Convert string to all upper case
        !!  Adapted from
        !!  Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
        !!  1995 Springer-Verlag, New York.
        function UpperCase(string)

            ! arguments
            character(len=*),   intent(in)  ::  string
            ! result
            character(len=len(string))  ::  UpperCase
            ! private variables
            integer ::  i,n
            ! start work
            UpperCase = string
            do i=1,len(string)
                n = index(lower_case_letters,UpperCase(i:i))
                if (n .ne. 0) UpperCase(i:i) = upper_case_letters(n:n)
            enddo
        end function UpperCase

        !>  \brief  Generate a random filename
        function RandomFilename(number_of_characters,prefix)
            integer,                intent(in)  ::  number_of_characters
            character(len=*),       intent(in)  ::  prefix
            character(len=number_of_characters+len(prefix)) ::  RandomFilename
            !
            integer     ::  character_counter
            real        ::  rnd
            integer     ::  irnd
            integer,    parameter   ::  ascii_code_of_first_digit       =   48
            integer,    parameter   ::  ascii_code_of_first_upper_alpha =   65
            integer,    parameter   ::  ascii_code_of_first_lower_alpha =   97
            !

            ! Start with the prefix
            RandomFilename(1:len(prefix)) = prefix

            call random_seed()
            do character_counter=len(prefix)+1,len(prefix)+number_of_characters
                call random_number(rnd)
                irnd = min(int(rnd * 62.0), 61) + 1
                if (irnd .le. 10) then
                    ! Digit
                    RandomFilename(character_counter:character_counter) = achar(irnd + ascii_code_of_first_digit       - 1)
                else if (irnd .le. 36) then
                    ! Upper case
                    RandomFilename(character_counter:character_counter) = achar(irnd + ascii_code_of_first_upper_alpha - 11)
                else
                    ! Lower case
                    RandomFilename(character_counter:character_counter) = achar(irnd + ascii_code_of_first_lower_alpha - 37)
                endif
            enddo
        end function RandomFilename

        !>  \brief  Return the 3-letter extension of a filename if present (without the period)
        pure function ExtensionFromFilename(filename)

            ! arguments
            character(len=*),   intent(in)  ::  filename
            ! return value
            character(len=3)    ::  ExtensionFromFilename
            ! private variables
            integer ::  length
            ! start work - we look for a simple extension, where the filename ends with .???
            length = len_trim(filename)
            if (scan(filename(1:length),'.',back=.true.) .eq. length-3) then
                ExtensionFromFilename = filename(length-2:length)
            else
                ExtensionFromFilename = '   '
            endif
        end function ExtensionFromFilename

        !>  \brief  Add a suffix to a filename, before its extension if it has one
        subroutine FilenameAddSuffix(filename,suffix)
            ! Arguments
            character(len=*),   intent(inout)   ::  filename    !<  Filename to be modified
            character(len=*),   intent(in)      ::  suffix          !<  a string to be added onto the filename (e.g. '_fsc1')
            ! Private variables
            character(len=:),   allocatable     ::  filename_ext, filename_tmp
            ! Start work
            filename_ext = ExtensionFromFilename(filename)
            call FilenameRemoveExtension(filename,filename_tmp)
            filename = trim(adjustl(filename_tmp)) // trim(adjustl(suffix)) // '.' // trim(adjustl(filename_ext))
        end subroutine FilenameAddSuffix

        !>  \brief  Replace the filename extension
        !!
        !!  Note, if the new extension is longer than the current extension, behaviour is not defined.
        subroutine FilenameReplaceExtension(filename, new_extension)

            implicit none
            ! arguments
            character(len=*),   intent(inout)   ::  filename        !<  filename to be modified
            character(len=*),   intent(in)      ::  new_extension   !<  the filename's extension will be replaced by this
            ! private variables
            character(len=:),   allocatable     ::  nnew_extension  ! cleaner version of the new extension, without the first dot, and without trailing blanks
            character(len=:),   allocatable     ::  temp
            character(len=:),   allocatable     ::  temp2
            ! start work

            !
            ! clean up the new extension
            !
            temp = trim(adjustl(new_extension))
            if (temp(1:1) .eq. '.') then
                nnew_extension = temp(2:len(temp))
            else
                nnew_extension = temp
            endif

            !
            ! if the new extension is longer than 3 characters, this won't work
            !
            if (len(nnew_extension) .gt. 3) then
                  call this_program%TerminateWithFatalError('UsefulFunctions::FilenameReplaceExtension', &
                                                 'can only deal with 3-character-long extensions for now')
            endif

            !
            ! remove the extension from the filename
            !
            call FilenameRemoveExtension(filename,temp2)

            !
            ! add the new extension to the filename
            !
            filename = trim(adjustl(temp2)) // '.' // nnew_extension

        end subroutine FilenameReplaceExtension

        !>  \brief  Return the filename without the .??? extension, if present
        subroutine FilenameRemoveExtension(filename_in, filename_out)

            ! arguments
            character(len=*),   intent(in)                ::  filename_in     !<  filename to be modified
            character(len=:),   allocatable,  intent(out) ::  filename_out    !<  modified filename

            ! private variables
            integer ::  length
            character(len=:),   allocatable               ::  buffer_string
            ! start work

            if (ExtensionFromFilename(filename_in) .eq. '   ') then
                ! input filename does not have an extension, nothing to do so just copy

                call CopyString(filename_in, filename_out)

                else

                call CopyString(trim(adjustl(filename_in)), buffer_string)

                length = len(buffer_string)
                call CopyString(buffer_string(1:length-4), filename_out)

                deallocate(buffer_string)
            endif
        end subroutine FilenameRemoveExtension

        logical function StringIsAComment(line)
            ! arguments
            character(len=*),   intent(in)  ::  line
            ! private variables
            integer ::  pos1
            ! find the first non-blank character
            pos1    =   FirstNonBlank(line)

            ! if we didn't find a non-blank character, then the line must be blank!
            if (pos1 .eq. 0) then
                StringIsAComment = .false.
            else if (scan(line(pos1:),'#!c;',back=.false.) .eq. 1) then
                ! the first non-blank character is a comment character. if that's a c, we need to check that the following character is a space.
                if (StringsAreEqual(line(pos1:pos1),'c')) then
                    if (StringsAreEqual(line(pos1+1:pos1+1),' ')) then
                        StringIsAComment = .true.
                    else
                        StringIsAComment = .false.
                    endif
                else
                    StringIsAComment = .true.
                endif
            else
            StringIsAComment = .false.
            endif
        end function StringIsAComment

            !>  \brief Works out whether a character string is blank
        pure logical function StringIsBlank(line)
            ! argument
            character(len=*),   intent(in)  ::  line
            ! start work
            if (trim(line) .eq. '') then
                 StringIsBlank = .true.
            else if (FirstNonBlank(line) .eq. 0) then
                 StringIsBlank = .true.
            else
                 StringIsBlank = .false.
            endif
        end function StringIsBlank

        !>  \brief  Find the first non-blank character in a string and return its position.
        pure integer function FirstNonBlank(string, back)
            implicit none
            ! arguments
            character(len=*),       intent(in)  ::  string
            logical,    optional,   intent(in)  ::  back        !<  if .true., we'll look for the last non-blank character in the string
            ! private variables
            logical ::  bback
            ! start work
            ! reverse?
            bback = .false.
            if (present(back)) bback = back
            FirstNonBlank =   verify(string,blanks_and_new_lines,back=bback)
        end function FirstNonBlank

            !>  \brief  find the first blank character in a string and return its position.
         pure integer function FirstBlank(string, back)
            ! arguments
            character(len=*),       intent(in)  ::  string
            logical,    optional,   intent(in)  ::  back        !<  if .true., we'll look for the last blank character in the string
            ! private variables
            logical ::  bback
            ! start work
            ! reverse?
            if (present(back)) then
                bback = back
            else
                bback = .false.
            endif

            FirstBlank =   scan(string,blanks_and_new_lines,back=bback)
        end function FirstBlank

        !>  \brief  Test whether two strings are equal, ignoring blank characters
        logical function StringsAreEqual(input_string1,input_string2,case_sensitive)

            ! arguments
            character(len=*),               intent(in)  ::  input_string1
            character(len=*),               intent(in)  ::  input_string2
            logical,            optional,   intent(in)  ::  case_sensitive
            ! private variables
            integer ::  first_non_blank_1, first_non_blank_2
            integer ::  last_non_blank_1,  last_non_blank_2
            logical ::  ccase_sensitive
            character(len=:), allocatable   ::  local_string1, local_string2
            ! start work

            ! Will we be case sensitive?
            ccase_sensitive = .true.
            if (present(case_sensitive)) ccase_sensitive = case_sensitive

            ! Copy to local string, with or without capitalization
            if (ccase_sensitive) then
                call CopyString(input_string1,local_string1)
                call CopyString(input_string2,local_string2)
            else
                local_string1 = UpperCase(input_string1)
                local_string2 = UpperCase(input_string2)
            endif


            ! Find positions of first and last non-blank characters
            first_non_blank_1 = verify(local_string1,blank_characters)
            last_non_blank_1  = verify(local_string1,blank_characters,back=.true.)
            first_non_blank_2 = verify(local_string2,blank_characters)
            last_non_blank_2  = verify(local_string2,blank_characters,back=.true.)


            if (first_non_blank_1 .eq. 0 .and. first_non_blank_2 .eq. 0) then
                ! both strings are blank
                StringsAreEqual = .true.
            else if (first_non_blank_1 .eq. 0 .or. first_non_blank_2 .eq. 0) then
                ! one string is blank, the other isn't
                StringsAreEqual = .false.
            else if (index(local_string1(first_non_blank_1:last_non_blank_1), &
                           local_string2(first_non_blank_1:last_non_blank_2)) .eq. 1  &
                    .and. last_non_blank_1-first_non_blank_1 .eq. last_non_blank_2-first_non_blank_2) then
                ! neither string is blank, and the strings match
                StringsAreEqual = .true.
            else
                ! all other cases
                StringsAreEqual = .false.
            endif
        end function StringsAreEqual

        !>  \brief  Turn integer variable into character variable
        function IntegerToString(intg) result(string)
            ! Argument
            integer,            intent(in)  ::  intg
            ! Result
            character(len=:),   allocatable ::  string
            ! Start work
            allocate(character(len=NumberOfDigitsInInteger(intg)) :: string)
            write(string,'(i0)') intg
        end function IntegerToString

        !>  \brief  Turn real variable into character variable
        function RealToString(self,number_of_decimals) result(string)
            real,               intent(in)  ::  self
            integer,            intent(in)  ::  number_of_decimals
            !
            character(len=:),   allocatable ::  string
            !
            allocate(character(len=NumberOfDigitsInInteger(int(self))+1+number_of_decimals) :: string)
            write(string,'(f0.'//IntegerToString(number_of_decimals)//')') self
        end function RealToString


        !>  \brief  Work out the number of digits in an integer
        pure integer function NumberOfDigitsInInteger(intg)
            ! Argument
            integer,   intent(in)  ::  intg
            ! Start work
            if (intg .eq. 0) then
                NumberOfDigitsInInteger = 1
            else
                NumberOfDigitsInInteger = int(log10(real(intg))) + 1
            endif
        end function NumberOfDigitsInInteger

        !> \brief Copies the Trimmed and Justified version of the input_string into output_string
        subroutine CopyString(input_string, output_string)
            !arguments
            character(len=*), intent(in)    :: input_string
            character(len=:), allocatable, intent(inout) :: output_string
            !private variables
            !integer                                      :: length_of_input_string

            output_string = trim(adjustl(input_string))
        end subroutine CopyString

        !>  \brief  Convert C string of characters to Fortran string of characters
        subroutine ConvertStringFromCToF(c_string, f_string)
            use iso_c_binding
            character(len=1,kind=c_char),   intent(in)  :: c_string(*)
            character(len=*),               intent(out) :: f_string
            integer :: i
            i=1
            do while(c_string(i)/=c_null_char .and. i<=len(f_string))
                f_string(i:i) = c_string(i)
                i=i+1
            end do
            if (i<=len(f_string)) f_string(i:) = ' '
        end subroutine ConvertStringFromCToF


        !>  \brief  Return a C string given an F string
        function GetCStringFromFString(f_string) result(c_string)
            character(len=*),               intent(in)  ::  f_string
            character(len=1,kind=c_char)                ::  c_string(len_trim(f_string)+1)
            !
            integer ::  i
            !
            do i=1,len_trim(f_string)
                c_string(i) = f_string(i:i)
            enddo
            c_string(len_trim(f_string)+1) = c_null_char
        end function GetCStringFromFString


        pure logical function StringStartsWithAQuestionMark(line)
            ! argument
            character(len=*),  intent(in)  ::  line
            ! local variables
            integer ::  first_non_blank
            ! start work
            first_non_blank = verify(line,blank_characters)
            if (line(first_non_blank:first_non_blank) .eq. '?') then
            StringStartsWithAQuestionMark = .true.
            else
            StringStartsWithAQuestionMark = .false.
        endif
        end function StringStartsWithAQuestionMark

        !>  \brief  works out whether a character string is a real
        pure logical function StringIsAReal(line)

            ! argument
            character(len=*), intent(in)  ::  line
            ! local variables
            integer ::  first_non_blank, last_non_blank
            ! start work
            first_non_blank = verify(line,blank_characters)
            last_non_blank  = verify(line,blank_characters,back=.true.)
            if (verify(line(first_non_blank:last_non_blank),'e+-0123456789.') .eq. 0) then
                StringIsAReal = .true.
            else
                StringIsAReal = .false.
            endif
        end function StringIsAReal

        !>  \brief  Remo
        function RemoveNonNumericCharacters(line) result (line_tmp)

            ! argument
            character(len=*), intent(in)  ::  line
            character(len=len(line)) :: line_tmp
            ! local variables
            integer ::  first_non_blank, last_non_blank
            integer     ::  character_counter
            integer     ::  output_character_counter
            ! start work
            line_tmp = ' '
            first_non_blank = verify(line,blank_characters)
            last_non_blank  = verify(line,blank_characters,back=.true.)
            output_character_counter = 0
            do character_counter=first_non_blank,last_non_blank
                if (verify(line(character_counter:character_counter),'e+-0123456789.') .eq. 0) then
                    output_character_counter = output_character_counter + 1
                    line_tmp(output_character_counter:output_character_counter) = line(character_counter:character_counter)
                endif
            enddo
        end function RemoveNonNumericCharacters

                !>  \brief  works out whether a character string is a real
        pure logical function StringIsAnInteger(line)

            ! argument
            character(len=*), intent(in)  ::  line
            ! local variables
            integer ::  first_non_blank, last_non_blank
            ! start work
            first_non_blank = verify(line,blank_characters)
            last_non_blank  = verify(line,blank_characters,back=.true.)
            if (verify(line(first_non_blank:last_non_blank),'+-0123456789') .eq. 0) then
                StringIsAnInteger = .true.
            else
                StringIsAnInteger = .false.
            endif
        end function StringIsAnInteger


        !>  \brief  emulates perl's split routine, which returns an array of strings
        subroutine Split(line,strings,separators)

            ! arguments
            character(len=*)                        ::  line        !<  line to be split
            character(len=*),           optional    ::  separators  !<  characters which separate words, if not present, default is blank characters (space, tabs...)
            ! results
            character(len=line_max_len),    allocatable ::  strings(:)
            ! private variables
            character(len=line_max_len)                 ::  buffer
            character(len=line_max_len),    allocatable ::  strings_tmp(:)
            integer                                 ::  pos1, pos2, word_count, ierr
            logical,                    parameter   ::  debug = .false.
            ! start work
            if (StringIsBlank(line)) then
                ! line is blank, crash out
                call this_program%TerminateWithFatalError('StringManipulations::Split', &
                                                 'Cannot split a blank line!')
            endif
            if (debug) print '(3a)', '**debug(split): line = |', trim(line), '|'
            ! start work in earnest
            buffer = trim(line)
            word_count = 0
            do
                if (present(separators)) then
                    ! find the first non-separator character
                    pos1 = verify(buffer,separators,back=.false.)
                    if (pos1 .ne. 0) then
                        ! find the last non-separator character after that
                        pos2 = pos1 + scan(buffer(pos1:),separators,back=.false.) - 2
                    endif
                else
                    ! find the first non-blank character
                    pos1 = FirstNonBlank(buffer)
                    if (pos1 .ne. 0) then
                        ! find the last non-blank character after that
                        pos2 = pos1 + FirstBlank(buffer(pos1:)) - 2
                        if (debug) print '(a,2(i0,1x))', '**debug(split): pos1, pos2 = ', pos1, pos2
                    endif
                endif
                if (pos1 .ne. 0 .and. pos2 .ne. 0) then ! if we found a word
                    if (word_count .gt. 0) strings_tmp = strings
                    word_count = word_count + 1
                    if (allocated(strings)) deallocate(strings)
                        allocate(strings(word_count),stat=ierr)
                    if (ierr .ne. 0) then
                    call this_program%TerminateWithFatalError('StringManipulations::Split', &
                                                 'could not allocate memory to output array!')
                endif
                if (word_count .gt. 1) strings(1:word_count-1) = strings_tmp(1:word_count-1)
                strings(word_count) = ''
                if (debug) print '(3a)', '**debug(split): adding this word: |', buffer(pos1:pos2), '|'
                write(strings(word_count),'(a)') buffer(pos1:pos2)
                if ( pos2+1 .ge. len_trim(buffer)) then ! if we reached end of the line
                    exit
                else
                    buffer = buffer(pos2+1:)
                endif
                else
                exit
            endif
        enddo

        if (debug) then
            print *, '**debug(split): result array = ', strings
        endif

    end subroutine Split

    !>  \brief  Count the number of records on a line
    function CountRecordsPerLine(line,separators) result(number_of_records)

            ! arguments
            character(len=*)                        ::  line        !<  line to be split
            character(len=*),           optional    ::  separators  !<  characters which separate words, if not present, default is blank characters (space, tabs...)

            ! private variables
            character(len=line_max_len)                 ::  buffer
            integer                                     ::  number_of_records
            integer                                     ::  pos1
            integer                                     ::  pos2

            ! start work

            if (StringIsBlank(line) .or. StringIsAComment(line)) then
                number_of_records = 0
            else
                buffer = trim(line)
                number_of_records = 0
                do
                    if (present(separators)) then
                        ! find the first non-separator character
                        pos1 = verify(buffer,separators,back=.false.)
                        if (pos1 .ne. 0) then
                            ! find the last non-separator character after that
                            pos2 = pos1 + scan(buffer(pos1:),separators,back=.false.) - 2
                        endif
                    else
                        ! find the first non-blank character
                        pos1 = FirstNonBlank(buffer)
                        if (pos1 .ne. 0) then
                            ! find the last non-blank character after that
                            pos2 = pos1 + FirstBlank(buffer(pos1:)) - 2
                        !if (debug) print '(a,2(i0,1x))', '**debug(split): pos1, pos2 = ', pos1, pos2
                        endif
                    endif

                    if (pos1 .ne. 0 .and. pos2 .ne. 0) then ! if we found a word
                        number_of_records = number_of_records + 1
                    endif

                    if ( pos2+1 .ge. len_trim(buffer)) then ! if we reached end of the line
                        exit
                    else
                        buffer = buffer(pos2+1:)
                    endif
                enddo
            endif



    end function CountRecordsPerLine

    !>  \brief  Parse flat-text file and retrieve a run-time parameter character string
    subroutine GetKeywordValueFromFile(filename, keyword, found_value, value_was_found,file_handle_unit_number)
        use UsefulFunctions, only : UnitIsOpen
        ! arguments
        character(len=*),               intent(in)      ::  filename        !<  name of file to be parsed
        character(len=*),               intent(in)      ::  keyword         !<  label of parameter to be read in
        character(len=:), allocatable,  intent(inout)   ::  found_value     !<  variable which will hold rtp value
        logical,            optional,   intent(out)     ::  value_was_found !<  indicates whether the parameter was found in the file
        integer,            optional,   intent(in)      ::  file_handle_unit_number !<  If the file is already open, supply the unit number
        ! private variables
        integer                                         ::  ios
        integer                                         ::  file_handle
        character(len=512)                              ::  io_err_msg
        ! initialise private variables
        ios     = 0

        ! debug
        !if (debug) write(*,'(6a)') '**debug(file_get_char): looking for ', trim(adjustl(parameter)),    &
        !                            ', which has current value of ', trim(adjustl(rtp_value)), ' in file ', &
        !                            trim(adjustl(filename))

        if (present(file_handle_unit_number)) then
            file_handle = file_handle_unit_number
        else
            ! open the file
            file_handle = this_program%GetAvailableUnit()
            open(unit=file_handle, file=filename,  iostat=ios, iomsg=io_err_msg)

            if (ios .ne. 0) then
                call this_program%TerminateWithFatalError('StringManipulations::GetKeywordValueFromFile', &
                                                        'Failed to open file: '//trim(adjustl(filename))//' ; '//trim(io_err_msg))
            endif
        endif

        ! use GetKeywordFromUnitNumber to do the actual work..

        if (present(value_was_found)) then
            call GetKeywordValueFromUnitNumber(file_handle, keyword, found_value, value_was_found)
        else
            call GetKeywordValueFromUnitNumber(file_handle, keyword, found_value)
        endif

        if (.not. present(file_handle_unit_number)) then
            close(file_handle)
            call this_program%ReleaseUnit(file_handle)
        endif

          ! debug
        !write(*,*) '**debug(file_get_char): was looking for ', trim(parameter), ', which now has value of ', trim(rtp_value), '.'
    end subroutine GetKeywordValueFromFile

        !>  \brief  Parse flat-text file and retrieve a run-time parameter character string
    subroutine GetKeywordValueFromUnitNumber(unit_number, keyword, found_value, value_was_found)

        use UsefulFunctions, only : UnitIsOpen
        ! arguments
        integer,                        intent(in)      ::  unit_number     !<  name of file to be parsed
        character(len=*),               intent(in)      ::  keyword         !<  label of parameter to be read in
        character(len=:), allocatable,  intent(inout)   ::  found_value     !<  variable which will hold rtp value
        logical,            optional,   intent(out)     ::  value_was_found !<  indicates whether the parameter was found in the file
        ! private variables

        character(len=line_max_len)                     ::  buffer      !< will hold a line from the file
        character(len=line_max_len)                     ::  label       !< will hold the label part of the line
        character(len=line_max_len)                     ::  value       !< will hold the value part of the line
        character(len=line_max_len),    allocatable     ::  words(:)
        integer                                         ::  ios         !< ios is negative if an end of record condition is encountered or if
                                                                        !! an endfile condition was detected.  it is positive if an error was
                                                                        !! detected.  ios is zero otherwise.
        character(len=line_max_len)                     ::  io_message
        integer                                         ::  line        !< line number
        logical                                         ::  ffound
        logical,                        parameter       ::  debug   =   .false.

        ! initialise private variables
        line    = 0
        ios     = 0
        ffound  = .false.

        if (UnitIsOpen(unit_number)) then
            ! rewind to first line of file
            rewind(unit_number, iostat=ios,iomsg=io_message)
            if (ios .ne. 0) then
                call this_program%TerminateWithFatalError('StringManipulations::GetKeywordValueFromFile', &
                                                 'Failed to rewind file: '//trim(io_message))
            endif

            ! iterate through the lines of the file (when eof, ios will be set to non-zero)
            ios = 0
            line = 0
            do while (ios == 0)
                read(unit_number, '(a)', iostat=ios) buffer
                if ( ios == 0 ) then
                    line    =   line + 1

                    !
                    ! if the line is a comment or is blank, we can cycle
                    !
                    if (StringIsAComment(buffer) .or. StringIsBlank(buffer)) then
                        !print *, '**debug(fil_get_char): cycling line: ', trim(buffer), '!'
                        cycle
                    endif

                    !
                    ! split the line into words, the first of which is the label, and the second the value
                    !
                    call split(buffer,words)
                    !if (debug) print '(a,i0,a)', '**debug(file_get_char): size, words array = ', size(words), words
                    label = words(1)
                    if (size(words) .gt. 1) then
                        value = words(2)
                    else
                        value = ''
                    endif

                    !if (debug) print '(a,i0,5a)', '**debug(file_get_char): line, label, value: |', line,    &
                    !                               '|', trim(label), '|', trim(value), '|'


                    !
                    ! is this the label we're looking for?
                    !
                    if ( StringsAreEqual(label, keyword,case_sensitive=.false.) ) then
                        ! we've found the label we were looking for - transfer it to output variable, and stop scanning through the file
                        call CopyString(value, found_value)
                    !   if (debug) write(*,'(4a)') '**debug(file_get_char): found ', trim(parameter), '=', trim(value)
                        ffound = .true.
                        exit
                    else
                        ! this line didn't contain the rtp we were looking for, do nothing
                    endif
                elseif ( ios .gt. 0) then
                    call this_program%TerminateWithFatalError('StringManipulations::GetKeywordValueFromFile', &
                                                          'io error encountered while parsing line')
                elseif ( ios .lt. 0) then
                    ! we reached the end of the file and still haven't found what we're looking for
                    ! if (debug) write(*,'(2a)') '**debug(file_get_char): reached end of file, did not find parameter ', trim(parameter)
                endif
            enddo

            ! rewind to first line of file
            rewind(unit_number, iostat=ios,iomsg=io_message)


        else !unit isn't open

            write(*,'(a)') '**Warning(StringManipulations::GetKeywordFromUnitNumber): Supplied unit is not open'
            write(*,*)

            ffound = .false.

        endif

        !cleanup
        if (allocated(words)) deallocate(words)

        !
        if (.not. ffound) then
            call CopyString(' ', found_value)
        endif

        if (present(value_was_found)) then
            value_was_found = ffound
        endif

            ! debug
        !   write(*,*) '**debug(file_get_char): was looking for ', trim(parameter), ', which now has value of ', trim(rtp_value), '.'
    end subroutine GetKeywordValueFromUnitNumber



end module StringManipulations
