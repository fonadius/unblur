!>  \brief  class & methods for gathering user inputs
module UserInputs

    use Globals
    use UserSuppliedParameters
    use UsefulFunctions
    use StringManipulations

    implicit none

    private

    type, public                :: UserInput

            integer                        :: defaults_file_unit = 9
            integer                        :: new_defaults_file_unit = 9

            character(len=:), allocatable  :: defaults_filename
            character(len=:), allocatable  :: new_defaults_filename

        contains
            procedure                   ::  Init
            !procedure                   ::  GetOptionFromUser
            procedure                   ::  GetLogicalFromUser
            procedure                   ::  GetIntegerFromUser
            procedure                   ::  GetRealFromUser
            procedure                   ::  GetFilenameFromUser
            procedure                   ::  UpdateDefaults
            procedure, private          ::  GetDefaultValue
            procedure, nopass, private  ::  AskQuestion
            procedure, nopass, private  ::  PrintHelp

            final               ::  Destroy
    end type

    contains

    !>  \brief  Initialise a user input object
    subroutine Init(self, wanted_program_name)
        use StringManipulations, only : RandomFilename
        ! arguments
        class(UserInput), intent(inout)        ::  self
        character(len=*), intent(in)           ::  wanted_program_name

        !variables
        integer                                ::  defaults_filename_length
        integer                                ::  ios



        integer,            parameter ::  keyword_column_width = 16
        character(len=4),   parameter ::  keyword_column_width_c   =   '16'
        integer,            parameter ::  value_column_width = 45
        character(len=4),   parameter ::  value_column_width_c   =   '45'
        integer,            parameter ::  setbyuser_column_width = 14
        character(len=4),   parameter ::  setbyuser_column_width_c   =   '14'
        integer,            parameter ::  description_column_width = 64
        character(len=4),   parameter ::  description_column_width_c   =   '64'



        defaults_filename_length = len_trim(adjustl(wanted_program_name))
        defaults_filename_length = defaults_filename_length + 5

        ! start work

        ! don't open defaults file this unless we're in interactive mode
        if (this_program%running_interactively) then
            self%defaults_filename = '.'//trim(adjustl(wanted_program_name))//'.dff'

            !new defaults filename is same as above but with current after the dot.

            defaults_filename_length = defaults_filename_length + 7
            self%new_defaults_filename = RandomFilename(16,prefix='.'//trim(this_program%program_name)//'_')

            self%defaults_file_unit = this_program%GetAvailableUnit()
            open(unit=self%defaults_file_unit, file=self%defaults_filename, iostat=ios)

            if (ios .ne. 0) then
                close(self%defaults_file_unit)
                call this_program%ReleaseUnit(self%defaults_file_unit)
            endif

            self%new_defaults_file_unit = this_program%GetAvailableUnit()
            ! Open a temporary file
            open(   unit=self%new_defaults_file_unit, &
                    file=self%new_defaults_filename, &
                    status='replace', &
                    iostat=ios)

            if (ios .ne. 0) then
                write(*,'(a)') '**Warning(UserInput::Init): Failed to write new defaults file'
                write(*,'(a)') '                            Is there write access to current dir?'
                write(*, *)
                close(self%new_defaults_file_unit)
                call this_program%ReleaseUnit(self%new_defaults_file_unit)
            endif
        endif

        ! if not running interactively, print headings

        if (.not. this_program%running_interactively) then
            write(*,*)
            write(*,*) 'Running non-interactively, with the following inputs :-'
            write(*,*)
            write(*,'(a'//keyword_column_width_c//' ,1x)',advance='no') 'Keyword'
            write(*,'(a'//value_column_width_c//' ,1x)',advance='no') 'Value'
            write(*,'(a'//setbyuser_column_width_c//' ,1x)') 'User Supplied?'
            write(*,'(a'//keyword_column_width_c//' ,1x)',advance='no') '-------'
            write(*,'(a'//value_column_width_c//' ,1x)',advance='no') '-----'
            write(*,'(a'//setbyuser_column_width_c//' ,1x)') '--------------'

            !write(*,'(a'//description_column_width_c//' ,1x)') 'Description'

        endif


    end subroutine Init


 !   function GetOptionFromUser(self, question_to_ask, help_text_and_description, associated_keyword, default_value, &





!    end function GetOptionFromUser

    function GetFilenameFromUser(self, question_to_ask, help_text_and_description, associated_keyword, default_value, &
                                                                                file_must_exist) result(returned_value)

        ! arguments
        class(UserInput), intent(inout)        ::  self
        type(UserSuppliedFilename)             ::  returned_value

        character(len=*), intent(in)           ::  question_to_ask
        character(len=*), intent(in)           ::  help_text_and_description
        character(len=*), intent(in)           ::  associated_keyword
        character(len=*), intent(in)           ::  default_value

        logical, intent(in)                    ::  file_must_exist

        ! private variables
        logical                                ::  value_was_found
        logical                                ::  default_value_was_found
        character(len=:), allocatable          ::  current_default_value

        character(len=:), allocatable          ::  user_supplied_entry

        ! Are we running in interactive mode?

        if (this_program%running_interactively) then

            ! Get the default value..

            call self%GetDefaultValue(associated_keyword, current_default_value, default_value_was_found)

            ! if we didn't find a default value, then we need to set it to the default-default supplied to the subroutine

            if (.not. default_value_was_found) then
                call CopyString(default_value, current_default_value)
            endif

            ! get an answer

            do

                call self%AskQuestion(question_to_ask, current_default_value, user_supplied_entry)

                ! is the returned value interpretable as a real..

                if (StringIsBlank(user_supplied_entry)) then != The user wants to accept the default value
                    call CopyString(current_default_value, user_supplied_entry)
                endif

                if (StringStartsWithAQuestionMark(user_supplied_entry)) then
                    call self%PrintHelp(help_text_and_description)
                    cycle
                endif

                if (file_must_exist) then

                    ! See if the file exists...
                    if (FileExists(user_supplied_entry)) then
                        exit

                    else

                        write(*,*)
                        write(*,*) 'No Such File - please enter a valid filename'
                        write(*,*)
                        cycle
                    endif
                else

                    ! check to see if we can open this filename.. (not doing this for now)



                    ! we can so exit

                    exit



                endif
            enddo

            returned_value%keyword = trim(adjustl(associated_keyword))
            returned_value%description = trim(adjustl(help_text_and_description))
            returned_value%set_by_user = .true.

            ! Copy in the actual filename
            returned_value%value = trim(adjustl(user_supplied_entry))


            ! Add the user_supplied_entry to the new defaults file..

            if (UnitIsOpen(self%new_defaults_file_unit)) then
                write(self%new_defaults_file_unit, '(a,a,a)') associated_keyword, ' ', user_supplied_entry
            endif


        else  ! we are not running interactively, we need to get the value from the command file..

            call GetKeywordValueFromFile(this_program%control_filename, associated_keyword, user_supplied_entry, value_was_found)

            ! check this value is ok, if it is then use it, if not take the default

            if (.not. value_was_found) then
                call CopyString(default_value, user_supplied_entry)
            endif

            ! check this value ok..

            if (file_must_exist) then

                !See if the file exists...
                 if (.not. FileExists(user_supplied_entry)) then
                    write(*,'(4a)') '**Error: Filename supplied for ', associated_keyword, &
                                    ' in invalid: ', trim(adjustl(user_supplied_entry))
                    call this_program%TerminateWithFatalError('UserInput::GetFilenameFromUser', &
                                    'Filename in control file is invalid')
                 endif
            else

                !check fo see if we can open the file (not doing this for now)




            endif



            ! all should be good, we should have our value, so update the new_defaults_file and return object

            call CopyString(associated_keyword, returned_value%keyword)
            call CopyString(help_text_and_description, returned_value%description)

            if (value_was_found) then
                returned_value%set_by_user = .true.
            else
                returned_value%set_by_user = .false.
            endif

            ! Copy in the actual filename
            returned_value%value = user_supplied_entry

            ! finally, print the info so there is some output when non-interactive

            call returned_value%PrintInfo()
        endif

    end function GetFilenameFromUser

    function GetLogicalFromUser(self, question_to_ask, help_text_and_description, associated_keyword, default_value) &
                                                                                                result(returned_value)

        ! arguments
        class(UserInput), intent(inout)        ::  self
        type(UserSuppliedLogical)              ::  returned_value

        character(len=*), intent(in)           ::  question_to_ask
        character(len=*), intent(in)           ::  help_text_and_description
        character(len=*), intent(in)           ::  associated_keyword
        character(len=*), intent(in)           ::  default_value

        ! private variables
        logical                                ::  value_was_found
        logical                                ::  default_value_was_found
        character(len=:), allocatable          ::  current_default_value

        character(len=:), allocatable          ::  user_supplied_entry
        integer                                ::  first_non_blank

        ! Are we running in interactive mode?

        if (this_program%running_interactively) then

            ! Get the default value..

            call self%GetDefaultValue(associated_keyword, current_default_value, default_value_was_found)

            ! if we didn't find a default value, then we need to set it to the default-default supplied to the subroutine

            if (.not. default_value_was_found) then
                call CopyString(default_value, current_default_value)
            endif

            ! get an answer

            do

                call self%AskQuestion(question_to_ask, current_default_value, user_supplied_entry)

                ! is the returned value interpretable as a real..

                if (StringIsBlank(user_supplied_entry)) then != The user wants to accept the default value
                    call CopyString(current_default_value, user_supplied_entry)
                endif

                if (StringStartsWithAQuestionMark(user_supplied_entry)) then
                    call self%PrintHelp(help_text_and_description)
                    cycle
                endif

                ! See if the value is a T/F/Y/N

                first_non_blank = FirstNonBlank(user_supplied_entry)

                if (verify(user_supplied_entry(first_non_blank:first_non_blank),'YyNnTtFf') .eq. 0) then

                    ! Value is valid
                    exit

                else

                write(*,*)
                write(*,*) 'Please enter Yes/No or True/False!'
                write(*,*)

                endif
            enddo

            ! now we should have a value as chars stored in returned value, we need to convert
            ! this to a real, then setup the UserSuppliedReal Object that is returned by the
            ! function.

            call CopyString(associated_keyword, returned_value%keyword)
            call CopyString(help_text_and_description, returned_value%description)
            returned_value%set_by_user = .true.

            select case (user_supplied_entry(first_non_blank:first_non_blank))
                case ('Y','y','T','t')
                    returned_value%value = .true.
                case default
                    returned_value%value = .false.
            end select
            ! Add the user_supplied_entry to the new defaults file..

            if (UnitIsOpen(self%new_defaults_file_unit)) then
                write(self%new_defaults_file_unit, '(a,a,a)') associated_keyword, ' ', user_supplied_entry
            endif


        else  ! we are not running interactively, we need to get the value from the command file..

            call GetKeywordValueFromFile(this_program%control_filename, &
                                         associated_keyword, user_supplied_entry, value_was_found)

            ! check this value is ok, if it is then use it, if not take the default

            if (.not. value_was_found) then
                call CopyString(default_value, user_supplied_entry)
            endif

            ! check this value ok..

            first_non_blank = FirstNonBlank(user_supplied_entry)

            if (verify(user_supplied_entry(first_non_blank:first_non_blank),'YyNnTtFf') .ne. 0) then
                 call this_program%TerminateWithFatalError('UserInput::GetLogicalFromUser', &
                                                            'Error in control file value!')
            endif

            ! all should be good, we should have our value, so update the new_defaults_file and return object

            call CopyString(associated_keyword, returned_value%keyword)
            call CopyString(help_text_and_description, returned_value%description)

            if (value_was_found) then
                returned_value%set_by_user = .true.
            else
                returned_value%set_by_user = .false.
            endif

            select case (user_supplied_entry(first_non_blank:first_non_blank))
                case ('Y','y','T','t')
                    returned_value%value = .true.
                case default
                    returned_value%value = .false.
            end select

            ! finally, print the info so there is some output when non-interactive

            call returned_value%PrintInfo()
        endif

    end function GetLogicalFromUser

    function GetIntegerFromUser(self, question_to_ask, help_text_and_description, associated_keyword, default_value, &
                             min_value, max_value,must_be_even,must_be_odd) result(returned_value)

        ! arguments
        class(UserInput), intent(inout)        ::  self
        type(UserSuppliedInteger)              ::  returned_value

        character(len=*), intent(in)           ::  question_to_ask
        character(len=*), intent(in)           ::  help_text_and_description
        character(len=*), intent(in)           ::  associated_keyword
        character(len=*), intent(in)           ::  default_value

        integer, optional, intent(in)          ::  min_value
        integer, optional, intent(in)          ::  max_value
        logical, optional, intent(in)          ::  must_be_even
        logical, optional, intent(in)          ::  must_be_odd

        ! private variables


        logical                                ::  value_was_found
        logical                                ::  default_value_was_found
        character(len=:), allocatable          ::  current_default_value
        character(len=:), allocatable          ::  user_supplied_entry
        integer                                ::  buffer_integer


        ! Are we running in interactive mode?

        if (this_program%running_interactively) then

            ! Get the default value..

            call self%GetDefaultValue(associated_keyword, current_default_value, default_value_was_found)

            ! if we didn't find a default value, then we need to set it to the default-default supplied to the subroutine

            if (.not. default_value_was_found) then
                call CopyString(default_value, current_default_value)
            endif

            ! get a useable integer..

            do

                call self%AskQuestion(question_to_ask, current_default_value, user_supplied_entry)

                ! is the returned value interpretable as a real..

                if (StringIsBlank(user_supplied_entry)) then != The user wants to accept the default value
                    call CopyString(current_default_value, user_supplied_entry)
                endif

                if (StringStartsWithAQuestionMark(user_supplied_entry)) then
                    call self%PrintHelp(help_text_and_description)
                    cycle
                endif

                ! See if the value is a real, if so, exit the loop

                if (StringIsAnInteger(user_supplied_entry)) then

                    ! if necessary check it is within the limits..

                    read(user_supplied_entry,*) buffer_integer

                    ! Check min/max
                    if (present(min_value) .and. present(max_value)) then
                        if (buffer_integer .lt. min_value .or. buffer_integer .gt. max_value) then
                            write(*,*)
                            write(*,'(a,i0,a,i0)') 'Please enter a number between ', min_value, ' and ', max_value
                            write(*,*)
                            cycle
                        endif
                    else if (present(min_value)) then
                        if (buffer_integer .lt. min_value) then
                            write(*,*)
                            write(*,'(a,f6.2)') 'Please enter a number greater than ', min_value
                            write(*,*)
                            cycle
                        endif
                    else if (present(max_value)) then
                        if (buffer_integer .gt. max_value) then
                            write(*,*)
                            write(*,'(a,f6.2)') 'Please enter a number less than ', max_value
                            write(*,*)
                            cycle
                        endif
                    endif

                    ! Check even/odd
                    if (present(must_be_odd)) then
                        if (must_be_odd .and. IsEven(buffer_integer)) then
                            write(*,*)
                            write(*,'(a)') 'Please enter an odd number '
                            write(*,*)
                            cycle
                        endif
                    else if (present(must_be_even)) then
                        if (must_be_even .and. IsOdd(buffer_integer)) then
                            write(*,*)
                            write(*,'(a)') 'Please enter an even number '
                            write(*,*)
                            cycle
                        endif
                    endif

                    exit

                else

                write(*,*)
                write(*,*) 'Please enter a valid Integer!'
                write(*,*)

                endif
            enddo

            ! now we should have a value as chars stored in returned value, we need to convert
            ! this to a real, then setup the UserSuppliedReal Object that is returned by the
            ! function.

            call CopyString(associated_keyword, returned_value%keyword)
            call CopyString(help_text_and_description, returned_value%description)
            returned_value%set_by_user = .true.
            read(user_supplied_entry,*) returned_value%value

            ! Add the user_supplied_entry to the new defaults file..
            if (UnitIsOpen(self%new_defaults_file_unit)) then
                write(self%new_defaults_file_unit, '(a,a,a)') associated_keyword, ' ', user_supplied_entry
            endif


        else  ! we are not running interactively, we need to get the value from the command file..

            call GetKeywordValueFromFile(this_program%control_filename, associated_keyword, user_supplied_entry, value_was_found)

            ! check this value is ok, if it is then use it, if not take the default

            if (.not. value_was_found) then
                call CopyString(default_value, user_supplied_entry)
            endif

            ! check this value is a real..

            if (.not. StringIsAnInteger(user_supplied_entry) .or. StringIsBlank(user_supplied_entry)) then
                 call this_program%TerminateWithFatalError('UserInput::GetIntegerFromUser', &
                                                           'Error in control file value!')
            endif

            ! all should be good, we should have our value, so update the new_defaults_file and return object

            call CopyString(associated_keyword, returned_value%keyword)
            call CopyString(help_text_and_description, returned_value%description)

            if (value_was_found) then
                returned_value%set_by_user = .true.
            else
                returned_value%set_by_user = .false.
            endif

            read(user_supplied_entry,*) returned_value%value

            ! finally, print the info so there is some output when non-interactive

            call returned_value%PrintInfo()
        endif

    end function GetIntegerFromUser

    function GetRealFromUser(self, question_to_ask, help_text_and_description, associated_keyword, default_value, &
                             min_value, max_value) result(returned_value)

        ! arguments
        class(UserInput), intent(inout)  ::  self
        type(UserSuppliedReal)           ::  returned_value

        character(len=*), intent(in)  ::  question_to_ask
        character(len=*), intent(in)  ::  help_text_and_description
        character(len=*), intent(in)  ::  associated_keyword
        character(len=*), intent(in)  ::  default_value

        real, optional, intent(in)             ::  min_value
        real, optional, intent(in)             ::  max_value

        ! private variables


        logical                                ::  value_was_found
        logical                                ::  default_value_was_found
        character(len=:), allocatable          ::  current_default_value
        character(len=:), allocatable          ::  user_supplied_entry
        real                                   ::  buffer_real


        ! Are we running in interactive mode?

        if (this_program%running_interactively) then

            ! Get the default value..

            call self%GetDefaultValue(associated_keyword, current_default_value, default_value_was_found)

            ! if we didn't find a default value, then we need to set it to the default-default supplied to the subroutine

            if (.not. default_value_was_found) then
                call CopyString(default_value, current_default_value)
            endif

            ! get a useable real..

            do

                call self%AskQuestion(question_to_ask, current_default_value, user_supplied_entry)

                ! is the returned value interpretable as a real..

                if (StringIsBlank(user_supplied_entry)) then != The user wants to accept the default value
                    call CopyString(current_default_value, user_supplied_entry)
                endif

                if (StringStartsWithAQuestionMark(user_supplied_entry)) then
                    call self%PrintHelp(help_text_and_description)
                    cycle
                endif

                ! See if the value is a real, if so, exit the loop

                if (StringIsAReal(user_supplied_entry)) then

                    ! if necessary check it is within the limits..

                    read(user_supplied_entry,*) buffer_real

                    if (present(min_value) .and. present(max_value)) then
                        if (buffer_real .lt. min_value .or. buffer_real .gt. max_value) then
                            write(*,*)
                            write(*,'(a,f6.2,a,f6.2)') 'Please enter a number between ', min_value, ' and ', max_value
                            write(*,*)
                            cycle
                        endif
                    else if (present(min_value)) then
                        if (buffer_real .lt. min_value) then
                            write(*,*)
                            write(*,'(a,f6.2)') 'Please enter a number greater than ', min_value
                            write(*,*)
                            cycle
                        endif
                    else if (present(max_value)) then
                        if (buffer_real .gt. max_value) then
                            write(*,*)
                            write(*,'(a,f6.2)') 'Please enter a number less than ', max_value
                            write(*,*)
                            cycle
                        endif

                    endif

                    exit

                else

                write(*,*)
                write(*,*) 'Please enter a valid number!'
                write(*,*)

                endif
            enddo

            ! now we should have a value as chars stored in returned value, we need to convert
            ! this to a real, then setup the UserSuppliedReal Object that is returned by the
            ! function.

            call CopyString(associated_keyword, returned_value%keyword)
            call CopyString(help_text_and_description, returned_value%description)
            returned_value%set_by_user = .true.
            read(user_supplied_entry,*) returned_value%value

            ! Add the user_supplied_entry to the new defaults file..

            if (UnitIsOpen(self%new_defaults_file_unit)) then
                write(self%new_defaults_file_unit, '(a,a,a)') associated_keyword, ' ', user_supplied_entry
            endif


        else  ! we are not running interactively, we need to get the value from the command file..

            call GetKeywordValueFromFile(this_program%control_filename, associated_keyword, user_supplied_entry, value_was_found)

            ! check this value is ok, if it is then use it, if not take the default

            if (.not. value_was_found) then
                call CopyString(default_value, user_supplied_entry)
            endif

            ! check this value is a real..

            if (.not. StringIsAReal(user_supplied_entry) .or. StringIsBlank(user_supplied_entry)) then
                 call this_program%TerminateWithFatalError('UserInput::GetRealFromUser', &
                                                 'Error in control file value!')
            endif

            ! all should be good, we should have our value, so update the new_defaults_file and return object

            call CopyString(associated_keyword, returned_value%keyword)
            call CopyString(help_text_and_description, returned_value%description)

            if (value_was_found) then
                returned_value%set_by_user = .true.
            else
                returned_value%set_by_user = .false.
            endif

            read(user_supplied_entry,*) returned_value%value

            ! finally, print the info so there is some output when non-interactive

            call returned_value%PrintInfo()
        endif

    end function GetRealFromUser




    subroutine AskQuestion(question_to_ask, current_default_value, user_supplied_entry)
        use iso_fortran_env
        !class(UserInput), intent(inout)   ::  self
        character(len=*), intent(in)      ::  question_to_ask        !<  Question to be put to the user
        character(len=*), intent(in)      ::  current_default_value  !<  Default option to be put to the user
        character(len=:), allocatable,  intent(inout)   ::  user_supplied_entry    !<  What the user typed

        character(len=line_max_len)       ::  user_answer
        integer                           ::  length_of_question
        integer                           ::  length_of_default_value
        integer                           ::  combined_length
        integer                           ::  number_of_spaces_needed
        integer                           ::  counter


        !start work


        !work out the length of the question, this is equal to len(question_to_ask) + (2 spaces)

        length_of_question = len(question_to_ask) + 2

        !work out the length of the default, this is equal to len(current_default_value) + 2 ([ and ])

        length_of_default_value = len(current_default_value) + 2

        combined_length = length_of_question + length_of_default_value

        if (combined_length < user_input_max_len) then ! no need for wrapping, reasonably simple

            number_of_spaces_needed = user_input_max_len - combined_length

            write(output_unit, '(a)', advance='no') trim(question_to_ask)

            ! write the correct number of spaces

            do counter = 0, number_of_spaces_needed
                write(output_unit, "(a)", advance="no")' '
            !write(*, advance='no') ' '
            enddo

            ! write the final colon

            write(output_unit,'(3a)', advance='no') ' [', trim(current_default_value), '] : '

            else

            ! we need some wrapping..


            write(*, '(a)') trim(question_to_ask)
            number_of_spaces_needed = user_input_max_len - length_of_default_value;

            if (number_of_spaces_needed < 0) number_of_spaces_needed = 0

            write(output_unit,'(3a)', advance='no') '[', trim(current_default_value), ']'

            ! write the correct number of spaces

            do counter = 0, number_of_spaces_needed
                write(output_unit, "(a)", advance="no")' '
            !write(*, advance='no') ' '
            enddo

            write(output_unit,'(a)', advance='no') ': '



        endif


        !ok, get an input from the user..
        read (unit=input_unit,fmt='(a)') user_answer

        ! if we aren't reading from a terminal, then print out this answer plus a carriage return..

        if (.not. this_program%reading_from_terminal) then
            write(output_unit,"(a)") trim(adjustl(user_answer))
        endif

        !copy this over to user_supplied_entry

        call CopyString(trim(adjustl(user_answer)), user_supplied_entry)

    end subroutine AskQuestion

    subroutine PrintHelp(text_to_print)

        !class(UserInput), intent(inout)   ::  self
        character(len=*), intent(in)      ::  text_to_print        !<  Question to be put to the user

        write (*, *)
        write (*, '(a)') text_to_print
        write (*, *)

    end subroutine PrintHelp


    subroutine GetDefaultValue(self, wanted_keyword, found_default_value, value_was_found)

        use StringManipulations

        ! arguments
        class(UserInput), intent(inout)  ::  self
        character(len=*),               intent(in)      ::  wanted_keyword        !<  label of parameter to be read in
        character(len=:), allocatable,  intent(inout)   ::  found_default_value   !<  the found value
        logical,            optional,   intent(out)     ::  value_was_found       !<  indicates whether the parameter was found in the file
        ! private variables

        logical                                         ::  ffound
        ! first_check if the default file is open, if it isn't, we can't get a default value

        if (UnitIsOpen(self%defaults_file_unit)) then

            call GetKeywordValueFromUnitNumber(self%defaults_file_unit, wanted_keyword, found_default_value, ffound)

        else

            ! set the default value to blank
            call CopyString(' ', found_default_value)
            ffound = .false.

        endif

        ! if we have the variable, set whether we found the value or not

        if (present(value_was_found)) then
            value_was_found = ffound
        endif

       end subroutine GetDefaultValue

    subroutine UpdateDefaults(self)

        ! arguments
        class(UserInput), intent(inout)  ::  self

        ! private variables
        character(len=line_max_len)          ::  buffer      !< will hold a line from the file
        character(len=:),  allocatable       ::  label       !< will hold the label part of the line
        character(len=:),  allocatable       ::  value       !< will hold the value part of the line
        character(len=:),  allocatable       ::  value2       !< will hold the value part of the line
        character(len=line_max_len),  allocatable       ::  words(:)
        integer                              ::  ios         !< ios is negative if an end of record condition is encountered or if
                                                                        !! an endfile condition was detected.  it is positive if an error was
                                                                        !! detected.  ios is zero otherwise.
        integer                              ::  line        !< line number
        logical                              ::  value_was_found

        !
        ! We will only update defaults if we are running interactively
        !
        if (this_program%running_interactively) then

            ! We need to copy all the values that are in the defaults file, and AREN'T in the new defaults
            ! file, to the new defaults file. Then we close the files, copy the new to the original, and then
            ! reopen them in case we aren't finished.
            if (FileIsOpen(self%defaults_filename) .and. FileIsOpen(self%new_defaults_filename)) then

                ! rewind defaults filename
                rewind(self%defaults_file_unit, iostat=ios)

                ! iterate through the lines of the file (when eof, ios will be set to non-zero)
                ios = 0
                line = 0

                do while (ios == 0)
                    read(self%defaults_file_unit, '(a)', iostat=ios) buffer
                    if ( ios == 0 ) then
                        line    =   line + 1

                        ! if the line is a comment or is blank, we can cycle
                        if (StringIsAComment(buffer) .or. StringIsBlank(buffer)) then
                            cycle
                        endif

                        ! split the line into words, the first of which is the label, and the second the value

                        call split(buffer,words)

                        label = trim(adjustl(words(1)))
                        if (size(words) .gt. 1) then
                            value = trim(adjustl(words(2)))
                        else
                            value = ''
                        endif

                        ! now we need to check if this label is already in the new_defaults file
                        call GetKeywordValueFromFile(   self%new_defaults_filename, label, value2, value_was_found)

                        ! If the label was not found in the new defaults file, write it there, along with the value
                        if (.not. value_was_found) then
                            !print *, '**', trim(label), trim(value)
                            write(self%new_defaults_file_unit, '(3a)') label, ' ', value
                        endif
                    endif
                enddo

                ! Now the new defaults files has everything up to date. We need to copy the new defaults file to the defaults file.
                if (FileIsOpen(self%defaults_filename)) then
                    close(self%defaults_file_unit)
                    call this_program%ReleaseUnit(self%defaults_file_unit)
                endif
                if (FileIsOpen(self%new_defaults_filename)) then
                    flush(self%new_defaults_file_unit)
                    close(self%new_defaults_file_unit)
                    call this_program%ReleaseUnit(self%new_defaults_file_unit)
                endif

                if (FileExists(self%new_defaults_filename)) call FileCopy(self%new_defaults_filename, self%defaults_filename)

                ! now delete the new_defaults_file, and reopen both files in case there is further use of user input
                if (FileIsOpen(self%new_defaults_filename)) then
                    close(self%new_defaults_file_unit,iostat=ios)
                    if (ios .ne. 0) then
                        write(*,'(a,i0)') '**debug(UserInput::UpdateDefaults): failed to close, status = ', ios
                    endif
                    call this_program%ReleaseUnit(self%new_defaults_file_unit)
                endif
                call FileDelete(self%new_defaults_filename)
            else
                ! The defaults files are not open, which suggests UpdateDefaults has already been run
                call this_program%TerminateWithFatalError('UserInput::UpdateDefaults','UpdateDefaults can only be called once')
            endif

            ! Cleanup
            if (allocated(label)) deallocate(label)
            if (allocated(value)) deallocate(value)
            if (allocated(value2)) deallocate(value2)
            if (allocated(words)) deallocate(words)

        endif ! end of test for running_interactively

    end subroutine UpdateDefaults

    subroutine Destroy(self)

        type(UserInput), intent(inout)   ::  self

        ! this destructor is getting called multiple times (stupid fortran) need some checks

        if (allocated(self%new_defaults_filename)) call self%UpdateDefaults()

        !call FileDelete(self%new_defaults_filename) ! This should not be necessary

        if (FileIsOpen(self%defaults_filename)) close(self%defaults_file_unit)
        if (FileIsOpen(self%new_defaults_filename)) close(self%new_defaults_file_unit)
        call FileDelete(self%new_defaults_filename)

        if (allocated(self%defaults_filename)) deallocate(self%defaults_filename)
        if (allocated(self%new_defaults_filename)) deallocate(self%new_defaults_filename)

    end subroutine Destroy

end module UserInputs
