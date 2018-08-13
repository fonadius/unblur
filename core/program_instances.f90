module ProgramInstances

    implicit none
    private

    integer,    parameter   ::  lowest_unit_number = 20
    integer,    parameter   ::  highest_unit_number = 200

    !>  \brief Hold data about a program instance
    type, public :: ProgramInstance
        private

        logical, public                         ::  running_interactively     !<  whether program is being run interactively
        logical, public                         ::  reading_from_terminal        !<  whether program is reading input from standard input

        character(len=:), public, allocatable   ::  program_name
        character(len=:), public, allocatable   ::  control_filename
        character(len=:), public, allocatable   ::  program_version

        ! start and finish time
        integer                                 ::  date_start(8)
        integer                                 ::  date_finish(8)


        ! Unit numbers - we keep track of the ones we're currently using
        logical                                 ::  unit_is_available(lowest_unit_number:highest_unit_number) = .true.

        contains
            procedure                           ::  Init
            procedure                           ::  SetStartDate
            procedure                           ::  ParseCommandLineArguments
            procedure                           ::  TerminateWithFatalError
            procedure                           ::  Terminate
            procedure                           ::  GetAvailableUnit
            procedure                           ::  ReleaseUnit
    end type


    contains


    !>  \brief  Find an available unit and reserve it for use
    integer function GetAvailableUnit(self)
        ! arguments
        class(ProgramInstance), intent(inout)   ::  self
        ! private variables
        logical     ::  success
        integer     ::  current_unit
        ! start work
        !$omp critical (ProgramInstance_UnitManagement)
        success = .false.
        do current_unit=lowest_unit_number,highest_unit_number
            if (self%unit_is_available(current_unit)) then
                success = .true.
                GetAvailableUnit = current_unit
                exit
            endif
        enddo
        if (success) then
            self%unit_is_available(GetAvailableUnit) = .false.
        else
            call self%TerminateWithFatalError('ProgramInstance::GetAvailableUnit','Could not find an unused io unit')
        endif
        !$omp end critical (ProgramInstance_UnitManagement)
    end function GetAvailableUnit

    !>  \brief  Let the world know that a particular unit is not being used any more
    subroutine ReleaseUnit(self,unit_number)
        class(ProgramInstance),     intent(inout)   ::  self
        integer,                    intent(in)      ::  unit_number
        ! start work
        !$omp critical (ProgramInstance_UnitManagement)
        if (unit_number .ge. lbound(self%unit_is_available,1) .and. &
            unit_number .le. ubound(self%unit_is_available,1)) then
            self%unit_is_available(unit_number) = .true.
        endif
        !$omp end critical (ProgramInstance_UnitManagement)
    end subroutine ReleaseUnit

    !>  \brief  Set the start date
    subroutine SetStartDate(self)
        use DatesAndTimes
        class(ProgramInstance), intent(inout)   ::  self
        type(DateAndTime)   ::  my_date_and_time
        self%date_start = my_date_and_time%DateAndTimeAsIntegers()
    end subroutine SetStartDate

    !>  \brief  Initialise a program instance
    subroutine Init(self, wanted_program_name, wanted_version, wanted_copyright_year, command_line_options)
        use DatesAndTimes
        use iso_fortran_env
#ifdef __INTEL_COMPILER
        use ifport         ! if we are using the fortran compiler include IFPORT to check if we are reading from tty.
#else
#ifdef NAGFOR_COMPILER
        use f90_unix_env
        use f90_unix_errno
#endif
        !
#endif
        ! arguments
        class(ProgramInstance),                         intent(inout)   ::  self
        character(len=*),                               intent(in)      ::  wanted_program_name
        character(len=*),                               intent(in)      ::  wanted_version
        character(len=*),                               intent(in)      ::  wanted_copyright_year
        character(len=:),   allocatable,    optional,   intent(inout)   ::  command_line_options(:)    !<  Command line options given by the user


        ! private variables
        type(DateAndTime)                      ::  my_date_and_time
        character(len=2),    parameter         ::  text_width    =    '30'
        integer                                ::  counter
        integer                                ::  length_of_program_name
        integer                                ::  number_of_needed_spaces
        integer                                ::  length_of_version_string
        integer                                ::   errno


#ifdef __INTEL_COMPILER
        if (this_image() .eq. 1) then
#endif


        ! start work

        self%date_start = my_date_and_time%DateAndTimeAsIntegers()
        self%program_name = trim(adjustl(wanted_program_name))
        self%program_version = trim(adjustl(wanted_version))

        !Set all unit slots to available
        self%unit_is_available = .true.

        ! Parse command line arguments

        write(*,*) ' '
        call self%ParseCommandLineArguments(command_line_options)

        ! Check to see if we are connected to a terminal
#ifdef NAGFOR_COMPILER
        call isatty(input_unit,self%reading_from_terminal,errno)
#else
        self%reading_from_terminal = isatty(input_unit)
#endif

        length_of_program_name = len(trim(adjustl(wanted_program_name)))
        length_of_program_name = length_of_program_name + 19
        number_of_needed_spaces = 30 - (length_of_program_name / 2)

        do counter = 0, number_of_needed_spaces
            write(output_unit, "(a)", advance="no")' '
        enddo

        write(*,'(3a)') '**  Welcome to ', trim(adjustl(wanted_program_name)), '  **'

        length_of_version_string = len(trim(adjustl(wanted_version)))
        length_of_version_string = length_of_version_string + 10
        number_of_needed_spaces = 30 - (length_of_version_string / 2)

        do counter = 0, number_of_needed_spaces
            write(output_unit, "(a)", advance="no")' '
        enddo


        write(*,'(2a)')'Version - ', trim(adjustl(wanted_version))
        write(*,'(a)')                         ' '

        if (self%running_interactively .and. self%reading_from_terminal) then
            write(output_unit,'(a'//text_width//',1x,a)')        'Input Mode:', 'Interactive'
        else if (.not. self%reading_from_terminal) then
            write(output_unit,'(a'//text_width//',1x,a)')        'Input Mode:', 'Batch'
        else
            write(output_unit,'(a'//text_width//',1x,a)')        'Input Mode:', 'Control File'
            write(output_unit,'(a'//text_width//',1x,a)')        'Control Filename:',    trim(adjustl(self%control_filename))
        endif

        write(output_unit,'(a'//text_width//',1x,a)')        'Date & Time:', my_date_and_time%DateAndTimeAsString()
        write(output_unit,'(a)')

        write(output_unit,'(/3a)') 'Copyright ', wanted_copyright_year, ' Howard Hughes Medical Institute. All rights reserved.'
        write(output_unit,'(a)') 'Use is subject to Janelia Farm Research Campus Software Copyright 1.1'
        write(output_unit,'(a/)') 'license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )'

#ifdef __INTEL_COMPILER
        endif
        flush(output_unit)
        sync all
#endif

    end subroutine Init

    !>  \brief  terminate the program
    subroutine Terminate(self)
        use DatesAndTimes

        ! arguments
        class(ProgramInstance), intent(inout)  ::  self
        ! private variables
        type(DateAndTime)   ::  my_date_and_time

        ! start work
#ifdef __INTEL_COMPILER
        flush(6)
        sync all
        if (this_image() .eq. 1) then
#endif
        self%date_finish = my_date_and_time%DateAndTimeAsIntegers()
        write(*,*)
        write(*,'(2a)') 'Total execution time : ', my_date_and_time%DurationAsString(  &
                                real(my_date_and_time%SecondsBetweenDates(self%date_start, self%date_finish)))
        write(*,'(4a)') my_date_and_time%DateAndTimeAsString(), '  : ', self%program_name, ' finished cleanly.'
        write(*,*)
#ifdef __INTEL_COMPILER
        endif
        sync all
#endif
        if (allocated(self%program_name)) deallocate(self%program_name)
        if (allocated(self%program_version)) deallocate(self%program_version)
        if (allocated(self%control_filename)) deallocate(self%control_filename)

        stop
    end subroutine Terminate



    subroutine TerminateWithFatalError(self,name_of_calling_function, wanted_error)
        use DatesAndTimes
#ifdef __INTEL_COMPILER
        use ifcore
#endif

        ! arguments
        class(ProgramInstance), intent(inout)   ::  self
        character(len=*),       intent(in)      ::  name_of_calling_function
        character(len=*),       intent(in)      ::  wanted_error

        ! private variables
        type(DateAndTime)                       ::  my_date_and_time

        ! start work
#if defined(__INTEL_COMPILER)
#if defined(WITH_DEBUG_SYMBOLS)
        call tracebackqq(user_exit_code=-1)
#endif
#endif
#ifdef __GFORTRAN__
#ifdef WITH_DEBUG_SYMBOLS
        call backtrace()
#endif
#endif
        self%date_finish = my_date_and_time%DateAndTimeAsIntegers()
        write(*,*)
        write(*,'(2a)') 'Total execution time : ', my_date_and_time%DurationAsString(  &
                                real(my_date_and_time%SecondsBetweenDates(self%date_start, self%date_finish)))
        write(*,'(/4a)',advance='no') my_date_and_time%DateAndTimeAsString(), ': Fatal error (', name_of_calling_function, '): '
        write(*,'(a//)') wanted_error

        ! Do some deallocation
        if (allocated(self%program_name)) deallocate(self%program_name)
        if (allocated(self%control_filename)) deallocate(self%control_filename)
        if (allocated(self%program_version)) deallocate(self%program_version)

#ifdef __INTEL_COMPILER
        error stop
#else
        stop
#endif

    end subroutine TerminateWithFatalError

    subroutine ParseCommandLineArguments(self,unprocessed_options)
#ifdef _OPENMP
        use omp_lib
#endif
        !

        class(ProgramInstance),                             intent(inout)   ::  self
        character(len=:),       allocatable,    optional,   intent(inout)   ::  unprocessed_options(:)

        ! local variables
        character(len=200)  ::  command_line    ! complete command line
        integer             ::  command_line_len    ! number of characters in command line
        integer             ::  arg_count   ! number of command-line arguments
        integer,    parameter   ::  max_arg_count   =   10   ! maximum number of arguments
        integer             ::  stat    ! status variable for error-checking
        character(len=200)  ::  current_arg ! argument currently being parsed
        integer             ::  i
        logical             ::  file_exists
        logical             ::  next_argument_should_be_number_of_threads
        integer             ::  number_of_threads
        !integer             ::  io_status
        !character(len=200)  ::  io_message
        character(len=200)  ::  number_of_threads_c
        logical             ::  omp_num_threads_set_by_command_line
        character(len=200)  ::  unprocessed_options_temp(max_arg_count)
        integer             ::  number_of_unprocessed_options


        self%control_filename = ''
        number_of_unprocessed_options = 0


#ifdef _OPENMP
        number_of_threads = omp_get_num_threads()
#endif
        omp_num_threads_set_by_command_line = .false.


        call get_command(command_line, command_line_len, stat)
        if (stat .eq. -1) then
            call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                              'command_line variable not long enough to hold the command line!')
        elseif (stat .gt. 0) then
            call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                              'command line could not be retrieved!')
        elseif (stat .lt. -1) then
            call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                              'fatal error in parse_command_line_arguments!')
        endif

        arg_count   =   command_argument_count()

        !if (arg_count .gt. max_arg_count) then
        !    write(*,'(a,i1,a,i1,a,i1,a)')   '**warning (parse_command_line_arguments): ', arg_count,                            &
        !                                    ' command-line arguments were supplied, but only expecting up to ', max_arg_count,  &
        !                                    '. ', arg_count - max_arg_count, ' arguments will be ignored.'
        !endif

        ! if there is more than one argument, we need to parse them
        if (arg_count .gt. 0) then
            next_argument_should_be_number_of_threads = .false.
            do i=1,arg_count

                call get_command_argument(i,current_arg,status=stat)
                if(stat .eq. -1) then
                    write(*,'(3a)')     '**error: argument supplied is longer', &
                                        ' than maximum allowable. either give a shorter argument, or recompile ',   &
                                        'program with larger character string.'
                    call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                                      'fatal error in parse_command_line_arguments!')
                elseif(stat .ne. 0) then
                    call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                                      'unknown fatal error in parse_command_line_arguments!')
                endif

                ! check if this is a flag
                if (current_arg(1:1) .eq. '-') then
                    if (next_argument_should_be_number_of_threads) then
                        call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                                'Malformed number of threads: '//trim(current_arg))
                    elseif (current_arg(1:17) .eq. '--omp-num-threads') then
                        if (len_trim(current_arg) .gt. 17) then
                            if (current_arg(18:18) .eq. '=') then
                                number_of_threads_c = current_arg(19:len_trim(current_arg))
                                if (StringIsAnInteger(number_of_threads_c)) then
                                    read(number_of_threads_c,*) number_of_threads
                                    omp_num_threads_set_by_command_line = .true.
                                else
                                    call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                                'Malformed number of threads: '//trim(number_of_threads_c))
                                endif
                            else
                                ! Malformed argument
                                call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                                                  'Malformed argument: '//trim(current_arg))
                            endif
                        else
                            ! We will expect the next argument to be the number of threads
                            next_argument_should_be_number_of_threads = .true.
                            if (i .eq. arg_count) then
                                call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                                                  'Did not find number of threads in command line arguments')
                            endif
                        endif
                    else
                        ! Unrecognised flag
                        number_of_unprocessed_options = number_of_unprocessed_options + 1
                        unprocessed_options_temp(number_of_unprocessed_options) = current_arg
                        !write(*,'(2a)') '**debug(parse_command_line_arguments): unrecognised flag: ', current_arg
                    endif
                else if (next_argument_should_be_number_of_threads) then
                    number_of_threads_c = current_arg
                    if (StringIsAnInteger(number_of_threads_c)) then
                        read(number_of_threads_c,*) number_of_threads
                        next_argument_should_be_number_of_threads = .false.
                        omp_num_threads_set_by_command_line = .true.
                    else
                        call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                                'Malformed number of threads: '//trim(number_of_threads_c))
                    endif
                else ! not a flag
                    inquire(file=trim(adjustl(current_arg)), exist=file_exists)
                    if (file_exists) then
                        self%control_filename = current_arg
                    else
                        write(*,'(2a)') '**warning(parse_command_line_arguments): this file does not exist: ', &
                                        trim(adjustl(current_arg))
                        write(*,'(3a)')  '**warning(parse_command_line_arguments): will run ', self%program_name, &
                                            ' in interactive mode'
                    endif
                endif
            enddo
        endif

        ! Deal with unprocessed options
        if (present(unprocessed_options)) then
            if (allocated(unprocessed_options)) deallocate(unprocessed_options)
            if (number_of_unprocessed_options .ge. 1) then
                allocate(character(len=200) :: unprocessed_options(number_of_unprocessed_options))
                unprocessed_options(1:number_of_unprocessed_options) = unprocessed_options_temp(1:number_of_unprocessed_options)
            endif
        endif

        ! Set number of threads
#ifdef _OPENMP
        if (omp_num_threads_set_by_command_line) then
            call omp_set_num_threads(number_of_threads)
        endif
#endif

        ! Determine whether we are running interactively
        self%running_interactively = self%control_filename .eq. ''


        contains


        !>  \brief  works out whether a character string is a real
        pure logical function StringIsAnInteger(line)
            use iso_c_binding
            ! argument
            character(len=*), intent(in)  ::  line
            ! local variables
            integer ::  first_non_blank, last_non_blank
            character(kind=c_char,len=*),   parameter :: BLANK_C_CHARACTERS = C_NULL_CHAR // C_HORIZONTAL_TAB
            character(len=*),               parameter :: BLANK_CHARACTERS =  ' '//BLANK_C_CHARACTERS
            ! start work
            first_non_blank = verify(line,blank_characters)
            last_non_blank  = verify(line,blank_characters,back=.true.)
            if (verify(line(first_non_blank:last_non_blank),'+-0123456789') .eq. 0) then
                StringIsAnInteger = .true.
            else
                StringIsAnInteger = .false.
            endif
        end function StringIsAnInteger

    end subroutine ParseCommandLineArguments

end module ProgramInstances
