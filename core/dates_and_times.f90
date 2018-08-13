module DatesAndTimes

    implicit none

    private

    integer,    parameter               ::  seconds_in_a_year        =   31536000
    integer,    parameter               ::  seconds_in_a_month       =   2592000
    integer,    parameter               ::  seconds_in_a_day         =   86400
    integer,    parameter               ::  seconds_in_an_hour       =   3600
    integer,    parameter               ::  seconds_in_a_minute      =   60
    integer,    parameter               ::  milliseconds_in_a_second =   1000

    type, public :: DateAndTime

        contains

        private

        procedure, nopass, public         :: DateAndTimeAsIntegers      !< Returns an integer array with date & time information (useful when you want to compute a duration)
        procedure, nopass, public         :: MilliSecondsBetweenDates   !< Return the number of milliseconds between two dates & times
        procedure, public                 :: SecondsBetweenDates        !< Return the number of seconds between two dates & times
        procedure, nopass, public         :: DateAndTimeAsString        !< Returns the date & time printed out in a nice format
        procedure, public                 :: DurationAsString           !< Convert a duration in seconds to a nice character string
        procedure, nopass, public         :: DurationAsIntegers         !< Convert a duration in seconds to a duration in days, hours, minutes and seconds


    end type

    contains
    !>    \brief    Returns an integer array with date & time information (useful when you want to compute a duration)
    function DateAndTimeAsIntegers() result(current_time)

        ! arguments

        integer :: current_time(8)
        call date_and_time(values=current_time)

    end function DateAndTimeAsIntegers

    !>  \brief  Return the number of milliseconds between two dates & times
    pure integer function MilliSecondsBetweenDates(date_first, date_second) result(milliseconds_between_dates)

        ! arguments

        integer,    intent(in)  ::  date_first(8), date_second(8)
        ! private variables
        integer ::  date_diff(8)
        ! start work

        date_diff = date_second - date_first

        milliseconds_between_dates = date_diff(8) + milliseconds_in_a_second * (date_diff(7) + date_diff(1) * seconds_in_a_year + &
                                     date_diff(2) * seconds_in_a_month + date_diff(3) * seconds_in_a_day + &
                                     date_diff(5) * seconds_in_an_hour + date_diff(6) * seconds_in_a_minute)

    end function MilliSecondsBetweenDates

    !>  \brief  Return the number of seconds between two dates & times
    pure integer function SecondsBetweenDates(self, date_first, date_second) result(seconds_between_dates)

        ! arguments
        class(DateAndTime),  intent(in) ::  self
        integer,    intent(in)  ::  date_first(8), date_second(8)

        ! private variables

        ! start work
        seconds_between_dates = int(self%MilliSecondsBetweenDates(date_first, date_second) / milliseconds_in_a_second)
    end function SecondsBetweenDates

    !>    \brief    Returns the date & time printed out in a nice format
    function DateAndTimeAsString() result(nice_date_and_time)

        ! arguments

        character(len=19)    ::    nice_date_and_time
        ! private variables
        character(len=8)     ::    current_date
        character(len=10)    ::    current_time

        call date_and_time(current_date, current_time)    ! this is an intrinsic function

        nice_date_and_time =   current_date(1:4) // '-' // current_date(5:6) // '-' // current_date(7:8) // ' ' // &
                               current_time(1:2) // ':' // current_time(3:4) // ':' // current_time(5:6) // ':' // current_time(7:8)

    end function DateAndTimeAsString

    !>    \brief    Convert a duration in seconds to a nice character string
    function DurationAsString(self, seconds) result(nice_time_from_seconds)

        ! arguments
        class(DateAndTime),  intent(in)   ::  self
        real                              ::    seconds
        ! result
        character(len=:),  allocatable    ::    nice_time_from_seconds
        ! private variables
        character(len=200)                ::    string
        integer                           ::    duration(4)

        ! start work

        duration = self%DurationAsIntegers(seconds)
        string = ''
        if (duration(1) .ne. 0) then
            write(string,'(i0,1x,a)') duration(1), 'day'
            if (duration(1) .gt. 1) then
                write(string,'(2a)') trim(adjustl(string)), 's,'
            else
                write(string,'(2a)') trim(adjustl(string)), ','
            endif
        endif
        if (duration(2) .ne. 0) then
            write(string,'(a,1x,i0,1x,a)') trim(adjustl(string)), duration(2), 'hour'
            if (duration(2) .gt. 1) then
                write(string,'(2a)') trim(adjustl(string)), 's,'
            else
                write(string,'(2a)') trim(adjustl(string)), ','
            endif
        endif
        if (duration(3) .ne. 0) then
            write(string,'(a,1x,i0,1x,a)') trim(adjustl(string)), duration(3), 'minute'
            if (duration(3) .gt. 1) then
                write(string,'(2a)') trim(adjustl(string)), 's and'
            else
                write(string,'(2a)') trim(adjustl(string)), ' and'
            endif
        endif
        write(string,'(a,1x,i0,1x,a)') trim(adjustl(string)), duration(4), 'second'
        if (duration(4) .ne. 1) then
            write(string,'(2a)') trim(adjustl(string)), 's'
        endif
        ! copy the result over
        string = adjustl(string)

		! i much prefer the first line, but gfortran has a bug which means
		! i have to use the second line instead, which i think is equivalent
        !allocate(character(len=len_trim(string)) :: nice_time_from_seconds)
        allocate(nice_time_from_seconds, source=(string(:len_trim(string))))

        nice_time_from_seconds = trim(string)

    end function DurationAsString

    !>    \brief    Convert a duration in seconds to a duration in days, hours, minutes and seconds
    function DurationAsIntegers(seconds) result(nice_duration)

        ! arguments
        real    ::    seconds
        ! result
        integer    ::    nice_duration(4)    !<    1. days, 2. hours, 3. minutes, 4. seconds
        ! private variables

        integer    ::    iseconds
        ! start work
        iseconds = nint(seconds)
        nice_duration(1) = iseconds/seconds_in_a_day
        nice_duration(2) = (iseconds-nice_duration(1)*seconds_in_a_day)/seconds_in_an_hour
        nice_duration(3) = (iseconds-nice_duration(1)*seconds_in_a_day-nice_duration(2)*seconds_in_an_hour)/seconds_in_a_minute
        nice_duration(4) = (iseconds-nice_duration(1)*seconds_in_a_day-nice_duration(2)*seconds_in_an_hour - &
                            nice_duration(3)*seconds_in_a_minute)
    end function DurationAsIntegers
end module DatesAndTimes
