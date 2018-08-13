!>  \brief Progress bar
module ProgressBars

    use Globals
    use DatesAndTimes

    implicit none

    private

    type,         public  ::  ProgressBar
        private
        integer             ::  total_number_of_ticks
        integer             ::  current_tick
        integer             ::  start_date(8)
        integer             ::  last_update_date(8)
        integer             ::  number_of_update_calls

    contains

        procedure,  public  ::  Begin
        procedure,  public  ::  Update
        procedure,  public  ::  Finish
        procedure,  public  ::  GetCurrentTick

    end type

    contains

    pure integer function GetCurrentTick(self)
        class(ProgressBar),     intent(in)  ::  self
        GetCurrentTick = self%current_tick
    end function GetCurrentTick


    !>  \brief  Start a new progressbar
    subroutine Begin(self, wanted_total_number_of_ticks)

        !Arguments
        class(ProgressBar),     intent(inout)   ::  self
        integer,                intent(in)      ::  wanted_total_number_of_ticks

        !Variables
        type(DateAndTime)                       ::  my_date_and_time

        self%total_number_of_ticks = wanted_total_number_of_ticks
        self%number_of_update_calls = 0
        self%current_tick = 0

        ! If there are less than 2 ticks, don't do anything, else we print the begin state

        if (wanted_total_number_of_ticks .gt. 1) then
            self%start_date = my_date_and_time%DateAndTimeAsIntegers()
            self%last_update_date = self%start_date
            write(*,*)
#ifdef __INTEL_COMPILER
            write(*, '(a,$,a)') '     0% [                              ] ???h:??m:??s     ',char(13)
#else
             write(*, '(a)',advance='no') '     0% [                              ] ???h:??m:??s     '
             write(*,'(a)')char(13)
#endif
        endif
	end subroutine Begin

	subroutine Update(self, current_tick)

        !Arguments
        class(ProgressBar),     intent(inout)   ::  self
        integer,                intent(in)      ::  current_tick

        !Variables
        integer                                 ::  seconds_remaining
		integer                                 ::  minutes_remaining
		integer                                 ::  hours_remaining

		integer                                 ::  percent_complete
		integer                                 ::  remaining_time
		integer                                 ::  filled_bar_size
		integer                                 ::  current_date(8)

		integer                                 ::  bar_position

		real                                    ::  current_seconds_per_tick

        type(DateAndTime)                       ::  my_date_and_time

        ! Start

        ! Zero or negative ticks are not allowed, and should indicate an error..

        if (current_tick .lt. 1) then
            call this_program%TerminateWithFatalError('ProgressBar::Update', 'Update called with less than 1 for current tick' )
        endif

        self%current_tick = current_tick

        ! if there are less than 2 ticks, we don't do anything


        if (self%total_number_of_ticks .gt. 1) then

            current_date = my_date_and_time%DateAndTimeAsIntegers()
            !self%number_of_update_calls = self%number_of_update_calls + 1

            ! Only proceed if we are more than 1 second from the previous update time, this
            ! is here to stop very fast things spending all their time drawing the bar.
            !
            ! if the process is fast enough, it will spend a lot of it's time asking for the
            ! current time, so the number of ticks that are called will need to be restricted


            if (my_date_and_time%SecondsBetweenDates(self%last_update_date, current_date) .ge. 1) then! .or. self%number_of_update_calls .le. 16) then

                self%last_update_date = current_date
                percent_complete = nint(100. / real(self%total_number_of_ticks) * real(current_tick))
                current_seconds_per_tick = real(my_date_and_time%SecondsBetweenDates(self%start_date, current_date))&
                                             / real(current_tick)
                remaining_time = nint((self%total_number_of_ticks - current_tick) * current_seconds_per_tick)
                filled_bar_size = nint(percent_complete * .3)

                if (remaining_time .gt. 3600) then
                    hours_remaining = remaining_time / 3600
                else
                    hours_remaining = 0
                endif
                if (remaining_time .gt. 60) then
                    minutes_remaining = (remaining_time / 60) - (hours_remaining * 60)
                else
                    minutes_remaining = 0
                endif

                seconds_remaining = remaining_time - ((hours_remaining * 60 + minutes_remaining) * 60)

                ! Sanity checking

                if (filled_bar_size > 30) filled_bar_size = 30

                ! draw the bar, starting with percent complete
#ifdef __INTEL_COMPILER
                write(*, '(a,$)') '   '
                write(*, '(i3,$)') percent_complete
                write(*, '(a,$)') '% ['
#else
                write(*, '(a)',advance='no') '   '
                write(*, '(i3)',advance='no') percent_complete
                write(*, '(a)',advance='no') '% ['
#endif

                do bar_position = 1, 30
                    if (bar_position .lt. filled_bar_size) then
#ifdef __INTEL_COMPILER
                        write(*, '(a,$)') '='
#else
                        write(*,'(a)',advance='no') '='
#endif
                    else
#ifdef __INTEL_COMPILER
                        write(*, '(a,$)') ' '
#else
                        write(*,'(a)',advance='no') ' '
#endif
                    endif
                enddo

#ifdef __INTEL_COMPILER
                write(*, '(a,$)') ']'
#else
                write(*,'(a)',advance='no') ']'
#endif

                ! ETA

                if (hours_remaining .gt. 999) then
#ifdef __INTEL_COMPILER
                    write(*, '(a,$)') '999h:99m:99s      '
#else
                    write(*,'(a)',advance='no') '999h:99m:99s      '
#endif
                else
#ifdef __INTEL_COMPILER
                    write(*, '(i3,a,$)') hours_remaining, 'h:'
#else
                    write(*,'(i3,a)',advance='no') hours_remaining, 'h:'
#endif
                    !now minutes, if less than 10 do a zero first

#ifdef __INTEL_COMPILER
                    if (minutes_remaining .lt. 10) write(*, '(a,$)') '0'
                    write(*, '(i0,a,$)') minutes_remaining, 'm:'
#else
                    if (minutes_remaining .lt. 10) write(*, '(a)',advance='no') '0'
                    write(*, '(i0,a)',advance='no') minutes_remaining, 'm:'
#endif

                    ! Same for seconds
#ifdef __INTEL_COMPILER
                    if (seconds_remaining .lt. 10) write(*, '(a,$)') '0'
                    write(*, '(i0,a,$)') seconds_remaining, 's      '
#else
                    if (seconds_remaining .lt. 10) write(*, '(a)',advance='no') '0'
                    write(*, '(i0,a)',advance='no') seconds_remaining, 's      '
#endif
                endif

                flush(6)

                ! Go back to the begining of the line
#ifdef __INTEL_COMPILER
                write(*, '(a,$)') char(13)
#else
                write(*,'(a)',advance='no') char(13)
#endif
            endif
        endif
    end subroutine Update


    subroutine Finish(self)

        !Arguments
        class(ProgressBar),     intent(inout)   ::  self

        ! Set the bar to 100% and end the line


        if (self%total_number_of_ticks .gt. 1) then
            write(*,'(a)')'   100% [==============================] done!          '
            write(*,*)
        endif

	end subroutine Finish


end module

