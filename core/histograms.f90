!>  \brief  Histogram class
!!
!!  In general, the intervals include the lower bound and exclude the upper bound
module Histograms
    implicit none
    private

    type, public :: Histogram
        private
        real                            ::  minimum_value               !<  The minimum value covered by the histogram
        real                            ::  maximum_value               !<  The maximum value covered by the histogram
        integer(kind=8),    allocatable ::  counts(:)                   !<  For each bin, the number of counts
        real,               allocatable ::  lower_bound(:)              !<  For each bin, its lower bound
        real                            ::  bin_width                   !<  The width of each bin
        logical                         ::  discard_extreme_values      !<  Do not keep track of values beyond the minimum and maximum
        integer(kind=8)                 ::  number_of_samples           !<  Total number of samples (counts)
        contains
            procedure, public   ::  Init
            procedure, public   ::  AddSampleValue
            procedure, public   ::  FirstValueWithCountAboveFractionOfMaxCount
            procedure, public   ::  LastValueWithCountAboveFractionOfMaxCount
            procedure, public   ::  QuantileFunction
            procedure           ::  UpdateBinBounds
            procedure           ::  UpdateBinWidth
            procedure           ::  Allocate
            procedure           ::  Deallocate
            procedure           ::  Reset
            procedure,  public  ::  PrintInfo
            procedure,  public  ::  GetBinWidth
    end type Histogram

    contains

    !>  \brief  Initialise a histogram
    subroutine Init(self,minimum_value,maximum_value,number_of_bins,discard_extreme_values)
        class(Histogram),           intent(inout)   ::  self
        real,                       intent(in)      ::  minimum_value
        real,                       intent(in)      ::  maximum_value
        integer,                    intent(in)      ::  number_of_bins
        logical,        optional,   intent(in)      ::  discard_extreme_values
        ! start work
        call self%Allocate(number_of_bins)
        call self%Reset()
        self%minimum_value = minimum_value
        self%maximum_value = maximum_value
        call self%UpdateBinWidth()
        call self%UpdateBinBounds()
        self%discard_extreme_values = .true.
        if (present(discard_extreme_values)) self%discard_extreme_values = discard_extreme_values
    end subroutine Init

    subroutine PrintInfo(self)
        class(Histogram),   intent(in)      ::  self
        !
        write(*,'(a)')              ' ** Histogram summary **'
        write(*,'(a,i0)')           'Number of bins = ', size(self%counts)
        write(*,'(a,f0.2)')         'Bin width = ', self%bin_width
        write(*,'(a,2(f0.3,1x))')   'Min, max = ', self%minimum_value, self%maximum_value
        write(*,'(a,2(l1,1x))')     'Allocated? ', allocated(self%counts), allocated(self%lower_bound)
        write(*,'(a)')              'Bin bounds: '
        write(*,*)                  self%lower_bound(:)
        write(*,'(a)')              'Counts: '
        write(*,*)                  self%counts(:)
    end subroutine PrintInfo

    !>  \brief  Add a sample value to the histogram
    subroutine AddSampleValue(self,sample_value)
        class(Histogram),   intent(inout)   ::  self
        real,               intent(in)      ::  sample_value
        ! private variable
        integer ::  which_bin
        ! start work
        if (sample_value .lt. self%minimum_value) then
            if (self%discard_extreme_values) then
                which_bin = 0
            else
                which_bin = 1
            endif
        else if (sample_value .gt. self%maximum_value) then
            if (self%discard_extreme_values) then
                which_bin = 0
            else
                which_bin = size(self%counts)
            endif
        else
            which_bin = min(int((sample_value - self%minimum_value)/self%bin_width)+1,size(self%counts))
        endif
        if (which_bin .ne. 0) self%counts(which_bin) = self%counts(which_bin) + 1
        !
        self%number_of_samples = self%number_of_samples + 1
    end subroutine AddSampleValue

    !>  \brief  Return the bin width
    real function GetBinWidth(self)
        class(Histogram),   intent(in)      ::  self
        GetBinWidth = self%bin_width
    end function

    !>  \brief  Estimate the threshold value at which the fraction of counts at or below that threshold is equal
    !!          to the given fraction
    !!  See http://mathworld.wolfram.com/QuantileFunction.html for a formal definition
    real function QuantileFunction(self,desired_fraction_at_or_below_threshold,conservative_low,conservative_high) result(threshold)
        class(Histogram),               intent(in)  ::  self
        real,                           intent(in)  ::  desired_fraction_at_or_below_threshold  !<  From 0.0 to 1.0
        logical,            optional,   intent(in)  ::  conservative_low                        !<  Return a conservative (on the low side) estimate of the threshold
        logical,            optional,   intent(in)  ::  conservative_high                       !<  Return a conservative (on the high side) estimate of the threshold
        ! private variables
        integer         ::  current_bin
        integer(kind=8) ::  total_counts, old_total_counts
        real            ::  desired_count_at_or_below_threshold
        logical         ::  cconservative_high, cconservative_low
        ! start work

        cconservative_high = .false.
        cconservative_low  = .false.
        if (present(conservative_high)) cconservative_high = conservative_high
        if (present(conservative_low )) cconservative_low  = conservative_low

        desired_count_at_or_below_threshold = desired_fraction_at_or_below_threshold * self%number_of_samples

        total_counts = 0
        do current_bin=1,size(self%counts)
            old_total_counts    = total_counts
            total_counts        = total_counts + self%counts(current_bin)
            if (total_counts .gt. desired_count_at_or_below_threshold) then
                ! By the end of the current bin, we have exceeded our target number of counts
                ! Let's work out where within that bin we met the target number of counts...
                if (cconservative_low) then
                    threshold = self%lower_bound(current_bin)
                else if (cconservative_high) then
                    threshold = self%lower_bound(min(current_bin+1,size(self%counts)))
                else
                    ! ... by linear interpolation
                    threshold =   self%lower_bound(current_bin) &
                                + self%bin_width * (desired_count_at_or_below_threshold - old_total_counts) &
                                                 / (total_counts - old_total_counts)
                endif
                exit
            endif
        enddo

    end function QuantileFunction

    !>  \brief  Return the first sample value at which the count is a given percentage of the maximum observed count
    real function FirstValueWithCountAboveFractionOfMaxCount(self,fraction_of_max_count) result(first_value)
        class(Histogram),   intent(in)      ::  self
        real,               intent(in)      ::  fraction_of_max_count
        ! private variables
        real    ::  threshold
        integer ::  current_bin
        ! start work
        threshold = maxval(self%counts) * fraction_of_max_count

        do current_bin=1,size(self%counts)
            if (self%counts(current_bin) .ge. threshold) exit
        enddo

        current_bin = min(current_bin,size(self%lower_bound))

        first_value = self%lower_bound(current_bin)
    end function FirstValueWithCountAboveFractionOfMaxCount

    !>  \brief  Return the last sample value at which the count is a given percentage of the maximum observed count
    real function LastValueWithCountAboveFractionOfMaxCount(self,fraction_of_max_count) result(last_value)
        class(Histogram),   intent(in)      ::  self
        real,               intent(in)      ::  fraction_of_max_count
        ! private variables
        real    ::  threshold
        integer ::  current_bin
        ! start work
        threshold = maxval(self%counts) * fraction_of_max_count

        do current_bin=size(self%counts),1,-1
            if (self%counts(current_bin) .ge. threshold) exit
        enddo

        current_bin = max(current_bin,1)

        if (current_bin .eq. size(self%counts)) then
            last_value = self%maximum_value
        else
            last_value = self%lower_bound(current_bin+1)
        endif
    end function LastValueWithCountAboveFractionOfMaxCount

    !>  \brief  Compute the bounds of each bin
    subroutine UpdateBinBounds(self)
        class(Histogram),   intent(inout)   ::  self
        ! private variables
        integer ::  current_bin
        ! start work
        do current_bin=1,size(self%lower_bound)
            self%lower_bound(current_bin) = self%minimum_value + self%bin_width * (current_bin-1)
        enddo
    end subroutine UpdateBinBounds

    !>  \brief  Return the width of the bins
    subroutine UpdateBinWidth(self)
        class(Histogram),   intent(inout)   ::  self
        self%bin_width = (self%maximum_value-self%minimum_value)/real(size(self%lower_bound),kind=8)
    end subroutine UpdateBinWidth

    !>  \brief  Deallocate memory
    subroutine Deallocate(self)
        class(Histogram),   intent(inout)   ::  self
        if (allocated(self%counts)) deallocate(self%counts)
        if (allocated(self%lower_bound)) deallocate(self%lower_bound)
    end subroutine Deallocate

    !>  \brief  Allocate memory
    subroutine Allocate(self,number_of_bins)
        class(Histogram),   intent(inout)   ::  self
        integer,            intent(in)      ::  number_of_bins
        call self%Deallocate()
        allocate(self%counts(number_of_bins),self%lower_bound(number_of_bins))
    end subroutine Allocate

    !>  \brief  Reset memory
    subroutine Reset(self)
        class(Histogram),   intent(inout)   ::  self
        self%counts = 0
        self%lower_bound = 0.0d0
        self%bin_width = 0.0d0
        self%discard_extreme_values = .true.
        self%number_of_samples = 0
    end subroutine Reset


end module Histograms
