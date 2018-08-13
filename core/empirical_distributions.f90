!>  \brief  Class to facilitate characterization of distributions of observed variables
module EmpiricalDistributions

    use Globals

    implicit none
    private

    type, public :: EmpiricalDistribution
        private
        real(kind=8)        ::  sum_of_samples
        real(kind=8)        ::  sum_of_squared_samples
        integer(kind=8)     ::  number_of_samples
        real                ::  minimum
        real                ::  maximum
        logical             ::  keep_sample_values = .false.
        real,   allocatable ::  sample_values(:)
        contains
            procedure, public   ::  Init
            procedure, public   ::  AddSampleValue
            procedure, public   ::  GetNumberOfSamples
            procedure, public   ::  GetSampleSum
            procedure, public   ::  GetSampleMean
            procedure, public   ::  GetSampleVariance
            procedure, public   ::  GetSampleSumOfSquares
            procedure, public   ::  GetUnbiasedEstimateOfPopulationVariance
            procedure, public   ::  PopulateHistogram
    end type EmpiricalDistribution


    contains

    !>  \brief  Initialisation
    pure subroutine Init(self,keep_sample_values)
        class(EmpiricalDistribution),               intent(inout)   ::  self
        logical,                        optional,   intent(in)      ::  keep_sample_values  !<  this is necessary when (for example) one wants to compute a histogram
        ! start work
        self%sum_of_samples = 0.0d0
        self%sum_of_squared_samples = 0.0d0
        self%number_of_samples = 0
        self%keep_sample_values = .false.
        self%minimum = huge(1.0e0)
        self%maximum = -huge(1.0e0)
        if (present(keep_sample_values)) self%keep_sample_values = keep_sample_values
        if (self%keep_sample_values) then
            if (allocated(self%sample_values)) deallocate(self%sample_values)
            allocate(self%sample_values(1024))
        endif
    end subroutine Init


    !>  \brief  Add an observed datum
    pure subroutine AddSampleValue(self,sample_value)
        class(EmpiricalDistribution),   intent(inout)   ::  self
        real,                           intent(in)      ::  sample_value
        ! private variables
        real, allocatable    :: temp(:)
        ! start work
        self%sum_of_samples = self%sum_of_samples + sample_value
        self%sum_of_squared_samples = self%sum_of_squared_samples + sample_value**2
        self%number_of_samples = self%number_of_samples + 1
        self%minimum = min(self%minimum,sample_value)
        self%maximum = max(self%maximum,sample_value)
        ! We may need to record the value
        if (self%keep_sample_values) then
            if (size(self%sample_values) .lt. self%number_of_samples) then
                ! We need to allocate more memory
                temp = self%sample_values
                deallocate(self%sample_values)
                allocate(self%sample_values((self%number_of_samples-1)*2))
                self%sample_values(1:size(temp)) = temp(:)
                deallocate(temp)
            endif
            self%sample_values(self%number_of_samples) = sample_value
        endif
    end subroutine AddSampleValue

    !>  \brief  Populate a histogram
    subroutine PopulateHistogram(self,my_histogram,number_of_bins)
        use Histograms
        class(EmpiricalDistribution),   intent(in)      ::  self
        type(Histogram),                intent(inout)   ::  my_histogram
        integer,                        intent(in)      ::  number_of_bins
        ! private variables
        integer(kind=8) ::  current_sample
        ! start work
        if (.not. self%keep_sample_values) then
            call this_program%TerminateWithFatalError('EmpiricalDistributions::PopulateHistogam','Sample values were not kept')
        endif
        call my_histogram%Init(self%minimum,self%maximum,number_of_bins)
        do current_sample=1,self%number_of_samples
            call my_histogram%AddSampleValue(self%sample_values(current_sample))
        enddo
    end subroutine PopulateHistogram


    !>  \brief  Return the sum of squares
    pure real function GetSampleSumOfSquares(self)
        class(EmpiricalDistribution),   intent(in)      ::  self
        GetSampleSumOfSquares = self%sum_of_squared_samples
    end function GetSampleSumOfSquares

    !>  \brief  Return the number of samples
    pure integer function GetNumberOfSamples(self)
        class(EmpiricalDistribution),   intent(in)      ::  self
        GetNumberOfSamples = self%number_of_samples
    end function GetNumberOfSamples

    !>  \brief  Return the sample mean
    pure real function GetSampleSum(self)
        class(EmpiricalDistribution),   intent(in)      ::  self
        ! start work
        GetSampleSum = self%sum_of_samples
    end function GetSampleSum

    !>  \brief  Return the sample mean
    pure real function GetSampleMean(self)
        class(EmpiricalDistribution),   intent(in)      ::  self
        ! start work
        if (self%number_of_samples .gt. 0) then
            GetSampleMean = self%sum_of_samples / self%number_of_samples
        else
            GetSampleMean = 0.0
        endif
    end function GetSampleMean

    !>  \brief  Return the sample variance
    pure real function GetSampleVariance(self)
        class(EmpiricalDistribution),   intent(in)      ::  self
        ! start work
        if (self%number_of_samples .gt. 0) then
            GetSampleVariance = (self%sum_of_squared_samples / self%number_of_samples) - &
                                (self%sum_of_samples / self%number_of_samples)**2
        else
            GetSampleVariance = 0.0
        endif
    end function GetSampleVariance

    !>  \brief  Return an unbiased estimate of the population variance
    pure real function GetUnbiasedEstimateOfPopulationVariance(self) result(variance_estimate)
        class(EmpiricalDistribution),   intent(in)      ::  self
        ! start work
        if (self%number_of_samples .gt. 0) then
            variance_estimate = self%GetSampleVariance() * self%number_of_samples / (self%number_of_samples - 1)
        else
            variance_estimate = 0.0
        endif
    end function GetUnbiasedEstimateOfPopulationVariance


end module EmpiricalDistributions
