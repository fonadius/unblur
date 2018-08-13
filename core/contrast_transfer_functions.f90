!>\brief    Contrast transfer functions for electron micrographs
module ContrastTransferFunctions
    use Globals
    use Units
    implicit none
    private
    public :: ContrastTransferFunction

    type :: ContrastTransferFunction
        private
        real    ::  spherical_aberration                =   0.0             !<  Spherical aberration
        integer ::  spherical_aberration_units          =   millimeters
        real    ::  wavelength                          =   0.0             !<  Wavelength of the electrons
        integer ::  wavelength_units                    =   angstroms
        real    ::  amplitude_contrast                  =   0.0             !<  Fraction of amplitude contrast (from 0. to 1., usually btw .07 and .14, see Mindell 03)
        real    ::  defocus_1                           =   0.0             !<  Underfocus along first axis  (positive values for underfocus; larger value = more underfocus)
        real    ::  defocus_2                           =   0.0             !<  Underfocus along second axis
        integer ::  defocus_1_units                     =   microns
        integer ::  defocus_2_units                     =   microns
        real    ::  defocus_half_range                  =   0.0             !<  Half-range of defocus variation during micrograph acquisition
        integer ::  defocus_half_range_units            =   microns
        real    ::  astigmatism_azimuth                 =   0.0             !<  Azimuth of first axis. 0.0 means axis is at 3 o'clock
        integer ::  astigmatism_azimuth_units           =   degrees
        real    ::  additional_phase_shift              =   0.0             !<  Additional phase shift of scattered electrons relative to unscattered electrons. This can be introduced e.g. by a Volta phase plate
        integer ::  additional_phase_shift_units        =   radians
        ! Fitting parameters
        real    ::  lowest_frequency_for_fitting        =   0.0
        integer ::  lowest_frequency_for_fitting_units  =   reciprocal_pixels
        real    ::  highest_frequency_for_fitting       =   0.5
        integer ::  highest_frequency_for_fitting_units =   reciprocal_pixels
        real    ::  astigmatism_tolerance               =   0.0                 !<  Expected (tolerated) astigmatism. During parameter search, astigmatism values much larger than this will be penalised against. Set to 0.0 to ignore this restraint.
        integer ::  astigmatism_tolerance_units         =   angstroms
        contains
            procedure,  public      ::  Init
            procedure,  public      ::  PrintInfo
            procedure,  public      ::  EvaluateAtSquaredSpatialFrequency
            procedure,  public      ::  CountNumberOfExtremaBeforeSquaredSpatialFrequency
            procedure,  public      ::  ComputeFrequencyOfAZero
            !
            procedure               ::  SetDefocusScalars
            procedure               ::  SetDefocusArray
            generic,    public      ::  SetDefocus => SetDefocusArray, SetDefocusScalars
            procedure,  public      ::  SetDefocusHalfRange
            procedure,  public      ::  SetAdditionalPhaseShift
            procedure,  public      ::  GetAdditionalPhaseShift
            procedure,  public      ::  HasPixelUnits
            procedure,  public      ::  SetLowestFrequencyForFitting
            procedure,  public      ::  GetLowestFrequencyForFitting
            procedure,  public      ::  GetHighestFrequencyForFitting
            procedure,  public      ::  GetAstigmatismTolerance
            procedure,  public      ::  GetAstigmatism
            procedure,  public      ::  GetDefocusParametersInAngstromsAndDegrees
            procedure,  public      ::  SetDefocusParametersInAngstromsAndDegrees
            procedure,  public      ::  GetDefocusParameters
            procedure               ::  ConvertToPixelsAndRadians
            procedure,  public      ::  GetDefocus1
            procedure,  public      ::  GetDefocus2
            procedure,  public      ::  GetDefocus1InAngstroms
            procedure,  public      ::  GetDefocus2InAngstroms
            procedure,  public      ::  GetAstigmatismAzimuthInDegrees
            procedure,  public      ::  GetAstigmatismAzimuthInRadians
    end type

    contains

    !>  \brief  Initialise a CTF object
    subroutine Init(self,acceleration_voltage,spherical_aberration,amplitude_contrast, &
                         defocus_1,defocus_2,astigmatism_azimuth, &
                         lowest_frequency_for_fitting,highest_frequency_for_fitting, &
                         astigmatism_tolerance, &
                         pixel_size,additional_phase_shift, &
                         defocus_half_range)
        ! arguments
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  acceleration_voltage            !<  Kilovolts
        real,                               intent(in)      ::  spherical_aberration            !<  mm
        real,                               intent(in)      ::  amplitude_contrast              !<  fraction
        real,                               intent(in)      ::  defocus_1                       !<  um
        real,                               intent(in)      ::  defocus_2                       !<  um
        real,                               intent(in)      ::  astigmatism_azimuth             !<  degrees
        real,                               intent(in)      ::  lowest_frequency_for_fitting    !<  1/A
        real,                               intent(in)      ::  highest_frequency_for_fitting   !<  1/A
        real,                               intent(in)      ::  astigmatism_tolerance           !<  A. Set to negative to indicate no restraint on astigmatism.
        real,                               intent(in)      ::  pixel_size                      !<  A
        real,                               intent(in)      ::  additional_phase_shift          !<  radians
        real,           optional,           intent(in)      ::  defocus_half_range              !<  um. Half-range of defocus variation within/during micrograph acquisition
        ! start work
        self%wavelength = akv_to_wl(acceleration_voltage)
        self%wavelength_units = angstroms
        self%spherical_aberration = spherical_aberration
        self%spherical_aberration_units = millimeters
        self%amplitude_contrast = amplitude_contrast
        self%defocus_1 = defocus_1
        self%defocus_2 = defocus_2
        self%defocus_1_units = microns
        self%defocus_2_units = microns
        self%astigmatism_azimuth = astigmatism_azimuth
        self%additional_phase_shift = additional_phase_shift
        self%additional_phase_shift_units = radians
        !if (present(lowest_frequency_for_fitting)) then
            self%lowest_frequency_for_fitting = lowest_frequency_for_fitting
            self%lowest_frequency_for_fitting_units = reciprocal_angstroms
        !endif
        !if (present(highest_frequency_for_fitting)) then
            self%highest_frequency_for_fitting = highest_frequency_for_fitting
            self%highest_frequency_for_fitting_units = reciprocal_angstroms
        !endif
        !if (present(astigmatism_tolerance)) then
        if (astigmatism_tolerance .lt. 0.0) then
            ! we will turn off the restraint on astigmatism
            self%astigmatism_tolerance = 0.0
        else if (astigmatism_tolerance .lt. 10.0) then
            self%astigmatism_tolerance = 10.0
        else
            self%astigmatism_tolerance = astigmatism_tolerance
        endif
        self%astigmatism_tolerance_units = angstroms
        if (present(defocus_half_range)) then
            self%defocus_half_range = defocus_half_range
        else
            self%defocus_half_range = 0.0
        endif
        self%defocus_half_range_units = microns
        !endif
        !if (present(pixel_size)) then
            call self%ConvertToPixelsAndRadians(pixel_size)
        !endif
    end subroutine Init

    subroutine SetLowestFrequencyForFitting(self,lowest_frequency_for_fitting)
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  lowest_frequency_for_fitting    !<  1/pixels

        self%lowest_frequency_for_fitting = lowest_frequency_for_fitting
        self%lowest_frequency_for_fitting_units = reciprocal_pixels
    end subroutine SetLowestFrequencyForFitting

    !>  \brief  Set the additional phase shift experienced by scattered electrons relative to unscattered electrons
    subroutine SetAdditionalPhaseShift(self,additional_phase_shift)
        ! arguments
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  additional_phase_shift  !<  radians
        ! private variablews

        ! start work
        self%additional_phase_shift = additional_phase_shift
        self%additional_phase_shift_units = radians
    end subroutine SetAdditionalPhaseShift

    !>  \brief  Get the additional phase shift experienced by scattered electrons relative to unscattered electrons (radians)
    pure real function GetAdditionalPhaseShift(self)
        ! arguments
        class(ContrastTransferFunction),    intent(in)   ::  self
        ! private variablews

        ! start work
        GetAdditionalPhaseShift = self%additional_phase_shift
    end function GetAdditionalPhaseShift

    !>  \brief  Set defocus half range (this affects the CTF envelope)
    subroutine SetDefocusHalfrange(self,defocus_half_range)
        ! arguments
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  defocus_half_range      !<  pixels
        !
        self%defocus_half_range = defocus_half_range
        self%defocus_half_range_units = pixels
    end subroutine SetDefocusHalfrange

    !>  \brief  Set defocus and astigmatism in pixels and radians
    subroutine SetDefocusScalars(self,defocus_1,defocus_2,astigmatism_azimuth)
        ! arguments
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  defocus_1               !<  pixels
        real,                               intent(in)      ::  defocus_2               !<  pixels
        real,                               intent(in)      ::  astigmatism_azimuth     !<  radians
        !
        self%defocus_1 = defocus_1
        self%defocus_1_units = pixels
        self%defocus_2 = defocus_2
        self%defocus_2_units = pixels
        self%astigmatism_azimuth = astigmatism_azimuth
        self%astigmatism_azimuth_units = radians
    end subroutine SetDefocusScalars
    subroutine SetDefocusArray(self,my_array)
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  my_array(:)
        call self%SetDefocusScalars(my_array(1),my_array(2),my_array(3))
    end subroutine SetDefocusArray
    subroutine SetDefocusParametersInAngstromsAndDegrees(self,defocus_parameters,pixel_size)
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  defocus_parameters(3)
        real,                               intent(in)      ::  pixel_size
        self%defocus_1 = defocus_parameters(1) / pixel_size
        self%defocus_1_units = pixels
        self%defocus_2 = defocus_parameters(2) / pixel_size
        self%defocus_2_units = pixels
        self%astigmatism_azimuth = convert(defocus_parameters(3),degrees,radians)
        self%astigmatism_azimuth_units = radians
    end subroutine

    function GetDefocusParameters(self) result(defocus)
        class(ContrastTransferFunction),    intent(in)      ::  self
        real                                                ::  defocus(3)
        defocus = [self%defocus_1,self%defocus_2,self%astigmatism_azimuth]
    end function
    function GetDefocusParametersInAngstromsAndDegrees(self,pixel_size) result (defocus)
        class(ContrastTransferFunction),    intent(in)      ::  self
        real                                                ::  defocus(3)
        real,                               intent(in)      ::  pixel_size
        defocus(1) = convert(self%defocus_1,self%defocus_1_units,angstroms,pixel_size)
        defocus(2) = convert(self%defocus_2,self%defocus_2_units,angstroms,pixel_size)
        defocus(3) = self%GetAstigmatismAzimuthInDegrees()
    end function
    real function GetDefocus1InAngstroms(self,pixel_size)
        class(ContrastTransferFunction),    intent(in)  ::  self
        real,   optional,                   intent(in)  ::  pixel_size
        GetDefocus1InAngstroms = convert(self%defocus_1,self%defocus_1_units,angstroms,pixel_size)
    end function
    real function GetDefocus2InAngstroms(self, pixel_size)
        class(ContrastTransferFunction),    intent(in)  ::  self
        real,   optional,                   intent(in)  ::  pixel_size
        GetDefocus2InAngstroms = convert(self%defocus_2,self%defocus_2_units,angstroms,pixel_size)
    end function
    real function GetDefocus1(self)
        class(ContrastTransferFunction),    intent(in)  ::  self
        GetDefocus1 = self%defocus_1
    end function
    real function GetDefocus2(self)
        class(ContrastTransferFunction),    intent(in)  ::  self
        GetDefocus2 = self%defocus_2
    end function
    real function GetAstigmatismAzimuthInDegrees(self) result(azimuth)
        class(ContrastTransferFunction),    intent(in)  ::  self
        azimuth = convert(self%astigmatism_azimuth,self%astigmatism_azimuth_units,degrees)
        azimuth = azimuth - 180.0e0 * nint(azimuth/180.0e0)
    end function
    real function GetAstigmatismAzimuthInRadians(self) result(azimuth)
        class(ContrastTransferFunction),    intent(in)  ::  self
        azimuth = convert(self%astigmatism_azimuth,self%astigmatism_azimuth_units,radians)
    end function

    subroutine PrintInfo(self)
        ! arguments
        class(ContrastTransferFunction),    intent(in)      ::  self
        ! start work
        write(*,'(a)')              '** ContrastTransferFunction **'
        write(*,'(a,f0.3,1x,a)')    'Wavelength = ', self%wavelength, unit_to_string(self%wavelength_units)
        write(*,'(a,f0.3,1x,a)')    'Spherical aberration = ', self%spherical_aberration, &
                                                                unit_to_string(self%spherical_aberration_units)
        write(*,'(a,f0.3)')         'Amplitude contrast = ', self%amplitude_contrast
        write(*,'(a,f0.3,1x,a)')    'Defocus 1 = ', self%defocus_1, unit_to_string(self%defocus_1_units)
        write(*,'(a,f0.3,1x,a)')    'Defocus 2 = ', self%defocus_2, unit_to_string(self%defocus_2_units)
        write(*,'(a,f0.3,1x,a)')    'Azimuth of astigmatism = ', self%astigmatism_azimuth, &
                                                                unit_to_string(self%astigmatism_azimuth_units)
        write(*,'(a,f0.3,1x,a)')    'Defocus half-range = ', self%defocus_half_range, &
                                                                unit_to_string(self%defocus_half_range_units)
        write(*,'(a,f0.3,1x,a)')    'Additional phase shift = ', self%additional_phase_shift, &
                                                                unit_to_string(self%additional_phase_shift_units)
        write(*,'(2(a,f0.3,a))')    'Frequencies for fitting: from ', self%lowest_frequency_for_fitting, &
                                    unit_to_string(self%lowest_frequency_for_fitting_units), ' to ', &
                                    self%highest_frequency_for_fitting, unit_to_string(self%highest_frequency_for_fitting_units)
        write(*,'(a,f0.1,1x,a)')    'Astigmatism tolerance = ', self%astigmatism_tolerance, &
                                                                unit_to_string(self%astigmatism_tolerance_units)
    end subroutine PrintInfo

    !>  \brief  Check that all units are in pixels / radians
    pure logical function HasPixelUnits(self)
        class(ContrastTransferFunction),    intent(in)      ::  self
        HasPixelUnits =         self%spherical_aberration_units == pixels &
                        .and.   self%wavelength_units == pixels &
                        .and.   self%defocus_1_units == pixels &
                        .and.   self%defocus_2_units == pixels &
                        .and.   self%astigmatism_azimuth_units == radians &
                        .and.   self%highest_frequency_for_fitting_units == reciprocal_pixels &
                        .and.   self%lowest_frequency_for_fitting_units == reciprocal_pixels &
                        .and.   self%astigmatism_tolerance_units == pixels &
                        .and.   self%additional_phase_shift_units == radians &
                        .and.   self%defocus_half_range_units == pixels
    end function HasPixelUnits

    pure real function GetLowestFrequencyForFitting(self) result(lowest_frequency)
        class(ContrastTransferFunction),    intent(in)  ::  self
        lowest_frequency = self%lowest_frequency_for_fitting
    end function

    pure real function GetHighestFrequencyForFitting(self) result(highest_frequency)
        class(ContrastTransferFunction),    intent(in)  ::  self
        highest_frequency = self%highest_frequency_for_fitting
    end function

    pure real function GetAstigmatismTolerance(self) result(astigmatism_tolerance)
        class(ContrastTransferFunction),    intent(in)  ::  self
        astigmatism_tolerance = self%astigmatism_tolerance
    end function

    pure real function GetAstigmatism(self) result(astigmatism)
        class(ContrastTransferFunction),    intent(in)  ::  self
        astigmatism = self%defocus_1 - self%defocus_2
    end function


    !>  \brief  Convert all values to pixel/radians units
    subroutine ConvertToPixelsAndRadians(self,pixel_size)
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  pixel_size
        ! start work
        call unit_conversion(self%spherical_aberration,self%spherical_aberration_units,pixels,pixel_size)
        call unit_conversion(self%wavelength,self%wavelength_units,pixels,pixel_size)
        call unit_conversion(self%defocus_1,self%defocus_1_units,pixels,pixel_size)
        call unit_conversion(self%defocus_2,self%defocus_2_units,pixels,pixel_size)
        call unit_conversion(self%astigmatism_azimuth,self%astigmatism_azimuth_units,radians)
        call unit_conversion(self%defocus_half_range,self%defocus_half_range_units,pixels,pixel_size)
        call unit_conversion(self%additional_phase_shift,self%additional_phase_shift_units,radians)
        call unit_conversion(self%lowest_frequency_for_fitting,self%lowest_frequency_for_fitting_units, &
                            reciprocal_pixels,pixel_size)
        call unit_conversion(self%highest_frequency_for_fitting,self%highest_frequency_for_fitting_units, &
                            reciprocal_pixels,pixel_size)
        call unit_conversion(self%astigmatism_tolerance,self%astigmatism_tolerance_units,pixels,pixel_size)
    end subroutine ConvertToPixelsAndRadians

    !>  \brief wrapper around eval_ctf which takes a ctf object as argument
    pure elemental function EvaluateAtSquaredSpatialFrequency(self,squared_spatial_frequency,azimuth, &
                                                                return_sign_only) result(ctf)
        ! arguments
        class(ContrastTransferFunction),intent(in)      ::  self
        real,                           intent(in)      ::  squared_spatial_frequency   !<  squared reciprocal pixels (0.25 is Nyquist)
        real,                           intent(in)      ::  azimuth                     !<  radians. azimuth at which to evalulate the ctf
        logical,        optional,       intent(in)      ::  return_sign_only            !<  return the sign of the ctf rather than its value
        !
        real(kind=4)                                    ::  ctf
        !
        ctf      = eval_ctf_slave(  self%spherical_aberration,self%wavelength,self%amplitude_contrast,    &
                                    self%defocus_1,self%defocus_2,self%astigmatism_azimuth, &
                                    self%defocus_half_range, &
                                    self%additional_phase_shift, &
                                    squared_spatial_frequency,azimuth,return_sign_only)
    end function EvaluateAtSquaredSpatialFrequency

    !>  \brief  Return the spatial frequency (reciprocal pixels) of the n-th zero, where n is given as which_zero
    real function ComputeFrequencyOfAZero(self,azimuth,which_zero)
        class(ContrastTransferFunction),    intent(in)  ::  self
        real,                               intent(in)  ::  azimuth         !<  Radians. Azimuth at which to evaluate the CTF.
        integer,                            intent(in)  ::  which_zero      !<  Number (starting at 1) of the zero you are interested in.
        !
        real,   allocatable ::  phase_shifts_of_zeroes(:)
        real,   allocatable ::  sq_sfs_of_zeroes(:)
        integer             ::  num_freqs
        ! start work

        ! allocate memory for temp results
        allocate(phase_shifts_of_zeroes(which_zero*8))
        allocate(sq_sfs_of_zeroes(2*size(phase_shifts_of_zeroes)))

        ! Solve the CTF for the phase shift - we will get lots of phase shifts, some of which will not be possible with current CTF
        call ctf_solve_for_phase_shift(self%amplitude_contrast,phase_shifts_of_zeroes,0.0e0)

        call ComputeSquaredFrequenciesGivenPhaseShifts(self,phase_shifts_of_zeroes,azimuth,sq_sfs_of_zeroes,num_freqs)

        ComputeFrequencyOfAZero = sqrt(sq_sfs_of_zeroes(which_zero))
    end function ComputeFrequencyOfAZero

    !>  \brief  Compute set of squared spatial frequencies at which the given phase shift is obtained by the CTF
    subroutine ComputeSquaredFrequenciesGivenPhaseShifts(self,phase_shifts,azimuth,squared_spatial_frequencies,number_of_solutions)
        use UsefulFunctions, only : QuickSort
        ! arguments
        type(ContrastTransferFunction), intent(in)      ::  self
        real,                           intent(in)      ::  phase_shifts(:)
        real,                           intent(in)      ::  azimuth
        real,                           intent(inout)   ::  squared_spatial_frequencies(:)  !<  Assumed to be 2x size of phase_shift array
        integer,                        intent(out)     ::  number_of_solutions
        ! private variables
        integer ::  i
        real    ::  current_sols(2)
        integer ::  current_num_sols
        ! start work

        number_of_solutions = 0
        do i=1,size(phase_shifts)
            call ctf_sq_sf_from_phase_shift(self%spherical_aberration,self%wavelength, &
                                            self%defocus_1,self%defocus_2,self%astigmatism_azimuth, &
                                            self%additional_phase_shift, &
                                            phase_shifts(i),azimuth, &
                                            current_sols,current_num_sols)
            if (current_num_sols .gt. 0) then
                squared_spatial_frequencies(number_of_solutions+1:number_of_solutions+current_num_sols) &
                                        = current_sols(1:current_num_sols)
                number_of_solutions = number_of_solutions + current_num_sols
            endif
        enddo

        ! Sort
        call QuickSort(squared_spatial_frequencies(1:number_of_solutions))
    end subroutine ComputeSquaredFrequenciesGivenPhaseShifts

    !>  \brief  Count how many extrema the CTF goes through before reaching the given spatial frequency at the given azimuth
    integer function CountNumberOfExtremaBeforeSquaredSpatialFrequency(self,squared_spatial_frequency,azimuth) &
                                                                result(number_of_extrema)
        class(ContrastTransferFunction),intent(in)      ::  self
        real,                           intent(in)      ::  squared_spatial_frequency   !<  squared reciprocal pixels (0.25 is Nyquist)
        real,                           intent(in)      ::  azimuth                     !<  radians. azimuth at which to evalulate the ctf
        ! private variable
        integer,    parameter   ::  maximum_number_of_rings = 128
        real,       allocatable ::  phase_shifts_of_minima(:) ! Phase shifts at which minima (-1.0) occur
        real,       allocatable ::  phase_shifts_of_maxima(:) ! Phase shifts at which maxima (+1.0) occur
        real,       allocatable ::  sq_sfs_of_minima(:) ! Squared spa freqs at which minima (-1.0) occur
        real,       allocatable ::  sq_sfs_of_maxima(:) ! Squared spa freqs at which maxima (+1.0) occur
        integer                 ::  num_minima
        integer                 ::  num_maxima
        !integer                 ::  i
        ! start work

        ! Allocate memory for temporary results
        allocate(phase_shifts_of_minima(maximum_number_of_rings))
        allocate(phase_shifts_of_maxima(maximum_number_of_rings))
        allocate(sq_sfs_of_maxima(maximum_number_of_rings*2))
        allocate(sq_sfs_of_minima(maximum_number_of_rings*2))

        ! Solve the CTF equation for phase shifts
        call ctf_solve_for_phase_shift(self%amplitude_contrast,phase_shifts_of_minima,-1.0e0)
        call ctf_solve_for_phase_shift(self%amplitude_contrast,phase_shifts_of_maxima,+1.0e0)


        ! Convert phase shifts to squared spatial frequencies
        call ComputeSquaredFrequenciesGivenPhaseShifts(self,phase_shifts_of_minima,azimuth,sq_sfs_of_minima,num_minima)
        call ComputeSquaredFrequenciesGivenPhaseShifts(self,phase_shifts_of_maxima,azimuth,sq_sfs_of_maxima,num_maxima)

        ! Let's count (note that now, phase_shifts_of_minima actually store squared spatial frequencies)
        number_of_extrema = count(sq_sfs_of_minima(1:num_minima) .le. squared_spatial_frequency &
                            .and. sq_sfs_of_minima(1:num_minima) .gt. 0.0e0)
        number_of_extrema = count(sq_sfs_of_maxima(1:num_maxima) .le. squared_spatial_frequency &
                            .and. sq_sfs_of_maxima(1:num_maxima) .gt. 0.0e0) &
                            + number_of_extrema

        ! We get pairs of solutions for sq sp freq, so let's divide by 2
        number_of_extrema = number_of_extrema / 2

    end function CountNumberOfExtremaBeforeSquaredSpatialFrequency

    !>  \brief returns the ctf, based on ctffind3 subroutine (see mindell 2003)
    pure elemental real function eval_ctf_slave(cs,wl,ampl_cont,dfmid1,dfmid2,angast,defocus_half_range,additional_phase_shift,&
                                                spa_freq_sq,ang,sign_only)
        real,                   intent(in)  ::  cs          !< spherical aberation (pixels)
        real,                   intent(in)  ::  wl          !< electron wavelength (pixels)
        real,                   intent(in)  ::  ampl_cont   !< fraction of amplitude contrast (from .07 to .14, see mindell 2003)
        real,                   intent(in)  ::  dfmid1      !< defocus along first axis (pixels)
        real,                   intent(in)  ::  dfmid2      !< defocus along second axis (for astigmatic ctf, dfmid1 .ne. dfmid2) (pixels)
        real,                   intent(in)  ::  angast      !< azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,                   intent(in)  ::  defocus_half_range      !<   half-range of defocus (pixels) during micrograph acquisition
        real,                   intent(in)  ::  additional_phase_shift  !<  In radians, additional phase shift experienced by scattered electrons, compared to scattered electrons
        real,                   intent(in)  ::  spa_freq_sq !< squared spatial frequency at which to compute the ctf (1/pixels^2)
        real,                   intent(in)  ::  ang         !< angle at which to compute the ctf (radians)
        logical,    optional,   intent(in)  ::  sign_only   !<  return 1.0 or -1.0

        ! private variables
        real    ::  wgh1    ! phase contrast
        real    ::  wgh2    ! amplitude contrast
        real    ::  phase_shift
        real    ::  b
        ! start work

        ! compute phase and amplitude contrast
        wgh1    =   sqrt(1-ampl_cont**2)
        wgh2    =   ampl_cont

        phase_shift = eval_ctf_phase_shift(cs,wl,dfmid1,dfmid2,angast,additional_phase_shift,spa_freq_sq,ang)
        eval_ctf_slave = - wgh1 * sin(phase_shift) - wgh2 * cos(phase_shift)

        ! Apply attenuation due to defocus uncertainty/sweep
        !! See notes of 13-Feb-2015 for derivation
        if (defocus_half_range .ne. 0.0 .and. spa_freq_sq .ne. 0.0) then
            ! b - the component of the phase shift which is defocus dependent
            b = pi * wl * spa_freq_sq
            eval_ctf_slave = eval_ctf_slave * sin(b*defocus_half_range)/(b*defocus_half_range)
        endif


        if (present(sign_only)) then
            if (sign_only) eval_ctf_slave = sign(1.0,eval_ctf_slave)
        endif
    end function eval_ctf_slave


    !>  \brief returns the argument (radians) to the sine and cosine terms of the ctf
    !!
    !!  We follow the convention, like the rest of the cryo-EM/3DEM field, that underfocusing the objective lens
    !!  gives rise to a positive phase shift of scattered electrons, whereas the spherical aberration gives a
    !!  negative phase shift of scattered electrons
    pure elemental real function eval_ctf_phase_shift(cs,wl,dfmid1,dfmid2,angast,additional_phase_shift,spa_freq_sq,ang) &
                                                    result(phase_shift)
        real,                   intent(in)  ::  cs          !< spherical aberation (pixels)
        real,                   intent(in)  ::  wl          !< electron wavelength (pixels)
        real,                   intent(in)  ::  dfmid1      !< defocus along first axis (pixels)
        real,                   intent(in)  ::  dfmid2      !< defocus along second axis (for astigmatic ctf, dfmid1 .ne. dfmid2) (pixels)
        real,                   intent(in)  ::  angast      !< azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,                   intent(in)  ::  additional_phase_shift  !<  In radians, additional phase shift experienced by scattered electrons, compared to scattered electrons
        real,                   intent(in)  ::  spa_freq_sq !< square of spatial frequency at which to compute the ctf (1/pixels^2)
        real,                   intent(in)  ::  ang         !< angle at which to compute the ctf (radians)
        ! private variables
        real    ::  df      ! defocus at point at which we're evaluating the ctf
        ! start work
        ! compute the defocus
        df = eval_df(ang,dfmid1,dfmid2,angast)
        ! compute the ctf argument
        phase_shift = pi * wl * spa_freq_sq * (df - 0.5 * wl**2 * spa_freq_sq * cs) + additional_phase_shift
    end function eval_ctf_phase_shift

    !>  \brief Returns the phase shift (in radians) experienced by scattered electrons relative to unscattered electrons
    pure subroutine eval_ctf_phase_shift_array(cs,wl,dfmid1,dfmid2,angast,additional_phase_shift,&
                                                spa_freq_sq,ang,phase_shift)
        real,                   intent(in)  ::  cs              !< spherical aberation (pixels)
        real,                   intent(in)  ::  wl              !< electron wavelength (pixels)
        real,                   intent(in)  ::  dfmid1          !< defocus along first axis (pixels)
        real,                   intent(in)  ::  dfmid2          !< defocus along second axis (for astigmatic ctf, dfmid1 .ne. dfmid2) (pixels)
        real,                   intent(in)  ::  angast          !< azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,                   intent(in)  ::  additional_phase_shift  !<  In radians, additional phase shift experienced by scattered electrons, compared to scattered electrons
        real,                   intent(in)  ::  spa_freq_sq(:)  !< square of spatial frequency at which to compute the ctf (1/pixels^2)
        real,                   intent(in)  ::  ang             !< angle at which to compute the ctf (radians)
        real,                   intent(out) ::  phase_shift(:)
        ! private variables
        real    ::  df      ! defocus at point at which we're evaluating the ctf
        !real    ::  ccos
        ! start work
        ! compute the defocus
        df = eval_df(ang,dfmid1,dfmid2,angast)
        ! compute the ctf argument
        phase_shift(:) = pi * wl * spa_freq_sq(:) * (df - 0.5 * wl**2 * spa_freq_sq(:) * cs) + additional_phase_shift
    end subroutine eval_ctf_phase_shift_array

    !>  \brief  Return the effective defocus given the astigmatism parameters and the azimuth of interest
    pure elemental real function eval_df(ang,dfmid1,dfmid2,angast) result (df)
        real,                   intent(in)  ::  ang         !< angle at which to compute the defocus (radians)
        real,                   intent(in)  ::  dfmid1      !< defocus along first axis (pixels)
        real,                   intent(in)  ::  dfmid2      !< defocus along second axis (for astigmatic ctf, dfmid1 .ne. dfmid2) (pixels)
        real,                   intent(in)  ::  angast      !< azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        !
        df = 0.5*(dfmid1+dfmid2+cos(2.0*(ang-angast))*(dfmid1-dfmid2))
    end function eval_df

    !>  \brief  Given the argument to the ctf, return the spatial frequency
    !pure elemental function ctf_sf_from_phase_shift(cs,wl,dfmid1,dfmid2,angast,additional_phase_shift,phase_shift,ang) &
    !                        result(sf)
    !    real,                   intent(in)  ::  cs          !< Spherical aberation (pixels)
    !    real,                   intent(in)  ::  wl          !< Electron wavelength (pixels)
    !    real,                   intent(in)  ::  dfmid1      !< Defocus along first axis (pixels)
    !    real,                   intent(in)  ::  dfmid2      !< Defocus along second axis (for astigmatic ctf, dfmid1 .NE. dfmid2) (pixels)
    !    real,                   intent(in)  ::  angast      !< Azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
    !    real,                   intent(in)  ::  additional_phase_shift  !<  In radians, additional phase shift of scattered electrons relative to non-scattered electrons, perhaps due to a phase plate
    !    real,                   intent(in)  ::  phase_shift !< Phase shift
    !    real,                   intent(in)  ::  ang         !< Angle at which to compute the ctf (radians)
    !    real                                ::  sf          !< Spatial frequency
    !    ! Private variables
    !    !real    ::  df      ! Defocus
    !    ! Start work
    !    !sf = ctf_sq_sf_from_phase_shift(cs,wl,dfmid1,dfmid2,angast,additional_phase_shift,phase_shift,ang)
    !    sf = sqrt(sf)
    !end function ctf_sf_from_phase_shift

    !>  \brief  Given the argument to the ctf (phase shift, in radians), return the squared spatial frequency
    !!
    !!  The solution below was found using Wolfram Alpha, see notes 150211
    subroutine ctf_sq_sf_from_phase_shift(cs,wl,dfmid1,dfmid2,angast,additional_phase_shift,phase_shift,ang,sq_sf,num_solutions)
        real,                   intent(in)  ::  cs              !< Spherical aberation (pixels)
        real,                   intent(in)  ::  wl              !< Electron wavelength (pixels)
        real,                   intent(in)  ::  dfmid1          !< Defocus along first axis (pixels)
        real,                   intent(in)  ::  dfmid2          !< Defocus along second axis (for astigmatic ctf, dfmid1 .NE. dfmid2) (pixels)
        real,                   intent(in)  ::  angast          !< Azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,                   intent(in)  ::  additional_phase_shift  !<  In radians, additional phase shift of scattered electrons relative to non-scattered electrons, perhaps due to a phase plate
        real,                   intent(in)  ::  phase_shift     !< Phase shift (radians)
        real,                   intent(in)  ::  ang             !< Angle at which to compute the ctf (radians)
        real,                   intent(out) ::  sq_sf(2)        !< Squared spatial frequency (pixels^-2)
        integer,                intent(out) ::  num_solutions   !<  There may be 0, 1 or 2 solutions
        ! Private variables
        real    ::  df      ! Defocus
        real    ::  phase_shift_without_phase_plate
        real    ::  a,b,c
        real    ::  det
        ! Start work

        ! Compute the defocus
        df = eval_df(ang,dfmid1,dfmid2,angast)

        phase_shift_without_phase_plate = phase_shift - additional_phase_shift

        !
        ! Solve
        !
        a = pi*wl*df
        b = 0.5*pi*wl**3*cs
        c = phase_shift_without_phase_plate

        det = a**2-4.0*b*c

        if ( det .ge. 0.0) then
            sq_sf(1) = (a+sqrt(det)) / (-2.0*b)
            sq_sf(2) = (a-sqrt(det)) / (-2.0*b)
            if (det .eq. 0.0) then
                num_solutions = 1
            else
                num_solutions = 2
            endif
        else
            ! No solutions
            sq_sf = 0.0
            num_solutions = 0
        endif

        ! Squared spatial frequencies must be >=0.0; let's remove bad solutions
        select case (num_solutions)
            case (1)
                if (sq_sf(1) .lt. 0.0) then
                    num_solutions = 0
                endif
            case (2)
                if (sq_sf(2) .lt. 0.0 .and. sq_sf(1) .ge. 0.0) then
                    num_solutions = 1
                else if (sq_sf(1) .lt. 0.0 .and. sq_sf(2) .ge. 0.0) then
                    num_solutions = 1
                    sq_sf(1) = sq_sf(2)
                else if (sq_sf(1) .lt. 0.0 .and. sq_sf(2) .lt. 0.0) then
                    num_solutions = 0
                endif
        end select

    end subroutine ctf_sq_sf_from_phase_shift

    !>  \brief  Given the argument to the ctf, return the defocus
    pure subroutine ctf_df_from_phase_shift(cs,wl,phase_shift_without_phase_plate,sf,df)
        real,                   intent(in)  ::  cs                              !< Spherical aberation (pixels)
        real,                   intent(in)  ::  wl                              !< Electron wavelength (pixels)
        real,                   intent(in)  ::  phase_shift_without_phase_plate !< Phase shift of scattered electrons due to non-phase-plate optics (radians)
        real,                   intent(in)  ::  sf                              !< Spatial frequency (1/pixels)
        real,                   intent(out) ::  df                              !< Defocus (pixels)
        ! Private variables

        ! Start work

        ! Compute the defocus
        df = (0.5 * cs * pi * sf**4 * wl**3 + phase_shift_without_phase_plate) / (pi * sf**2 * wl)
    end subroutine ctf_df_from_phase_shift

    !>  \brief  Find the ctf argument (phase shift) values at which CTF=ctf_value
    !!
    !!  According to Wolfram Alpha, the solutions to a*cos(t) + b*sin(t) = c are:
    !!  1)  t = 2 (atan((b-sqrt(a^2+b^2-c^2))/(a+c))+pi*n),   with n an integer
    !!  2)  t = 2 (atan((b+sqrt(a^2+b^2-c^2))/(a+c))+pi*n),   with n an integer
    !!
    !!  The above two solutions only hold if:
    !!                  a+c != 0
    !!                 -b*sqrt(a**2+b**2-c**2)+a**2+ac+b**2 != 0   [for solution 1]
    !!                  b*sqrt(a**2+b**2-c**2)+a**2+ac+b**2 != 0   [for solution 2]
    !!
    !!  In our case, a = - ampl_const
    !!           and b = - sqrt(1-ampl_const**2)
    !!               c = ctf_value
    !!  and t is the "argument" to the ctf (i.e. the phase shift), which can be "converted" to a spatial frequency
    subroutine ctf_solve_for_phase_shift(ampl_cont,sols,ctf_value)
        real,                   intent(in)      ::  ampl_cont
        real,   allocatable,    intent(inout)   ::  sols(:)     !<  Solutions
        real,   optional,       intent(in)      ::  ctf_value   !<  Value of the ctf at the solutions
        ! Private variables
        integer ::  i,j
        real    ::  a,b
        real    ::  c
        real    ::  sqrt_arg
        real    ::  sqrt_val
        real    ::  sol_first_term_1,sol_first_term_2
        ! Start work
        if (.not. allocated(sols)) then
            call this_program%TerminateWithFatalError('ctf_solve_for_arg','array of solutions is not allocated')
        endif
        if (size(sols) .lt. 1) then
            write(*,'(a,i0)') '**error(ctf_solve_for_arg): array of solutions has unexpected size ', size(sols)
            call this_program%TerminateWithFatalError('ctf_solve_for_arg','array of solutions has unexpected size')
        endif
        c = 0.0e0
        if (present(ctf_value)) c = ctf_value
        !
        a = -ampl_cont
        b = -sqrt(1-ampl_cont**2)
        ! Note that since 0.0 <= ampl_cont <= 1.0:
        ! a**2+b**2 = 1.0
        !sqrt_arg = a**2+b**2-c**2
        sqrt_arg = max(1.0-c**2,0.0)
        sqrt_val = sqrt(sqrt_arg)
        sol_first_term_1 = atan((b-sqrt_val)/(a+c))
        sol_first_term_2 = atan((b+sqrt_val)/(a+c))

        ! Loop over the zeroes
        do i=1,size(sols)/2
            j = (i-1)*2+1
            sols(j) = 2.0 * (sol_first_term_1 + pi*real(i-1-size(sols)/4))
            j = (i-1)*2+2
            sols(j) = 2.0 * (sol_first_term_2 + pi*real(i  -size(sols)/4))
        enddo
    end subroutine ctf_solve_for_phase_shift

    !>  \brief  convert acceleration voltage in kv into electron wavelength in angstroms
    pure real function akv_to_wl(akv)
        real, intent(in) :: akv !<  kV
        akv_to_wl = 12.26/sqrt(1000.0*akv+0.9784*(1000.0*akv)**2/(10.0**6.0))
    end function akv_to_wl


end module
