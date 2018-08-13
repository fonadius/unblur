!
! Copyright 2014 Howard Hughes Medical Institute
! All rights reserved
! Use is subject to Janelia Farm Research Campus Software Copyright 1.1
! license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
!
!
module ImageCTFComparisons
    use Images
    use ContrastTransferFunctions
    type ImageCTFComparison
        type(Image),    allocatable     ::  spectra(:)
        type(ContrastTransferFunction)  ::  ctf
        real                            ::  pixel_size  !<  Angstroms
        logical                         ::  find_phase_shift = .false.
        logical                         ::  fit_defocus_sweep = .false.
    end type

    logical,    parameter               ::  fit_squared_ctf = .false. ! whether to fit squared ctf or absolute ctf; Amplitude spectrum has the advantage that it does not accentuate strong artefacts (e.g. central cross)
    logical,    parameter               ::  apply_defocus_uncertainty_envelope = .true. !<  whether to take into account the change in defocus during frame acquisiton when computing CTFs to fit to the spectra

    contains

    !>  \brief  Scoring function to be minimised
    function ctffind_objective_function(array_of_values,scoring_parameters) bind(c)
        use iso_c_binding
        use Units
        use ContrastTransferFunctions
        use StringManipulations, only : IntegerToString
        ! arguments
        type(c_ptr), value ::  array_of_values
        type(c_ptr), value ::  scoring_parameters
        ! result
        real(kind=c_double) ::  ctffind_objective_function
        ! private variables
        real(kind=8),               pointer ::  values_to_try(:)    !<  1. Average defocus (first frame; Angstroms); 2. Astigmatism (defocus1-defocus2; Angstroms) (first frame); 3. Astigmatism azimuth (degrees); 4. Additional phase shift; 5(4 if no phase shift). Average defocus (final frame, relative to first frame)
        type(ImageCTFComparison),   pointer ::  self
        type(ContrastTransferFunction)      ::  my_ctf
        logical,                parameter   ::  debug = .false.
        integer                             ::  number_of_dimensions
        integer                             ::  current_spectrum
        real                                ::  additional_defocus
        real                                ::  start_defocus1, start_defocus2
        real                                ::  astigmatism_radians
        real                                ::  current_correlation_terms(3), correlation_terms(3)
        real                                ::  astig_penalty
        integer                             ::  final_defocus_dimension_number
        real                                ::  defocus_change_per_frame
        ! start work

        ! Get at the comparison object
        call c_f_pointer(scoring_parameters,self)

        ! Work out the number of search dimensions
        number_of_dimensions = 3
        if (self%find_phase_shift)  number_of_dimensions = number_of_dimensions + 1
        if (self%fit_defocus_sweep) number_of_dimensions = number_of_dimensions + 1

        ! Get to the array of values
        call c_f_pointer(array_of_values,values_to_try,[number_of_dimensions])

        ! Convert defocus values from Angstroms to pixels
        start_defocus1 = convert(real(values_to_try(1)+values_to_try(2)),angstroms,pixels,self%pixel_size)
        start_defocus2 = convert(real(values_to_try(1)-values_to_try(2)),angstroms,pixels,self%pixel_size)
        astigmatism_radians = convert(real(values_to_try(3)),degrees,radians)

        ! Set defocus values & astigmatism
        my_ctf = self%ctf
        call my_ctf%SetDefocus( start_defocus1,start_defocus2, astigmatism_radians )

        ! Set phase shift
        if (self%find_phase_shift) then
            call my_ctf%SetAdditionalPhaseShift(real(values_to_try(4)))
        endif

        ! Set defocus half-range
        if (.not. apply_defocus_uncertainty_envelope) then
            call my_ctf%SetDefocusHalfRange(0.0e0)
        endif

        ! Compute score
        correlation_terms = 0.0d0
        if (self%fit_defocus_sweep) then
            ! Which dimension is the additional defocus parameter?
            if (self%find_phase_shift) then
                final_defocus_dimension_number = 5
            else
                final_defocus_dimension_number = 4
            endif
            !
            defocus_change_per_frame = convert( &
                                    real((values_to_try(final_defocus_dimension_number)-values_to_try(1)) &
                                    /(size(self%spectra)-1),kind=4), &! Linear sweep is assumed
                                    angstroms,pixels,self%pixel_size)
            if (apply_defocus_uncertainty_envelope) call my_ctf%SetDefocusHalfRange(0.5*defocus_change_per_frame)
            !
            do current_spectrum=1,size(self%spectra)
                additional_defocus = (current_spectrum-1)*defocus_change_per_frame
                call my_ctf%SetDefocus( start_defocus1 + additional_defocus, &
                                        start_defocus2 + additional_defocus, &
                                        astigmatism_radians )
                call self%spectra(current_spectrum)%CTFOperation(my_ctf,correlation_terms=current_correlation_terms,&
                                                                 astigmatism_penalty=astig_penalty)
                correlation_terms = correlation_terms + current_correlation_terms
            enddo
            ! Finalise the score
            if (correlation_terms(2) .lt. 0.0 .or. correlation_terms(3) .lt. 0.0) then
                ctffind_objective_function = 0.0
            else
                ctffind_objective_function = - correlation_terms(1) / sqrt(correlation_terms(2)*correlation_terms(3))
            endif
            ctffind_objective_function = ctffind_objective_function + astig_penalty
        else
            ctffind_objective_function = - self%spectra(1)%GetCorrelationWithCTF(my_ctf,fit_squared_ctf)
        endif

        if (debug) then
            !$omp critical (ctffind_objective_function_debug_print)
            write(*,'(a,'//IntegerToString(size(values_to_try))//'(f0.2,1x),f0.3)') &
                                            '**debug(ctffind_objective_function): values, score: ', &
                                            values_to_try,ctffind_objective_function
            !$omp end critical (ctffind_objective_function_debug_print)
        endif
    end function ctffind_objective_function

end module ImageCTFComparisons


!>  \brief  A container for ctffind routines
module ctffind_routines

    implicit none

    contains

    !>  \brief  Given the filename of the output summary from CTFFind, populate an array of CTF objects
    subroutine ParseCTFFindResultsFile( results_filename, pixel_size, &
                                        ctf )
        use ContrastTransferFunctions
        use NumericTextFiles
        ! arguments
        character(len=*),                               intent(in)      ::  results_filename
        real,                                           intent(in)      ::  pixel_size          !<  The CTF objects will be returned with defocus values in pixels, but are read in from the file in Angstroms
        type(ContrastTransferFunction), allocatable,    intent(inout)   ::  ctf(:)
        ! private variables
        type(NumericTextFile)                       ::  numeric_results_file
        integer                                     ::  number_of_images
        integer                                     ::  image_counter
        real                                        ::  acceleration_voltage
        real                                        ::  spherical_aberration
        real                                        ::  amplitude_contrast
        real                                        ::  current_line_results(7)
        ! start work

        if (allocated(ctf)) deallocate(ctf)
        call numeric_results_file%Init(results_filename,OPEN_TO_READ)

        ! How many images?
        number_of_images = numeric_results_file%GetNumberOfDataLines()
        allocate(ctf(number_of_images))

        ! Get illumination parameters
        call GetParameterValueFromCTFFindResultsFile(results_filename,'voltage:',acceleration_voltage)
        call GetParameterValueFromCTFFindResultsFile(results_filename,'aberration:',spherical_aberration)
        call GetParameterValueFromCTFFindResultsFile(results_filename,'contrast:',amplitude_contrast)



        ! Loop over images to
        do image_counter=1,number_of_images
            call numeric_results_file%ReadNextDataLine(current_line_results)
            call ctf(image_counter)%Init(acceleration_voltage,spherical_aberration,amplitude_contrast, &
                                         current_line_results(2)/10000.0,current_line_results(3)/10000.0,&
                                         current_line_results(4), &
                                         0.0,1.0,100.0,pixel_size,current_line_results(5))
        enddo

    end subroutine ParseCTFFindResultsFile

    subroutine GetParameterValueFromCTFFindResultsFile(results_filename,keyword,value_found)
        use StringManipulations
        use UsefulFunctions, only : grep
        !
        character(len=*),               intent(in)      ::  results_filename
        character(len=*),               intent(in)      ::  keyword
        real,                           intent(out)     ::  value_found
        !
        character(len=line_max_len)                     ::  current_line
        character(len=line_max_len),    allocatable     ::  words(:)
        integer                                         ::  word_counter
        integer                                         ::  keyword_position
        character(len=:),               allocatable     ::  word
        !
        current_line = grep(results_filename,keyword)
        call Split(current_line,words)
        do word_counter = 1, size(words)
            if (StringsAreEqual(words(word_counter),keyword)) keyword_position = word_counter
        enddo
        word = RemoveNonNumericCharacters(words(keyword_position+1))
        read(word,*) value_found
    end subroutine GetParameterValueFromCTFFindResultsFile



    subroutine PrepareDiagnosticImage(  spectrum,ctf, &
                                        rotational_average,spatial_frequency, &
                                        rotational_average_astig,rotational_average_astig_fit, &
                                        frc_of_fit,frc_of_fit_sigma)
        use Images
        use ImageCTFComparisons, only : fit_squared_ctf
        use ContrastTransferFunctions
        ! Arguments
        type(Image),                    intent(inout)   ::  spectrum            !<  Amplitude spectrum to modify with an overlay
        type(ContrastTransferFunction), intent(inout)   ::  ctf                 !<  CTF which was fit to the spectrum
        real(kind=8),   intent(inout)   ::  rotational_average(:)
        real(kind=8),   intent(inout)   ::  spatial_frequency(:)
        real(kind=8),   intent(inout)   ::  rotational_average_astig(:)
        real(kind=8),   intent(inout)   ::  rotational_average_astig_fit(:)
        real(kind=8),   intent(inout)   ::  frc_of_fit(:)
        real(kind=8),   intent(inout)   ::  frc_of_fit_sigma(:)
        ! Private variables
        real                ::  min_rad_for_spec_stats
        real                ::  max_rad_for_spec_stats
        real                ::  average
        real                ::  sigma
        logical,        parameter   ::  debug   =   .false.
        logical             ::  spectrum_is_blank
        logical,        parameter   ::  skip_rescaling = .false.
        ! Start work

        spectrum_is_blank = spectrum%IsConstant()

        call spectrum%AddConstant((-1.0)*spectrum%GetAverageOfValuesOnEdges())

        ! Work out radii over which we will compute statistics for prettyfying the spectrum
        min_rad_for_spec_stats = ctf%ComputeFrequencyOfAZero(0.0,2)
        max_rad_for_spec_stats = min(max(ctf%ComputeFrequencyOfAZero(0.0,3),&
                                        ctf%GetHighestFrequencyForFitting()),0.5)
        if (min_rad_for_spec_stats .gt. 0.4) then
            min_rad_for_spec_stats = ctf%GetLowestFrequencyForFitting()
        endif
        ! Prettify the spectrum
        if (.not. spectrum_is_blank) then
            call spectrum%ComputeAverageAndSigmaOfValuesInSpectrum(                                         &
                                                                min_rad_for_spec_stats          &
                                                                *spectrum%GetLogicalDimension(1),           &
                                                                max_rad_for_spec_stats         &
                                                                *spectrum%GetLogicalDimension(1),           &
                                                                average,sigma)
            call spectrum%ApplyCircularMask(5.0,inverse=.true.)
            call spectrum%SetMaximumValueOnCentralCross(average)
            if (sigma .gt. 0.0) then
                call spectrum%SetMinimumAndMaximumValue(average-4.0*sigma,average+4.0*sigma)
            endif
            call spectrum%ComputeAverageAndSigmaOfValuesInSpectrum(                                         &
                                                                min_rad_for_spec_stats          &
                                                                *spectrum%GetLogicalDimension(1),           &
                                                                max_rad_for_spec_stats         &
                                                                *spectrum%GetLogicalDimension(1),           &
                                                                average,sigma)
            call spectrum%AddConstant(-1.0*average)
            if (sigma .gt. 0.0) then
                call spectrum%MultiplyByConstant(1.0e0/sigma)
            endif
            call spectrum%AddConstant(average)
        endif
        call spectrum%Compute1DRotationalAverage(rotational_average)
        if (debug) call spectrum%WriteToDisk('dbg_average_spectrum_before_rescaling.mrc')
        call spectrum%ComputeRotationalAverageOfPowerSpectrum(ctf,spatial_frequency,            &
                                                                      rotational_average_astig,                 &
                                                                      rotational_average_astig_fit,             &
                                                                      frc_of_fit,frc_of_fit_sigma,              &
                                                                      rescale_input=.not. skip_rescaling,       &
                                                                      squared_ctf_was_fit=fit_squared_ctf)


        call spectrum%ComputeAverageAndSigmaOfValuesInSpectrum(                                         &
                                                            min_rad_for_spec_stats          &
                                                            *spectrum%GetLogicalDimension(1),           &
                                                            max_rad_for_spec_stats         &
                                                            *spectrum%GetLogicalDimension(1),           &
                                                            average,sigma)

        call spectrum%SetMinimumAndMaximumValue(average-1.0*sigma,average+2.0*sigma)
        if (debug) call spectrum%WriteToDisk('dbg_average_spectrum_before_overlay.mrc')
        call spectrum%OverlayCTF(ctf,squared_ctf=fit_squared_ctf)


    end subroutine PrepareDiagnosticImage

    subroutine PrepareAmplitudeSpectrumForFitting(spectrum,pixel_size,minimum_resolution,ctf)
        use Images
        use ContrastTransferFunctions
        use UsefulFunctions, only : IsEven
        ! Arguments
        type(Image),                    intent(inout)   ::  spectrum            !<  Amplitude spectrum to be filtered/background-subtracted
        real,                           intent(in)      ::  pixel_size          !<  In Angstroms, pixel size of micrograph from which amplitude spectrum was computed
        real,                           intent(in)      ::  minimum_resolution  !<  In reciprocal Angstroms, lowest resolution to be used for CTF fitting
        type(ContrastTransferFunction), intent(in)      ::  ctf                 !<  CTF object set up to correspond to the lowest defocus which will searched for
        ! Private variables
        logical,    parameter   ::  debug = .false.
        real                    ::  distance_between_zeroes
        integer                 ::  convolution_box_size
        type(Image)             ::  convoluted_spectrum
        real                    ::  average
        real                    ::  sigma
        real                    ::  threshold
        ! Start work

        ! Try to weaken artefacts by doing some mild thresholding
        if (debug) call spectrum%WriteToDisk('dbg_average_spectrum_before_cross_removal.mrc')
        call spectrum%ComputeAverageAndSigmaOfValuesInSpectrum(                                         &
                                spectrum%GetLogicalDimension(1)*pixel_size/minimum_resolution, &
                                spectrum%GetLogicalDimension(1)*1.0,           &
                                average,sigma, &
                                cross_half_width=12.0)
        call spectrum%SetMinimumAndMaximumValue(average-1.0*sigma,average+10.0*sigma)


        ! Compute low-pass filtered version of spectrum
        distance_between_zeroes = ctf%ComputeFrequencyOfAZero(0.0,3)-ctf%ComputeFrequencyOfAZero(0.0,2)
        convolution_box_size = max( &
                                ! the ctffind4 way of doing things, which sometimes is too small
                                distance_between_zeroes * spectrum%GetLogicalDimension(1) * 3.0, &
                                ! the ctffind3 way of doing things
                                2.0 * pixel_size/minimum_resolution* &
                                spectrum%GetLogicalDimension(1) * sqrt(2.0)&
                                )
        if (IsEven(convolution_box_size)) convolution_box_size = convolution_box_size + 1
        if (debug) call spectrum%WriteToDisk('dbg_average_spectrum_before_conv.mrc')
        call spectrum%SpectrumBoxConvolution(   convolution_box_size,convoluted_spectrum, &
                                                minimum_radius=spectrum%logical_dimensions(1) &
                                                                *pixel_size/minimum_resolution)

        if (debug) call convoluted_spectrum%WriteToDisk('dbg_background.mrc')

        ! Subtract low-pass filtered power spectrum from power spectrum. This should remove background slope.
        spectrum%real_values(1:spectrum%GetLogicalDimension(1),:,:) &
        =   spectrum%real_values(1:spectrum%GetLogicalDimension(1),:,:)**2 &
        -   convoluted_spectrum%real_values(1:spectrum%GetLogicalDimension(1),:,:)**2

        if (debug) call spectrum%WriteToDisk('dbg_average_spectrum_background_subtracted.mrc')


        ! Threshold high values
        threshold = spectrum%GetMaximumValue(minimum_distance_from_center=3,minimum_distance_from_edge=3)
        call spectrum%SetMaximumValue(threshold)

        ! Debug dump
        if (debug) call spectrum%WriteToDisk('dbg_av_spec.mrc')

    end subroutine PrepareAmplitudeSpectrumForFitting

end module ctffind_routines
