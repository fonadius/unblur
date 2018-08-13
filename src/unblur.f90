program unblur

    use Globals
    use UserSuppliedParameters
    use UserInputs
    use ImageFiles
    use Images
    use Peaks
    use Curves
    use ProgressBars
    use NumericTextFiles
    use UsefulFunctions
    use StringManipulations
    use ElectronDoses
    use ContrastTransferFunctions
    use unblur_functions, only : unblur_refine_alignment
    !use ctffind_routines, only : ParseCTFFindResultsFile

    implicit none

    ! Variables associated with user input

    type(UserInput)                 ::  my_user_input

    type(UserSuppliedInteger)       ::  number_of_frames_per_movie

    type(UserSuppliedInteger)       ::  maximum_number_of_iterations


    type(UserSuppliedInteger)       ::  width_of_vertical_line
    type(UserSuppliedInteger)       ::  width_of_horizontal_line

    type(UserSuppliedReal)          ::  pixel_size
    type(UserSuppliedReal)          ::  inner_radius
    type(UserSuppliedReal)          ::  outer_radius
    type(UserSuppliedReal)          ::  bfactor
    type(UserSuppliedReal)          ::  threshold_maximum_shift

    type(UserSuppliedFilename)      ::  input_filename
    type(UserSuppliedFilename)      ::  output_filename
    type(UserSuppliedFilename)      ::  aligned_frames_filename
    type(UserSuppliedFilename)      ::  shifts_output_filename
    type(UserSuppliedFilename)      ::  frc_output_filename

    type(UserSuppliedLogical)       ::  apply_drift_filter
    type(UserSuppliedLogical)       ::  mask_central_cross
    type(UserSuppliedLogical)       ::  set_expert_options
    type(UserSuppliedLogical)       ::  verbose_output

    type(UserSuppliedLogical)       ::  apply_dose_filter
    type(UserSuppliedLogical)       ::  divide_dose_filter
    type(UserSuppliedLogical)       ::  save_aligned_frames
    type(UserSuppliedReal)          ::  acceleration_voltage
    type(UserSuppliedReal)          ::  exposure_per_frame
    type(UserSuppliedReal)          ::  pre_exposure_amount

    !type(UserSuppliedLogical)       ::  apply_ctf
    !type(UserSuppliedLogical)       ::  phase_flip_only
    !type(UserSuppliedFilename)      ::  ctffind_results

    ! Variable associated with Images

    type(ImageFile)                 ::  input_imagefile

    type(Image),    allocatable     ::  image_stack(:)
    type(Image),    allocatable     ::  unbinned_image_stack(:)
    type(Image)                     ::  sum_image
    type(Image)                     ::  sum_image_2


    real                            ::  maximum_shift
    real                            ::  maximum_x_shift
    real                            ::  maximum_y_shift
    real                            ::  current_drift_x
    real                            ::  current_drift_y

    real                            ::  final_score
    real                            ::  final_score_numerator
    real                            ::  final_score_denominator

    real                            ::  average_x_shift
    real                            ::  average_y_shift

    integer                         ::  number_of_input_images
    integer                         ::  number_of_movies

    integer                         ::  image_counter
    integer                         ::  bin_counter
    integer                         ::  movie_counter

    integer                         ::  iterations_performed

    integer                         ::  order_of_polynomial
    integer                         ::  number_of_middle_image


    type(Curve)                     ::  x_shifts
    type(Curve)                     ::  y_shifts
    type(Curve)                     ::  fsc_curve

    type(ProgressBar)               ::  my_progress_bar

    type(NumericTextFile)           ::  my_shifts_file
    type(NumericTextFile)           ::  my_frc_file

    real, allocatable               ::  current_x_shifts(:)
    real, allocatable               ::  current_y_shifts(:)
    real, allocatable               ::  old_x_shifts(:)
    real, allocatable               ::  old_y_shifts(:)
    real, allocatable               ::  total_x_shifts(:)
    real, allocatable               ::  total_y_shifts(:)

    real, allocatable               ::  frc_x_data(:)
    real, allocatable               ::  frc_y_data(:)

    !real                            ::  temp_real(2) ! for writing shifts to file

    ! Option for pre-binning
    logical,        parameter       ::  allow_pre_binning = .true.
    integer                         ::  pre_binning_factor
    integer                         ::  binned_dimensions(3)
    real                            ::  binned_pixel_size
    real                            ::  inner_radius_after_binning
    real                            ::  outer_radius_after_binning

    ! Dose filtering
    type(ElectronDose)              ::  my_electron_dose
    !type(Image),    allocatable     ::  dose_filters_squared(:)
    !type(Image)                     ::  dose_filters_squared_sum

    integer                         ::  i,j
    real                            ::  x,y, x_sq, y_sq
    real                            ::  azimuth
    real                            ::  current_denominator
    real                            ::  current_critical_dose
    real                            ::  current_optimal_dose
    real                            ::  current_dose
    real                            ::  current_ctf
    real,           parameter       ::  minimum_absolute_ctf_value = 0.001

    type(ContrastTransferFunction), allocatable ::  frame_ctf(:)



    ! Initialise this_program with the program name and version number
    call this_program%Init('UnBlur', '1.0.2','2015')


    ! Do the user input
    call my_user_input%Init(this_program%program_name)

    input_filename   = my_user_input%GetFilenameFromUser('Input stack filename', &
                                   'The input file, containing your raw movie frames', 'input_filename', 'my_movie.mrc', .true.)

    ! Open the input file
    call input_imagefile%Init(input_filename%value)
    number_of_input_images = input_imagefile%GetStackSize()

    number_of_frames_per_movie      =   my_user_input%GetIntegerFromUser('Number of frames per movie', &
                                        'How many frames per micrograph?', &
                                        'number_of_frames_per_movie',IntegerToString(number_of_input_images), &
                                        min_value=1,max_value=number_of_input_images)

    output_filename                 =   my_user_input%GetFilenameFromUser('Output aligned sum file', &
                                        'The output file, containing a weighted sum of the aligned input frames', &
                                        'output_filename', 'my_aligned_sum.mrc', .false.)

    shifts_output_filename          =   my_user_input%GetFilenameFromUser('Output shifts file', &
                                        'This file will contain the computed X/Y shifts for each frame.', &
                                        'shifts_filename', 'my_shifts.txt', .false.)

    pixel_size                      =   my_user_input%GetRealFromUser('Pixel size of images (A)', &
                                        'Pixel size of input images in Angstroms', 'pixel_size', '1', min_value = 0.0e0)

    apply_dose_filter               =   my_user_input%GetLogicalFromUser('Apply Dose filter?', &
                                        'Apply a dose-dependent filter to frames before summing them', &
                                        'apply_dose_filter','NO')

    if (apply_dose_filter%value) then
         exposure_per_frame         =   my_user_input%GetRealFromUser('Exposure per frame (e/A^2)', &
                                        'Exposure per frame, in electrons per square Angstrom','exposure_per_frame', &
                                        '1.0',min_value=0.00001e0)
         acceleration_voltage       =   my_user_input%GetRealFromUser('Acceleration voltage (kV)', &
                                        'Acceleration voltage during imaging','acceleration_voltage','300.0', &
                                        min_value=200.0e0,max_value=300.0e0)

         pre_exposure_amount         =   my_user_input%GetRealFromUser('Pre-exposure amount(e/A^2)', &
                                        'Amount of pre-exposure prior to the first frame, in electrons per square Angstrom','pre_exposure_amount', &
                                        '0.0',min_value=0.0000)

    endif

    save_aligned_frames             =   my_user_input%GetLogicalFromUser('Save Aligned Frames?', &
                                        'If yes, the aligned frames will be saved to the supplied filename', &
                                        'save_aligned_frames','NO')

    if (save_aligned_frames%value) then

    aligned_frames_filename          =   my_user_input%GetFilenameFromUser('Aligned frames output filename', &
                                        'The output file, containing the aligned frames', &
                                        'aligned_frames_filename', 'my_aligned_frames.mrc', .false.)


    endif



    set_expert_options              =   my_user_input%GetLogicalFromUser('Set Expert Options?', &
                                        'Set these for more control, hopefully not needed', 'set_expert_options', 'NO')



    if (set_expert_options%value) then

       ! apply_ctf                       =   my_user_input%GetLogicalFromUser('Apply CTF to frames?', &
       !                                     'Apply CTF to frames before aligning them','apply_ctf','no')

       ! if (apply_ctf%value) then
       !     phase_flip_only             =   my_user_input%GetLogicalFromUser('Phase flip only?', &
       !                                     'Multiply by the sign of the CTF rather than its value','phase_flip_only','no')

       !     ctffind_results             =   my_user_input%GetFilenameFromUser('Results from CTFfind', &
       !                                     'Filename of CTFfind output','ctffind_results','my_ctf_results.txt', &
       !                                     file_must_exist=.true.)
       ! endif

        frc_output_filename         =   my_user_input%GetFilenameFromUser('Output FRC file', &
                                        'This file will contain the computed FSC of the two half sums', &
                                        'frc_filename', 'my_frc.txt', .false.)


        inner_radius                =   my_user_input%GetRealFromUser('Minimum shift for initial search (Angstroms)', &
                                        'Initial search will be limited to between the inner and outer radii.  ', &
                                        'inner_radius', '2.0', min_value = 0.0e0)

        outer_radius                =   my_user_input%GetRealFromUser('Outer radius shift limit (Angstroms)', &
                                        'The maximum shift of each alignment step will be limited to this value. ', &
                                        'outer_radius', '200.0', min_value = inner_radius%value)

        bfactor                     =   my_user_input%GetRealFromUser('B-factor to apply to images (A^2)', &
                                        'Initial search will be limited to between the inner and outer radii', &
                                        'bfactor', '1500',min_value=0.0)

        !apply_drift_filter          =   my_user_input%GetLogicalFromUser('Apply Drift Filter?', &
        !                                'If yes, a corrective drift filter will be  applied', 'use_drift_filter', 'NO')
        apply_drift_filter%value = .false.


        !mask_central_cross          =   my_user_input%GetLogicalFromUser('Mask central cross in Fourier space?', &
        !                                'This can help reduce problems caused by line artifacts from the detector', &
        !                                'mask_with_cross', 'YES')
        mask_central_cross%value    =   .true.
        mask_central_cross%set_by_user = .false.

        if (mask_central_cross%value) then

            width_of_vertical_line  =   my_user_input%GetIntegerFromUser('Half-width of central vertical line of Fourier mask', &
                                        'The vertical line mask will be twice this size. The central cross mask helps'//&
                                        ' reduce problems by line artefacts from the detector', 'vertical_width', '1', 1)

            width_of_horizontal_line=   my_user_input%GetIntegerFromUser('Half-width of central horizontal line of Fourier mask', &
                                        'The horizontal line mask will be twice this size. The central cross mask helps'//&
                                        ' reduce problems by line artefacts from the detector', 'horizontal_width', '1', 1)
        endif

        threshold_maximum_shift     =   my_user_input%GetRealFromUser('Termination shift threshold', &
                                        'Alignment will iterate until the maximum shift is below this value', &
                                        'threshold_shift', '0.1')

        maximum_number_of_iterations=   my_user_input%GetIntegerFromUser('Maximum number of iterations', &
                                        'Alignment will stop at this number, even if the threshold shift is not reached', &
                                        'max_iterations', '10', min_value = 1)

        if (apply_dose_filter%value) then

            divide_dose_filter      =   my_user_input%GetLogicalFromUser('Restore Noise Power?', &
                                        'divide by dose filter squared', &
                                        'RESTORE_NOISE_POWER','YES')
        endif

        verbose_output              =   my_user_input%GetLogicalFromUser('Verbose Output?', &
                                        'All the output!', 'verbose_output', 'NO')
    else

        inner_radius%value = pixel_size%value * 2e0
        outer_radius%value = 80e0
        bfactor%value = 1500
        apply_drift_filter%value = .false.
        mask_central_cross%value = .true.
        width_of_vertical_line%value = 1e0
        width_of_horizontal_line%value = 1e0
        threshold_maximum_shift%value = pixel_size%value * 0.5
        maximum_number_of_iterations%value = 20
        verbose_output%value = .false.
        divide_dose_filter%value = .true.
        !apply_ctf%value = .false.
    endif

    call my_user_input%UpdateDefaults()
    write(*,*)

    ! Open the output file for the shifts
    call my_shifts_file%Init(shifts_output_filename%value, OPEN_TO_WRITE, number_of_frames_per_movie%value)

    ! Do some basic checks..
    number_of_movies = number_of_input_images / number_of_frames_per_movie%value

    ! Initialise the curves
    call x_shifts%Init(number_of_frames_per_movie%value)
    call y_shifts%Init(number_of_frames_per_movie%value)

    ! Loop over movies
    do movie_counter=1,number_of_movies
        ! If we're going to bin, let's find a good binning factor
        pre_binning_factor = 1
        ! Bin the images
        binned_pixel_size = pixel_size%value

        ! Perform the initial alignment...
        inner_radius_after_binning = inner_radius%value / binned_pixel_size
        outer_radius_after_binning = outer_radius%value / binned_pixel_size

        if (inner_radius%value .ge. 1.0e0) then
            inner_radius_after_binning = max(inner_radius_after_binning,1.0e0)
        endif

        call unblur_refine_alignment(   image_stack, &
                                        1, & ! one iteration only
                                        bfactor%value / binned_pixel_size**2, &
                                        inner_radius_after_binning, &
                                        outer_radius_after_binning, &
                                        max(threshold_maximum_shift%value,&
                                            inner_radius_after_binning+0.01), &
                                        current_x_shifts,current_y_shifts)

        total_x_shifts = current_x_shifts
        total_y_shifts = current_y_shifts

        call unblur_refine_alignment(   image_stack, &
                                        maximum_number_of_iterations%value, &
                                        bfactor%value / binned_pixel_size**2, &
                                        0.0e0, & ! no inner radius limit for peak search
                                        outer_radius_after_binning, &
                                        1., &
                                        current_x_shifts,current_y_shifts)

        total_x_shifts = total_x_shifts + current_x_shifts
        total_y_shifts = total_y_shifts + current_y_shifts

        ! create the final sum and save the aligned images if requird.
        call sum_image%Allocate(mould=image_stack(1),in_real_space=.false.)
        sum_image = (0e0, 0e0)


        do image_counter=1, number_of_frames_per_movie%value
            call sum_image%AddImage(image_stack(image_counter))

            if (save_aligned_frames%value) then
                ! Write out the aligned frames
                !call image_stack(image_counter)%WriteToImageFile(output_aligned_stack_file, (movie_counter-1)*number_of_frames_per_movie%value+image_counter)
                call image_stack(image_counter)%WriteToDisk(aligned_frames_filename%value, (movie_counter-1)*number_of_frames_per_movie%value+image_counter, delete_if_already_exists=.true.)
            endif
        enddo
        ! Write out the sum

        call sum_image%WriteToDisk(output_filename%value, movie_counter, delete_if_already_exists=.true.)

        ! Write out shifts file
        if (movie_counter .eq. 1) then
            call my_shifts_file%WriteCommentLine('Unblur shifts file for input stack : '//input_filename%value)
            call my_shifts_file%WriteCommentLine('Shifts below are given in Angstroms')
            call my_shifts_file%WriteCommentLine('Number of micrographs: '//IntegerToString(number_of_movies))
            call my_shifts_file%WriteCommentLine('Number of frames per movie: '//&
                                                  IntegerToString(number_of_frames_per_movie%value))
            call my_shifts_file%WriteCommentLine('Pixel size (A): '//RealToString(pixel_size%value,4))
            call my_shifts_file%WriteCommentLine('2 lines per micrograph. 1: X-Shift (A); 2: Y-Shift (A)')
            call my_shifts_file%WriteCommentLine('-------------------------')
        endif

        call my_shifts_file%WriteCommentLine('Micrograph '//IntegerToString(movie_counter)//&
                                             ' of '//IntegerToString(number_of_movies))
        call my_shifts_file%WriteDataLine(total_x_shifts * pixel_size%value)
        call my_shifts_file%WriteDataLine(total_y_shifts * pixel_size%value)
    enddo ! end of loop over movies

    ! Deallocation

    if (allocated(image_stack)) deallocate(image_stack)
    if (allocated(current_x_shifts)) deallocate(current_x_shifts)
    if (allocated(current_y_shifts)) deallocate(current_y_shifts)
    if (allocated(total_x_shifts)) deallocate(total_x_shifts)
    if (allocated(total_y_shifts)) deallocate(total_y_shifts)

    call input_imagefile%Close()

end program unblur
