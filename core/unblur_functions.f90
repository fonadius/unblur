module unblur_functions

    implicit none

    contains

    subroutine unblur_refine_alignment( stack_of_images, &
                                        maximum_number_of_iterations,unitless_bfactor, &
                                        inner_radius_for_peak_search,outer_radius_for_peak_search, &
                                        maximum_shift_convergence_threshold, &
                                        additional_shifts_x,additional_shifts_y)
        use Images
        use Peaks
        use ProgressBars
        use Curves
        use UsefulFunctions, only : IsOdd
        ! arguments
        type(Image),                    intent(inout)   ::  stack_of_images(:)                  !<  Images to align to each other. Should already have approximate shifts applied if they're known. The images are returned with additional shifts applied.
        integer,                        intent(in)      ::  maximum_number_of_iterations        !<  For the alignment procedure
        real,                           intent(in)      ::  unitless_bfactor                    !<  B factor (unitless) to apply during alignment
        logical,                        intent(in)      ::  mask_central_cross                  !<  Whether to mask the central cross in Fourier space during alignment
        integer,                        intent(in)      ::  width_of_vertical_line              !<  Width of vertical line of the cross
        integer,                        intent(in)      ::  width_of_horizontal_line            !<  Width of horizontal line of the central cross
        logical,                        intent(in)      ::  apply_circular_mask_to_reference    !<  During alignment iterations
        real,                           intent(in)      ::  circular_mask_radius                !<  In pixels
        real,                           intent(in)      ::  inner_radius_for_peak_search        !<  Pixels. Set to 0.0 to search for peaks even at the origin of the cross-correlation function.
        real,                           intent(in)      ::  outer_radius_for_peak_search        !<  pixels
        real,                           intent(in)      ::  maximum_shift_convergence_threshold !<  Pixels. Alignment is deemed to have converged when the maximum image shift in either X or Y is at or below this
        real,           allocatable,    intent(inout)   ::  additional_shifts_x(:)              !<  Output (i.e. wiped on input). Shifts in X which had to be applied to each frame to bring it into register with the others
        real,           allocatable,    intent(inout)   ::  additional_shifts_y(:)              !<  Output (i.e. wiped on input). Shifts in Y which had to be applied to each frame to bring it into register with the others
        real,                           intent(in)      ::  pixel_size                          !<  In Angstroms. This is only used for print out information.
        logical,        optional,       intent(in)      ::  indicate_progress                   !<  Print out progress bars etc. False by default.
        integer,        optional,       intent(in)      ::  running_average_width               !<  How many frames should be averaged together? Defaults to 1. Must be odd.
        ! private variables
        type(Image)                 ::  sum_of_images
        type(Image)                 ::  sum_of_images_minus_current
        type(Image)                 ::  current_frame_image
        integer                     ::  current_iteration
        integer                     ::  image_counter
        integer                     ::  inner_image_counter
        type(Peak)                  ::  my_peak
        !real                        ::  maximum_x_shift, maximum_y_shift
        real                        ::  maximum_shift
        real,       allocatable     ::  old_shifts_x(:), old_shifts_y(:), current_shifts_x(:), current_shifts_y(:)
        type(Curve)                 ::  x_shifts
        type(Curve)                 ::  y_shifts
        integer                     ::  number_of_middle_image
        integer                     ::  order_of_polynomial
        logical                     ::  iindicate_progress
        type(ProgressBar)           ::  my_progress_bar
        integer                     ::  number_of_images_processed
        integer                     ::  rrunning_average_width
        integer                     ::  first_frame_to_align, last_frame_to_align
        real                        ::  acceleration_x, acceleration_y
        ! Options for smoothing of shifts
        integer,    parameter           ::  no_smoothing = 0
        integer,    parameter           ::  polynomial = 1
        integer,    parameter           ::  spline = 2
        integer,    parameter           ::  smoothing_type = spline
        ! Debug
        logical,    parameter           ::  debug = .false.
        ! start work

        iindicate_progress = .false.
        rrunning_average_width = 1
        first_frame_to_align = 1 + (rrunning_average_width-1)/2
        last_frame_to_align  = size(stack_of_images) - (rrunning_average_width-1)/2

        ! Prepare arrays to hold results
        additional_shifts_x = 0.0e0
        additional_shifts_y = 0.0e0
        old_shifts_x = additional_shifts_x
        old_shifts_y = additional_shifts_y

        ! Which is the middle image?
        if (IsOdd(size(stack_of_images))) then
            number_of_middle_image = (size(stack_of_images) + 1) / 2
        else
            number_of_middle_image = size(stack_of_images) / 2 + 1
        endif

        ! Prepare the sum of images
        call sum_of_images%Allocate(mould=stack_of_images(1),in_real_space=.false.)
        sum_of_images = (0.0e0,0.0e0)
        do image_counter = 1, size(stack_of_images)
            call sum_of_images%AddImage(stack_of_images(image_counter))
        enddo

        ! Prepare the smoothing curve
        call x_shifts%Init(size(stack_of_images))
        call y_shifts%Init(size(stack_of_images))

        ! Perform the main alignment loop, loop until the maximum shift is less than 0.1 pixels
        do current_iteration = 1, maximum_number_of_iterations

            !$omp parallel default(shared) private(image_counter, sum_of_images_minus_current, my_peak,inner_image_counter,current_frame_image)
            call sum_of_images_minus_current%Reset()
            call current_frame_image%Reset()

            !$omp do
            do image_counter=first_frame_to_align,last_frame_to_align

                ! Prepare the current image to be aligned
                current_frame_image = stack_of_images(image_counter)

                ! Prepare the reference image (sum of all images, minus current image)
                sum_of_images_minus_current = sum_of_images
                call sum_of_images_minus_current%SubtractImage(current_frame_image)
                call sum_of_images_minus_current%ApplyBFactor(unitless_bfactor)

                ! Compute the cross-correlation function, and find a peak in there
                call sum_of_images_minus_current%CalculateCrossCorrelationImageWith(current_frame_image)
                my_peak = sum_of_images_minus_current%FindPeakWithParabolaFit(min_radius=inner_radius_for_peak_search, &
                                                                              max_radius=outer_radius_for_peak_search)
                ! Update shifts for this frame
                additional_shifts_x(image_counter) = additional_shifts_x(image_counter) + my_peak%CoOrdinates(1)
                additional_shifts_y(image_counter) = additional_shifts_y(image_counter) + my_peak%CoOrdinates(2)

            enddo ! end of loop over images

            !$omp enddo

            !$omp end parallel

            ! Apply [spline|polynomial|no] smoothing to the shifts
            call x_shifts%ClearData()
            call y_shifts%ClearData()
            do image_counter=1, size(stack_of_images)
                call x_shifts%AddPoint(real(image_counter), additional_shifts_x(image_counter))
                call y_shifts%AddPoint(real(image_counter), additional_shifts_y(image_counter))
            enddo

            select case (smoothing_type)
                case (spline)
                    call x_shifts%FitSplineToData()
                    call x_shifts%CopySplineModel(additional_shifts_x)
                    call y_shifts%FitSplineToData()
                    call y_shifts%CopySplineModel(additional_shifts_y)
                case (no_smoothing)
                    call x_shifts%CopyYData(additional_shifts_x)
                    call y_shifts%CopyYData(additional_shifts_y)
                case (polynomial)
                    order_of_polynomial = min(6,size(stack_of_images)-1)
                    call x_shifts%FitPolynomialToData(order_of_polynomial)
                    call x_shifts%CopyPolynomialModel(additional_shifts_x)
                    call y_shifts%FitPolynomialToData(order_of_polynomial)
                    call y_shifts%CopyPolynomialModel(additional_shifts_y)
            end select


            ! Subtract off the shift from the middle image, so that the alignment is kept such that that image has
            ! a zero shift, we also want to work out the maximum shift based on before we do this so as to find the true answer
            !maximum_x_shift = maxval(abs(additional_shifts_x - old_shifts_x))
            !maximum_y_shift = maxval(abs(additional_shifts_y - old_shifts_y))

            !if (maximum_x_shift .gt. maximum_y_shift) then
            !    maximum_shift = maximum_x_shift
            !else
            !    maximum_shift = maximum_y_shift
            !endif

            ! Subtract middle images shift
            additional_shifts_x = additional_shifts_x - additional_shifts_x(number_of_middle_image)
            additional_shifts_y = additional_shifts_y - additional_shifts_y(number_of_middle_image)

            ! Work out the actual shift to apply for this iteration
            current_shifts_x = additional_shifts_x - old_shifts_x
            current_shifts_y = additional_shifts_y - old_shifts_y

            ! What is the largest absolute shift we'll be applying for this iteration?
            maximum_shift = max(maxval(abs(current_shifts_x)),maxval(abs(current_shifts_y)))

            ! Update old_shifts (the old_shifts for the middle frame will be 0.0)
            old_shifts_x = additional_shifts_x
            old_shifts_y = additional_shifts_y


            ! Shift the images
            !$omp parallel default(shared) private(image_counter)
            !$omp do
            do image_counter=1, size(stack_of_images)

                call stack_of_images(image_counter)%PhaseShift( current_shifts_x(image_counter), &
                                                                current_shifts_y(image_counter), &
                                                                0e0)
            enddo
            !$omp enddo
            !$omp end parallel

            ! If either of the convergence criterion are met then break
            if (     current_iteration .eq. maximum_number_of_iterations &
                .or. maximum_shift     .le. maximum_shift_convergence_threshold) exit

            ! Need to make the sum image for the next iteration
            sum_of_images = (0e0, 0e0)
            do image_counter=1, size(stack_of_images)
                call sum_of_images%AddImage(stack_of_images(image_counter))
            enddo

        enddo ! end of loop over alignment iterations

    end subroutine unblur_refine_alignment

end module unblur_functions
