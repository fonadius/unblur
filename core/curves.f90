!>  \brief  Class to deal with Curves
module Curves

    use Globals

    implicit none

    integer,        parameter               ::  linear_spline  = 1
    integer,        parameter               ::  cubic_spline   = 2
    integer,        parameter               ::  quintic_spline = 3
    integer,        parameter               ::  heptic_spline  = 4

    type,         public  ::  Curve

        integer                             ::  number_of_points = 0            !<  How many points we currently have for this object
        integer                             ::  allocated_space_for_points = 0  !<  How much space we have allocated for points
        real,   allocatable                 ::  data_x(:)                       !<  X data of the curve - the methods should autotrack that there is enough memory allocated
        real,   allocatable                 ::  data_y(:)                       !<  Y data of the curve - the methods should autotrack that there is enough memory allocated

        ! Arrays to hold fit data
        real,   allocatable                 ::  polynomial_values_at_data_x(:)  !<  Holds the y values for the polynomial fit of the input_data, it is only used if the polynomial routines are called
        real,   allocatable                 ::  spline_values_at_data_x(:)      !<  Y values for the spline fit

        ! Polynomial specifics..
        integer                             ::  polynomial_order                !<  If a polynomial is to be / has been fit, what order was it
        real, allocatable                   ::  polynomial_coefficients(:)      !<  If a polynomial has been fit, what are it's coefficients

        ! For convenience and pretty printing
        character(len=32)                   ::  title_x
        character(len=32)                   ::  title_y

        logical                             ::  have_been_initialised = .false.

   contains
        ! Memory, initialisation etc
        procedure                           ::  CurveInit
        procedure, public                   ::  Init   => CurveInit                         !<  Initialise the curve object.
        procedure                           ::  Destroy
        procedure                           ::  CopyFrom
        generic, public                     ::  assignment(=) => CopyFrom

        ! Data in/out

        procedure, public                   ::  ClearData                       !<  Reset the curve object
        procedure, public                   ::  AddPoint                        !<  Add a point to the curve data
        procedure, public                   ::  AddPoints                       !<  Add an array of points to the curve data
        procedure, public                   ::  PrintCurve                      !<  Write Out the Curve info - useful mainly for debugging.
        procedure, public                   ::  PrintPolynomial                 !<  Write Out the Polynomial fit to the data - useful mainly for debugging.


        procedure, public                   ::  CopyPolynomialModel
        procedure, public                   ::  CopySplineModel
        procedure, public                   ::  CopyXData
        procedure, public                   ::  CopyYData

        procedure, public                   ::  WriteToDisk
        procedure, public                   ::  ReadFromDisk

        ! Get / set
        procedure, public                   ::  GetNumberOfDataPoints
        procedure, public                   ::  GetX
        procedure, public                   ::  GetY

        ! Fit methods
        procedure                           ::  FitPolynomialToData
        procedure                           ::  FitSplineToData

    end type

    contains

    !>  \brief  Initialise the object
    subroutine CurveInit(self, number_of_points_to_allocate,title_x,title_y)

        ! Arguments
        class(Curve),                       intent(inout)   ::  self
        integer,                optional,   intent(in)      ::  number_of_points_to_allocate
        character(len=*),       optional,   intent(in)      ::  title_x, title_y    !<  Titles for the X and Y data series

        ! Variables
        integer ::  nnumber_of_points_to_allocate

        ! Begin

        ! We want to pre-initialise and blank for a set of points, this can be specified as an argument, if not we'll do 256

        if (present(number_of_points_to_allocate)) then
            nnumber_of_points_to_allocate = number_of_points_to_allocate
        else
            nnumber_of_points_to_allocate = 256
        endif

        ! check data isn't already allocated..

        if (allocated(self%data_x)) deallocate(self%data_x)
        if (allocated(self%data_y)) deallocate(self%data_y)

        ! allocate

        allocate(self%data_x(nnumber_of_points_to_allocate))
        allocate(self%data_y(nnumber_of_points_to_allocate))

        self%allocated_space_for_points = nnumber_of_points_to_allocate

        ! There should be no points currently..

        self%number_of_points = 0

        ! Set the titles
        self%title_x = 'X'
        self%title_y = 'Y'
        if (present(title_x)) self%title_x = title_x
        if (present(title_y)) self%title_y = title_y

        ! We've been initialised
        self%have_been_initialised = .true.

    end subroutine CurveInit

    !>  \brief  Clear up memory
    subroutine Destroy(self)
        class(Curve),   intent(inout)   ::  self
        if (allocated(self%data_x)) deallocate(self%data_x)
        if (allocated(self%data_y)) deallocate(self%data_y)
        if (allocated(self%polynomial_values_at_data_x)) deallocate(self%polynomial_values_at_data_x)
        if (allocated(self%spline_values_at_data_x)) deallocate(self%spline_values_at_data_x)
        if (allocated(self%polynomial_coefficients)) deallocate(self%polynomial_coefficients)
        self%have_been_initialised = .false.
    end subroutine Destroy

    !>  \brief  Assignment
    subroutine CopyFrom(self,other)
        class(Curve),   intent(inout)   ::  self        !<  LHS of assignment
        class(Curve),   intent(in)      ::  other       !<  RHS of assignment
        !

        !
        call self%Destroy()
        call self%Init(other%number_of_points)
        self%title_x                      = other%title_x
        self%title_y                      = other%title_y
        self%number_of_points             = other%number_of_points
        self%allocated_space_for_points   = other%allocated_space_for_points
        if (allocated(other%data_x)) then
        self%data_x                       = other%data_x
        endif
        if (allocated(other%data_y)) then
        self%data_y                       = other%data_y
        endif
        if (allocated(other%polynomial_values_at_data_x)) then
        self%polynomial_values_at_data_x  = other%polynomial_values_at_data_x
        endif
        if (allocated(other%spline_values_at_data_x)) then
        self%spline_values_at_data_x      = other%spline_values_at_data_x
        endif
        self%polynomial_order             = other%polynomial_order
        if (allocated(other%polynomial_coefficients)) then
        self%polynomial_coefficients      = other%polynomial_coefficients
        endif
        self%have_been_initialised        = other%have_been_initialised
    end subroutine CopyFrom

    pure integer function GetNumberOfDataPoints(self)
        class(Curve),   intent(in)  ::  self
        GetNumberOfDataPoints = self%number_of_points
    end function GetNumberOfDataPoints

    pure real function GetX(self,point_number)
        class(Curve),   intent(in)  ::  self
        integer,        intent(in)  ::  point_number
        GetX = self%data_x(point_number)
    end function GetX

    pure real function GetY(self,point_number)
        class(Curve),   intent(in)  ::  self
        integer,        intent(in)  ::  point_number
        GetY = self%data_y(point_number)
    end function GetY

    subroutine PrintCurve(self)
        ! Arguments
        class(Curve),           intent(inout)   ::  self
        ! Variables
        integer     ::  point_counter

        ! Begin
        if (.not. self%have_been_initialised) call self%Init()

        ! Print header row
        write(*,'(/a32,1x,a32)') self%title_x, self%title_y
        ! Print data rows
        do point_counter = 1, self%number_of_points
            write(*, '(f32.4, 1x, f32.4)') self%data_x(point_counter), self%data_y(point_counter)
        enddo
        write(*,'(a)') ' '

    end subroutine PrintCurve

    !>  \brief  Write out curve data to a file on disk
    subroutine WriteToDisk(self,filename)
        use NumericTextFiles
        ! Arguments
        class(Curve),           intent(in)      ::  self
        character(len=*),       intent(in)      ::  filename
        ! Private variables
        real                        ::  current_values(2)
        type(NumericTextFile)       ::  my_text_file
        character(len=line_max_len) ::  header_line
        integer                     ::  i
        ! Start work
        call my_text_file%Init(filename,OPEN_TO_WRITE,2)
        write(header_line,'(a1,1x,2(a8,1x))') '!', self%title_x, self%title_y
        call my_text_file%WriteCommentLine(header_line)
        do i=1,self%number_of_points
            current_values(1) = self%data_x(i)
            current_values(2) = self%data_y(i)
            call my_text_file%WriteDataLine(current_values)
        enddo
    end subroutine WriteToDisk


    !>  \brief  Read curve data from a file on disk (current data are overwritten)
    subroutine ReadFromDisk(self,filename)
        use NumericTextFiles
        ! Arguments
        class(Curve),           intent(inout)   ::  self
        character(len=*),       intent(in)      ::  filename
        ! Private variables
        real                        ::  current_values(2)
        type(NumericTextFile)       ::  my_text_file
        integer                     ::  i
        ! Start work
        call my_text_file%Init(filename,OPEN_TO_READ,2)
        call self%ClearData()
        do i=1,my_text_file%number_of_data_lines
            call my_text_file%ReadNextDataLine(current_values)
            call self%AddPoint(current_values(1),current_values(2))
        enddo
    end subroutine ReadFromDisk


    subroutine CopyPolynomialModel(self, array_to_copy_to)

        !Arguments
        class(Curve),           intent(inout)   ::  self
        real, allocatable,      intent(inout)   ::  array_to_copy_to(:)


        if (.not. self%have_been_initialised) call self%Init()

        if (allocated(array_to_copy_to)) then

            if (size(array_to_copy_to) .lt. self%number_of_points) then
                deallocate(array_to_copy_to)
                allocate(array_to_copy_to(self%number_of_points))
            endif
        else
            allocate(array_to_copy_to(self%number_of_points))
        endif


        array_to_copy_to(1:self%number_of_points) = self%polynomial_values_at_data_x(1:self%number_of_points)


    end subroutine CopyPolynomialModel

    subroutine CopySplineModel(self, array_to_copy_to)

        !Arguments
        class(Curve),           intent(inout)   ::  self
        real, allocatable,      intent(inout)   ::  array_to_copy_to(:)


        if (.not. self%have_been_initialised) call self%Init()

        if (allocated(array_to_copy_to)) then

            if (size(array_to_copy_to) .lt. self%number_of_points) then
                deallocate(array_to_copy_to)
                allocate(array_to_copy_to(self%number_of_points))
            endif
        else
            allocate(array_to_copy_to(self%number_of_points))
        endif


        array_to_copy_to(1:self%number_of_points) = self%spline_values_at_data_x(1:self%number_of_points)


    end subroutine CopySplineModel


    subroutine CopyXData(self, array_to_copy_to)

        !Arguments
        class(Curve),           intent(inout)   ::  self
        real, allocatable,      intent(inout)   ::  array_to_copy_to(:)


        if (.not. self%have_been_initialised) call self%Init()

        if (allocated(array_to_copy_to)) then

            if (size(array_to_copy_to) .lt. self%number_of_points) then
                deallocate(array_to_copy_to)
                allocate(array_to_copy_to(self%number_of_points))
            endif
        else
            allocate(array_to_copy_to(self%number_of_points))
        endif


        array_to_copy_to(1:self%number_of_points) = self%data_x(1:self%number_of_points)

    end subroutine CopyXData


    subroutine CopyYData(self, array_to_copy_to)

        !Arguments
        class(Curve),           intent(inout)   ::  self
        real, allocatable,      intent(inout)   ::  array_to_copy_to(:)


        if (.not. self%have_been_initialised) call self%Init()

        if (allocated(array_to_copy_to)) then

            if (size(array_to_copy_to) .lt. self%number_of_points) then
                deallocate(array_to_copy_to)
                allocate(array_to_copy_to(self%number_of_points))
            endif
        else
            allocate(array_to_copy_to(self%number_of_points))
        endif


        array_to_copy_to(1:self%number_of_points) = self%data_y(1:self%number_of_points)

    end subroutine CopyYData


    subroutine PrintPolynomial(self)

        !Arguments
        class(Curve),           intent(inout)   ::  self

        !Variables

        integer                                 ::  point_counter

        !Begin

        if (.not. self%have_been_initialised) call self%Init()

        if (allocated(self%polynomial_values_at_data_x)) then

            do point_counter = 1, self%number_of_points
                write(*, '(a, f8.2, a, f8.2)')  'X = ', self%data_x(point_counter), &
                                                ' Y = ', self%polynomial_values_at_data_x(point_counter)
            enddo
        endif

    end subroutine PrintPolynomial

    !>  \brief Basically set the number of points to 0, so newly added values will be the only values..
    subroutine ClearData(self)
        !Arguments
        class(Curve),           intent(inout)   ::  self

        !Begin

        self%number_of_points = 0

    end subroutine ClearData

    !>  \brief Add a new single point to the data
    subroutine AddPoint(self, x_value, y_value)
        !Arguments
        class(Curve),           intent(inout)   ::  self
        real,                   intent(in)      ::  x_value
        real,                   intent(in)      ::  y_value

        !Variables
        integer                                ::  new_allocation_size
        real, allocatable                       ::  temporary_array(:)

        ! check initialisation
        if (.not. self%have_been_initialised) call self%Init()

        ! check initialisation

        if (.not. self%have_been_initialised) call self%Init()

        !Begin - check we have enough space in memory

        if (self%allocated_space_for_points .le. self%number_of_points) then

            ! We need to allocate more space.. it doesn't make sense to only allocate 4 more bytes, so allocate double the current allocation (up to a limit of 1MB)

            if (self%allocated_space_for_points .le. 131072) then
                new_allocation_size = self%allocated_space_for_points * 2
            else
                new_allocation_size = self%allocated_space_for_points + 262144
            endif

            ! Resize X

            allocate(temporary_array(new_allocation_size))
            temporary_array(1:self%allocated_space_for_points) = self%data_x
            call move_alloc(temporary_array, self%data_x)


            ! Resize Y

            allocate(temporary_array(new_allocation_size))
            temporary_array(1:self%allocated_space_for_points) = self%data_y
            call move_alloc(temporary_array, self%data_y)

            self%allocated_space_for_points = new_allocation_size

        endif

        ! The allocation should be fine now - copy the data

        self%number_of_points = self%number_of_points + 1
        self%data_x(self%number_of_points) = x_value
        self%data_y(self%number_of_points) = y_value


    end subroutine AddPoint

    !>  \brief Add a new array of points to the data
    subroutine AddPoints(self, x_values, y_values)
        !Arguments
        class(Curve),           intent(inout)   ::  self
        real,      intent(in)                   ::  x_values(:)
        real,      intent(in)                   ::  y_values(:)

        !Variables
        integer                                 ::  new_allocation_size
        real, allocatable                       ::  temporary_array(:)
        integer                                 ::  size_of_input_arrays


        ! check initialisation
        if (.not. self%have_been_initialised) call self%Init()

        !Begin - work out how many values we are adding, and whether both arrays have the same number

        !if (.not. allocated(x_values) .or. .not. allocated(y_values)) then
        !    call this_program%TerminateWithFatalError('Curve::AddPoints', 'An input array is not allocated')
        !endif

        size_of_input_arrays = size(x_values)

        if (size_of_input_arrays .ne. size(y_values)) then
            call this_program%TerminateWithFatalError('Curve::AddPoints', 'The input arrays are different sizes')
        endif

        ! Check the memory allocation..

        if (self%allocated_space_for_points .le. self%number_of_points + size_of_input_arrays) then

            ! We need to allocate more space.. it doesn't make sense to only allocate 4 more bytes, so allocate double the current allocation (up to a limit of 1MB)
            ! or if that isn't enough, then add enough for these points.

            if (self%allocated_space_for_points .le. 131072) then
                new_allocation_size = self%allocated_space_for_points * 2
            else
                new_allocation_size = self%allocated_space_for_points + 262144
            endif

            if (new_allocation_size .lt. self%allocated_space_for_points + size_of_input_arrays) then
                new_allocation_size = self%allocated_space_for_points + size_of_input_arrays
            endif

            ! Resize X

            allocate(temporary_array(new_allocation_size))
            temporary_array(1:self%allocated_space_for_points) = self%data_x
            call move_alloc(temporary_array, self%data_x)

            ! Resize Y

            allocate(temporary_array(new_allocation_size))
            temporary_array(1:self%allocated_space_for_points) = self%data_y
            call move_alloc(temporary_array, self%data_y)

            self%allocated_space_for_points = new_allocation_size

        endif

        ! The allocation should be fine now - copy the data

        self%data_x(self%number_of_points + 1:self%number_of_points + size_of_input_arrays) = x_values
        self%data_y(self%number_of_points + 1:self%number_of_points + size_of_input_arrays) = y_values
        self%number_of_points = self%number_of_points + size_of_input_arrays

    end subroutine AddPoints

    subroutine FitSplineToData(self,order_of_spline,interpolate,print_result_stats)
        use GCV_SPLINES
        ! Arguments
        class(Curve),               intent(inout)   ::  self
        integer,        optional,   intent(in)      ::  order_of_spline     !<  Cubic by default. Can be linear_spline, cubic_spline, quintic_spline or heptic_spline.
        logical,        optional,   intent(in)      ::  interpolate         !<  False by default. Force interpolation: the curve must go exactly through every data point.
        logical,        optional,   intent(in)      ::  print_result_stats  !<  Print out result statistics
        ! Private variables
        integer,        parameter   ::  number_of_curves_to_fit = 1
        real(kind=8)                ::  data_x(self%number_of_points), data_y(self%number_of_points)
        real(kind=8)                ::  weights_x(self%number_of_points)
        real(kind=8)                ::  weights_y(number_of_curves_to_fit)
        integer                     ::  oorder_of_spline
        integer                     ::  optimization_mode
        real(kind=8)                ::  optimization_mode_value
        logical                     ::  iinterpolate
        real(kind=8)                ::  spline_coefficients(self%number_of_points*number_of_curves_to_fit)
        !integer                     ::  number_of_coefficients
        real(kind=8)                ::  fitting_work_array(6*(self%number_of_points*heptic_spline+1)+heptic_spline)
        real(kind=8),   allocatable ::  evaluation_work_array(:)
        integer                     ::  error_number
        integer                     ::  current_point
        integer                     ::  current_segment
        logical                     ::  pprint_result_stats
        ! Start work

        ! Check that the x values are strictly increasing
        do current_point=2,size(self%data_x)
            if (self%data_x(current_point) .le. self%data_x(current_point-1)) then
                write(*,'(a)') 'Error: x values must be strictly increasing'
                print *, self%data_x
                call this_program%TerminateWithFatalError('Curves::FitSplineToData','Bad input data')
            endif
        enddo


        ! Check the value of order of spline makes sense
        if (present(order_of_spline)) then
            if (order_of_spline .lt. linear_spline .or. order_of_spline .gt. heptic_spline) then
                write(*,'(a,i0)') 'Error: bad value for order of spline: ', order_of_spline
                call this_program%TerminateWithFatalError('Curves::FitSplineToData','Bad parameter')
            endif
        endif


        ! We need the data in 64bit floats
        data_x = self%data_x
        data_y = self%data_y

        ! Weights
        weights_x = 1.0d0
        !weights_x = (/ (real(size(weights_x))/real(current_point),current_point=1,size(weights_x)) /)
        weights_y = 1.0d0

        ! Order of the spline
        oorder_of_spline = cubic_spline
        if (present(order_of_spline)) oorder_of_spline = order_of_spline

        ! Do we want to interpolate?
        iinterpolate = .false.
        if (present(interpolate)) iinterpolate = interpolate

        ! Print out some stats when done?
        pprint_result_stats = .false.
        if (present(print_result_stats)) pprint_result_stats = print_result_stats

        ! Optimization type
        if (iinterpolate) then
            optimization_mode = 1 ! known p
            optimization_mode_value = 0.0d0
        else
            optimization_mode = 2 ! generalized cross validation (GCV)
            optimization_mode_value = 0.1d0 ! not used with GCV, I think
        endif

        ! Coefficients computed by gcvspl
        !allocate(spline_coefficients(self%number_of_points,number_of_curves_to_fit))

        ! Work array - will also contain results
        !allocate(fitting_work_array(6*(self%number_of_points*oorder_of_spline+1)+oorder_of_spline))

        !call GCV_EXAMPLE()
        !call this_program%Terminate()

        ! Generalised Cross Validation SPLine fitting routine from netlib

        call gcvspl(data_x, &
                    data_y, &
                    self%number_of_points, &
                    weights_x, &
                    weights_y, &
                    oorder_of_spline, &
                    self%number_of_points, &
                    number_of_curves_to_fit, &
                    optimization_mode, &
                    optimization_mode_value, &
                    spline_coefficients, &
                    self%number_of_points, &
                    fitting_work_array, &
                    error_number &
                    )


        if (error_number .ne. 0) then
            write(*,'(a,i0,a)') 'Error ', error_number, ' when performing GCV spline fit'
            print *, 'data_x = ', data_x
            print *, 'data_y = ', data_y
            print *, self%number_of_points
            call this_program%TerminateWithFatalError('Curves::FitSplineToData','Error when fitting')
        endif

        if (pprint_result_stats) then
            write(*,'(a,f0.3)') 'Generalized cross validation value: ', fitting_work_array(1)
            write(*,'(a,f0.3)') 'Mean squared residual: ', fitting_work_array(2)
            write(*,'(a,f0.3)') 'Estimate of # of deg of freedom of residual sum of squares: ', fitting_work_array(3)
            write(*,'(a,f0.3)') 'Smoothing parameter: ', fitting_work_array(4)
            write(*,'(a,f0.3)') 'Estimate of true mean squared error: ', fitting_work_array(5)
            write(*,'(a,f0.3)') 'Gauss-Markov error variance: ', fitting_work_array(6)
        endif

        ! Let's now evaluate the spline that was fit at the points
        if (allocated(self%spline_values_at_data_x)) deallocate(self%spline_values_at_data_x)
        allocate(self%spline_values_at_data_x(self%number_of_points * 2))
        allocate(evaluation_work_array(2*oorder_of_spline * 2))

        do current_point=1,self%number_of_points

            current_segment = current_point
            self%spline_values_at_data_x(current_point) &
                                         = real( splder(0,oorder_of_spline,self%number_of_points, &
                                                        data_x(current_point), &
                                                        data_x, &
                                                        spline_coefficients, &
                                                        current_segment, &
                                                        evaluation_work_array) &
                                                )

        enddo





    end subroutine FitSplineToData


    subroutine FitPolynomialToData(self, order_of_polynomial)
        !Arguments
        class(Curve),           intent(inout)   ::  self
        integer, optional,       intent(in)     ::  order_of_polynomial

        !variables

        integer                                 ::  used_order_of_polynomial
        !integer                                 ::  point_counter
        integer                                 ::  dummy
        real                                    ::  standard_deviation

        !Begin

        ! check initialisation
        if (.not. self%have_been_initialised) call self%Init()

        if (present(order_of_polynomial)) then
            used_order_of_polynomial = order_of_polynomial
        else
            used_order_of_polynomial = 6
        end if

        if (self%number_of_points .le. 2) then
            call this_program%TerminateWithFatalError('Curve::FitPolynomialtoData','Cannot fit a polynomial to less than 3 points.')
        endif

        if (self%number_of_points - 1 .lt. used_order_of_polynomial) then
            write(*,*)
            write(*,'(a)') '**Warning - Curve:FitPolynomialtoData'
            write(*,'(a)') '            Requested order of polynomial is greater than number of points - 2.'
            write(*,'(a,i0)') '            Setting polynomial order : ', self%number_of_points - 2
            write(*, *)

            used_order_of_polynomial = self%number_of_points - 2
        endif

        ! need to make sure the polynomial areas are allocated appropriately.

        if (allocated(self%polynomial_values_at_data_x)) then

            if (size(self%polynomial_values_at_data_x) .ne. self%number_of_points) then
                deallocate(self%polynomial_values_at_data_x)
                allocate(self%polynomial_values_at_data_x(self%number_of_points))
            endif
        else
            allocate(self%polynomial_values_at_data_x(self%number_of_points))
        endif

        if (allocated(self%polynomial_coefficients)) then

            if (size(self%polynomial_coefficients) .ne. used_order_of_polynomial) then
                deallocate(self%polynomial_coefficients)
                allocate(self%polynomial_coefficients(used_order_of_polynomial))
            endif
        else
            allocate(self%polynomial_coefficients(used_order_of_polynomial))
        endif

        ! Just call the LS_POLY routine

        call LS_POLY(used_order_of_polynomial, 0.0e0,self%number_of_points, dummy, &
                    self%data_x, self%data_y, self%polynomial_coefficients, standard_deviation, &
                    self%polynomial_values_at_data_x)
    end subroutine FitPolynomialToData



!*****************************************************************
!*         LEAST SQUARES POLYNOMIAL FITTING PROCEDURE            *
!* ------------------------------------------------------------- *
!* This program least squares fits a polynomial to input data.   *
!* forsythe orthogonal polynomials are used in the fitting.      *
!* The number of data points is n.                               *
!* The data is input to the subroutine in x[i], y[i] pairs.      *
!* The coefficients are returned in c[i],                        *
!* the smoothed data is returned in v[i],                        *
!* the order of the fit is specified by m.                       *
!* The standard deviation of the fit is returned in d.           *
!* There are two options available by use of the parameter e:    *
!*  1. if e = 0, the fit is to order m,                          *
!*  2. if e > 0, the order of fit increases towards m, but will  *
!*     stop if the relative standard deviation does not decrease *
!*     by more than e between successive fits.                   *
!* The order of the fit then obtained is l.                      *
!*****************************************************************
!* Reference: BASIC Scientific Subroutines, Vol. II *
!* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
!*                                                  *
!*                F90 version by J-P Moreau, Paris  *
!*                        (www.jpmoreau.fr)         *
Subroutine LS_POLY(m,e1,n,l,x,y,c_out,dd,v)
  implicit none

  !Labels: 10,15,20,30,50

  !Arguments

  integer,            intent(inout) :: m
  integer,            intent(in)    :: n
  integer,            intent(inout) :: l
  real,               intent(in)    :: e1
  real,               intent(inout) :: dd
  real,               intent(in)    :: x(:)
  real,               intent(in)    :: y(:)
  real,               intent(inout) :: c_out(:)
  real,               intent(inout) :: v(:)

  ! Variables

  real, allocatable                 :: a(:)
  real, allocatable                 :: b(:)
  real, allocatable                 :: c(:)
  real, allocatable                 :: c2(:)
  real, allocatable                 :: d(:)
  real, allocatable                 :: e(:)
  real, allocatable                 :: f(:)

  integer i,l2,n1
  real a1,a2,b1,b2,c1,d1,f1,f2,v1,v2,vv,w


  ! Allocate

  allocate(a(m+1))
  allocate(b(m+1))
  allocate(c(0:m+1))
  allocate(c2(m+1))
  allocate(d(n))
  allocate(e(n))
  allocate(f(m+1))

  n1 = m + 1; l=0
  v1 = 1.e7
  ! Initialize the arrays
  do i=1, n1
    a(i) = 0.e0; b(i) = 0.e0; f(i) = 0.e0
  end do
  do i=1, n
    v(i) = 0.e0; d(i) = 0.e0
  end do
  d1 = sqrt(float(n)); w = d1;
  do i=1, n
    e(i) = 1.e0 / w
  end do
  f1 = d1; a1 = 0.e0
  do i=1, n
    a1 = a1 + x(i) * e(i) * e(i)
  end do
  c1 = 0.e0
  do i=1, n
    c1 = c1 + y(i) * e(i)
  end do
  b(1) = 1.e0 / f1; f(1) = b(1) * c1
  do i=1, n
    v(i) = v(i) + e(i) * c1
  end do
  m = 1
! Save latest results
10 do i=1, l
    c2(i) = c(i)
  end do
  l2 = l; v2 = v1; f2 = f1; a2 = a1; f1 = 0.e0
  do i=1, n
    b1 = e(i)
    e(i) = (x(i) - a2) * e(i) - f2 * d(i)
    d(i) = b1
    f1 = f1 + e(i) * e(i)
  end do
  f1 = sqrt(f1)
  do i=1, n
    e(i) = e(i) / f1
  end do
  a1 = 0.e0
  do i=1, n
    a1 = a1 + x(i) * e(i) * e(i)
  end do
  c1 = 0.e0
  do i=1, n
    c1 = c1 + e(i) * y(i)
  end do
  m = m + 1; i = 0
15 l = m - i; b2 = b(l); d1 = 0.e0
  if (l > 1)  d1 = b(l - 1)
  d1 = d1 - a2 * b(l) - f2 * a(l)
  b(l) = d1 / f1; a(l) = b2; i = i + 1
  if (i.ne.m) goto 15
  do i=1, n
    v(i) = v(i) + e(i) * c1
  end do
  do i=1, n1
    f(i) = f(i) + b(i) * c1
    c(i) = f(i)
  end do
  vv = 0.e0
  do i=1, n
    vv = vv + (v(i) - y(i)) * (v(i) - y(i))
  end do
  !Note the division is by the number of degrees of freedom
  vv = sqrt(vv / float(n - l - 1)); l = m
  if (e1.eq.0.e0) goto 20
  !Test for minimal improvement
  if (abs(v1 - vv) / vv < e1) goto 50
  !if error is larger, quit
  if (e1 * vv > e1 * v1) goto 50
  v1 = vv
20 if (m.eq.n1) goto 30
  goto 10
!Shift the c[i] down, so c(0) is the constant term
30 do i=1, l
    c(i - 1) = c(i)
  end do
  c(l) = 0.e0
  ! l is the order of the polynomial fitted
  l = l - 1; dd = vv
  return
! Aborted sequence, recover last values
50 l = l2; vv = v2
  do i=1, l
    c(i) = c2(i)
    c_out = c(i)
  end do
  goto 30
end subroutine




end module
