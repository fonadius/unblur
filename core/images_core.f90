!
! Copyright 2014 Howard Hughes Medical Institute
! All rights reserved
! Use is subject to Janelia Farm Research Campus Software Copyright 1.1
! license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
!
!

!>  \brief  Image class
!!
!!  Volumes are images with 3rd logical dimension > 1.
!!
!!  Some definitions:
!!  - logical address: set of 3 integer index values which define a voxel, such that the origin has logical adress (0,0,0). Positive and negative.
!!  - physical address: set of 3 integer index values which define the position within the Fortran 3D array of a voxel (i,j,k). Allowed values for the nth dimension are from 1 to self%dim(n).
!!  - coordinates (shortened to coos): set of 3 real values which define the position of a point within the volume, relative to the origin. In units of pixels or 1/pixels
!!
!!
!!  By convention, we never set the logical_dimensions directly. Instead, we always use the Allocate method, which ensures that relevant properties related to looping & addressing are updated.
!!
!!  Below is a description of the layout of the information in memory for a 4x4 image. The layout is such because of the way FFTW achieves in-place Fourier transforms.
!!
!!
!!  Real space physical address
!!  ===========================
!!
!!  ( 1, 4) ( 2, 4) ( 3, 4) ( 4, 4) ( 5, 4) ( 6, 4)
!!
!!  ( 1, 3) ( 2, 3) ( 3, 3) ( 4, 3) ( 5, 3) ( 6, 3)
!!
!!  ( 1, 2) ( 2, 2) ( 3, 2) ( 4, 2) ( 5, 2) ( 6, 2)
!!
!!  ( 1, 1) ( 2, 1) ( 3, 1) ( 4, 1) ( 5, 1) ( 6, 1)
!!
!!
!!
!!  Real space order of pixel storage in memory
!!  ===========================================
!!
!!  (   19) (   20) (   21) (   22) (   23) (   24)
!!
!!  (   13) (   14) (   15) (   16) (   17) (   18)
!!
!!  (    7) (    8) (    9) (   10) (   11) (   12)
!!
!!  (    1) (    2) (    3) (    4) (    5) (    6)
!!
!!
!!
!!  Real space logical address
!!  ==========================
!!
!!  (-2, 1) (-1, 1) ( 0, 1) ( 1, 1) ( n/a ) ( n/a )
!!
!!  (-2, 0) (-1, 0) ( 0, 0) ( 1, 0) ( n/a ) ( n/a )
!!
!!  (-2,-1) (-1,-1) ( 0,-1) ( 1,-1) ( n/a ) ( n/a )
!!
!!  (-2,-2) (-1,-2) ( 0,-2) ( 1,-2) ( n/a ) ( n/a )
!!
!!
!!
!!  Fourier space physical address
!!  ==============================
!!
!!  (     1,4     ) (     2,4     ) (     3,4     )
!!
!!  (     1,3     ) (     2,3     ) (     3,3     )
!!
!!  (     1,2     ) (     2,2     ) (     3,2     )
!!
!!  (     1,1     ) (     2,1     ) (     3,1     )
!!
!!
!!
!!  Fourier space logical address
!!  =============================
!!
!!  (     0,-1    ) (     1,-1    ) (     2,-1    )
!!
!!  (     0, 2    ) (     1, 2    ) (     2, 2    )
!!
!!  (     0, 1    ) (     1, 1    ) (     2, 1    )
!!
!!  (     0, 0    ) (     1, 0    ) (     2, 0    )
!!
!!
!!
module Images
    use Globals
    use, intrinsic :: iso_c_binding
    implicit none

    private

    ! Stuff needed by Frealix. Could/should be moved elsewhere
    public  ::  IndexOfCentralPixelGivenLogicalDimension

    ! Unit tests
    public  ::  image_fft_unit_tests, image_base_unit_tests

    ! This cannot be type-bound because it works on an array
    public  ::  SearchForRotationalAlignmentWith

    type, public :: Image
        !private
        integer                                             ::  logical_dimensions(3)                           =   (/0,0,0/)   !<  Logical (X,Y,Z) dimensions of the image.
                                                                                                                                !!Note that this does not necessarily correspond to memory allocation dimensions (ie physical dimensions).
        logical                                             ::  is_in_real_space                                =   .true.      !<  Whether the image is in real or Fourier space
        logical                                             ::  object_is_centered_in_box                       =   .true.      !<  Whether the object or region of interest is near the center of the box (as opposed to near the corners and wrapped around). This
                                                                                                                                !!  refers to real space and is meaningless in Fourier space.
        ! Looping in real and Fourier space may require these
        integer                                             ::  physical_upper_bound_complex(3)                 =   (/0,0,0/)   !<  In each dimension, the upper bound of the complex image's physical addresses
        integer                                             ::  physical_address_of_box_center(3)               =   (/0,0,0/)   !<  In each dimension, the address of the pixel at the origin
        integer                                             ::  physical_index_of_first_negative_frequency(3)   =   (/0,0,0/)   !<  In each dimension, the physical index of the first pixel which stores negative frequencies
        real                                                ::  fourier_voxel_size(3)                          =   (/0.0,0.0,0.0/) !<  Distance from Fourier voxel to Fourier voxel, expressed in reciprocal pixels
        integer                                             ::  logical_upper_bound_complex(3)                  =   (/0,0,0/)   !<  In each dimension, the upper bound of the complex image's logical addresses
        integer                                             ::  logical_lower_bound_complex(3)                  =   (/0,0,0/)   !<  In each dimension, the lower bound of the complex image's logical addresses
        integer                                             ::  logical_upper_bound_real(3)                     =   (/0,0,0/)   !<  In each dimension, the upper bound of the real image's logical addresses
        integer                                             ::  logical_lower_bound_real(3)                     =   (/0,0,0/)   !<  In each dimension, the lower bound of the real image's logical addresses


        ! Arrays to hold voxel values
        type(c_ptr)                                         ::  p                                               =   c_null_ptr  !<  C pointer to the image data array, which will be allocated by FFTW
        real(kind=c_float),             &
#ifdef __INTEL_COMPILER
                                        contiguous, &                                                                           !   The contiguous attribute should help with optimization, vectorization etc.
#endif
                                        pointer, public     ::  real_values(:,:,:)                              =>  null()      !<  Real array to hold values for REAL images.
        complex(kind=c_float_complex),  &
#ifdef __INTEL_COMPILER
                                        contiguous, &                                                                           !   The contiguous attribute should help with optimization, vectorization etc.
#endif
                                        pointer, public     ::  complex_values(:,:,:)                           =>  null()      !<  Complex array to hold values for COMP images.
        logical                                             ::  is_in_memory                                    =   .false.     !<  Whether image values are in-memory, in other words whether the image has memory space allocated to its data array. Default = .FALSE.

        ! FFTW-specific
        type(c_ptr)                                         ::  plan_fwd                                        =   c_null_ptr  !<  FFTW plan for the image (fwd)
        type(c_ptr)                                         ::  plan_bwd                                        =   c_null_ptr  !<  FFTW plan for the image (bwd)
        logical                                             ::  planned                                         =   .false.     !<  Whether the plan has been setup by/for FFTW
        contains
            procedure,  public  ::  PrintInfo
            procedure,  public  ::  PrintFourierStatistics
            procedure,  public  ::  PrintValues
            ! Fast Fourier transforms
            procedure,  public  ::  ForwardFFT
            procedure,  public  ::  BackwardFFT
            ! Disk I/O
            procedure           ::  ReadFromImageFile
            procedure           ::  ReadFromDiskGivenFilename
            procedure           ::  WriteToImageFile
            procedure           ::  WriteToDiskGivenFilename
            generic,    public  ::  ReadFromDisk => ReadFromImageFile, ReadFromDiskGivenFilename
            generic,    public  ::  WriteToDisk  => WriteToImageFile,  WriteToDiskGivenFilename
            ! Queries
            procedure,  public  ::  IsInRealSpace
            procedure,  public  ::  IsInFourierSpace
            procedure           ::  IsSquare
            procedure           ::  IsAVolume
            procedure           ::  IsBinarised
            procedure           ::  IsConstant
            procedure           ::  HasNan
            procedure,  public  ::  HasEvenDimensions
            procedure,  public  ::  IsAllocated
            procedure,  public  ::  ObjectIsCenteredInBox
            ! Comparisons
            procedure           ::  HasSameDimensionsAs
            procedure           ::  IsInSameSpaceAs
            generic,    public  ::  operator(.SameDimensions.) => HasSameDimensionsAs
            generic,    public  ::  operator(.SameSpace.) => IsInSameSpaceAs
            ! Correlations
            procedure,  public  ::  CalculateCrossCorrelationImageWith
            procedure,  public  ::  GetFSCWith
            procedure,  public  ::  GetCorrelationWithCTF
            procedure,  public  ::  ComputeFrealign9SimilarityMeasure
            ! CTF
            procedure           ::  CTFOperation
            procedure,  public  ::  OverlayCTF
            procedure,  public  ::  ApplyCTF
            ! Pixel/voxel value statistics
            procedure           ::  GetMinimumValue
            procedure,  public  ::  GetMaximumValue
            procedure           ::  GetAverageOfValues
            procedure           ::  GetAverageOfValuesOnEdges
            procedure           ::  GetAverageOfValuesInRing
            procedure           ::  GetSigmaOfValues
            procedure           ::  GetVarianceOfValues
            procedure           ::  GetSumOfValues
            procedure           ::  GetSumOfValuesDouble
            procedure,  public  ::  ComputeHistogramOfValues
            procedure,  public  ::  ComputeAverageAndSigmaOfValuesInSpectrum
            ! Arithmetics
            procedure,  public  ::  MultiplyByConstant
            procedure,  public  ::  SetMaximumValue
            procedure,  public  ::  SetMinimumAndMaximumValue
            procedure,  public  ::  SetMaximumValueOnCentralCross
            procedure,  public  ::  AddConstant

            procedure,  public  ::  TaperEdges
            procedure,  public  ::  ApplyMirror
            procedure,  public  ::  ApplyRotationalAverage
            procedure,  public  ::  ApplyCylindricalAverage
            procedure,  public  ::  Compute1DRotationalAverage
            procedure,  public  ::  ComputeRotationalAverageOfPowerSpectrum
            procedure,  public  ::  Compute1DProjection
            ! Normalisations
            procedure,  public  ::  ZeroFloat
            procedure,  public  ::  ZeroFloatAndNormalise
            procedure,  public  ::  NormaliseVariance
            procedure           ::  ComputeSigmaOfResolutionZones
            procedure           ::  NormaliseSigmaOfResolutionZones
            ! Noise (the Gaussian and other noises are in a separate module because of the GSL dependency)
            procedure,  public  ::  FillWithWhiteNoise
            ! Two images (todo: should ClipInto, AddInto and CutFromOtherImage be consolidated? Their differences should be clarified anyway - perhaps better naming?)
            procedure,  public  ::  AddImage
            procedure,  public  ::  MultiplyImage
            procedure,  public  ::  DivideImage
            procedure,  public  ::  SubtractImage
            procedure,  public  ::  ClipInto
            procedure,  public  ::  AddInto
            procedure,  public  ::  CutFromOtherImage
            ! Masking
            procedure,  public  ::  ApplyCircularMask
            procedure,  public  ::  ApplySoftCircularMask
            procedure,  public  ::  ApplyRectangularMask
            procedure,  public  ::  ApplySoftRectangularMask
            procedure,  public  ::  MaskCentralCross
            ! Filtering
            procedure,  public  ::  ApplyBFactor
            procedure,  public  ::  ApplyBandPassFilter
            procedure,  public  ::  ApplyDriftFilter
            procedure,  public  ::  ComputeDriftFilter
            procedure,  public  ::  SpectrumBoxConvolution
            procedure,  public  ::  RemoveLinearGradient
            procedure,  public  ::  ReplaceWithLinearGradient
            procedure,  public  ::  ApplyMedianFilter
            ! Amplitudes & phases
            procedure           ::  GetAmplitudeSpectrum
            procedure           ::  GetPhaseSpectrum
            procedure           ::  ComputeAmplitudeSpectrum
            procedure           ::  ComputePhaseSpectrum
            procedure           ::  SpectrumFromHalfToFull
            procedure,  public  ::  RandomisePhases
            ! Shifts, rotations, transformations
            procedure,  public  ::  PhaseShift
            procedure           ::  SwapRealSpaceQuadrants
            procedure,  public  ::  ApplyAffineTransformation
            procedure,  public  ::  CorrectMagnificationDistortion
            ! Peaks
            procedure,  public  ::  FindPeakWithIntegerCoordinates
            procedure,  public  ::  FindPeakWithParabolaFit
            ! Get pixel values
            procedure           ::  GetSlice
            procedure           ::  GetRealValueByLinearInterpolationNoBoundsCheckVolume
            procedure           ::  GetRealValueByLinearInterpolationNoBoundsCheckImage
            procedure           ::  GetComplexValueByLinearInterpolationNoBoundsCheckVolume
            procedure           ::  GetComplexValueByLinearInterpolationNoBoundsCheckVolumeVector
            procedure           ::  GetComplexValueByNearestNeighborNoBoundsCheckVolume
            ! Insert & extract slices
            procedure,  public  ::  ExtractComplexSlice
            procedure,  public  ::  ExtractComplexSliceVector
            procedure,  public  ::  InsertComplexSlice
            procedure,  public  ::  InsertRealSlice
            ! Put pixel values
            procedure           ::  PutRealValueByLinearInterpolationNoBoundsCheckImage
            procedure           ::  PutRealValueByLinearInterpolationNoBoundsCheckVolume
            procedure           ::  PutComplexValueByLinearInterpolationNoBoundsCheckVolume
            ! Get & set dimensions (set should be kept private - Allocate should be used instead)
            procedure           ::  SetLogicalDimensions_1
            procedure           ::  SetLogicalDimensions_2
            generic             ::  SetLogicalDimensions => SetLogicalDimensions_1,SetLogicalDimensions_2
            procedure,  public  ::  GetLogicalDimensions
            procedure,  public  ::  GetLogicalDimension
            procedure,  public  ::  GetPhysicalAddressOfBoxCenter
            procedure,  public  ::  GetPhysicalIndexOfBoxCenter
            procedure,  public  ::  GetFourierVoxelSize
            procedure,  public  ::  GetPhysicalUpperBoundComplex
            procedure           ::  UpdateLoopingAndAddressing
            procedure           ::  UpdatePhysicalAddressOfBoxCenter
            procedure,  public  ::  Resize

            procedure           ::  GetMaximumRadius
            procedure           ::  GetMaximumDiagonalRadius
            procedure           ::  GetRadiiGivenFraction
            ! Convenience functions when looping
            procedure           ::  PhysicalAddressGivenLogicalAddressInFourierSpace
            procedure           ::  LogicalIndexGivenPhysicalIndexInFourierSpace
            procedure           ::  FourierComponentIsRedundant
            procedure           ::  FourierComponentHasExplicitHermitianMate
            procedure           ::  NumberOfResolutionZonesFromTheirWidth
            procedure,  nopass  ::  ResolutionZoneFromSpatialFrequency
            procedure           ::  PhysicalAddressOfHermitianMate
            procedure           ::  PhysicalIndicesOfHermitianMate
            ! Memory management
            procedure,  public  ::  Allocate
            procedure,  public  ::  Deallocate
            procedure,  public  ::  Reset
            procedure           ::  AppendToArray
            final               ::  Destruct
            final               ::  DestructArray
            ! Assignments
            procedure           ::  AssignComplexToImage
            procedure           ::  AssignRealToImage
            procedure           ::  AssignImageToImage
            generic,    public  ::  assignment(=) => AssignComplexToImage,AssignRealToImage,AssignImageToImage
            ! Processing
            procedure           ::  LookForCrossBetaLayerLine

    end type Image

    contains


    !>  \brief  Add a random number to each voxel, drawn from a uniform distribution between 0.0 and 1.0
    subroutine FillWithWhiteNoise(self)
        class(Image),   intent(inout)   ::  self
        call random_seed()
        call random_number(self%real_values)
    end subroutine FillWithWhiteNoise

    subroutine AddConstant(self,constant)
        class(Image),           intent(inout)       ::  self
        real,                   intent(in)          ::  constant
        if (self%IsInRealSpace()) then
            self%real_values = self%real_values + constant
        else
            self%complex_values = self%complex_values + constant
        endif
    end subroutine AddConstant

    !>  \brief  Multiply pixel values in an image by a constant
    subroutine MultiplyByConstant(self,constant)
        class(Image),   intent(inout)   ::  self    !<  Image object whose pixel/voxel values are to be multiplied
        real,           intent(in)      ::  constant
        ! Start work
        if (self%IsInRealSpace()) then
            self%real_values = self%real_values * constant
        else
            self%complex_values = self%complex_values * constant
        endif
    end subroutine MultiplyByConstant


    subroutine ZeroFloat(self)
        class(Image),           intent(inout)       ::  self

        self%real_values = self%real_values - self%GetAverageOfValues()
    end subroutine ZeroFloat


    !>  \brief  On output, image densities will be 0.0 on average and will have the desired standard deviation
    !>  \todo   Could be made faster
    subroutine ZeroFloatAndNormalise(self,new_standard_deviation)
        class(Image),           intent(inout)       ::  self
        real,                   intent(in)          ::  new_standard_deviation
        ! private variables

        ! start work
        self%real_values = self%real_values - self%GetAverageOfValues()
        self%real_values = self%real_values / self%GetSigmaOfValues() * new_standard_deviation
    end subroutine ZeroFloatAndNormalise

    !>  \brief  Scale the density values such that their variance on output is 1.0 (leave the average unchanged)
    subroutine NormaliseVariance(self)
        class(Image),           intent(inout)       ::  self
        ! private variables
        real    ::  average_value
        ! start work
        average_value = self%GetAverageOfValues()
        self%real_values = (self%real_values - average_value) / self%GetVarianceOfValues() + average_value
    end subroutine NormaliseVariance


    !>  \brief  Average 3d image cylindrically
    subroutine ApplyCylindricalAverage(self)
        class(Image),    intent(inout)   ::  self
        ! private variables
        real,       allocatable ::  average_value(:)            !  average value as a function of radius
        integer,    allocatable ::  num_of_values(:)
        integer,    parameter   ::  num_central_slices  =   23  !   number of central slices over which the average should be measured
        integer,    parameter   ::  radial_oversampling =   4
        integer ::  start_slice, finish_slice
        integer ::  i,j,k,ori(3), addr_logi(3)
        integer ::  num_of_shells, current_shell
        integer ::  ierr
        logical,    parameter   ::  debug   =   .false.

        ! start work

        ! only meant for real images
        if (.not. self%is_in_real_space) then
            call this_program%TerminateWithFatalError('Image::ApplyCylindricalAverage','Image is not in real space')
        endif

        ! get origin address
        if (self%is_in_real_space .and. self%object_is_centered_in_box) then
            ori = self%physical_address_of_box_center
        else
            ori = [1,1,1]
        endif

        ! work out the start and finish slices
        if (self%logical_dimensions(3) .le. num_central_slices) then
            start_slice = 1
            finish_slice = self%logical_dimensions(3)
        else
            start_slice  = ori(3) - (num_central_slices-1)/2
            finish_slice = ori(3) + (num_central_slices-1)/2
        endif

        ! work out dimensions for allocatable array
        num_of_shells = maxval(self%logical_dimensions(1:2)-ori(1:2)+1)*radial_oversampling

        ! allocate arrays
        allocate(average_value(num_of_shells),num_of_values(num_of_shells),stat=ierr)
        if (ierr .ne. 0) call this_program%TerminateWithFatalError('Image::ApplyCylindricalAverage','Image is not in real space')

        ! initialise arrays
        average_value = 0.0
        num_of_values = 0

        if (debug) then
            write(*,'(a)') '**debug(ApplyCylindricalAverage): first loop over image...'
        endif

        ! loop over image to compute average value as a function of radius
        do k=start_slice,finish_slice
            do j=1,self%logical_dimensions(2)
                addr_logi(2) = (j - ori(2))**2
                do i=1,self%logical_dimensions(1)
                    addr_logi(1) = (i - ori(1))**2
                    ! what shell are in?
                    current_shell = nint(sqrt(real(sum(addr_logi(1:2)))+1.0)*real(radial_oversampling))
                    ! if we are within our limits, accumulate statistics
                    if (current_shell .le. size(num_of_values)) then
                        num_of_values(current_shell) = num_of_values(current_shell) + 1
                        average_value(current_shell) = average_value(current_shell) + self%real_values(i,j,k)
                    endif
                enddo
            enddo
        enddo

        if (debug) then
            write(*,'(a)') '**debug(ApplyCylindricalAverage): computing average as a funciton of radius...'
        endif

        ! do the actual averaging
        where (num_of_values == 0)
            average_value = 0.0
        else where
            average_value = average_value / real(num_of_values)
        end where

        if (debug) then
            write(*,'(a)') '**debug(ApplyCylindricalAverage): replacing values in image...'
        endif

        ! loop over the image to replace values with cylindrical average
        do k=1,self%logical_dimensions(3)
            do j=1,self%logical_dimensions(2)
                addr_logi(2) = (j - ori(2))**2
                do i=1,self%logical_dimensions(1)
                    addr_logi(1) = (i - ori(1))**2
                    ! what shell are in?
                    current_shell = nint(sqrt(real(sum(addr_logi(1:2)))+1.0)*real(radial_oversampling))
                    ! replace value
                    if (current_shell .le. size(num_of_values)) then
                        self%real_values(i,j,k) = average_value(current_shell)
                    else
                        ! self%real_values   =   0.0
                    endif
                enddo
            enddo
        enddo

        if (debug) then
            write(*,'(a)') '**debug(ApplyCylindricalAverage): all done'
        endif

    end subroutine ApplyCylindricalAverage

    !>  \brief  Rotational average
    subroutine ApplyRotationalAverage(self)
        class(Image),   intent(inout)   ::  self
        ! private variables
        real(kind=8),   allocatable ::  average(:)
        integer                     ::  i,j,k
        integer                     ::  i_logi,j_logi,k_logi
        real                        ::  rad
        ! start work
        if (.not. self%IsInRealSpace()) then
            call this_program%TerminateWithFatalError('Image::ApplyRotationalAverage','Not implemented for Fourier space')
        endif

        call self%Compute1DRotationalAverage(average)

        !
        do k=1,self%logical_dimensions(3)
            k_logi = (k-self%physical_address_of_box_center(3))**2
            do j=1,self%logical_dimensions(2)
                j_logi = (j-self%physical_address_of_box_center(2))**2 + k_logi
                do i=1,self%logical_dimensions(1)
                    i_logi = (i-self%physical_address_of_box_center(1))**2 + j_logi
                    !
                    rad = sqrt(real(i_logi)) + 1.0
                    !
                    self%real_values(i,j,k) = average(int(rad)  ) * (rad-int(rad)    )
                    self%real_values(i,j,k) = average(int(rad)+1) * (int(rad)-rad+1.0) + self%real_values(i,j,k)
                enddo
            enddo
        enddo


    end subroutine ApplyRotationalAverage

    !>  \brief  Compute the 1D rotational average
    !!          Each bin is 1 pixel wide and there are as many bins as fit in the diagonal of the image
    subroutine Compute1DRotationalAverage(self,average)
        class(Image),                   intent(in)      ::  self
        real(kind=8),                   intent(inout)   ::  average(:)
        ! private variables
        integer ::  number_of_bins
        integer ::  i,j,k
        integer ::  i_logi,j_logi,k_logi
        real    ::  rad
        real(kind=8)    ::  number_of_values(size(average))
        ! start work

        if (.not. self%IsInRealSpace()) then
            call this_program%TerminateWithFatalError('Image::Compute1DRotationalAverage','Not implemented for Fourier space')
        endif

        ! How many bins?
        !!  \todo   Use GetMaximumDiagonalRadius instead of this
        !number_of_bins = ceiling(sqrt(sum((1.0-real(self%physical_address_of_box_center))**2))+1)
        number_of_bins = size(average)

        ! Allocate memory for the result
        !if (allocated(average)) deallocate(average)
        !allocate(average(number_of_bins),number_of_values(number_of_bins))
        average = 0.0
        number_of_values = 0.0

        !
        do k=1,self%logical_dimensions(3)
            k_logi = (k-self%physical_address_of_box_center(3))**2
            do j=1,self%logical_dimensions(2)
                j_logi = (j-self%physical_address_of_box_center(2))**2 + k_logi
                do i=1,self%logical_dimensions(1)
                    i_logi = (i-self%physical_address_of_box_center(1))**2 + j_logi
                    !
                    rad = sqrt(real(i_logi))+1.0 ! Because the first element, indexed at 1, is the origin
                    !
                    average(int(rad)  ) = average(int(rad)  ) + (rad-int(rad)    )*self%real_values(i,j,k)
                    average(int(rad)+1) = average(int(rad)+1) + (int(rad)-rad+1.0)*self%real_values(i,j,k)
                    !
                    !number_of_values(int(rad)  ) = number_of_values(int(rad)  ) + 1
                    !number_of_values(int(rad)+1) = number_of_values(int(rad)+1) + 1
                    !
                    number_of_values(int(rad)  ) = number_of_values(int(rad)  ) + (rad-int(rad)    )
                    number_of_values(int(rad)+1) = number_of_values(int(rad)+1) + (int(rad)-rad+1.0)
                enddo
            enddo
        enddo

        where (number_of_values .ne. 0.0)
            average = average / number_of_values
        elsewhere
            average = 0.0
        endwhere

    end subroutine Compute1DRotationalAverage

    !>  \brief  Compute the 1D projection along one of the axes
    subroutine Compute1DProjection(self,along_this_axis,projection)
        class(Image),                   intent(inout)   ::  self
        integer,                        intent(in)      ::  along_this_axis !<  1,2 or 3. For example, if 1, the projection will be computed along the Y and Z axes, such that a 1D profile along X is obtained
        real,           allocatable,    intent(inout)   ::  projection(:)   !<  Output projection profile
        !
        integer     ::  i
        !
        if (allocated(projection)) deallocate(projection)
        allocate(projection(self%logical_dimensions(along_this_axis)))

        select case (along_this_axis)
            case(1)
                do i=1,self%logical_dimensions(1)
                    projection(i) = sum(self%real_values(i,:,:))
                enddo
            case(2)
                do i=1,self%logical_dimensions(2)
                    projection(i) = sum(self%real_values(1:self%logical_dimensions(1),i,:))
                enddo
            case(3)
                do i=1,self%logical_dimensions(3)
                    projection(i) = sum(self%real_values(1:self%logical_dimensions(1),:,i))
                enddo
            case default
                call this_program%TerminateWithFatalError('Image::Compute1DProjection','Bad value for axis number')
        end select

    end subroutine Compute1DProjection

    !>  \brief  Mirror the image along the given axis
    subroutine ApplyMirror(self,along_this_axis)
        class(Image),   intent(inout)   ::  self
        integer,        intent(in)      ::  along_this_axis !<  1, 2 or 3. For example, if 1, the mirror will be about the YZ plane, along the X axis.
        !
        integer             ::  i
        real,   allocatable ::  temp(:,:)
        !
        if (.not. self%IsInRealSpace()) then
            call this_program%TerminateWithFatalError('Image::ApplyMirror','Not implemented for Fourier space')
        endif

        select case(along_this_axis)
            case (1)
                allocate(temp(self%logical_dimensions(2),self%logical_dimensions(3)))
                do i=1,self%physical_address_of_box_center(1)-1
                    temp(:,:) = self%real_values(i,:,:)
                    self%real_values(i,:,:) = self%real_values(self%logical_dimensions(1)-i+1,:,:)
                    self%real_values(self%logical_dimensions(1)-i+1,:,:) = temp(:,:)
                enddo
            case (2)
                allocate(temp(self%logical_dimensions(1),self%logical_dimensions(3)))
                do i=1,self%physical_address_of_box_center(2)-1
                    temp(:,:) = self%real_values(1:self%logical_dimensions(1),i,:)
                    self%real_values(:,i,:) = self%real_values(:,self%logical_dimensions(2)-i+1,:)
                    self%real_values(1:self%logical_dimensions(1),self%logical_dimensions(2)-i+1,:) = temp(:,:)
                enddo
            case (3)
                allocate(temp(self%logical_dimensions(1),self%logical_dimensions(2)))
                do i=1,self%physical_address_of_box_center(3)-1
                    temp(:,:) = self%real_values(1:self%logical_dimensions(1),:,i)
                    self%real_values(:,:,i) = self%real_values(:,:,self%logical_dimensions(3)-i+1)
                    self%real_values(:,:,self%logical_dimensions(3)-i+1) = temp(:,:)
                enddo
            case default
                call this_program%TerminateWithFatalError('Image::ApplyMirror','Bad value for axis number')
        end select
    end subroutine ApplyMirror

    !>  \brief  Taper edges of image so that there are no sharp discontinuities in real space
    !!  This is a re-implementation of the MRC program taperedgek.for (Richard Henderson, 1987)
    subroutine TaperEdges(self)
        ! Arguments
        class(Image),               intent(inout)   ::  self
        ! Private variables
        integer,    parameter       ::  averaging_strip_width(3)    =   100
        integer,    parameter       ::  tapering_strip_width(3)     =   500
        integer,    parameter       ::  smoothing_half_width(3)     =   1
        integer                     ::  current_dimension, number_of_dimensions
        integer                     ::  second_dimension, third_dimension
        integer                     ::  i,j,k
        integer                     ::  j_shift,k_shift
        integer                     ::  jj,kk
        integer                     ::  number_of_values_in_running_average
        real,       allocatable     ::  average_for_current_edge_start  (:,:)
        real,       allocatable     ::  average_for_current_edge_finish (:,:)
        real,       allocatable     ::  average_for_current_edge_average(:,:)
        real,       allocatable     ::  smooth_average_for_current_edge_start  (:,:)
        real,       allocatable     ::  smooth_average_for_current_edge_finish (:,:)
        ! Start work

        number_of_dimensions = 2
        if (self%IsAVolume()) number_of_dimensions = 3

        do current_dimension=1,number_of_dimensions
            select case (current_dimension)
                case (1)
                    second_dimension = 2
                    third_dimension  = 3
                case (2)
                    second_dimension = 1
                    third_dimension  = 3
                case (3)
                    second_dimension = 1
                    third_dimension  = 2
            end select


            ! Deallocate memory
            if (allocated(average_for_current_edge_start))      deallocate(average_for_current_edge_start)
            if (allocated(average_for_current_edge_finish))     deallocate(average_for_current_edge_finish)
            if (allocated(average_for_current_edge_average))    deallocate(average_for_current_edge_average)
            if (allocated(smooth_average_for_current_edge_start))   deallocate(smooth_average_for_current_edge_start)
            if (allocated(smooth_average_for_current_edge_finish))  deallocate(smooth_average_for_current_edge_finish)

            ! Allocate memory
            allocate(   average_for_current_edge_start  (self%logical_dimensions(second_dimension), &
                                                         self%logical_dimensions(third_dimension)), &
                        average_for_current_edge_finish (self%logical_dimensions(second_dimension), &
                                                         self%logical_dimensions(third_dimension)), &
                        average_for_current_edge_average(self%logical_dimensions(second_dimension), &
                                                         self%logical_dimensions(third_dimension)), &
                        smooth_average_for_current_edge_start  (self%logical_dimensions(second_dimension), &
                                                                self%logical_dimensions(third_dimension)), &
                        smooth_average_for_current_edge_finish (self%logical_dimensions(second_dimension), &
                                                                self%logical_dimensions(third_dimension)) &
                    )


            ! Initialise memory
            average_for_current_edge_start         = 0.0e0
            average_for_current_edge_finish        = 0.0e0
            average_for_current_edge_average       = 0.0e0
            smooth_average_for_current_edge_start  = 0.0e0
            smooth_average_for_current_edge_finish = 0.0e0

            !
            ! Deal with X=0 and X=self%logical_dimensions(1) edges
            !

            i=1
            do k=1,self%logical_dimensions(third_dimension)
                do j=1,self%logical_dimensions(second_dimension)
                    select case (current_dimension)
                        case (1)
                            average_for_current_edge_start(j,k)                         &
                                = sum(self%real_values(1:averaging_strip_width(current_dimension),j,k)) &
                                    / averaging_strip_width(current_dimension)
                            average_for_current_edge_finish(j,k)                                                         &
                                = sum(self%real_values(self%logical_dimensions(current_dimension) &
                                                       -averaging_strip_width(1)+1:self%logical_dimensions(current_dimension) &
                                                       ,j,k)) &
                                    / averaging_strip_width(current_dimension)
                        case (2)
                            average_for_current_edge_start(j,k)                         &
                                = sum(self%real_values(j,1:averaging_strip_width(current_dimension),k)) &
                                    / averaging_strip_width(current_dimension)
                            average_for_current_edge_finish(j,k)                                                         &
                                = sum(self%real_values(j,self%logical_dimensions(current_dimension) &
                                                       -averaging_strip_width(1)+1:self%logical_dimensions(current_dimension) &
                                                       ,k)) &
                                    / averaging_strip_width(current_dimension)
                        case (3)
                            average_for_current_edge_start(j,k)                         &
                                = sum(self%real_values(j,k,1:averaging_strip_width(current_dimension))) &
                                    / averaging_strip_width(current_dimension)
                            average_for_current_edge_finish(j,k)                                                         &
                                = sum(self%real_values(j,k,self%logical_dimensions(current_dimension) &
                                                       -averaging_strip_width(1)+1:self%logical_dimensions(current_dimension) &
                                                       )) &
                                    / averaging_strip_width(current_dimension)
                    end select
                enddo
            enddo

            average_for_current_edge_average = 0.5e0 * (  average_for_current_edge_finish &
                                                        + average_for_current_edge_start )
            average_for_current_edge_start   =    average_for_current_edge_start &
                                                - average_for_current_edge_average
            average_for_current_edge_finish  =    average_for_current_edge_finish &
                                                - average_for_current_edge_average


            ! Apply smoothing parallel to edge in the form of a running average
            do k=1,self%logical_dimensions(third_dimension)
                do j=1,self%logical_dimensions(second_dimension)
                    number_of_values_in_running_average = 0
                    ! Loop over neighbourhood of non-smooth arrays
                    do k_shift=-smoothing_half_width(third_dimension),smoothing_half_width(third_dimension)
                        kk = k+k_shift
                        if (kk .lt. 1 .or. kk .gt. self%logical_dimensions(third_dimension)) cycle
                        do j_shift=-smoothing_half_width(second_dimension),smoothing_half_width(second_dimension)
                            jj = j+j_shift
                            if (jj .lt. 1 .or. jj .gt. self%logical_dimensions(second_dimension)) cycle
                            number_of_values_in_running_average = number_of_values_in_running_average + 1
                            smooth_average_for_current_edge_start (j,k) = smooth_average_for_current_edge_start(j,k) &
                                                                        + average_for_current_edge_start(jj,kk)
                            smooth_average_for_current_edge_finish(j,k) = smooth_average_for_current_edge_finish(j,k) &
                                                                        + average_for_current_edge_finish(jj,kk)
                        enddo
                    enddo
                    ! Now we can compute the average
                    smooth_average_for_current_edge_start (j,k) = smooth_average_for_current_edge_start(j,k)  &
                                                                / number_of_values_in_running_average
                    smooth_average_for_current_edge_finish(j,k) = smooth_average_for_current_edge_finish(j,k) &
                                                                / number_of_values_in_running_average
                enddo
            enddo

            ! Taper the image
            do i=1,self%logical_dimensions(current_dimension)
                if (i .le. tapering_strip_width(current_dimension)) then
                    select case (current_dimension)
                        case (1)
                            self%real_values(i,:,:) = self%real_values(i,:,:) &
                                                    - smooth_average_for_current_edge_start (:,:) &
                                                    * (tapering_strip_width(current_dimension)-i+1) &
                                                    /tapering_strip_width(current_dimension)
                        case (2)
                            self%real_values(1:self%logical_dimensions(1),i,:)  &
                                                    = self%real_values(1:self%logical_dimensions(1),i,:) &
                                                    - smooth_average_for_current_edge_start (:,:) &
                                                    * (tapering_strip_width(current_dimension)-i+1) &
                                                    /tapering_strip_width(current_dimension)
                        case (3)
                            self%real_values(1:self%logical_dimensions(1),:,i) &
                                                    = self%real_values(1:self%logical_dimensions(1),:,i) &
                                                    - smooth_average_for_current_edge_start (:,:) &
                                                    * (tapering_strip_width(current_dimension)-i+1) &
                                                    /tapering_strip_width(current_dimension)
                    end select
                else if (i .ge. self%logical_dimensions(current_dimension)-tapering_strip_width(current_dimension)+1) then
                    select case (current_dimension)
                        case (1)
                            self%real_values(i,:,:) = self%real_values(i,:,:) &
                                                    - smooth_average_for_current_edge_finish(:,:) &
                                                    * (tapering_strip_width(current_dimension)+i &
                                                    -self%logical_dimensions(current_dimension)) &
                                                    /tapering_strip_width(current_dimension)
                        case (2)
                            self%real_values(1:self%logical_dimensions(1),i,:) &
                                                    = self%real_values(1:self%logical_dimensions(1),i,:) &
                                                    - smooth_average_for_current_edge_finish(:,:) &
                                                    * (tapering_strip_width(current_dimension)+i &
                                                    -self%logical_dimensions(current_dimension)) &
                                                    /tapering_strip_width(current_dimension)
                        case (3)
                            self%real_values(1:self%logical_dimensions(1),:,i) &
                                                    = self%real_values(1:self%logical_dimensions(1),:,i) &
                                                    - smooth_average_for_current_edge_finish(:,:) &
                                                    * (tapering_strip_width(current_dimension)+i &
                                                    -self%logical_dimensions(current_dimension)) &
                                                    /tapering_strip_width(current_dimension)
                    end select
                endif
            enddo
        enddo ! End of loop over dimensions



    end subroutine TaperEdges

    !>  \brief  Add values into the other image
    subroutine AddInto(self,other_image,address_of_center,overwrite)
        ! Arguments
        class(Image),               intent(in)      ::  self
        type(Image),                intent(inout)   ::  other_image
        integer,        optional,   intent(in)      ::  address_of_center(3)    !<  Address of pixel in other_image at which self will be centered
        logical,        optional,   intent(in)      ::  overwrite               !<  Overwrite existing densities in other_image rather than adding to them
        ! Private variables
        integer ::  aaddress_of_center(3)
        integer ::  i,j,k
        integer ::  i_logi,j_logi,k_logi
        integer ::  ii,jj,kk
        logical ::  ooverwrite
        ! Start work

        ! Optional arguments
        aaddress_of_center = other_image%physical_address_of_box_center
        if (present(address_of_center)) aaddress_of_center = address_of_center

        ooverwrite = .false.
        if (present(overwrite)) ooverwrite = overwrite

        ! The images should be in the same space
        if (.not. other_image%IsInSameSpaceAs(self)) then
            call this_program%TerminateWithFatalError('Image::PasteInto','Images are not in same space')
        endif

        !
        if (self%IsInRealSpace()) then
            if (ooverwrite) then
                ! Loop over input image
                do k=1,self%logical_dimensions(3)
                    k_logi = k-self%physical_address_of_box_center(3)
                    kk = aaddress_of_center(3) + k_logi
                    if (kk .lt. 1 .or. kk .gt. other_image%logical_dimensions(3)) cycle
                    do j=1,self%logical_dimensions(2)
                        j_logi = j-self%physical_address_of_box_center(2)
                        jj = aaddress_of_center(2) + j_logi
                        if (jj .lt. 1 .or. jj .gt. other_image%logical_dimensions(2)) cycle
                        do i=1,self%logical_dimensions(1)
                            i_logi = i-self%physical_address_of_box_center(1)
                            ii = aaddress_of_center(1) + i_logi
                            if (ii .lt. 1 .or. ii .gt. other_image%logical_dimensions(1)) cycle
                            !
                            other_image%real_values(ii,jj,kk) = self%real_values(i,j,k)
                        enddo
                    enddo
                enddo
            else
                ! Loop over input image
                do k=1,self%logical_dimensions(3)
                    k_logi = k-self%physical_address_of_box_center(3)
                    kk = aaddress_of_center(3) + k_logi
                    if (kk .lt. 1 .or. kk .gt. other_image%logical_dimensions(3)) cycle
                    do j=1,self%logical_dimensions(2)
                        j_logi = j-self%physical_address_of_box_center(2)
                        jj = aaddress_of_center(2) + j_logi
                        if (jj .lt. 1 .or. jj .gt. other_image%logical_dimensions(2)) cycle
                        do i=1,self%logical_dimensions(1)
                            i_logi = i-self%physical_address_of_box_center(1)
                            ii = aaddress_of_center(1) + i_logi
                            if (ii .lt. 1 .or. ii .gt. other_image%logical_dimensions(1)) cycle
                            !
                            other_image%real_values(ii,jj,kk) = other_image%real_values(ii,jj,kk) + self%real_values(i,j,k)
                        enddo
                    enddo
                enddo
            endif ! end of test for overwrite
        else
            call this_program%TerminateWithFatalError('Image::PasteInto','Not implemented for Fourier space')
        endif


    end subroutine AddInto

    !>  \brief  Copy values into the other image. The centering of the image or of its Fourier transform is preserved, but its logical dimensions may change
    subroutine ClipInto(self,other_image,address_of_center,padding_value_real)
        ! Arguments
        class(Image),               intent(in)      ::  self
        type(Image),                intent(inout)   ::  other_image
        integer,        optional,   intent(in)      ::  address_of_center(3)    !<  Address of pixel in self at which other_image will be centered
        real,           optional,   intent(in)      ::  padding_value_real
        ! Private variables
        integer             ::  i,j,k
        integer             ::  ii,jj,kk
        integer             ::  ii_logi,jj_logi,kk_logi
        integer             ::  logical_address(3)
        integer             ::  physical_address(3)
        integer             ::  aaddress_of_center(3)
        real                ::  ppadding_value_real
        complex, parameter  ::  ppadding_value_complex = (0.0e0,0.0e0)
        ! Start work

        ! Optional arguments
        if (present(address_of_center)) then
            aaddress_of_center = address_of_center
            if (.not. self%is_in_real_space) then
                write(*,'(a)') '**error(ClipInto): clipping off-center only supported in real space'
                call this_program%TerminateWithFatalError('Image::ClipInto','Feature is not available')
            endif
        else
            aaddress_of_center = self%physical_address_of_box_center
        endif

        ppadding_value_real = 0.0e0
        if (present(padding_value_real)) ppadding_value_real = padding_value_real

        ! The other image takes the following attributes from self
        other_image%is_in_real_space = self%is_in_real_space
        other_image%object_is_centered_in_box = self%object_is_centered_in_box

        ! We may need to allocate the other_image
        if (.not. other_image%IsAllocated()) call other_image%Allocate()

        !
        if (self%is_in_real_space) then
            if (self%object_is_centered_in_box) then
                ! In real space, with object centered in the box
                ! Loop over the other image
                do kk=1,other_image%logical_dimensions(3)
                    kk_logi = kk-other_image%physical_address_of_box_center(3)
                    k = aaddress_of_center(3) + kk_logi
                    if (k .lt. 1 .or. k .gt. self%logical_dimensions(3)) then
                        other_image%real_values(:,:,kk) = ppadding_value_real
                        cycle
                    endif

                    do jj=1,other_image%logical_dimensions(2)
                        jj_logi = jj-other_image%physical_address_of_box_center(2)
                        j = aaddress_of_center(2) + jj_logi
                        if (j .lt. 1 .or. j .gt. self%logical_dimensions(2)) then
                            other_image%real_values(:,jj,kk) = ppadding_value_real
                            cycle
                        endif

                        do ii=1,other_image%logical_dimensions(1)
                            ii_logi = ii-other_image%physical_address_of_box_center(1)
                            i = aaddress_of_center(1) + ii_logi
                            if (i .lt. 1 .or. i .gt. self%logical_dimensions(1)) then
                                other_image%real_values(ii,jj,kk) = ppadding_value_real
                            else
                                other_image%real_values(ii,jj,kk) = self%real_values(i,j,k)
                            endif

                        enddo
                    enddo
                enddo

            else
                ! In real space, object not centered in box
                call this_program%TerminateWithFatalError('Image::ClipInto','Object is not centered in box')
            endif
        else
            ! In Fourier space

            ! Loop over the other image
            do kk=1,other_image%physical_upper_bound_complex(3)
                logical_address(3) = other_image%LogicalIndexGivenPhysicalIndexInFourierSpace(kk,3)
                if (     logical_address(3) .gt. self%logical_upper_bound_complex(3) &
                    .or. logical_address(3) .lt. self%logical_lower_bound_complex(3)) cycle
                do jj=1,other_image%physical_upper_bound_complex(2)
                    logical_address(2) = other_image%LogicalIndexGivenPhysicalIndexInFourierSpace(jj,2)
                    if (     logical_address(2) .gt. self%logical_upper_bound_complex(2) &
                        .or. logical_address(2) .lt. self%logical_lower_bound_complex(2)) cycle
                    do ii=1,other_image%physical_upper_bound_complex(1)
                        logical_address(1) = ii-1
                        if (     logical_address(1) .gt. self%logical_upper_bound_complex(1) &
                            .or. logical_address(1) .lt. self%logical_lower_bound_complex(1)) cycle
                        call self%PhysicalAddressGivenLogicalAddressInFourierSpace(logical_address,physical_address)
                        if (any(physical_address .lt. 1) .or. physical_address(1) .gt. size(self%complex_values,1) &
                                                         .or. physical_address(2) .gt. size(self%complex_values,2) &
                                                         .or. physical_address(3) .gt. size(self%complex_values,3)  ) then
                            other_image%complex_values(ii,jj,kk) = ppadding_value_complex
                        else
                            other_image%complex_values(ii,jj,kk) = self%complex_values( physical_address(1), &
                                                                                        physical_address(2), &
                                                                                        physical_address(3))
                        endif
                    enddo
                enddo
            enddo
        endif

    end subroutine ClipInto

    !>  \brief  Resize the image so it has the new logical dimensions given as argument
    subroutine Resize(self,new_logical_dimensions,padding_value_real)
        class(Image),               intent(inout)   ::  self
        integer,                    intent(in)      ::  new_logical_dimensions(3)
        real,           optional,   intent(in)      ::  padding_value_real
        ! Private variable
        type(Image) ::  temp_image
        ! Start work
        call temp_image%Allocate(dims=new_logical_dimensions)
        call self%ClipInto(temp_image,padding_value_real=padding_value_real)
        self = temp_image
    end subroutine Resize

    !>  \brief  Find the rotation angle which, when applied to self around the 3rd axis, brings it in register with other_image. Angle is left-handed positive.
    subroutine SearchForRotationalAlignmentWith(self,other_image,search_half_range,search_step_size,best_rotation, &
                                                print_progress, &
                                                minimum_radius,maximum_radius)
        use EmpiricalDistributions
        use ProgressBars
        ! Arguments
        class(Image),           intent(inout)   ::  self(:)
        type(Image),            intent(in)      ::  other_image(:)
        real,                   intent(in)      ::  search_half_range       !<  degrees
        real,                   intent(in)      ::  search_step_size        !<  degrees
        logical,                intent(in)      ::  print_progress
        real,       optional,   intent(in)      ::  minimum_radius          !<  As a fraction of the image size (i.e. 0.5 is the edge of the box)
        real,       optional,   intent(in)      ::  maximum_radius          !<  As a fraction of image size (0.5 is edge)
        ! Private variables
        integer                     ::  i,j
        real                        ::  i_logi,j_logi
        real                        ::  i_logi_frac,j_logi_frac
        real                        ::  ii_phys, jj_phys
        real                        ::  best_rotation, current_rotation_rad
        real                        ::  current_interpolated_value
        type(EmpiricalDistribution) ::  cc_numerator_dist, cc_denom_self_dist, cc_denom_other_dist
        real                        ::  current_cc, best_cc
        real(kind=8)                ::  current_rotation
        real                        ::  mminimum_radius, mmaximum_radius
        real                        ::  inverse_logical_dimensions(3)
        integer                     ::  current_image
        type(ProgressBar)           ::  my_progress_bar
        logical,    parameter       ::  debug = .false.
        ! Start work

        ! Both images must be in real space and 2D
        if (.not. self(1)%IsInRealSpace() .or. .not. other_image(1)%IsInRealSpace()) then
            call this_program%TerminateWithFatalError('Image::SearchWithRotationalAlignmentWith','Image must be in real space')
        endif
        if (self(1)%IsAVolume()) then
            call this_program%TerminateWithFatalError('Image::SearchWithRotationalAlignmentWith','Image must be 2D')
        endif

        ! Both sets need to have the same number of images
        if (size(self) .ne. size(other_image)) then
            call this_program%TerminateWithFatalError('Image::SearchWithRotationalAlignmentWith','Image stacks not same size')
        endif

        !
        mminimum_radius = 0.0
        if (present(minimum_radius)) mminimum_radius = minimum_radius**2
        mmaximum_radius = 1.0
        if (present(maximum_radius)) mmaximum_radius = maximum_radius**2

        !
        inverse_logical_dimensions = 1.0e0 / other_image(1)%logical_dimensions

        ! Set up progress bar
        if (print_progress) then
            call my_progress_bar%Begin(nint((2.0*search_half_range+search_step_size)/search_step_size))
        endif


        best_cc = -1.0e0
        best_rotation = 0.0e0
        current_rotation = - search_half_range
        do while (current_rotation .lt. search_half_range + search_step_size)
            if (print_progress) then
                call my_progress_bar%Update(my_progress_bar%GetCurrentTick()+1)
            endif
            current_rotation_rad = current_rotation / 180.0e0 * pi
            call cc_numerator_dist%Init()
            call cc_denom_self_dist%Init()
            call cc_denom_other_dist%Init()
            ! Loop over the array of images
            do current_image=1,size(self)
                ! Loop over the other image
                do j=1,other_image(current_image)%logical_dimensions(2)
                    j_logi=j-other_image(current_image)%physical_address_of_box_center(2)
                    j_logi_frac = (j_logi * inverse_logical_dimensions(2))**2
                    do i=1,other_image(current_image)%logical_dimensions(1)
                        i_logi=i-other_image(current_image)%physical_address_of_box_center(1)
                        i_logi_frac = (i_logi * inverse_logical_dimensions(1))**2 + j_logi_frac
                        if (i_logi_frac .ge. mminimum_radius .and. i_logi_frac .le. mmaximum_radius) then
                            ! We do ccw rotation to go from other_image (reference) to self (input image)
                            ii_phys = i_logi * cos(current_rotation_rad) - j_logi * sin(current_rotation_rad) &
                                        + self(current_image)%physical_address_of_box_center(1)
                            jj_phys = i_logi * sin(current_rotation_rad) + j_logi * cos(current_rotation_rad) &
                                        + self(current_image)%physical_address_of_box_center(2)
                            !
                            if (        ii_phys .gt. 1.0 .and. ii_phys .lt. real(self(current_image)%logical_dimensions(1)) &
                                .and.   jj_phys .gt. 1.0 .and. jj_phys .lt. real(self(current_image)%logical_dimensions(2)) ) then
                                call self(current_image)%GetRealValueByLinearInterpolationNoBoundsCheckImage(&
                                                                                                current_interpolated_value, &
                                                                                                ii_phys,jj_phys)
                                call cc_numerator_dist%AddSampleValue(current_interpolated_value&
                                                                        *other_image(current_image)%real_values(i,j,1))
                                call cc_denom_other_dist%AddSampleValue(other_image(current_image)%real_values(i,j,1)**2)
                                call cc_denom_self_dist%AddSampleValue(current_interpolated_value**2)
                            endif
                        endif
                    enddo ! i
                enddo ! j
            enddo   ! end of loop over image stack
            current_cc = cc_numerator_dist%GetSampleSum() &
                         / sqrt(cc_denom_other_dist%GetSampleSum()*cc_denom_self_dist%GetSampleSum())
            if (debug) write(*,'(a,f0.1,a,f0.4)') '**debug(Image::SearchForRotationalAlignmentWith): rot: ', current_rotation, &
                                                    ' score: ', current_cc
            if (current_cc .gt. best_cc) then
                 best_cc = current_cc
                 best_rotation = current_rotation
            endif
            ! Increment the rotation
            current_rotation = current_rotation + search_step_size
        enddo

        if (print_progress) then
            call my_progress_bar%Finish()
        endif


    end subroutine SearchForRotationalAlignmentWith


    !>  \brief  Add another image to this one
    subroutine AddImage(self,other_image,set_phases_to_zero,subtract)
        ! Arguments
        class(Image),           intent(inout)   ::  self
        type(Image),            intent(in)      ::  other_image
        logical,    optional,   intent(in)      ::  set_phases_to_zero
        logical,    optional,   intent(in)      ::  subtract            !<  Subtract other_image from self, rather than add
        ! Private variables
        logical ::  sset_phases_to_zero
        logical ::  ssubtract
        integer ::  i,j,k
        ! Start work
        sset_phases_to_zero = .false.
        if (present(set_phases_to_zero)) sset_phases_to_zero = set_phases_to_zero
        ssubtract = .false.
        if (present(subtract)) ssubtract = subtract
        if (.not. self%is_in_memory) then
            call self%Allocate(mould=other_image)
            self%real_values = 0.0
            if (sset_phases_to_zero) then
                self%is_in_real_space = .false.
            endif
        endif
        if (.not. self%HasSameDimensionsAs(other_image)) then
            call this_program%TerminateWithFatalError('Image::AddImage','Images do not have the same dimensions')
        endif
        if (.not. self%IsInSameSpaceAs(other_image)) then
            call this_program%TerminateWithFatalError('Image::AddImage','Images are not in the same space')
        endif
        if (sset_phases_to_zero) then
            if (self%is_in_real_space .or. other_image%is_in_real_space) then
                call this_program%TerminateWithFatalError('Image::AddImage','Image is in real space, but trying to drop phases.')
            endif
            if (ssubtract) then
                self%complex_values = self%complex_values - cmplx(cabs(other_image%complex_values),0.0)
            else
                self%complex_values = self%complex_values + cmplx(cabs(other_image%complex_values),0.0)
            endif
        else
            if (ssubtract) then
                self%real_values = self%real_values - other_image%real_values
            else
                do k=1,self%logical_dimensions(3)
                    do j=1,self%logical_dimensions(2)
                        do i=1,self%logical_dimensions(1)
                            self%real_values(i,j,k) = self%real_values(i,j,k) + other_image%real_values(i,j,k)
                        enddo
                    enddo
                enddo
            endif
        endif
    end subroutine AddImage


    !>  \brief  Multiply By Another Image
    subroutine MultiplyImage(self,other_image)
        ! Arguments
        class(Image),           intent(inout)   ::  self
        type(Image),            intent(in)      ::  other_image
        ! Private variables

        ! Start work

        if (.not. self%HasSameDimensionsAs(other_image)) then
            call this_program%TerminateWithFatalError('Image::MultiplyImage','Images do not have the same dimensions')
        endif
        !if (.not. self%IsInSameSpaceAs(other_image)) then
        !    call this_program%TerminateWithFatalError('Image::MultiplyImage','Images are not in the same space')
        !endif

        if (self%IsInRealSpace() .and. other_image%IsInRealSpace()) then
            self%real_values = self%real_values * other_image%real_values
        else if (self%IsInFourierSpace() .and. other_image%IsInFourierSpace()) then
            self%complex_values = self%complex_values * other_image%complex_values
        else if (self%IsInFourierSpace() .and. other_image%IsInRealSpace()) then
            self%complex_values = self%complex_values * other_image%real_values
        else
             call this_program%TerminateWithFatalError('Image::MultiplyImage','Images don not have compatible types')
        endif


    end subroutine MultiplyImage

        !>  \brief  Divide By Another Image
    subroutine DivideImage(self,other_image)
        ! Arguments
        class(Image),           intent(inout)   ::  self
        type(Image),            intent(in)      ::  other_image
        ! Private variables

        ! Start work

        if (.not. self%HasSameDimensionsAs(other_image)) then
            call this_program%TerminateWithFatalError('Image::MultiplyImage','Images do not have the same dimensions')
        endif
        !if (.not. self%IsInSameSpaceAs(other_image)) then
        !    call this_program%TerminateWithFatalError('Image::MultiplyImage','Images are not in the same space')
        !endif

        if (self%IsInRealSpace() .and. other_image%IsInRealSpace()) then
            self%real_values = self%real_values / other_image%real_values
        else if (self%IsInFourierSpace() .and. other_image%IsInFourierSpace()) then
            self%complex_values = self%complex_values / other_image%complex_values
        else if (self%IsInFourierSpace() .and. other_image%IsInRealSpace()) then
            self%complex_values = self%complex_values / other_image%real_values
        else
             call this_program%TerminateWithFatalError('Image::DivideImage','Images don not have compatible types')
        endif


    end subroutine DivideImage


    !>  \brief  Cut an image the size and shape of the current image, out of another image at the
    !!          co-ordinates specified.  the resulting cut image is stored in this image
    !!
    !!  \todo   This could be merged with or call ClipInto, I think
    subroutine CutFromOtherImage(self,other_image, x_position, y_position, z_position)
        ! Arguments
        class(Image),           intent(inout)   ::  self
        type(Image),            intent(in)      ::  other_image
        integer,                intent(in)      ::  x_position
        integer,                intent(in)      ::  y_position
        integer, optional,      intent(in)      ::  z_position

        !variables

        integer                                 ::  start_x
        integer                                 ::  start_y
        integer                                 ::  start_z

        integer                                 ::  y_counter
        integer                                 ::  z_counter
        integer                                 ::  x_counter

        !do some checks

        if (.not. self%IsInRealSpace() .or. .not. other_image%IsInRealSpace()) then
            call this_program%TerminateWithFatalError('Image::CutFromOtherImage', 'Images must be real' )
        endif

        start_x = x_position - self%logical_dimensions(1) / 2
        start_y = y_position - self%logical_dimensions(2) / 2

        if (present(z_position)) then
            start_z = z_position - self%logical_dimensions(3) / 2
        else
            start_z = 0
        endif

        do z_counter = 1, self%logical_dimensions(3)
            do y_counter = 1, self%logical_dimensions(2)
                do x_counter = 1, self%logical_dimensions(1)

                    if (start_x + x_counter .lt. 1 .or. &
                        start_x + x_counter .gt. other_image%logical_dimensions(1) .or. &
                        start_y + y_counter .lt. 1 .or. &
                        start_y + y_counter .gt. other_image%logical_dimensions(2) .or. &
                        start_z + z_counter .lt. 1 .or. &
                        start_z + z_counter .gt. other_image%logical_dimensions(3)) then
                        self%real_values(x_counter, y_counter, z_counter) = 0.0e0
                    else
                        self%real_values(x_counter, y_counter, z_counter) = &
                            other_image%real_values(start_x + x_counter, start_y + y_counter, start_z + z_counter)
                    endif

                enddo
            enddo
        enddo

    end subroutine CutFromOtherImage

    !>  \brief  Apply a median filter (naive algorithm)
    subroutine ApplyMedianFilter(self,output_image)
        !use napack_sort2, only : sort2
        use UsefulFunctions, only : QuickSort
        ! arguments
        class(Image),                       intent(in)      ::  self
        type(Image),                        intent(inout)   ::  output_image
        ! private variables
        integer,        parameter   ::  window_size = 3
        real(kind=8)                ::  window_values(window_size**2)
        !real(kind=8)                ::  window_values_sorted(window_size**2)
        !real(kind=8)                ::  work_array(window_size**2)
        integer                     ::  i,j
        ! start work

        if (self%IsAVolume()) call this_program%TerminateWithFatalError('Image::ApplyMedianFilter','Volumes not supported')

        output_image = self

        do j=2,self%GetLogicalDimension(2)-1
            do i=2,self%GetLogicalDimension(1)-1
                window_values(1:3) = self%real_values(i-1:i+1,j-1,1)
                window_values(4:6) = self%real_values(i-1:i+1,j  ,1)
                window_values(7:9) = self%real_values(i-1:i+1,j+1,1)
                !
                !call sort2(window_values,window_values_sorted,work_array,window_size**2)
                call QuickSort(window_values)
                !
                output_image%real_values(i,j,1) = window_values(5)
            enddo
        enddo
    end subroutine ApplyMedianFilter


    !>  \brief  When objects drift during image acquisition, information is lost in the direction of the drift.
    !!          This method applies the corresponding filter enveloppe to the image
    subroutine ApplyDriftFilter(self, drift_x, drift_y, drift_z_optional )
        ! arguments
        class(image),                       intent(inout)   ::  self
        real,                               intent(in)      ::  drift_x             !<  In pixels, how far the object moved during image acquisition
        real,                               intent(in)      ::  drift_y             !<  In pixels, how far the object moved during image acquisition
        real,                   optional,   intent(in)      ::  drift_z_optional    !<  In pixels, how far the object moved during image acquisition
        ! private variables
        integer ::  i,j,k
        real    ::  spatial_frequency(3)
        real    ::  drift_along_current_direction
        real    ::  attenuation
        real    ::  drift_z
        real,   parameter   ::  tim_fudge_factor = 1.0e0 ! Tim suggested using a fudge factor to amplify the effect of the drift filter. He suggested 3.
        ! start work

        !
        drift_z = 0.0e0
        if (present(drift_z_optional)) drift_z = drift_z_optional

        if (self%is_in_real_space) then
            do k= self%logical_lower_bound_real(3),self%logical_upper_bound_real(3)
                spatial_frequency(3) = k * self%fourier_voxel_size(3)
                do j=self%logical_lower_bound_real(2),self%logical_upper_bound_real(2)
                    spatial_frequency(2) = j * self%fourier_voxel_size(2)
                    do i=self%logical_lower_bound_real(1),self%logical_upper_bound_real(1)
                        spatial_frequency(1) = i *self%fourier_voxel_size(1)

                        ! Drift in radians (= Tv in my notes of 9-Jan-2014, 10-Jul-2013, 13-Oct-2014)
                        drift_along_current_direction = pi * (                       &
                                                        spatial_frequency(1) * drift_x * tim_fudge_factor +     &
                                                        spatial_frequency(2) * drift_y * tim_fudge_factor +     &
                                                        spatial_frequency(3) * drift_z * tim_fudge_factor      &
                                                                    )
                        ! The attenuation (see derivation in 9-Jan-2014 or 10-Jul-2013 notes) is a sinc function
                        ! See also the more elegant Ruben Diaz derivation, March 2015
                        if (abs(drift_along_current_direction) .gt. 0.00001e0) then
                            attenuation = sin(drift_along_current_direction) / drift_along_current_direction
                        else
                            attenuation = 1.0e0
                        endif
                        ! Apply the attenuation
                        self%real_values(i + self%physical_address_of_box_center(1), &
                                         j + self%physical_address_of_box_center(2), &
                                         k + self%physical_address_of_box_center(3)) = &
                            self%real_values(i + self%physical_address_of_box_center(1), &
                                             j + self%physical_address_of_box_center(2), &
                                             k + self%physical_address_of_box_center(3)) &
                                          * attenuation
                    enddo
                enddo
            enddo
        else
            do k=1,self%logical_dimensions(3)
                spatial_frequency(3) = self%LogicalIndexGivenPhysicalIndexInFourierSpace(k,3) * self%fourier_voxel_size(3)
                do j=1,self%logical_dimensions(2)
                    spatial_frequency(2) = self%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2) * self%fourier_voxel_size(2)
                    do i=1,self%physical_upper_bound_complex(1)
                        spatial_frequency(1) = (i-1)*self%fourier_voxel_size(1)
                        ! Drift in radians (= Tv in my notes of 9-Jan-2014 and 10-Jul-2013)
                        drift_along_current_direction = pi * (                       &
                                                        spatial_frequency(1) * drift_x * tim_fudge_factor +     &
                                                        spatial_frequency(2) * drift_y * tim_fudge_factor +     &
                                                        spatial_frequency(3) * drift_z * tim_fudge_factor       &
                                                                    )
                        ! The attenuation (see derivation in 9-Jan-2014 or 10-Jul-2013 notes)
                        !if (abs(drift_along_current_direction) .gt. 0.00001e0) then
                        !    attenuation = sqrt(   sin(drift_along_current_direction)**2         &
                        !                        + (1-cos(drift_along_current_direction))**2     &
                        !                        ) &
                        !                        / abs(drift_along_current_direction)
                        !else
                        !    attenuation = 1.0e0
                        !endif
                        ! See notes from Ruben Diaz, March 2015
                        if(abs(drift_along_current_direction) .gt. 0.00001e0) then
                            attenuation = sin(drift_along_current_direction) / drift_along_current_direction
                        else
                            attenuation = 1.0e0
                        endif

                        ! Apply the attenuation
                        self%complex_values(i,j,k) = self%complex_values(i,j,k) * attenuation
                    enddo
                enddo
            enddo
        endif
    end subroutine ApplyDriftFilter

    !>  \brief  When objects drift during image acquisition, information is lost in the direction of the drift.
    !!          This method applies the corresponding filter enveloppe to the image
    subroutine ComputeDriftFilter(self, drift_x, drift_y, drift_z )
        ! arguments
        class(image),                       intent(inout)   ::  self
        real,                               intent(in)      ::  drift_x             !<  In pixels, how far the object moved during image acquisition
        real,                               intent(in)      ::  drift_y             !<  In pixels, how far the object moved during image acquisition
        real,                   optional,   intent(in)      ::  drift_z             !<  In pixels, how far the object moved during image acquisition
        ! start work
        if (.not. self%IsInRealSpace()) then
            call this_program%TerminateWithFatalError('Image::ComputeDriftFilter', 'Image must be real' )
        endif

        self = 1.0e0
        call self%ApplyDriftFilter(drift_x,drift_y,drift_z)

    end subroutine ComputeDriftFilter

    !>  \brief  Subtract another image from this one
    subroutine SubtractImage(self,other_image,set_phases_to_zero)
        ! Arguments
        class(Image),           intent(inout)   ::  self
        type(Image),            intent(in)      ::  other_image
        logical,    optional,   intent(in)      ::  set_phases_to_zero
        ! Start work
        call self%AddImage(other_image,set_phases_to_zero,subtract=.true.)
    end subroutine SubtractImage


    !>  \brief  apply a bfactor to the image / volume, in place
    subroutine ApplyBFactor(    self,bfactor,   &
                                lp_cosine,lp_cosine_freq,lp_cosine_falloff, &
                                hp_cosine,hp_cosine_freq,hp_cosine_falloff, &
                                print_filter_shape,squared,filter_shape)
        ! arguments
        class(image),                       intent(inout)   ::  self
        real,                               intent(in)      ::  bfactor             !<  b factor, unitless (to do this, take the b factor in a^2 and divide by the square of the pixel size). using the crystallographic convention, a la henderson & grigorieff
        logical,                optional,   intent(in)      ::  lp_cosine           !<  whether also to apply a cosine-edge low-pass filter. false by default.
        real,                   optional,   intent(in)      ::  lp_cosine_freq      !<  spatial frequency (in 1/pixel; 0.5 is nyquist) at which the cosine falloff ends. in other words, there is a guarantee that the falloff is complete by that point.
        real,                   optional,   intent(in)      ::  lp_cosine_falloff   !<  width (in 1/pixels) of the cosine falloff
        logical,                optional,   intent(in)      ::  hp_cosine           !<  whether also to apply a cosine-edge high-pass filter. false by default.
        real,                   optional,   intent(in)      ::  hp_cosine_freq      !<  spatial frequency (in 1/pixel; 0.5 is nyquist) at which the cosine falloff ends. in other words, the filter reaches 1.0 by that point.
        real,                   optional,   intent(in)      ::  hp_cosine_falloff   !<  width (in 1/pixels) of the cosine falloff
        logical,                optional,   intent(in)      ::  print_filter_shape  !<  print the values of the filter to output. defaults to false.
        logical,                optional,   intent(in)      ::  squared             !<  whether the weighting should be squared (useful for e.g. equation 1 of grigorieff 1998)
        real,   allocatable,    optional,   intent(inout)   ::  filter_shape(:)     !<  if present, the filter values will be returned in this array
        ! private variables
        integer ::  i,j,k
        logical ::  llp_cosine
        real    ::  llp_cosine_freq, llp_cosine_freq_min, llp_cosine_falloff, lp_cosine_attenuation
        logical ::  hhp_cosine
        real    ::  hhp_cosine_freq, hhp_cosine_freq_min, hhp_cosine_falloff, hp_cosine_attenuation
        logical ::  pprint_filter_shape
        real    ::  squaring_factor
        logical,    parameter   ::  debug   =   .false.
        real    ::  filter_value
        integer ::  box_center(3)
        real    ::  coos(3), freq_sq
        real    ::  freq
        real,   allocatable ::  ffilter_shape(:)
        ! start work

        ! will we be squaring?
        squaring_factor = 1.0
        if (present(squared)) then
            if (squared) then
                squaring_factor = 2.0
            endif
        endif
        if (debug) write(*,'(a,f0.1)') '**debug(ApplyBFactor): squaring factor = ', squaring_factor

        ! will we be doing a cosine-edge low-pass filter too?
        llp_cosine = .false.
        if (present(lp_cosine)) llp_cosine = lp_cosine
        hhp_cosine = .false.
        if (present(hp_cosine)) hhp_cosine = hp_cosine


        ! if we are doing a cosine-edge low-pass filter, check we have all the information we need
        if (llp_cosine) then
            if (present(lp_cosine_freq) .and. present(lp_cosine_falloff)) then
                llp_cosine_freq = lp_cosine_freq
                llp_cosine_falloff = lp_cosine_falloff
            else
                call this_program%TerminateWithFatalError('Image::ApplyBFactor', &
                                    'You asked for a cosine-edged low-pass filter but did not supply the required arguments')
            endif
        else
            llp_cosine_freq = maxval(self%logical_dimensions(1:3)) + 2.
            llp_cosine_falloff = 1.
        endif

        ! if we are doing a cosine-edge low-pass filter, check we have all the information we need
        if (hhp_cosine) then
            if (present(hp_cosine_freq) .and. present(hp_cosine_falloff)) then
                hhp_cosine_freq = hp_cosine_freq
                hhp_cosine_falloff = hp_cosine_falloff
            else
                call this_program%TerminateWithFatalError('Image::ApplyBFactor', &
                                    'You asked for a cosine-edged low-pass filter but did not supply the required arguments')
            endif
        else
            hhp_cosine_freq = 0.0
            hhp_cosine_falloff = 1.
        endif

        ! will we be printing out the shape of the filter?
        if (present(print_filter_shape)) then
            pprint_filter_shape = print_filter_shape
        else
            pprint_filter_shape = .false.
        endif

        ! check the image is in memory
        if (.not. self%is_in_memory) then
            call this_program%TerminateWithFatalError('Image::ApplyBFactor', 'Image is not allocated')
        endif

        ! the square of the dimension of the image, times 4 because we follow the crystallography / henderson / grigorieff convention
!       nsam_sq = (self%logical_dimensions(1)*self%logical_dimensions(2)) * 4

        ! work out parameters for the cosine filters
        llp_cosine_freq_min =   max(llp_cosine_freq - llp_cosine_falloff,0.0)
        hhp_cosine_freq_min =   max(hhp_cosine_freq - hhp_cosine_falloff,0.0)

        if (debug) then
            write(*,'(a,3(f0.3,1x))') '**debug(ApplyBFactor): lp: ', llp_cosine_freq_min, llp_cosine_freq, llp_cosine_falloff
            write(*,'(a,3(f0.3,1x))') '**debug(ApplyBFactor): hp: ', hhp_cosine_freq_min, hhp_cosine_freq, hhp_cosine_falloff
        endif

        ! origin
        box_center = self%physical_address_of_box_center
        if (debug) write(*,'(a,3(i0,1x))') '**debug(ApplyBFactor): origin = ', box_center

        ! allocate filter_shape
        allocate(ffilter_shape(self%physical_address_of_box_center(1)))


        if (debug) write(*,'(a)',advance='no') '**debug(ApplyBFactor): bfactor (not squared) = '
        if (self%is_in_real_space) then
            if (self%object_is_centered_in_box) then

                do k=1,self%logical_dimensions(3)
                    coos(3) = ((k - box_center(3)) * self%fourier_voxel_size(3))**2
                    do j=1,self%logical_dimensions(2)
                        coos(2) = ((j - box_center(2)) * self%fourier_voxel_size(2))**2
                        do i=1,self%logical_dimensions(1)
                            coos(1) = ((i - box_center(1)) * self%fourier_voxel_size(1))**2
                            ! compute squared radius, in units of reciprocal pixels
                            freq_sq = sum(coos)
                            freq = sqrt(freq_sq)
                            ! apply the low-pass filter
                            if (freq .lt. llp_cosine_freq_min) then
                                lp_cosine_attenuation = 1.0
                            else if (freq .lt. llp_cosine_freq) then
                                lp_cosine_attenuation = (cos(((freq-llp_cosine_freq_min)/llp_cosine_falloff)*pi)+1.)*0.5e0
                            else
                                lp_cosine_attenuation = 0.0
                            endif
                            ! apply the high-pass filter
                            if (freq .lt. hhp_cosine_freq_min) then
                                hp_cosine_attenuation = 0.0
                            else if (freq .lt. hhp_cosine_freq) then
                                hp_cosine_attenuation = (sin(((freq-hhp_cosine_freq_min)/hhp_cosine_falloff)*pi)+1.)*0.5e0
                            else
                                hp_cosine_attenuation = 1.0
                            endif
                            filter_value =      exp(-bfactor * squaring_factor * freq_sq * 0.25e0)  &
                                            *   lp_cosine_attenuation * hp_cosine_attenuation! divided by 4.0 because we follow the crystallography / henderson / grigorieff convention
                            self%real_values(i,j,k) = self%real_values(i,j,k) * filter_value
                            if (j .eq. box_center(2) .and. k .eq. box_center(3) .and. i .ge. box_center(1)) then
                                if (debug) write(*,'(g15.3,1x)',advance='no') filter_value
                                ffilter_shape(i-box_center(1)+1) = filter_value
                            endif
                        enddo
                    enddo
                enddo

            else

                do k=1,self%logical_dimensions(3)
                    ! the distance from the origin to the current pixel, in reciprocal pixel units
                    coos(3) = (self%LogicalIndexGivenPhysicalIndexInFourierSpace(k,3) * self%fourier_voxel_size(3))**2
                    do j=1,self%logical_dimensions(2)
                        coos(2) = (self%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2) * self%fourier_voxel_size(2))**2
                        do i=1,self%logical_dimensions(1)
                            coos(1) = ((i-1) * self%fourier_voxel_size(1))**2

                            ! compute squared radius, in units of reciprocal pixels
                            freq_sq = sum(coos)
                            freq = sqrt(freq_sq)
                            ! apply the low-pass filter
                            if (freq .lt. llp_cosine_freq_min) then
                                lp_cosine_attenuation = 1.0
                            else if (freq .lt. llp_cosine_freq) then
                                lp_cosine_attenuation = (cos(((freq-llp_cosine_freq_min)/llp_cosine_falloff)*pi)+1.)*0.5e0
                            else
                                lp_cosine_attenuation = 0.0
                            endif
                            ! apply the high-pass filter
                            if (freq .lt. hhp_cosine_freq_min) then
                                hp_cosine_attenuation = 0.0
                            else if (freq .lt. hhp_cosine_freq) then
                                hp_cosine_attenuation = (sin(((freq-hhp_cosine_freq_min)/hhp_cosine_falloff)*pi)+1.)*0.5e0
                            else
                                hp_cosine_attenuation = 1.0
                            endif
                            filter_value =      exp(-bfactor * squaring_factor * freq_sq * 0.25e0)  &
                                            *   lp_cosine_attenuation * hp_cosine_attenuation ! divide by 4.0 because we follow the crystallography / henderson / grigorieff convention
                            self%real_values(i,j,k) = self%real_values(i,j,k) * filter_value
                            if (j .eq. 1 .and. k .eq. 1) then
                                if (debug) write(*,'(g15.3,1x)',advance='no') filter_value
                                ffilter_shape(i-box_center(1)+1) = filter_value
                            endif
                        enddo
                    enddo
                enddo

            endif
        else
            ! image is in Fourier space

            do k=1,self%logical_dimensions(3)
                ! the distance from the origin to the current pixel, in reciprocal pixel units
                coos(3) = (self%LogicalIndexGivenPhysicalIndexInFourierSpace(k,3) * self%fourier_voxel_size(3))**2
                do j=1,self%logical_dimensions(2)
                    coos(2) = (self%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2) * self%fourier_voxel_size(2))**2
                    do i=1,self%physical_upper_bound_complex(1)
                        coos(1) = ((i-1) * self%fourier_voxel_size(1))**2
                        ! compute squared radius, in units of reciprocal pixels
                        freq_sq = sum(coos)
                        freq = sqrt(freq_sq)
                        ! apply the low-pass filter
                        if (freq .lt. llp_cosine_freq_min) then
                            lp_cosine_attenuation = 1.0
                        else if (freq .lt. llp_cosine_freq) then
                            lp_cosine_attenuation = (cos(((freq-llp_cosine_freq_min)/llp_cosine_falloff)*pi)+1.)*0.5e0
                        else
                            lp_cosine_attenuation = 0.0
                        endif
                        ! apply the high-pass filter
                        if (freq .lt. hhp_cosine_freq_min) then
                            hp_cosine_attenuation = 0.0
                        else if (freq .lt. hhp_cosine_freq) then
                            hp_cosine_attenuation = (sin(((freq-hhp_cosine_freq_min)/hhp_cosine_falloff)*pi)+1.)*0.5e0
                        else
                            hp_cosine_attenuation = 1.0
                        endif
                        filter_value =      exp(-bfactor * squaring_factor * freq_sq * 0.25e0)  &
                                        *   lp_cosine_attenuation * hp_cosine_attenuation! 4.0 because we follow the crystallography / henderson / grigorieff convention
                        self%complex_values(i,j,k) = self%complex_values(i,j,k) * filter_value
                        if (j .eq.1 .and. k .eq. 1) then
                            if (debug) write(*,'(g15.3,1x)',advance='no') filter_value
                            ffilter_shape(i) = filter_value
                        endif
                    enddo
                enddo
            enddo

        endif
        if (debug) write(*,'(a)') '.'

        ! print out shape of the filter if user asked for it
        if (pprint_filter_shape) then
            write(*,'(/a)') 'weight of combined b-factor & low-pass filters'
            write(*,'(a16,5x,a10)') 'freq (1/pixel)', 'weight'
            do i=1,size(ffilter_shape)
                freq = real(i-1) / real(self%logical_dimensions(2))
                write(*,'(f16.3,5x,f10.4)') freq, ffilter_shape(i)
            enddo
        endif

        !
        if (present(filter_shape)) filter_shape = ffilter_shape

    end subroutine ApplyBFactor


    !>  \brief  apply a cosine edged band-pass filter to the image, and return it in the space it was given.
    !!          note: set parameters to 0.0 to get a low-pass- or high-pass-only filter.
    !!  \todo   Consider removing the automatic FFT
    !!  \todo   Optimise: work with squared frequency inside the loop?
#ifdef skip_runtime_checks
    pure &
#endif
    subroutine ApplyBandPassFilter(self,hp_freq,lp_freq,width_of_cosine_falloff)
        class(image),       intent(inout)   ::  self                    !<  Input image. can be real or complex, either way filter is applied in reciprocal space and image is returned in same space as given.
        real,               intent(in)      ::  hp_freq                 !<  Frequencies at this radius and below will be 0.0 in output. In reciprocal pixel units (Nyquist is 0.5). Ignored if set to 0.0.
        real,               intent(in)      ::  lp_freq                 !<  Frequencies at this radius and above will be 0.0 in output. In reciprocal pixel units (Nyquist is 0.5). Ignored if set to 0.0.
        real,   optional,   intent(in)      ::  width_of_cosine_falloff !<  In reciprocal pixel units (frequencies at hp_freq + cosine_falloff will be left intact)
        ! private variables
        logical                 ::  do_ft
        real                    ::  ccosine_falloff
        !logical,    parameter   ::  debug   =   .false.
        integer                 ::  i,j,k
        real                    ::  freq, coos(3), falloff_factor
        ! start work

#ifndef skip_runtime_checks
        ! check filter parameters
        if (hp_freq .ge. 0.5) then
            write(*,'(a,f5.1,a)') '**error(ApplyBandPassFilter): high-pass frequency is higher than nyquist (', &
                                    hp_freq, '). if you meant to do a low-pass-only filter, set high-pass parameter to 0.0'
            call this_program%TerminateWithFatalError('Image::ApplyBandPassFilter','Bad filter parameter')
        endif
        if (lp_freq .lt. 0.0 .or. hp_freq .lt. 0.0) then
            write(*,'(a,2(f0.4,1x))') '**error(ApplyBandPassFilter): negative filter parameters do not make sense: ', &
                                        lp_freq, hp_freq
            call this_program%TerminateWithFatalError('Image::ApplyBandPassFilter','Bad filter parameter')
        endif
        if (hp_freq .gt. lp_freq .and. hp_freq .ne. 0.0 .and. lp_freq .ne. 0.0) then
            write(*,'(a,2(f0.4,1x))') '**error(ApplyBandPassFilter): high-pass frequency cannot be '//&
                                        'lower than low-pass frequency: ', lp_freq, hp_freq
            call this_program%TerminateWithFatalError('Image::ApplyBandPassFilter','Bad filter parameter')
        endif
#endif

        ! cosine falloff is 0.05 unless otherwise specified
        ccosine_falloff = 0.05
        if (present(width_of_cosine_falloff)) ccosine_falloff = width_of_cosine_falloff

        ! if necessary, ft the input image
        do_ft = self%IsInRealSpace()
        if (do_ft) call self%ForwardFFT()

        ! Apply the filter, loop over image
        do k=1,self%logical_dimensions(3)
            ! the distance from the origin to the current pixel, in reciprocal pixel units
            coos(3) = (real(self%LogicalIndexGivenPhysicalIndexInFourierSpace(k,3))*self%fourier_voxel_size(3))**2
            do j=1,self%logical_dimensions(2)
                coos(2) = (real(self%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2))*self%fourier_voxel_size(2))**2
                do i=1,self%physical_upper_bound_complex(1)
                    coos(1) = (real((i-1))*self%fourier_voxel_size(1))**2
                    ! compute radius, in units of reciprocal pixels
                    freq = sqrt(sum(coos))
                    ! apply the high-pass
                    if (hp_freq .ne. 0.0) then
                        if (freq .le. hp_freq) then
                            self%complex_values(i,j,k) = (0.0,0.0)
                        else if (freq .le. hp_freq + ccosine_falloff) then
                            falloff_factor = ((freq - hp_freq)/ccosine_falloff)*pi  ! from 0 to pi
                            falloff_factor = cos(falloff_factor) ! from 1 to -1
                            falloff_factor = -falloff_factor ! from -1 to 1
                            falloff_factor = (falloff_factor + 1.0e0)*0.5e0  ! from 0 to 1
                            self%complex_values(i,j,k) = self%complex_values(i,j,k) * falloff_factor
                        endif
                    endif
                    ! apply the low-pass
                    if (lp_freq .ne. 0.0) then
                        if (freq .ge. lp_freq) then
                            self%complex_values(i,j,k) = 0.0
                        else if (freq .ge. lp_freq - ccosine_falloff) then
                            falloff_factor = ((freq - (lp_freq - ccosine_falloff))/ccosine_falloff)*pi ! from 0 to pi
                            falloff_factor = cos(falloff_factor) ! from 1 to -1
                            falloff_factor = (falloff_factor + 1.0e0)*0.5e0 ! from 1 to 0
                            self%complex_values(i,j,k) = self%complex_values(i,j,k) * falloff_factor
                        endif
                    endif
                enddo
            enddo
        enddo

        ! if the image was given as a real image, return it as a real image
        if (do_ft) call self%BackwardFFT()
    end subroutine ApplyBandPassFilter

    !>  \brief  In real space, replace image density values with values of planar gradient across the image
    subroutine ReplaceWithLinearGradient(self)
        class(Image),               intent(inout)   ::  self
        !
        call self%RemoveLinearGradient(replace_with_gradient=.true.)
    end subroutine ReplaceWithLinearGradient

    !>  \brief  In real space, remove any linear gradient across the image
    subroutine RemoveLinearGradient(self,replace_with_gradient)
        class(Image),               intent(inout)   ::  self
        logical,        optional,   intent(in)      ::  replace_with_gradient
        ! private variable
        real(kind=8)    ::  average_at_first_pixel(3)
        real(kind=8)    ::  average_at_last_pixel(3)
        real            ::  linear_gradient_value(3)
        real            ::  average_value
        integer         ::  i,j,k
        logical         ::  rreplace_with_gradient
        real            ::  inv_num_dimensions
        ! start work

        ! optionally, we could replace the image with its gradient rather than remove the gradient
        rreplace_with_gradient = .false.
        if (present(replace_with_gradient)) rreplace_with_gradient = replace_with_gradient

        !
        average_at_first_pixel(1) = sum(self%real_values(1                           ,:,:)) &
                                    /(self%logical_dimensions(2)*self%logical_dimensions(3))
        average_at_first_pixel(2) = sum(self%real_values(1:self%logical_dimensions(1),1,:)) &
                                    /(self%logical_dimensions(1)*self%logical_dimensions(3))
        average_at_first_pixel(3) = sum(self%real_values(1:self%logical_dimensions(1),:,1)) &
                                    /(self%logical_dimensions(1)*self%logical_dimensions(2))

        average_at_last_pixel(1)  = sum(self%real_values(self%logical_dimensions(1),:,:))&
                                    /(self%logical_dimensions(2)*self%logical_dimensions(3))
        average_at_last_pixel(2)  = sum(self%real_values(1:self%logical_dimensions(1),self%logical_dimensions(2),:))&
                                    /(self%logical_dimensions(1)*self%logical_dimensions(3))
        average_at_last_pixel(3)  = sum(self%real_values(1:self%logical_dimensions(1),:,self%logical_dimensions(3)))&
                                    /(self%logical_dimensions(1)*self%logical_dimensions(2))

        average_value = self%GetAverageOfValues()

        ! this line is not done in CTFFIND's BOXIMG, which I think is a harmless bug
        average_value = average_value * count(self%logical_dimensions .gt. 1)
        inv_num_dimensions = 1.0 / count(self%logical_dimensions .gt. 1)

        ! Remove gradient
        do k=1,self%logical_dimensions(3)
            if (self%logical_dimensions(3) .gt. 1) then
                linear_gradient_value(3) = average_at_first_pixel(3)+&
                                          (average_at_last_pixel(3)-average_at_first_pixel(3))/(self%logical_dimensions(3)-1)*(k-1)
            else
                linear_gradient_value(3) = 0.0
            endif
            do j=1,self%logical_dimensions(2)
                linear_gradient_value(2) = average_at_first_pixel(2)+(average_at_last_pixel(2)-average_at_first_pixel(2))&
                                                                        /(self%logical_dimensions(2)-1)*(j-1)
                do i=1,self%logical_dimensions(1)
                    linear_gradient_value(1) = average_at_first_pixel(1)+(average_at_last_pixel(1)-average_at_first_pixel(1))&
                                                                            /(self%logical_dimensions(1)-1)*(i-1)
                    !
                    if (rreplace_with_gradient) then
                        self%real_values(i,j,k) = sum(linear_gradient_value) * inv_num_dimensions
                    else
                        self%real_values(i,j,k) = self%real_values(i,j,k) &
                            !- sum(linear_gradient_value) + average_value
                            - sum(linear_gradient_value) * inv_num_dimensions
                    endif
                enddo
            enddo
        enddo

    end subroutine RemoveLinearGradient

    !>  \brief  Real-space box convolution meant for 2D power spectra
    !!
    !!  This is adapted from the MSMOOTH routine from CTFFIND, with a different wrap-around behaviour
    subroutine SpectrumBoxConvolution(self,box_size,output_image,minimum_radius)
        use UsefulFunctions, only : IsOdd
        class(Image),   intent(in)      ::  self
        integer,        intent(in)      ::  box_size            !<  Must be odd (pixels). Edge dimension of square convolution kernel.
        type(Image),    intent(inout)   ::  output_image
        real,           intent(in)      ::  minimum_radius      !<  Below this radius (pixels), the input image is left unchanged
        ! private variables
        integer     ::  num_voxels
        integer     ::  half_box_size
        integer     ::  i,j
        integer     ::  i_sq, j_sq
        integer     ::  ii_sq,jj_sq
        real        ::  radius
        integer     ::  ii,jj
        integer     ::  l,m
        !integer     ::  i_logi,j_logi
        integer     ::  i_friedel,j_friedel
        integer     ::  ii_friedel, jj_friedel
        integer     ::  iii,jjj
        real        ::  average_just_outside_minimum_radius
        integer     ::  cross_half_width_to_ignore
        ! start work

        if (.not. IsOdd(box_size)) call this_program%TerminateWithFatalError('Image::BoxConvolution','box_size must be odd')
        half_box_size = (box_size-1)/8
        if (self%IsAVolume()) call this_program%TerminateWithFatalError('Image::BoxConvolution','2D images only')

        ! allocate output
        call output_image%Allocate(self)

        ! Cross half width
        cross_half_width_to_ignore = max(min(6,int(self%GetLogicalDimension(1)/80.0)),2)

        !
        average_just_outside_minimum_radius = self%GetAverageOfValuesInRing(minimum_radius,minimum_radius+1.0)

        !$omp   parallel default(shared) &
        !$omp   private(i,j,ii,jj,m,l,i_friedel,j_friedel,num_voxels,radius,i_sq,j_sq,ii_friedel,jj_friedel,iii,jjj,ii_sq,jj_sq)

        ! loop over output image. To save time, we only loop over one half of the image.
        !$omp   do
        do j=1,self%logical_dimensions(2)
            j_friedel=2*self%physical_address_of_box_center(2)-j
            j_sq = (j-self%physical_address_of_box_center(2))**2
            do i=1,self%physical_address_of_box_center(1)
                i_friedel=2*self%physical_address_of_box_center(1)-i
                i_sq = (i-self%physical_address_of_box_center(1))**2
                !
                radius = sqrt(real(i_sq+j_sq))
                !
                if (radius .le. minimum_radius) then
                    output_image%real_values(i,j,1) = self%real_values(i,j,1)
                else
                    output_image%real_values(i,j,1) = 0.0e0
                    num_voxels = 0
                    do m=-half_box_size,half_box_size
                        jj = j+m
                        if (jj .lt. 1) jj = jj + self%logical_dimensions(2)
                        if (jj .gt. self%logical_dimensions(2)) jj = jj - self%logical_dimensions(2)
                        jj_friedel = 2*self%physical_address_of_box_center(2)-jj
                        jj_sq = (jj-self%physical_address_of_box_center(2))**2
                        do l=-half_box_size,half_box_size
                            ii = i+l
                            if (ii .lt. 1) ii = ii + self%logical_dimensions(1)
                            if (ii .gt. self%logical_dimensions(1)) ii = ii - self%logical_dimensions(1)
                            ii_friedel = 2*self%physical_address_of_box_center(1)-ii
                            ii_sq = (ii-self%physical_address_of_box_center(1))**2
                            !
                            radius = sqrt(real(ii_sq+jj_sq))
                            !
                            !if (radius .ge. minimum_radius) then
                                !
                                if (ii .gt. self%physical_address_of_box_center(1)) then
                                    iii = ii_friedel
                                    jjj = jj_friedel
                                    if (jjj .gt. self%logical_dimensions(2) .or. iii .gt. self%logical_dimensions(1)) cycle
                                else
                                    iii = ii
                                    jjj = jj
                                endif
                                if (abs(iii - self%physical_address_of_box_center(1)) .le. cross_half_width_to_ignore .or. &
                                    abs(jjj - self%physical_address_of_box_center(2)) .le. cross_half_width_to_ignore) then
                                    cycle
                                endif
                                if (iii .le. self%logical_dimensions(1) .and. jjj .le. self%logical_dimensions(2)) then ! it sometimes happens that we end up on Nyquist Friedel mates that we don't have (perhaps this can be fixed)
                                    output_image%real_values(i,j,1) = output_image%real_values(i,j,1) &
                                                                        + self%real_values(iii,jjj,1)
                                endif
                                num_voxels = num_voxels + 1
                            !endif
                        enddo
                    enddo
                    if (num_voxels .eq. 0) then
                        output_image%real_values(i,j,1) = self%real_values(i,j,1)
                    else
                        output_image%real_values(i,j,1) = output_image%real_values(i,j,1) / real(num_voxels)
                    endif
                endif
                if (j_friedel .lt. output_image%logical_dimensions(2)) then
                    output_image%real_values(i_friedel,j_friedel,1) = output_image%real_values(i,j,1)
                endif
            enddo
        enddo
        !$omp   end do


        !$omp   end parallel

        ! There are a few pixels that are not set by the logic above
        output_image%real_values(output_image%physical_address_of_box_center(1)+1:output_image%logical_dimensions(1),1,1) = &
        output_image%real_values(output_image%physical_address_of_box_center(1)-1:2:-1,1,1)
        !
        output_image%real_values(output_image%physical_address_of_box_center(1)+1:output_image%logical_dimensions(1), &
                                 output_image%logical_dimensions(2),1) = &
        output_image%real_values(output_image%physical_address_of_box_center(1)-1:2:-1, &
                                 output_image%logical_dimensions(2),1)

    end subroutine SpectrumBoxConvolution


    !>  \brief  Impose a sharp circular mask on image. Pixels at a distance .GE. rad from the origin will be set to 0.0
#ifdef skip_runtime_checks
    pure &
#endif
    subroutine ApplyCircularMask(self, rad,inverse)
        class(Image),           intent(inout)   ::  self
        real,                   intent(in)      ::  rad !< Any voxels at this distance or greater from origin will be zeroed. In pixels (for real images) or reciprocal pixels (for complex images; so that 0.5 is Nyquist)
        logical,    optional,   intent(in)  ::  inverse !<  Invert the mask (center is masked out
        ! Private variables
        integer ::  i,j,k
        real    ::  coos(3)
        real    ::  freq
        integer ::  origin(3)
        integer ::  x,y,z
        real    ::  rad_tmp, rad_sq
        logical :: iinverse
        ! Start work
#ifndef skip_runtime_checks
        if (.not. self%IsInRealSpace()) then
            if (rad .gt. 0.5) then
                write(*,'(a,F0.3)') '**ERROR(circ_mask): masking complex image with radius > 0.5, i.e. beyond Nyquist: ', rad
                call this_program%TerminateWithFatalError('Image::ApplyCircularMask','Feature not fully tested/implemented')  ! Remove this eventually!
            endif
        endif
#endif

        iinverse = .false.
        if (present(inverse)) iinverse = inverse

        if (self%IsInRealSpace() .and. self%ObjectIsCenteredInBox()) then
            origin = self%GetPhysicalAddressOfBoxCenter()
        else
            origin = [1,1,1]
        endif

        rad_sq = rad**2

        if (self%IsInRealSpace()) then
            do k=1,self%GetLogicalDimension(3)
                z = (k-origin(3))**2
                do j=1,self%GetLogicalDimension(2)
                    y = (j-origin(2))**2
                    do i=1,self%GetLogicalDimension(1)
                        x = (i-origin(1))**2
                        rad_tmp = real(x+y+z)
                        if (iinverse) then
                            if (rad_tmp .le. rad_sq) then
                                self%real_values(i,j,k) = 0.0
                            endif
                        else
                            if (rad_tmp .ge. rad_sq) then
                                self%real_values(i,j,k) = 0.0
                            endif
                        endif
                    enddo
                enddo
            enddo
        else
            do k=1,self%GetLogicalDimension(3)
                ! The distance from the origin to the current pixel, in reciprocal pixel units
                coos(3) = (self%LogicalIndexGivenPhysicalIndexInFourierSpace(k,3)*self%GetFourierVoxelSize(3))**2
                do j=1,self%GetLogicalDimension(2)
                    coos(2) = (self%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2)*self%GetFourierVoxelSize(2))**2
                    do i=1,self%physical_upper_bound_complex(1)
                        coos(1) = (real(i-1)*self%GetFourierVoxelSize(1))**2
                        ! Compute frequency, in units of reciprocal pixels
                        freq = sqrt(sum(coos))
                        if (iinverse) then
                            if (freq .le. rad) then
                                self%complex_values(i,j,k) = (0.0,0.0)
                            endif
                        else
                            if (freq .ge. rad) then
                                self%complex_values(i,j,k) = (0.0,0.0)
                            endif
                        endif
                    enddo
                enddo
            enddo
        endif

    end subroutine ApplyCircularMask


    !>  \brief  Apply a cosine-edge circular mask to an image in real space. Mask has 1s near the center, and 0s near the edge
#ifdef skip_runtime_checks
    pure &
#endif
    subroutine ApplySoftCircularMask(self,radius_after_falloff,falloff,wipe_input,mask_with_average,cylinder,inverse)
        ! arguments
        class(Image),           intent(inout)   ::  self
        real,                   intent(in)      ::  radius_after_falloff    !<  radius at which the falloff reaches 0.0 (in pixels, or fraction of max radius)
        real,                   intent(in)      ::  falloff                 !<  distance from the start of the falloff (1.0) to the end of the falloff(0.0). the mask reaches 0.0 at radius_after_falloff. in pixels or fraction.
        logical,    optional,   intent(in)      ::  wipe_input              !<  whether to wipe / initialise the input to 1.0. .false. by default
        logical,    optional,   intent(in)      ::  mask_with_average       !<  if .true., the mask will go to the average value within the mask, rather than 0.0. if the mask is inverted, this refers to the average value in the area that will be masked out.
        logical,    optional,   intent(in)      ::  cylinder                !<  cosine-edge cylinder rather than sphere
        logical,    optional,   intent(in)      ::  inverse                 !<  invert the filter.
        ! private variables
        integer         ::  i,j,k
        integer         ::  smallest_dim
        real            ::  dist, radius_before_falloff, rradius_after_falloff
        logical         ::  wwipe_input, mmask_with_average
        logical         ::  radius_given_as_fraction
        real(kind=8)    ::  sum_in_mask
        real(kind=8)    ::  count_in_mask
        real            ::  mask_value
        real            ::  falloff_factor
        integer         ::  origin(3)
        real            ::  x,y,z
        integer         ::  ix,iy,iz
        real            ::  inv_max_rad_sq(3)
        logical         ::  ccylinder
        logical         ::  iinverse
        ! start work

#ifndef skip_runtime_checks
        ! built-in assumption the image is square/cubic
        !> \todo Get rid of square/cubic assumption
        if (.not. self%IsSquare()) then
            !call this_program%TerminateWithFatalError('Image::ApplySoftCircularMask','Input image is not square or cubic')
        endif
#endif

        ! wipe input?
        wwipe_input = .false.
        if (present(wipe_input)) wwipe_input = wipe_input

        ! mask with average?
        mmask_with_average = .false.
        if (present(mask_with_average)) mmask_with_average = mask_with_average

        ! cylinder?
        ccylinder = .false.
        if (present(cylinder)) ccylinder = cylinder

        ! invert?
        iinverse = .false.
        if (present(inverse)) iinverse = inverse

#ifndef skip_runtime_checks
        ! mask with average is not compatible with wipe input
        if (mmask_with_average .and. wwipe_input) then
            call this_program%TerminateWithFatalError('Image::ApplySoftCircularMask','Cannot mask with average and wipe input')
        endif

        ! check type
        !>  \todo Make this work in Fourier space so it can be used to filter images
        if (.not. self%is_in_real_space) then
            call this_program%TerminateWithFatalError('Image::ApplySoftCircularMask','Image must be in real space')
        endif

        ! check the radius makes sense
        if (self%IsAVolume()) then
            smallest_dim = minval(self%logical_dimensions(1:3))
        else
            smallest_dim = minval(self%logical_dimensions(1:2))
        endif

        if (radius_after_falloff .le. 0.0) then
            write(*,'(a,f12.3)') '**error(ApplySoftCircularMask): radius is too small: ', radius_after_falloff
            call this_program%TerminateWithFatalError('Image::ApplySoftCircularMask','Radius is too small')
        else if (radius_after_falloff .gt. smallest_dim) then
            write(*,'(a,f12.3,1x,3(i0,1x))') '**error(ApplySoftCircularMask): radius is too large: ', &
                                        radius_after_falloff, self%logical_dimensions
            call this_program%TerminateWithFatalError('Image::ApplySoftCircularMask','Radius is too large')
        endif
#endif

        ! are we dealing with a fraction? if so and the image isn't square/cubic, we'll get an elipse rather than disc
        if (radius_after_falloff .gt. 0.0 .and. radius_after_falloff .le. 1.0) then
#ifndef skip_runtime_checks
            if (falloff .gt. 1.0) then
                write(*,'(a,f0.1,a,f8.2,a)')    '**error(ApplySoftCircularMask): radius was specified as a fraction (', &
                                                radius_after_falloff,     &
                                                '), but the falloff was not (', falloff, ')'
                call this_program%TerminateWithFatalError(  'Image::ApplySoftCircularMask',&
                                                            'Radius specified as fraction, but falloff was not')
            endif
#endif
            radius_given_as_fraction = .true.
        else
            radius_given_as_fraction = .false.
        endif

        rradius_after_falloff = radius_after_falloff
        radius_before_falloff = rradius_after_falloff-falloff

        ! if required, wipe the input
        if (wwipe_input) then
            self%real_values = 1.0
        endif

        ! origin
        origin = self%physical_address_of_box_center

        ! maximum radius (squared) in each dimension
        if (radius_given_as_fraction) then
            if (self%object_is_centered_in_box .and. self%is_in_real_space) then
                if (self%IsAVolume()) then
                    inv_max_rad_sq = 1.0/(origin-1)**2
                else
                    inv_max_rad_sq(1:2) = 1.0/((origin(1:2)-1)**2)
                    inv_max_rad_sq(3) = 0.0
                endif
            else
                inv_max_rad_sq = 1.0/real(int(self%logical_dimensions/2))**2
            endif
        endif

        if (mmask_with_average) then
            ! loop over pixels to determine average value within mask
            sum_in_mask = 0.0d0
            count_in_mask = 0
            if (radius_given_as_fraction) then
                do k=1,self%logical_dimensions(3)
                    z = (k - origin(3))**2 * inv_max_rad_sq(3)
                    do j=1,self%logical_dimensions(2)
                        y = (j - origin(2))**2 * inv_max_rad_sq(2)
                        do i=1,self%logical_dimensions(1)
                            x = (i - origin(1))**2 * inv_max_rad_sq(1)
                            ! work out distance from origin (fraction of max radius)
                            if (ccylinder) then
                                dist = sqrt(x+y)
                            else
                                dist = sqrt(x+y+z)
                            endif
                            if (dist .ge. radius_before_falloff .and. dist .le. rradius_after_falloff) then
                                falloff_factor = (cos(((dist-radius_before_falloff)/falloff)*pi)+1.)/2.0
                                sum_in_mask = sum_in_mask + dble(self%real_values(i,j,k) * falloff_factor)
                                count_in_mask = count_in_mask + dble(falloff_factor)
                            else if (dist .lt. radius_before_falloff) then
                                sum_in_mask = sum_in_mask + dble(self%real_values(i,j,k))
                                count_in_mask = count_in_mask + 1.0d0
                            endif
                        enddo
                    enddo
                enddo
            else
                do k=1,self%logical_dimensions(3)
                    iz = (k - origin(3))**2
                    do j=1,self%logical_dimensions(2)
                        iy = (j - origin(2))**2
                        do i=1,self%logical_dimensions(1)
                            ix = (i - origin(1))**2
                            ! work out distance from origin
                            if (ccylinder) then
                                dist = sqrt(real(ix+iy))
                            else
                                dist = sqrt(real(ix+iy+iz))
                            endif
                            ! accumulate statistics
                            if (dist .ge. radius_before_falloff .and. dist .le. rradius_after_falloff) then
                                falloff_factor = (cos(((dist-radius_before_falloff)/falloff)*pi)+1.)/2.0
                                sum_in_mask = sum_in_mask + dble(self%real_values(i,j,k) * falloff_factor)
                                count_in_mask = count_in_mask + dble(falloff_factor)
                            else if (dist .lt. radius_before_falloff) then
                                sum_in_mask = sum_in_mask + dble(self%real_values(i,j,k))
                                count_in_mask = count_in_mask + 1.0d0
                            endif
                        enddo
                    enddo
                enddo
            endif   ! end of test for radius_given_as_fraction
            mask_value = real(sum_in_mask / count_in_mask)
        else
            mask_value = 0.0
        endif

        ! loop over pixels to apply the mask
        if (radius_given_as_fraction) then
            do k=1,self%logical_dimensions(3)
                z = (k - origin(3))**2 * inv_max_rad_sq(3)
                do j=1,self%logical_dimensions(2)
                    y = (j - origin(2))**2 * inv_max_rad_sq(2)
                    do i=1,self%logical_dimensions(1)
                        x = (i - origin(1))**2 * inv_max_rad_sq(1)
                        ! work out distance from origin (fraction of max radius)
                        if (ccylinder) then
                            dist = sqrt(x+y)
                        else
                            dist = sqrt(x+y+z)
                        endif
                        if (dist .ge. radius_before_falloff .and. dist .le. rradius_after_falloff) then
                            falloff_factor = (cos(((dist-radius_before_falloff)/falloff)*pi)+1.)/2.0
                            if (iinverse) falloff_factor = (falloff_factor-1.0)*(-1.0)
                            self%real_values(i,j,k) = self%real_values(i,j,k) * falloff_factor + (1.-falloff_factor)*mask_value
                        else if (dist .gt. rradius_after_falloff .and. .not. iinverse) then
                            self%real_values(i,j,k) = mask_value
                        else if (dist .le. rradius_after_falloff .and. iinverse) then
                            self%real_values(i,j,k) = mask_value
                        endif
                    enddo
                enddo
            enddo
        else
            do k=1,self%logical_dimensions(3)
                z = (k - origin(3))**2
                do j=1,self%logical_dimensions(2)
                    y = (j - origin(2))**2
                    do i=1,self%logical_dimensions(1)
                        x = (i - origin(1))**2
                        ! work out distance from origin
                        if (ccylinder) then
                            dist = sqrt(x+y)
                        else
                            dist = sqrt(x+y+z)
                        endif
                        if (dist .ge. radius_before_falloff .and. dist .le. rradius_after_falloff) then
                            falloff_factor = (cos(((dist-radius_before_falloff)/falloff)*pi)+1.)/2.0
                            if (iinverse) falloff_factor = (falloff_factor-1.0)*(-1.0)
                            self%real_values(i,j,k) = self%real_values(i,j,k) * falloff_factor + (1.-falloff_factor)*mask_value
                        else if (dist .gt. rradius_after_falloff .and. .not. iinverse) then
                            self%real_values(i,j,k) = mask_value
                        else if (dist .le. rradius_after_falloff .and. iinverse) then
                            self%real_values(i,j,k) = mask_value
                        endif
                    enddo
                enddo
            enddo
        endif

    end subroutine ApplySoftCircularMask

    !>  \brief  mask real image with a cosine-falloff-edged rectangle. rectangle is filled with 1s near the center, and 0s near the edge.
    pure subroutine ApplySoftRectangularMask(img,rect_dim,azimuth,falloff,wipe_input,mask_with_average,shift_x,shift_y, &
                            mask_with_average_edge_value)
        ! arguments
        class(image),           intent(inout)   ::  img                 !<  image to mask
        real,                   intent(in)      ::  rect_dim(2)         !<  dimensions of the rectangle, in pixels
        real,                   intent(in)      ::  azimuth             !<  in degrees, the clockwise angle that specifies the first axis of the rectangle. 0.0 indicates the x axis.
        real,                   intent(in)      ::  falloff             !<  distance from the start of the falloff (1.0) to the end of the falloff(0.0). the mask reaches 0.0 at rect_dim/2. in pixels.
        logical,    optional,   intent(in)      ::  wipe_input          !<  wipe / initialise the input to 1.0. .true. by default
        logical,    optional,   intent(in)      ::  mask_with_average   !<  if .true., the mask will go to the average value within the mask, rather than 0.0.
        real,       optional,   intent(in)      ::  shift_x,shift_y     !<  the center of the rectangle should be shifted by this much in x and y
        logical,    optional,   intent(in)      ::  mask_with_average_edge_value    !<  Mask with average grey value of edge pixels
        ! private variables
        logical         ::  wwipe_input, mmask_with_average, mmask_with_average_edge_value
        real            ::  mask_value
        !logical,    parameter   ::  debug   =   .false.
        ! start work

        ! wipe input?
        if (present(wipe_input)) then
            wwipe_input = wipe_input
        else
            wwipe_input = .true.
        endif

        ! mask with average?
        mmask_with_average = .false.
        if (present(mask_with_average)) mmask_with_average = mask_with_average
        mmask_with_average_edge_value = .false.
        if (present(mask_with_average_edge_value)) mmask_with_average_edge_value = mask_with_average_edge_value

        ! if required, wipe the input
        if (wwipe_input) then
            img%real_values = 1.0
        endif

        ! Work out the mask value
        if (mmask_with_average) then
            ! loop over image to compute the average value within the mask
            call cosine_rect_slave( img,rect_dim,azimuth,falloff,get_stats=.true.,mask_value=mask_value, &
                                    shift_x=shift_x,shift_y=shift_y)
        else if (mmask_with_average_edge_value) then
            mask_value = img%GetAverageOfValuesOnEdges()
        else
            mask_value = 0.0
        endif

        ! loop over pixels to apply the mask
        call cosine_rect_slave(img,rect_dim,azimuth,falloff,get_stats=.false.,mask_value=mask_value,shift_x=shift_x,shift_y=shift_y)

    end subroutine ApplySoftRectangularMask

    !>  \brief  convenience routine containing the looping structure to define a smooth-edged rectangle. can be called
    !!          either to compile statistics to determine average value within mask, or to apply mask
    pure subroutine cosine_rect_slave(img,rect_dim,azimuth,falloff,get_stats,mask_value,shift_x,shift_y)
        use UsefulFunctions, only : RadiansSingle
        implicit none
        ! arguments
        type(image),            intent(inout)   ::  img                 !<  image to mask
        real,                   intent(in)      ::  rect_dim(2)         !<  dimensions of the rectangle, in pixels
        real,                   intent(in)      ::  azimuth             !<  in degrees, the clockwise angle that specifies the first axis of the rectangle. 0.0 indicates the x axis.
        real,                   intent(in)      ::  falloff             !<  distance from the start of the falloff (1.0) to the end of the falloff(0.0). the mask reaches 0.0 at rect_dim/2. in pixels.
        logical,                intent(in)      ::  get_stats           !<  if true, we'll compute average value within mask rather than apply a mask
        real,                   intent(inout)   ::  mask_value          !<  greyvalue for the mask. if get_stats, this will be an output.
        real,       optional,   intent(in)      ::  shift_x             !<  center of the rectangle should be shifted by this much in x
        real,       optional,   intent(in)      ::  shift_y             !<  center of the rectangle should be shifted by this much in y
        ! private variables
        real    ::  origin(3)
        real    ::  cos_az, sin_az
        real    ::  rect_min(2), rect_max(2), edge_min(2), edge_max(2)
        real(kind=8)    ::  sum_in_mask
        real(kind=8)    ::  count_in_mask
        integer         ::  i,j
        real            ::  x,y,xx,yy,xx_tmp,yy_tmp
        real            ::  rad
        real            ::  falloff_factor
        ! start work

        ! values we'll need when looping
        origin      = [1.0,1.0,1.0]
        if (img%object_is_centered_in_box .and. img%is_in_real_space) origin = real(img%physical_address_of_box_center)
        if (present(shift_x)) origin(1) = origin(1) + shift_x
        if (present(shift_y)) origin(2) = origin(2) + shift_y
        cos_az      =  cos(RadiansSingle(azimuth))
        sin_az      =  sin(RadiansSingle(azimuth))
        rect_min    = -rect_dim/2.
        rect_max    =  rect_dim/2.
        edge_min    =  rect_min + falloff
        edge_max    =  rect_max - falloff

        ! loop over pixels and work out whether we are within the rectangle
        sum_in_mask = 0.0d0
        count_in_mask = 0
        do j=1,img%logical_dimensions(2)
            y = real(j - origin(2))
            xx_tmp = - sin_az * y
            yy_tmp = + cos_az * y
            do i=1,img%logical_dimensions(1)
                x = real(i - origin(1))
                ! work out coordinates before rotation
                xx = cos_az * x + xx_tmp
                yy = sin_az * x + yy_tmp
                ! is this point within the rectangle?
                if (        xx .ge. rect_min(1) .and. xx .le. rect_max(1)   &
                    .and.   yy .ge. rect_min(2) .and. yy .le. rect_max(2))  then
                    ! we are within the rectangle
                    ! now, check whether we are within the cosine-edge region of the rectangle
                    if (        xx .le. edge_min(1) .or. xx .ge. edge_max(1)    &
                        .or.    yy .le. edge_min(2) .or. yy .ge. edge_max(2)) then
                        ! we are within the falloff edge of the rectangle
                        ! we need to work out the falloff factor
                        ! we follow different logic depending on whether we are near the corners of the rectangle
                        if (xx .le. edge_min(1)) then
                            if (yy .le. edge_min(2)) then
                                !  we are in a corner. work out radius to corner
                                rad = sqrt((xx-edge_min(1))**2+(yy-edge_min(2))**2)
                            else if (yy .ge. edge_max(2)) then
                                ! we are in a corner
                                rad = sqrt((xx-edge_min(1))**2+(yy-edge_max(2))**2)
                            else
                                ! we are along an edge
                                rad = edge_min(1) - xx
                            endif
                        else if (xx .ge. edge_max(1)) then
                            if (yy .le. edge_min(2)) then
                                !  we are in a corner. work out radius to corner
                                rad = sqrt((xx-edge_max(1))**2+(yy-edge_min(2))**2)
                            else if (yy .ge. edge_max(2)) then
                                !  we are in a corner. work out radius to corner
                                rad = sqrt((xx-edge_max(1))**2+(yy-edge_max(2))**2)
                            else
                                ! we are along an edge
                                rad = xx-edge_max(1)
                            endif
                        else if (yy .le. edge_min(2)) then
                            ! we are along an edge
                            rad = edge_min(2)-yy
                        else if (yy .ge. edge_max(2)) then
                            ! we are along an edge
                            rad = yy-edge_max(2)
                        else
                            ! we must have got lost!
                            !write(*,'(a)') '**error(cosine_rect_slave): programming error.'
                            !call terminate('fatal error in cosine_rect_slave')
                        endif
                        ! now we've wored out how far we are from the start of the falloff, we can work out the falloff factor
                        if (rad .le. falloff) then
                            falloff_factor = (cos((rad/falloff)*pi)+1.)/2.0
                        else
                            falloff_factor = 0.0d0
                        endif
                    else
                        ! we are within the rectangle, not on its edge
                        falloff_factor = 1.0d0
                    endif ! end of test for whether we're in the rectangle and within its border region

                    ! now we've worked out the falloff factor for where we are, let's act on it
                    if (get_stats) then
                        ! update stats
                        sum_in_mask = sum_in_mask + dble(img%real_values(i,j,1) * falloff_factor)
                        count_in_mask = count_in_mask + dble(falloff_factor)
                    else
                        ! apply mask
                        img%real_values(i,j,1) = img%real_values(i,j,1) * falloff_factor + (1.-falloff_factor)*mask_value
                    endif
                else
                    ! the point is not within the rectangle at all
                    if (get_stats) then
                        ! nothing to do
                    else
                        img%real_values(i,j,1) = mask_value
                    endif
                endif
            enddo   ! end of loop over first dimension
        enddo ! end of loop over second dimension

        if (get_stats) then
            mask_value = real(sum_in_mask / count_in_mask)
        endif

    end subroutine cosine_rect_slave


    !>  \brief  mask real image with a rectangle.
    pure subroutine ApplyRectangularMask(img,rect_dim,azimuth,wipe_input,mask_with_average,shift_x,shift_y)
        use UsefulFunctions, only : RadiansSingle
        ! arguments
        class(image),           intent(inout)   ::  img                 !<  image to mask
        real,                   intent(in)      ::  rect_dim(2)         !<  dimensions of the rectangle, in pixels. length (along azimuth) first.
        real,                   intent(in)      ::  azimuth             !<  in degrees, the clockwise angle that specifies the first axis of the rectangle. 0.0 indicates the x axis.
        logical,    optional,   intent(in)      ::  wipe_input          !<  wipe / initialise the input to 1.0. .true. by default
        logical,    optional,   intent(in)      ::  mask_with_average   !<  if .true., the mask will go to the average value within the mask, rather than 0.0.
        real,       optional,   intent(in)      ::  shift_x,shift_y     !<  the center of the rectangle should be shifted by this much in x and y
        ! private variables
        logical         ::  wwipe_input, mmask_with_average
        real            ::  mask_value
        !logical,    parameter   ::  debug   =   .false.
        ! start work

        ! wipe input?
        if (present(wipe_input)) then
            wwipe_input = wipe_input
        else
            wwipe_input = .true.
        endif

        ! mask with average?
        if (present(mask_with_average)) then
            mmask_with_average = mask_with_average
        else
            mmask_with_average = .false.
        endif

        ! if required, wipe the input
        if (wwipe_input) then
            img%real_values = 1.0
        endif

        if (mmask_with_average) then
            ! loop over image to compute the average value within the mask
            call rect_mask_slave(img,rect_dim,azimuth,get_stats=.true.,mask_value=mask_value,shift_x=shift_x,shift_y=shift_y)
        else
            mask_value = 0.0
        endif

        ! loop over pixels to apply the mask
        call rect_mask_slave(img,rect_dim,azimuth,get_stats=.false.,mask_value=mask_value,shift_x=shift_x,shift_y=shift_y)

    end subroutine ApplyRectangularMask

    !>  \brief  Zero pixels outside specified rectangle
    pure subroutine rect_mask_slave(img,rect_dim,azimuth,get_stats,mask_value,shift_x,shift_y)
        use UsefulFunctions, only : RadiansSingle
        ! arguments
        type(image),            intent(inout)   ::  img                 !<  image to mask
        real,                   intent(in)      ::  rect_dim(2)         !<  dimensions of the rectangle, in pixels. first is length (along azimuth), second is width.
        real,                   intent(in)      ::  azimuth             !<  in degrees, the clockwise angle that specifies the first axis of the rectangle. 0.0 indicates the x axis.
        logical,                intent(in)      ::  get_stats           !<  if true, we'll compute average value within mask rather than apply a mask
        real,                   intent(inout)   ::  mask_value          !<  greyvalue for the mask. if get_stats, this will be an output.
        real,       optional,   intent(in)      ::  shift_x             !<  center of the rectangle should be shifted by this much in x
        real,       optional,   intent(in)      ::  shift_y             !<  center of the rectangle should be shifted by this much in y
        ! private variables
        real            ::  origin(3)
        real            ::  cos_az, sin_az
        real            ::  rect_min(2), rect_max(2)
        real(kind=8)    ::  sum_in_mask
        real(kind=8)    ::  count_in_mask
        integer         ::  i,j
        real            ::  x,y,xx,yy,xx_tmp,yy_tmp
        real            ::  falloff_factor
        ! start work

        ! values we'll need when looping
        origin      = [1.0,1.0,1.0]
        if (img%object_is_centered_in_box .and. img%is_in_real_space) origin = real(img%physical_address_of_box_center)
        if (present(shift_x)) origin(1) = origin(1) + shift_x
        if (present(shift_y)) origin(2) = origin(2) + shift_y
        cos_az      =  cos(RadiansSingle(azimuth))
        sin_az      =  sin(RadiansSingle(azimuth))
        rect_min    = -rect_dim/2.
        rect_max    =  rect_dim/2.

        ! loop over pixels and work out whether we are within the rectangle
        sum_in_mask = 0.0d0
        count_in_mask = 0
        do j=1,img%logical_dimensions(2)
            y = real(j - origin(2))
            xx_tmp = - sin_az * y
            yy_tmp = + cos_az * y
            do i=1,img%logical_dimensions(1)
                x = real(i - origin(1))
                ! work out coordinates before rotation
                xx = cos_az * x + xx_tmp
                yy = sin_az * x + yy_tmp
                ! is this point within the rectangle?
                if (        xx .ge. rect_min(1) .and. xx .le. rect_max(1)   &
                    .and.   yy .ge. rect_min(2) .and. yy .le. rect_max(2))  then
                    ! we are within the rectangle
                    falloff_factor = 1.0
                    if (get_stats) then
                        ! update stats
                        sum_in_mask = sum_in_mask + dble(img%real_values(i,j,1) * falloff_factor)
                        count_in_mask = count_in_mask + dble(falloff_factor)
                    endif
                else
                    falloff_factor = 0.0
                    if (get_stats) then
                        ! nothing to do
                    else
                        img%real_values(i,j,1) = mask_value
                    endif
                endif
            enddo   ! end of loop over first dimension
        enddo ! end of loop over second dimension

        if (get_stats) then
            mask_value = real(sum_in_mask / count_in_mask)
        endif

    end subroutine rect_mask_slave



     !>  \brief  Find a peak with integer accuracy
    function FindPeakWithIntegerCoordinates(self, min_radius, max_radius) result(found_peak)

        use Peaks

        !Arguments

        class(Image),           intent(in)   ::  self
        real, optional,         intent(in)   ::  min_radius
        real, optional,         intent(in)   ::  max_radius

        !Variables
        type(Peak)                           ::  found_peak
        integer                              ::  location_of_maximum_value(3)
        real                                 ::  max_found_value

        real                                 ::  applied_max_radius
        real                                 ::  applied_min_radius
        real                                 ::  distance_from_origin
        real                                 ::  inv_max_rad_sq(3)

        integer                              ::  i
        integer                              ::  j
        integer                              ::  k


        real                                 ::  x
        real                                 ::  y
        real                                 ::  z

        logical                              ::  radii_are_fractional

        !Start

        max_found_value = -huge(1.0e0)

        !
        ! Parse the radii
        !
        radii_are_fractional = .false.
        applied_min_radius   = 0.0e0
        applied_max_radius   = self%GetMaximumDiagonalRadius()**2
        if (present(min_radius)) applied_min_radius = min_radius**2
        if (present(max_radius)) applied_max_radius = max_radius**2
        ! We will understand the radii to be fractional only if both are >0.0 and <1.0
        if (        applied_min_radius .gt. 0.0e0 .and. applied_min_radius .lt. 1.0 &
            .and.   applied_max_radius .gt. 0.0e0 .and. applied_max_radius .lt. 1.0 ) then
            radii_are_fractional = .true.
        endif

        !write(*,'(a,1x,2(f0.2,1x),l1)') '**debug(FindPeakWithIntegerCoordinates): min, max radii, fractional?', sqrt(applied_min_radius), sqrt(applied_max_radius), radii_are_fractional



       if (self%IsAVolume()) then
            inv_max_rad_sq = 1.0 / (self%physical_address_of_box_center - 1)**2
            else
            inv_max_rad_sq(1:2) = 1.0 / ((self%physical_address_of_box_center(1:2)-1)**2)
        end if

        ! loop over the data to find the max value and location

        if (radii_are_fractional) then
            do k=1,self%logical_dimensions(3)
                z = (k - self%physical_address_of_box_center(3))**2 * inv_max_rad_sq(3)
                do j=1,self%logical_dimensions(2)
                    y = (j - self%physical_address_of_box_center(2))**2 * inv_max_rad_sq(2)
                    do i=1,self%logical_dimensions(1)
                        x = (i - self%physical_address_of_box_center(1))**2 * inv_max_rad_sq(1)

                        distance_from_origin = x + y + z

                        if (distance_from_origin .ge. applied_min_radius .and. distance_from_origin .le. applied_max_radius) then

                            if (self%real_values(i, j, k) .gt. max_found_value) then
                                max_found_value = self%real_values(i, j, k)
                                location_of_maximum_value = [i, j, k]
                            endif

                        endif

                    enddo
                enddo
            enddo
        else
            do k=1,self%logical_dimensions(3)
                z = (k - self%physical_address_of_box_center(3))**2
                do j=1,self%logical_dimensions(2)
                    y = (j - self%physical_address_of_box_center(2))**2
                    do i=1,self%logical_dimensions(1)
                        x = (i - self%physical_address_of_box_center(1))**2
                        !work out distance from origin

                        distance_from_origin = x + y + z

                        if (distance_from_origin .ge. applied_min_radius .and. distance_from_origin .le. applied_max_radius) then

                            if (self%real_values(i, j, k) .gt. max_found_value) then
                                max_found_value = self%real_values(i, j, k)
                                location_of_maximum_value = [i, j, k]
                            endif

                        endif

                    enddo
                enddo
            enddo

        endif

        found_peak%CoOrdinates = location_of_maximum_value - self%physical_address_of_box_center
        found_peak%PeakHeight = max_found_value

        !write(*,'(a,2(f0.2,1x))') '**debug(FindPeakWithIntegerCoordinates): peak coos ', found_peak%CoOrdinates(1:2)
    end function FindPeakWithIntegerCoordinates

    !>  \brief  Mask out the central cross in Fourier space to try to reduce line artifacts
    subroutine MaskCentralCross(self, vertical_half_width, horizontal_half_width)

        !Arguments

        class(Image),           intent(inout)   ::  self
        integer, optional,      intent(in)      ::  vertical_half_width
        integer, optional,      intent(in)      ::  horizontal_half_width

        !Variables

        logical                              ::  must_fft
        integer                              ::  pixel_counter
        integer                              ::  width_counter
        integer                              ::  used_vertical_half_width
        integer                              ::  used_horizontal_half_width


        !start

        if (present(vertical_half_width)) then
            used_vertical_half_width = vertical_half_width
        else
            used_vertical_half_width = 1
        endif

        if (present(horizontal_half_width)) then
            used_horizontal_half_width = horizontal_half_width
        else
            used_horizontal_half_width = 1
        endif

        must_fft = .false.
        if (self%is_in_real_space) must_fft = .true.


        if (must_fft) call self%ForwardFFT()

        ! do the masking..

        do pixel_counter=1, self%GetLogicalDimension(2)
            do width_counter = 1, used_vertical_half_width
                self%complex_values(width_counter, pixel_counter, 1) = (0e0, 0e0)
            enddo
        enddo

        do pixel_counter=1, self%physical_upper_bound_complex(1)
            do width_counter = 1, used_horizontal_half_width

                self%complex_values(pixel_counter, width_counter, 1) = (0e0, 0e0)

                if (width_counter .gt. 1) then
                    self%complex_values(pixel_counter, self%logical_dimensions(2) + 2 - width_counter, 1) = (0e0, 0e0)
                endif

            enddo
        enddo

        if (must_fft) call self%BackwardFFT()

    end subroutine MaskCentralCross

    !>  \brief  Find a peak with sub pixel accuracy, by fitting a parabola (2d only at the moment)
    function FindPeakWithParabolaFit(self, min_radius, max_radius) result(found_peak)  ! This code is based on code written by Ardan Patwardhan

        use Peaks

        !Arguments

        class(Image),           intent(in)   ::  self
        real, optional,         intent(in)   ::  min_radius
        real, optional,         intent(in)   ::  max_radius

        !Variables
        type(Peak)                           ::  found_peak
        type(Peak)                           ::  integer_peak

        real                                 ::  average_of_square
        real                                 ::  scale_factor

        real                                 ::  coefficient_one
        real                                 ::  coefficient_two
        real                                 ::  coefficient_three
        real                                 ::  coefficient_four
        real                                 ::  coefficient_five
        real                                 ::  coefficient_six

        real                                 ::  denominator

        real                                 ::  scaled_square(3,3)

        real                                 ::  y_max
        real                                 ::  x_max


        integer                              ::  best_i
        integer                              ::  best_j
        integer                              ::  current_i
        integer                              ::  current_j
        integer                              ::  address_i
        integer                              ::  address_j


        !Start - first of all, find the integer peak..

        ! The 3D version of this is pretty complex, and for now I don't want to implement it - however i do have it.
        ! So, if you are reading this, and want a 3D version of this function - ask me to do it - Tim.

        if (self%IsAVolume()) then
            call this_program%TerminateWithFatalError('Image::FindPeakWithParabolaFit', '3D volumes not currently supported')
        endif

        if (present(min_radius) .and. present(max_radius)) then
            integer_peak = self%FindPeakWithIntegerCoordinates(min_radius=min_radius, max_radius=max_radius)
        else if (present(min_radius)) then
            integer_peak = self%FindPeakWithIntegerCoordinates(min_radius=min_radius)
        else if (present(max_radius)) then
            integer_peak = self%FindPeakWithIntegerCoordinates(max_radius=max_radius)
        else
            integer_peak = self%FindPeakWithIntegerCoordinates()
        endif

        best_i = integer_peak%CoOrdinates(1) + self%physical_address_of_box_center(1)
        best_j = integer_peak%CoOrdinates(2) + self%physical_address_of_box_center(2)

        average_of_square = 0.0e0
        average_of_square = average_of_square + sum(self%real_values(best_i-1:best_i+1, best_j-1:best_j+1, 1))
        average_of_square = average_of_square / 9.0e0

        if (average_of_square .ne. 0.0e0) then
            scale_factor = 1.0e0 / average_of_square
        else
            scale_factor = 1.0e0
        endif

        if (best_i .gt. 1                          .and. best_j .gt. 1 .and. &
            best_i .lt. self%logical_dimensions(1) .and. best_j .lt. self%logical_dimensions(2)) then
            scaled_square =  self%real_values(best_i-1:best_i+1, best_j-1:best_j+1, 1)
        else
            do current_j = 1, 3
                do current_i = 1, 3
                    address_i = best_i - 2 + current_i
                    address_j = best_j - 2 + current_j
                    if (address_i .ge. 1                          .and. address_j .ge. 1 .and. &
                        address_i .le. self%logical_dimensions(1) .and. address_j .le. self%logical_dimensions(2)) then
                        scaled_square(current_i, current_j) = self%real_values(address_i, address_j, 1)
                    else
                        scaled_square(current_i, current_j) = 0.
                    endif
                enddo
            enddo
        endif

        scaled_square =  scaled_square * scale_factor

        ! now fit a parabola to a 3x3 box

        coefficient_one = (26.0e0 * scaled_square(1,1) - scaled_square(1,2) + &
                            2.0e0 * scaled_square(1,3) - scaled_square(2,1) - &
                           19.0e0 * scaled_square(2,2) - &
                            7.0e0 * scaled_square(2,3) + &
                            2.0e0 * scaled_square(3,1) - &
                            7.0e0 * scaled_square(3,2) + &
                           14.0e0 * scaled_square(3,3)) / 9.0e0

        coefficient_two = ( 8.0e0 * scaled_square(1,1) - &
                            8.0e0 * scaled_square(1,2) + &
                            5.0e0 * scaled_square(2,1) - &
                            8.0e0 * scaled_square(2,2) + &
                            3.0e0 * scaled_square(2,3) + &
                            2.0e0 * scaled_square(3,1) - &
                            8.0e0 * scaled_square(3,2) + &
                            6.0e0 * scaled_square(3,3)) / (-1.0e0 * 6.0e0)

        coefficient_three = (scaled_square(1,1) - &
                            2.0e0 * scaled_square(1,2) + scaled_square(1,3) + scaled_square(2,1) - &
                            2.0e0 * scaled_square(2,2) + scaled_square(2,3) + scaled_square(3,1) - &
                            2.0e0 * scaled_square(3,2) + scaled_square(3,3)) / 6.0e0;

        coefficient_four = (8.0e0 * scaled_square(1,1) + &
                            5.0e0 * scaled_square(1,2) + &
                            2.0e0 * scaled_square(1,3) - &
                            8.0e0 * scaled_square(2,1) - &
                            8.0e0 * scaled_square(2,2) - &
                            8.0e0 * scaled_square(2,3) + &
                            3.0e0 * scaled_square(3,2) + &
                            6.0e0 * scaled_square(3,3)) / (-1.0e0 * 6.0e0)

        coefficient_five = (scaled_square(1,1) - scaled_square(1,3) - scaled_square(3,1) + scaled_square(3,3)) / 4.0e0;

        coefficient_six = (scaled_square(1,1) + scaled_square(1,2) + scaled_square(1,3) - &
                            2.0e0 * scaled_square(2,1) - &
                            2.0e0 * scaled_square(2,2) - &
                            2.0e0 * scaled_square(2,3) + &
                            scaled_square(3,1) + scaled_square(3,2) + scaled_square(3,3)) / 6.0e0;

        denominator = 4.0e0 * coefficient_three* coefficient_six - coefficient_five**2

        if (denominator .eq. 0) then ! something has gone wrong, just return the integer peak
            found_peak = integer_peak
        else
            y_max = (coefficient_four * coefficient_five - 2.0e0 * coefficient_two * coefficient_six) / denominator
            x_max = (coefficient_two * coefficient_five - 2.0e0 * coefficient_four * coefficient_three) / denominator

            y_max = y_max - 2.0e0
            x_max = x_max - 2.0e0

            if (y_max .gt. 1.05e0 .or. y_max .lt. -1.05e0) then ! we are probably strangely shaped and the fit is diverging.. just take the integer
                y_max = 0.0e0
            endif

            if (x_max .gt. 1.05e0 .or. x_max .lt. -1.05e0) then
                x_max = 0.0e0
            endif

            found_peak = integer_peak
            found_peak%CoOrdinates(1) = found_peak%CoOrdinates(1) + x_max
            found_peak%CoOrdinates(2) = found_peak%CoOrdinates(2) + y_max

            found_peak%PeakHeight = 4.0e0 * coefficient_one * coefficient_three * coefficient_six - &
                                            coefficient_one * coefficient_five**2 - coefficient_two**2 * coefficient_six + &
                                            coefficient_two * coefficient_four*coefficient_five - &
                                            coefficient_four**2 * coefficient_three
            found_peak%PeakHeight = found_peak%PeakHeight * average_of_square / denominator

            if (abs((found_peak%PeakHeight - integer_peak%PeakHeight) / &
                    (found_peak%PeakHeight + integer_peak%PeakHeight)) > 0.15) then !Probably this value has gone screwy, use the integer value
                found_peak%PeakHeight = integer_peak%PeakHeight
            endif
        endif

     end function FindPeakWithParabolaFit

    !>  \brief  Apply a phase shift to an image
    subroutine PhaseShift(self, x_shift, y_shift, z_shift)

        use UsefulFunctions

        class(Image),           intent(inout)   ::  self
        real,                   intent(in)      ::  x_shift
        real,                   intent(in)      ::  y_shift
        real,                   intent(in)      ::  z_shift

        integer                                 ::  i
        integer                                 ::  j
        integer                                 ::  k

        integer                                 ::  i_logical
        integer                                 ::  j_logical
        integer                                 ::  k_logical

        real                                    ::  phase_z
        real                                    ::  phase_y
        real                                    ::  phase_x

        complex                                 ::  phase_shift_to_apply

        logical                                 ::  need_to_fft

        ! Start work

        need_to_fft = .false.

        if (self%IsInRealSpace()) then
            need_to_fft = .true.
            call self%ForwardFFT()
        endif


        ! Loop over voxels
        do k=1, self%physical_upper_bound_complex(3)
            k_logical = self%LogicalIndexGivenPhysicalIndexInFourierSpace(k, 3)
            phase_z = ReturnPhaseFromShift(z_shift, k_logical, self%logical_dimensions(3))

            do j=1,self%physical_upper_bound_complex(2)
                j_logical = self%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2)
                 phase_y = ReturnPhaseFromShift(y_shift, j_logical, self%logical_dimensions(2))
                 do i=1,self%physical_upper_bound_complex(1)
                    i_logical = i-1
                    phase_x = ReturnPhaseFromShift(x_shift, i_logical, self%logical_dimensions(1))

                    phase_shift_to_apply = Return3DPhaseFromIndividualDimensions(phase_x, phase_y, phase_z)

                    self%complex_values(i, j, k) = self%complex_values(i, j, k) * phase_shift_to_apply
                enddo
            enddo
        enddo

        if (need_to_fft) call self%BackwardFFT()


    end subroutine PhaseShift

    !>  \brief Swap the qudrants of a real-space image. the "centre" pixel is moved to the first pixel (corner with address 1,1,1), and vice versa.
    !!  The image is always returned in the space it was given
    subroutine SwapRealSpaceQuadrants(self)
        use UsefulFunctions, only : IsEven
        ! Arguments
        class(Image),   intent(inout)   ::  self

        ! Private Variables
        real                            ::  shifts_to_apply(3)
        logical                         ::  need_to_fft
        integer                         ::  current_dim

        ! Start work

        ! Check image type on entry
        need_to_fft = self%IsInRealSpace()

        ! Work out shift to apply
        if (.not. self%object_is_centered_in_box) then
            shifts_to_apply(1:3) = real(self%physical_address_of_box_center) - 1.0e0
        else
            do current_dim=1,3
                if (IsEven(self%logical_dimensions(current_dim))) then
                    shifts_to_apply(current_dim) = real(self%physical_address_of_box_center(current_dim)) - 1.0e0
                else
                    shifts_to_apply(current_dim) = real(self%physical_address_of_box_center(current_dim))
                endif
            enddo
        endif


        if (self%logical_dimensions(3) .eq. 1) shifts_to_apply(3) = 0.0e0

        ! Apply appropriate phase shift

        if (need_to_fft) call self%ForwardFFT()

        call self%PhaseShift(shifts_to_apply(1), shifts_to_apply(2), shifts_to_apply(3))

        if (need_to_fft) call self%BackwardFFT()

        ! Keep track of where the real-space origin is now
        self%object_is_centered_in_box = .not. self%object_is_centered_in_box

   end subroutine SwapRealSpaceQuadrants

    !>  \brief  Correct for magnification distortion in the image
    subroutine CorrectMagnificationDistortion(self,output_image,distortion_angle,x_scale,y_scale)
        use Units
        ! Arguments
        class(Image),   intent(inout)   ::  self
        type(Image),    intent(inout)   ::  output_image        !<  Output, corrected image
        real,           intent(in)      ::  distortion_angle    !<  Angle of the distortion (degrees)
        real,           intent(in)      ::  x_scale             !<  Scale factor in X (1.0 = no distortion)
        real,           intent(in)      ::  y_scale             !<  Scale factor in Y (1.0 = no distortion)
        ! Private variables
        real                            ::  angle_in_radians
        real                            ::  x_scale_factor
        real                            ::  y_scale_factor
        real                            ::  average_edge_value
        integer                         ::  x
        integer                         ::  y
        real                            ::  new_x
        real                            ::  new_y
        real                            ::  final_x
        real                            ::  final_y
        real                            ::  interpolated_value
        ! Optionally, can pad the input images in Fourier space before doing real-space
        ! interpolation, to reduce artefacts in the output images
        integer,    parameter           ::  fourier_padding_factor  =   1
        integer                         ::  original_dimensions(3)
        integer                         ::  padded_dimensions(3)

        ! Start work
        angle_in_radians = convert(distortion_angle,degrees,radians)
        x_scale_factor = 1.0e0 / x_scale
        y_scale_factor = 1.0e0 / y_scale

        if (fourier_padding_factor .gt. 1) then
            original_dimensions = self%GetLogicalDimensions()
            padded_dimensions = original_dimensions * fourier_padding_factor
            padded_dimensions(3) = 1
            call self%ForwardFFT()
            call self%Resize(padded_dimensions)
            call self%BackwardFFT()
        endif

        average_edge_value = self%GetAverageOfValuesOnEdges()


        ! Loop over the output image
        !$omp parallel default(shared) private(y, x, new_x, new_y, final_x, final_y, interpolated_value)
        !$omp do
        do y = 1, self%GetLogicalDimension(2)
            do x = 1, self%GetLogicalDimension(1)

            ! Rotation

            new_x =  real(y - self%GetPhysicalIndexOfBoxCenter(2))*sin(-angle_in_radians) &
                    +real(x - self%GetPhysicalIndexOfBoxCenter(1))*cos(-angle_in_radians)
            new_y =  real(y - self%GetPhysicalIndexOfBoxCenter(2))*cos(-angle_in_radians) &
                    -real(x - self%GetPhysicalIndexOfBoxCenter(1))*sin(-angle_in_radians)

            ! Scaling

            new_x = new_x * x_scale_factor
            new_y = new_y * y_scale_factor

            new_x = new_x + self%GetPhysicalIndexOfBoxCenter(1)
            new_y = new_y + self%GetPhysicalIndexOfBoxCenter(2)


            ! Rotate back

            final_x = real(new_y - self%GetPhysicalIndexOfBoxCenter(2))*sin(angle_in_radians) &
                     +real(new_x - self%GetPhysicalIndexOfBoxCenter(1))*cos(angle_in_radians)
            final_y = real(new_y - self%GetPhysicalIndexOfBoxCenter(2))*cos(angle_in_radians) &
                     -real(new_x - self%GetPhysicalIndexOfBoxCenter(1))*sin(angle_in_radians)

            final_x = final_x + self%GetPhysicalIndexOfBoxCenter(1)
            final_y = final_y + self%GetPhysicalIndexOfBoxCenter(2)

            ! Current value

            if (     int(final_x) .lt. 1 .or. int(final_x) + 1 .gt. self%GetLogicalDimension(1) &
                .or. int(final_y) .lt. 1 .or. int(final_y) + 1 .gt. self%GetLogicalDimension(2)) then
                output_image%real_values(x, y, 1) = average_edge_value
            else
                call self%GetRealValueByLinearInterpolationNoBoundsCheckImage(interpolated_value, final_x, final_y)
                output_image%real_values(x, y, 1) = interpolated_value
            endif

            enddo
        enddo
        !$omp enddo
        !$omp end parallel

        if (fourier_padding_factor .gt. 1) then
            call output_image%ForwardFFT()
            call output_image%Resize(original_dimensions)
            call output_image%BackwardFFT()
        endif

    end subroutine CorrectMagnificationDistortion


    !>  \brief  Apply an affine transformation (e.g. some combination of rotation and shift). The transform matrix
    !!          must be (4,4) and must describe the transformation from output to input.
    !!  \todo   Change the convention, so that the matrix defines the transformation from input to output (requires that every single call be checked carefully)
    subroutine ApplyAffineTransformation(self,output_image,transformation_matrix)
        class(Image),   intent(in)      ::  self
        type(Image),    intent(inout)   ::  output_image
        real,           intent(in)      ::  transformation_matrix(4,4)
        ! private variables
        real    ::  output_coordinates_homogeneous(4)
        real    ::  input_coordinates_homogeneous(4)
        integer ::  i,j,k
        real    ::  background_value
        ! start work

        ! Make sure the output image is correctly allocated
        call output_image%Allocate(mould=self)

        !
        output_coordinates_homogeneous(1:3) = 0.0
        output_coordinates_homogeneous(4) = 1.0
        input_coordinates_homogeneous(1:3) = 0.0
        input_coordinates_homogeneous(4) = 1.0

        ! Are we in real space?
        if (self%IsInRealSpace()) then
            ! We are in real space
            ! choose a background value
            background_value = self%real_values(1,1,1)
            if (self%IsAVolume()) then
                ! Loop over output voxels
                do k=1,output_image%logical_dimensions(3)
                    output_coordinates_homogeneous(3) = real(k-output_image%physical_address_of_box_center(3))
                    do j=1,output_image%logical_dimensions(2)
                        output_coordinates_homogeneous(2) = real(j-output_image%physical_address_of_box_center(2))
                        do i=1,output_image%logical_dimensions(1)
                            output_coordinates_homogeneous(1) = real(i-output_image%physical_address_of_box_center(1))
                            !
                            output_image%real_values(i,j,k) = background_value
                            !
                            input_coordinates_homogeneous = matmul(transformation_matrix,output_coordinates_homogeneous)
                            ! Check we are within the input volume
                            if (        input_coordinates_homogeneous(1) .lt. self%logical_lower_bound_real(1) &
                                .or.    input_coordinates_homogeneous(2) .lt. self%logical_lower_bound_real(2) &
                                .or.    input_coordinates_homogeneous(3) .lt. self%logical_lower_bound_real(3) &
                                .or.    input_coordinates_homogeneous(1) .gt. self%logical_upper_bound_real(1) &
                                .or.    input_coordinates_homogeneous(2) .gt. self%logical_upper_bound_real(2) &
                                .or.    input_coordinates_homogeneous(3) .gt. self%logical_upper_bound_real(3) ) cycle
                            ! Interpolate
                            call self%GetRealValueByLinearInterpolationNoBoundsCheckVolume(output_image%real_values(i,j,k), &
                                                                                           input_coordinates_homogeneous(1), &
                                                                                           input_coordinates_homogeneous(2), &
                                                                                           input_coordinates_homogeneous(3) )
                        enddo
                    enddo
                enddo
            else
                ! Loop over output voxels
                do j=1,output_image%logical_dimensions(2)
                    output_coordinates_homogeneous(2) = real(j-output_image%physical_address_of_box_center(2))
                    do i=1,output_image%logical_dimensions(1)
                        output_coordinates_homogeneous(1) = real(i-output_image%physical_address_of_box_center(1))
                        !
                        output_image%real_values(i,j,1) = background_value
                        !
                        input_coordinates_homogeneous = matmul(transformation_matrix,output_coordinates_homogeneous)
                        input_coordinates_homogeneous(1:2) = input_coordinates_homogeneous(1:2) + &
                                                                self%physical_address_of_box_center(1:2)
                        ! Check we are within the input volume
                        if (        input_coordinates_homogeneous(1) .lt. 1.0 &
                            .or.    input_coordinates_homogeneous(2) .lt. 1.0 &
                            .or.    input_coordinates_homogeneous(1) .gt. real(self%logical_dimensions(1))-0.01 &
                            .or.    input_coordinates_homogeneous(2) .gt. real(self%logical_dimensions(2))-0.01 ) cycle
                        ! Interpolate
                        call self%GetRealValueByLinearInterpolationNoBoundsCheckImage(output_image%real_values(i,j,1), &
                                           input_coordinates_homogeneous(1), &
                                           input_coordinates_homogeneous(2))
                    enddo
                enddo
            endif
        else
            ! We are in Fourier space
            call this_program%TerminateWithFatalError('Image::ApplyAffineTransformation','Not implemented')
        endif ! end of test for real space

    end subroutine ApplyAffineTransformation


    !>  \brief  Compute rotational average of a power spectrum, taking into account astigmatism
    subroutine ComputeRotationalAverageOfPowerSpectrum( self,ctf,spatial_frequency,average,average_fit, &
                                                        frc,frc_sigma,rescale_input,squared_ctf_was_fit)
        use ContrastTransferFunctions
        use UsefulFunctions, only : NormalizedCrossCorrelation
        use Units, only : convert, radians, degrees
        class(Image),                   intent(inout)   ::  self                    !<  Power spectrum. Only modified if subtract_background_from_input
        type(ContrastTransferFunction), intent(in)      ::  ctf                     !<  CTF parameters which fit the power spectrum
        real(kind=8),                   intent(inout)   ::  spatial_frequency(:)    !<  Spatial frequency (1/pixel; Nyquist is 0.5)
        real(kind=8),                   intent(inout)   ::  average(:)              !<  Rotational average
        real(kind=8),                   intent(inout)   ::  average_fit(:)          !<  Average CTF, which should be comparable to the rotational average
        real(kind=8),                   intent(inout)   ::  frc(:)                  !<  Correlation between the power spectrum and the fit CTF
        real(kind=8),                   intent(inout)   ::  frc_sigma(:)            !<  3-sigma threshold for correlation of power spectrum with fit CTF
        logical,                        intent(in)      ::  rescale_input           !<  Modify input power spectrum
        logical,                        intent(in)      ::  squared_ctf_was_fit     !<  Whether the square of the CTF or just the absolute values were fit
        ! private variables
        type(Image) ::  number_of_extrema, ctf_values
        integer     ::  i,j
        real        ::  i_logi,j_logi
        real        ::  i_logi_sq,j_logi_sq
        real        ::  current_spatial_frequency_squared
        real        ::  current_azimuth
        integer     ::  number_of_bins
        integer     ::  current_bin
        integer     ::  chosen_bin
        real        ::  azimuth_of_mid_defocus
        real        ::  ctf_diff_from_current_bin
        real        ::  number_of_extrema_profile(size(spatial_frequency)), ctf_values_profile(size(spatial_frequency))
        integer     ::  number_of_values(size(spatial_frequency))
        logical,    parameter   ::  rescale_average = .true.
        integer,    parameter   ::  rescale_based_on_maximum_number = 2
        integer                 ::  current_maximum_number
        logical                 ::  at_a_minimum, at_a_maximum
        logical                 ::  minimum_at_previous_bin, maximum_at_previous_bin
        real                    ::  background(size(spatial_frequency)), peak(size(spatial_frequency))
        integer                 ::  location_of_previous_minimum, location_of_previous_maximum
        integer                 ::  normalisation_bin_number
        real                    ::  angular_distance_from_axes
        real,       parameter   ::  min_angular_distance_from_axes_degrees = 10.0
        ! debug
        type(Image)             ::  background_image, peak_image
        real                    ::  rad
        logical                 ::  actually_do_rescaling
        logical                 ::  self_is_blank
        ! start work

        self_is_blank = self%IsConstant()

        ! For each pixel, count how many extrema precede it
        if (.not. self_is_blank) then
            call number_of_extrema%Allocate(mould=self)
            call ctf_values%Allocate(mould=self)
            do j=1,self%logical_dimensions(2)
                j_logi = real(j-self%physical_address_of_box_center(2))*self%fourier_voxel_size(2)
                j_logi_sq = j_logi**2
                do i=1,self%logical_dimensions(1)
                    i_logi = real(i-self%physical_address_of_box_center(1))*self%fourier_voxel_size(1)
                    i_logi_sq = i_logi**2
                    ! Where are we?
                    current_spatial_frequency_squared = j_logi_sq+i_logi_sq
                    if (current_spatial_frequency_squared .gt. 0.0) then
                        current_azimuth = atan2(j_logi,i_logi)
                    else
                        current_azimuth = 0.0
                    endif
                    !
                    ctf_values%real_values(i,j,1) = ctf%EvaluateAtSquaredSpatialFrequency(current_spatial_frequency_squared, &
                                                                                          current_azimuth)
                    !
                    number_of_extrema%real_values(i,j,1) = ctf%CountNumberOfExtremaBeforeSquaredSpatialFrequency( &
                                                                            current_spatial_frequency_squared,current_azimuth)
                enddo
            enddo
        endif

        ! How many bins?
        !number_of_bins = ceiling(sqrt(sum((1.0-real(self%physical_address_of_box_center))**2))+1)
        number_of_bins = size(average)
        !if (allocated(average)) deallocate(average)
        !if (allocated(spatial_frequency)) deallocate(spatial_frequency)
        !if (allocated(average_fit)) deallocate(average_fit)
        !allocate(spatial_frequency(number_of_bins),average(number_of_bins))
        !allocate(number_of_extrema_profile(number_of_bins),ctf_values_profile(number_of_bins),average_fit(number_of_bins))
        !allocate(number_of_values(number_of_bins))

        average = 0.0d0
        spatial_frequency = 0.0d0
        average_fit = 0.0d0
        number_of_extrema_profile = 0.0
        ctf_values_profile = 0.0
        average_fit = 0.0
        number_of_values = 0

        if (.not. self_is_blank) then
            ! For each bin of our 1D profile, we compute the CTF. We choose the azimuth to be mid way between the two defoci of the astigmatic CTF
            azimuth_of_mid_defocus = ctf%GetAstigmatismAzimuthInRadians() + pi*0.25e0
            ! We don't want the azimuth too close to the axes, which may have been blanked by the central-cross-artefact-suppression-system (tm)
            angular_distance_from_axes = convert(mod(azimuth_of_mid_defocus,pi/2.0),radians,degrees)
            if (abs(angular_distance_from_axes) .lt. min_angular_distance_from_axes_degrees) then
                if (angular_distance_from_axes .gt. 0.0) then
                    azimuth_of_mid_defocus =  convert(min_angular_distance_from_axes_degrees,degrees,radians)
                else
                    azimuth_of_mid_defocus = -convert(min_angular_distance_from_axes_degrees,degrees,radians)
                endif
            endif
            if (abs(angular_distance_from_axes) .gt. 90.0-min_angular_distance_from_axes_degrees) then
                if (angular_distance_from_axes .gt. 0.0) then
                    azimuth_of_mid_defocus =  pi*0.5 - convert(min_angular_distance_from_axes_degrees,degrees,radians)
                else
                    azimuth_of_mid_defocus = -pi*0.5 + convert(min_angular_distance_from_axes_degrees,degrees,radians)
                endif
            endif
            ! Now that we've chosen an azimuth, we can compute the CTF for each bin of our 1D profile.
            do current_bin=1,number_of_bins
                current_spatial_frequency_squared = (real(current_bin-1)*self%fourier_voxel_size(2))**2
                spatial_frequency(current_bin) = sqrt(current_spatial_frequency_squared)
                ctf_values_profile(current_bin) = ctf%EvaluateAtSquaredSpatialFrequency(current_spatial_frequency_squared, &
                                                                                        azimuth_of_mid_defocus)
                number_of_extrema_profile(current_bin) = ctf%CountNumberOfExtremaBeforeSquaredSpatialFrequency( &
                                                                        current_spatial_frequency_squared,azimuth_of_mid_defocus)
            enddo

            ! Now we can loop over the power spectrum again and decide to which bin to add it
            do j=1,self%logical_dimensions(2)
                do i=1,self%logical_dimensions(1)
                    ctf_diff_from_current_bin = huge(1.0e0)
                    chosen_bin = ComputePowerSpectrumBinNumber(number_of_bins,number_of_extrema_profile,number_of_extrema, &
                                                                i,j,ctf_values,ctf_values_profile)
                    average(chosen_bin) = average(chosen_bin) + self%real_values(i,j,1)
                    number_of_values(chosen_bin) = number_of_values(chosen_bin) + 1
                enddo
            enddo

            ! Do the actual averaging
            where (number_of_values .eq. 0)
                average = 0.0e0
            elsewhere
                average = average / real(number_of_values)
            endwhere
            if (squared_ctf_was_fit) then
                average_fit = ctf_values_profile**2
            else
                average_fit = abs(ctf_values_profile)
            endif
        endif

        ! Do the rescaling so that the peaks and troughs are at 0.0 and 1.0
        ! and compute CC and threshold to help evaluate quality of the fit
        ! \todo Move this to seperate routine
        !if (allocated(frc)) deallocate(frc)
        !if (allocated(frc_sigma)) deallocate(frc_sigma)
        !allocate(frc(number_of_bins),frc_sigma(number_of_bins))
        !allocate(background(number_of_bins),peak(number_of_bins))
        background = 0.0e0
        peak = 1.0e0
        frc = 1.0e0
        frc_sigma = 0.0e0
        location_of_previous_maximum = 1
        location_of_previous_minimum = 1
        current_maximum_number= 0
        at_a_maximum = .false.
        at_a_minimum = .true.  ! Note, this may not be true if we have the perfect phase plate
        normalisation_bin_number = 1

        if (.not. self_is_blank) then
            !
            do current_bin=2,number_of_bins-1
                ! Remember where we were before - minimum, maximum or neither
                maximum_at_previous_bin = at_a_maximum
                minimum_at_previous_bin = at_a_minimum
                ! Are we at a CTF min or max?
                at_a_minimum =      average_fit(current_bin) .le. average_fit(current_bin-1) &
                            .and.   average_fit(current_bin) .le. average_fit(current_bin+1)
                at_a_maximum =      average_fit(current_bin) .ge. average_fit(current_bin-1) &
                            .and.   average_fit(current_bin) .ge. average_fit(current_bin+1)
                ! It could be that the CTF is constant in this region, in which case we stay at a minimum if we were there
                if (at_a_maximum .and. at_a_minimum) then
                    at_a_minimum = minimum_at_previous_bin
                    at_a_maximum = maximum_at_previous_bin
                endif
                ! Fill in values for the background or peak by linear interpolation
                if (at_a_minimum) then
                    if (rescale_average .and. .not. self_is_blank) then
                        do i=location_of_previous_minimum+1,current_bin
                            ! linear interpolation of average values at the peaks and troughs of the CTF
                            background(i) =   average(location_of_previous_minimum) &
                                            * real(current_bin-i)/real(current_bin-location_of_previous_minimum)&
                                            + average(current_bin) &
                                            * real(i-location_of_previous_minimum)/real(current_bin-location_of_previous_minimum)
                        enddo
                    endif
                    location_of_previous_minimum = current_bin
                endif
                if (at_a_maximum) then
                    if (.not. maximum_at_previous_bin .and. average_fit(current_bin) .gt. 0.7) then
                        current_maximum_number = current_maximum_number + 1
                    endif
                    if (rescale_average .and. .not. self_is_blank) then
                        do i=location_of_previous_maximum+1,current_bin
                            ! linear interpolation of average values at the peaks and troughs of the CTF
                            peak(i)       =   average(location_of_previous_maximum) &
                                            * real(current_bin-i)/real(current_bin-location_of_previous_maximum)&
                                            + average(current_bin) &
                                            * real(i-location_of_previous_maximum)/real(current_bin-location_of_previous_maximum)
                            if (current_maximum_number .eq. rescale_based_on_maximum_number) then
                                normalisation_bin_number = current_bin
                            endif
                        enddo
                    endif
                    ! Let's do some scoring
                    frc(location_of_previous_maximum+1:current_bin) =  &
                                    NormalizedCrossCorrelation(average    (location_of_previous_maximum+1:current_bin), &
                                                               average_fit(location_of_previous_maximum+1:current_bin))
                    frc_sigma(location_of_previous_maximum+1:current_bin) = 2.0d0 &
                                                            /sqrt(real(current_bin-location_of_previous_maximum,kind=8))
                    ! Keep track of previous maximum
                    location_of_previous_maximum = current_bin
                endif
                if (at_a_maximum .and. at_a_minimum) then
                   call this_program%TerminateWithFatalError('ComputeRotationalAverageOfPowerSpectrum', &
                                                             'At a minimum and a maximum simultaneously')
                endif
            enddo
        endif

        if (rescale_input .and. .not. self_is_blank) then
            call background_image%Allocate(mould=self)
            call peak_image%Allocate(mould=self)
            actually_do_rescaling = (peak(normalisation_bin_number) - background(normalisation_bin_number)) .gt. 0.0 ! if not, could be a sign something went wrong
            ! Now we can loop over the power spectrum again and decide to which bin to add it
            do j=1,self%logical_dimensions(2)
                do i=1,self%logical_dimensions(1)
                    chosen_bin = ComputePowerSpectrumBinNumber(number_of_bins,number_of_extrema_profile,number_of_extrema, &
                                                                i,j,ctf_values,ctf_values_profile)

                    if (actually_do_rescaling) then
                        self%real_values(i,j,1) = (self%real_values(i,j,1) - background(chosen_bin)) &
                                                / min(max(peak(chosen_bin)-background(chosen_bin),0.2),5.0)
                    endif
                    ! We want the peaks to reach at least 0.1
                    ! Todo: only do this up to frequency where we are fitting
                    !if ((peak(chosen_bin)-background(chosen_bin))/(peak(normalisation_bin_number)-background(normalisation_bin_number)) .lt. 0.1e0 &
                    !                .and. abs(peak(chosen_bin)-background(chosen_bin)) .gt. 0.01) then
                    !    self%real_values(i,j,1) = self%real_values(i,j,1) / (peak(chosen_bin)-background(chosen_bin)) * (peak(normalisation_bin_number) &
                    !             - background(normalisation_bin_number)) * 0.1e0
                    !endif
                    rad = sqrt( real(i-self%GetPhysicalIndexOfBoxCenter(1))**2 &
                             +real(j-self%GetPhysicalIndexOfBoxCenter(2))**2)
                    if ( rad .gt. 20.0 .and. rad .lt. 255.0) then
                        background_image%real_values(i,j,1) = background(chosen_bin)
                        peak_image%real_values(i,j,1) = peak(chosen_bin)
                    else
                        background_image%real_values(i,j,1) = 0.0
                        peak_image%real_values(i,j,1) = 0.0
                    endif
                enddo
            enddo
        endif


        !if (rescale_average) average = (average - background) / (peak - background) * 0.95e0 ! all peaks go to 1.0
        !if (rescale_average) average = (average - background) ! only do background subtraction
        if (rescale_average .and. .not. self_is_blank) then
            if ((peak(normalisation_bin_number) - background(normalisation_bin_number)) .gt. 0.0) then
                average = (average - background) &
                        / (peak(normalisation_bin_number) - background(normalisation_bin_number) ) &
                        * 0.95d0
                ! We want the peaks to reach at least 0.1
                where ((peak-background)/(peak(normalisation_bin_number)-background(normalisation_bin_number)) .lt. 0.1e0 &
                                .and. abs(peak-background) .gt. 0.000001)
                    average = average / (peak-background) * (peak(normalisation_bin_number) &
                             - background(normalisation_bin_number)) * 0.1e0
                endwhere
            else
                write(*,'(a)') '**warning(Image::ComputeRotationalAverageOfPowerSpectrum): '//&
                                'unable to rescale average experimental power'
            endif
        endif

    end subroutine ComputeRotationalAverageOfPowerSpectrum

    integer function ComputePowerSpectrumBinNumber(number_of_bins,number_of_extrema_profile,number_of_extrema, &
                                                    i,j,ctf_values,ctf_values_profile) result(chosen_bin)
        integer,        intent(in)  ::  number_of_bins
        real,           intent(in)  ::  number_of_extrema_profile(:)
        type(Image),    intent(in)  ::  number_of_extrema
        integer,        intent(in)  ::  i,j                         !<  Physical memory address
        type(Image),    intent(in)  ::  ctf_values
        real,           intent(in)  ::  ctf_values_profile(:)
        ! private variables
        integer     ::  current_bin
        real        ::  diff_number_of_extrema, diff_number_of_extrema_previous, diff_number_of_extrema_next
        real        ::  ctf_diff_from_current_bin
        real        ::  ctf_diff_from_current_bin_old
        ! Let's find the bin which has the same number of preceding extrema and the most similar ctf value
        ctf_diff_from_current_bin = huge(1.0e0)
        chosen_bin = 0
        do current_bin=1,number_of_bins
            diff_number_of_extrema  = abs(number_of_extrema%real_values(i,j,1) - number_of_extrema_profile(current_bin))
            if (current_bin .gt. 1) then
                diff_number_of_extrema_previous = abs(number_of_extrema%real_values(i,j,1) &
                                                    - number_of_extrema_profile(current_bin-1))
            else
                diff_number_of_extrema_previous = huge(1.0e0)
            endif
            if (current_bin .lt. number_of_bins) then
                diff_number_of_extrema_next     = abs(number_of_extrema%real_values(i,j,1) &
                                                    - number_of_extrema_profile(current_bin+1))
            else
                diff_number_of_extrema_next = huge(1.0e0)
            endif
            if (number_of_extrema%real_values(i,j,1) .gt. number_of_extrema_profile(number_of_bins)) then
                chosen_bin = number_of_bins
            else
                if (        diff_number_of_extrema .le. 0.01 &
                    .or.    (     diff_number_of_extrema .lt. diff_number_of_extrema_previous &
                            .and. diff_number_of_extrema .le. diff_number_of_extrema_next &
                            .and. number_of_extrema_profile(max(current_bin-1,1)) &
                                .ne. number_of_extrema_profile(min(current_bin+1,number_of_bins))) &
                    ) then
                    ! We're nearly there
                    ! Let's look for the position of the nearest CTF value
                    ctf_diff_from_current_bin_old = ctf_diff_from_current_bin
                    ctf_diff_from_current_bin = abs(ctf_values%real_values(i,j,1) - ctf_values_profile(current_bin))
                    if (ctf_diff_from_current_bin .lt. ctf_diff_from_current_bin_old) then
                        chosen_bin = current_bin
                    endif
                endif
            endif
        enddo
        if (chosen_bin .eq. 0) then
            print *, number_of_extrema_profile
            print *, i, j, number_of_extrema%real_values(i,j,1), ctf_values%real_values(i,j,1)
            print *, diff_number_of_extrema, diff_number_of_extrema_previous, diff_number_of_extrema_next
            call this_program%TerminateWithFatalError('ComputeRotationalAverageOfPowerSpectrum','Could not find bin')
        endif
    end function ComputePowerSpectrumBinNumber

    !>  \brief  Apply CTF to an image
    subroutine ApplyCTF(self,ctf,phase_flip_only,minimum_absolute_ctf_value)
        use ContrastTransferFunctions
        ! arguments
        class(Image),                   intent(inout)   ::  self
        type(ContrastTransferFunction), intent(in)      ::  ctf
        logical,                        intent(in)      ::  phase_flip_only
        real,                           intent(in)      ::  minimum_absolute_ctf_value  !<  Give 0.0 to get actual CTF. If greater than 0.0, do not allow CTF values lower than this.
        ! private variables
        integer ::  i,j
        real    ::  i_logi,j_logi
        real    ::  i_logi_sq,j_logi_sq
        real    ::  current_spatial_frequency_squared
        real    ::  current_azimuth
        real    ::  current_ctf_value
        ! start work

        ! Meant for images only
        if (self%IsAVolume()) then
            call this_program%TerminateWithFatalError('Image::ApplyCTF','Image must be two-dimensional')
        endif

        ! Make sure in Fourier space
        if (self%IsInRealSpace()) call self%ForwardFFT()

        !
        ! Let's loop over half of the image (ignore Friedel mates)
        do j=1,self%physical_upper_bound_complex(2)
            j_logi = self%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2) * self%fourier_voxel_size(2)
            j_logi_sq = j_logi**2
            do i=1,self%physical_upper_bound_complex(1)
                i_logi = real(i-1) * self%fourier_voxel_size(1)
                i_logi_sq = i_logi**2
                ! Where are we?
                current_spatial_frequency_squared = j_logi_sq+i_logi_sq
                if (current_spatial_frequency_squared .gt. 0.0) then
                    current_azimuth = atan2(j_logi,i_logi)
                else
                    current_azimuth = 0.0
                endif
                current_ctf_value = ctf%EvaluateAtSquaredSpatialFrequency(current_spatial_frequency_squared,current_azimuth)
                if (phase_flip_only) then
                    current_ctf_value = sign(1.0,current_ctf_value)
                else if (minimum_absolute_ctf_value .gt. 0.0) then
                    if (abs(current_ctf_value) .lt. minimum_absolute_ctf_value) then
                        current_azimuth = sign(minimum_absolute_ctf_value,current_ctf_value)
                    endif
                endif
                self%complex_values(i,j,1) = self%complex_values(i,j,1) * current_ctf_value
            enddo
        enddo

    end subroutine ApplyCTF


    !>  \brief  Overlay a CTF function on half of the image, a la CTFFIND
    subroutine OverlayCTF(self,ctf,squared_ctf)
        use ContrastTransferFunctions
        ! arguments
        class(Image),                   intent(inout)   ::  self        !<  Filtered power spectrum
        type(ContrastTransferFunction), intent(inout)   ::  ctf
        logical,                        intent(in)      ::  squared_ctf !<  Whether to correlated with the squared CTF (or just its absolute)
        ! start work
        call self%CTFOperation(ctf,overlay_ctf=.true.,squared_ctf=squared_ctf)
    end subroutine OverlayCTF

    !>  \brief   Compute the cross-correlation between an image (assumed to be a filtered power spectrum) and a CTF
    real function GetCorrelationWithCTF(self,ctf,squared_ctf) result(correlation)
        use ContrastTransferFunctions
        ! arguments
        class(Image),                   intent(inout)   ::  self        !<  Filtered power spectrum
        type(ContrastTransferFunction), intent(inout)   ::  ctf
        logical,                        intent(in)      ::  squared_ctf !<  Whether to correlated with the squared CTF (or just its absolute)
        ! start work
        call self%CTFOperation(ctf,correlation=correlation,squared_ctf=squared_ctf)
    end function GetCorrelationWithCTF


    !>  \brief  Compute the cross-correlation between an image (assumed to be a filtered power spectrum) and a CTF
    !!  \todo   Find a better name for this method
    subroutine CTFOperation(self,ctf,correlation,overlay_ctf,squared_ctf,correlation_terms,astigmatism_penalty)
        use ContrastTransferFunctions
        use EmpiricalDistributions
        ! arguments
        class(Image),                   intent(inout)   ::  self                    !<  Filtered power spectrum
        type(ContrastTransferFunction), intent(in)      ::  ctf
        real,           optional,       intent(out)     ::  correlation             !<  Normalised cross-correlation between the spectrum and the CTF
        logical,        optional,       intent(in)      ::  overlay_ctf             !<  Overwrite one half of the image with the CTF
        logical,        optional,       intent(in)      ::  squared_ctf             !<  Correlate with, or overlay, the squared CTF rather than the absolute of the CTF
        real,           optional,       intent(out)     ::  correlation_terms(3)    !<  Terms for the computation of the normalised cross-correlation. 1- numerator (cross-product); 2- denom1 (variance of spectrum); 3- denom2 (variance of fit ctf)
        real,           optional,       intent(out)     ::  astigmatism_penalty     !<  Astigmatism restraint penalty
        ! private variables
        integer         ::  i,j
        real(kind=4)    ::  i_logi,j_logi
        real(kind=4)    ::  i_logi_sq,j_logi_sq
        real(kind=4)    ::  current_spatial_frequency_squared, current_azimuth
        real(kind=4)    ::  current_ctf_squared, current_real_value
        real(kind=4)    ::  cross_product, norm_image, norm_ctf
        integer         ::  number_of_values
        logical         ::  ctffind_overlay
        real(kind=4)    ::  inverse_logical_dimensions(3)
        real(kind=4)    ::  lowest_freq, highest_freq
        type(EmpiricalDistribution) ::  values_in_rings, values_in_fitting_range
        real            ::  target_sigma
        logical         ::  ignore_central_cross_when_scoring = .true.
        integer,    parameter   ::  central_cross_half_width = 10
        logical         ::  use_ctf_squared
        real            ::  astig_penalty
        ! start work

        ! The image must be real and two-dimensional
        if (.not. self%is_in_real_space) then
            call this_program%TerminateWithFatalError('Image::CTFOperation','Image must be in real space')
        endif
        if (self%IsAVolume()) then
            call this_program%TerminateWithFatalError('Image::CTFOperation','Image must be two-dimensional')
        endif

        ! The CTF must have the correct units
        if (.not. ctf%HasPixelUnits()) then
            call this_program%TerminateWithFatalError('Image::CTFOperation','CTF must be in pixels and radians')
        endif

        !
        ctffind_overlay = .false.
        if (present(overlay_ctf)) ctffind_overlay = overlay_ctf

        !
        use_ctf_squared = .false.
        if (present(squared_ctf)) use_ctf_squared = squared_ctf

        ! Init
        cross_product       = 0.0d0
        norm_image          = 0.0d0
        norm_ctf            = 0.0d0
        number_of_values    = 0
        inverse_logical_dimensions = 1.0 / real(self%logical_dimensions)
        lowest_freq         =   ctf%GetLowestFrequencyForFitting()**2
        highest_freq        =   ctf%GetHighestFrequencyForFitting()**2
        if(present(correlation_terms)) correlation_terms = 0.0e0
        !if (ctffind_overlay) then
        !    ! We don't want to mask the overlay
        !    lowest_freq = 0.0
        !    highest_freq = 0.5
        !endif
        if (ctffind_overlay) then
            self%real_values = self%real_values - self%GetAverageOfValuesOnEdges()
            call values_in_rings%Init()
            call values_in_fitting_range%Init()
        endif


        ! Let's loop over half of the image (ignore Friedel mates)
        do j=1,self%logical_dimensions(2)
            j_logi = real(j-self%physical_address_of_box_center(2))*inverse_logical_dimensions(2)
            j_logi_sq = j_logi**2
            do i=1,self%physical_address_of_box_center(1)
                i_logi = real(i-self%physical_address_of_box_center(1))*inverse_logical_dimensions(1)
                i_logi_sq = i_logi**2
                ! Where are we?
                current_spatial_frequency_squared = j_logi_sq+i_logi_sq

                if (current_spatial_frequency_squared .gt. lowest_freq .and. &
                    current_spatial_frequency_squared .le. highest_freq) then

                    ! Here is the current value - we cast it from 32 bit to 64 bit
                    current_real_value = self%real_values(i,j,1)

                    if (current_spatial_frequency_squared .gt. 0.0) then
                        current_azimuth = atan2(j_logi,i_logi)
                    else
                        current_azimuth = 0.0
                    endif
                    ! Evaluate the CTF at this point
                    if (use_ctf_squared) then
                        current_ctf_squared &
                                  = ctf%EvaluateAtSquaredSpatialFrequency(current_spatial_frequency_squared,current_azimuth)**2
                    else
                        current_ctf_squared &
                                  = abs(ctf%EvaluateAtSquaredSpatialFrequency(current_spatial_frequency_squared,current_azimuth))
                    endif

                    ! Accumulate results
                    if (.not. ignore_central_cross_when_scoring .or. &
                        (      i .lt. self%physical_address_of_box_center(1)-central_cross_half_width &
                       .and.  (j .lt. self%physical_address_of_box_center(2)-central_cross_half_width &
                       .or.    j .gt. self%physical_address_of_box_center(2)+central_cross_half_width)) ) then
                        number_of_values = number_of_values + 1
                        cross_product = cross_product + current_real_value * current_ctf_squared
                        norm_image    = norm_image    + current_real_value**2
                        norm_ctf      = norm_ctf      + current_ctf_squared**2
                    endif
                    ! Overlay image
                    if (ctffind_overlay) then
                        if (use_ctf_squared) then
                            if (current_ctf_squared .gt. 0.25) then
                                call values_in_rings%AddSampleValue(self%real_values(i,j,1))
                            endif
                        else
                            if (current_ctf_squared .gt. 0.5) then
                                call values_in_rings%AddSampleValue(self%real_values(i,j,1))
                            endif
                        endif
                        call values_in_fitting_range%AddSampleValue(self%real_values(i,j,1))
                        if (j .lt. self%physical_address_of_box_center(2)) then
                            self%real_values(i,j,1) = current_ctf_squared
                        endif
                    endif
                !else
                !    if (ctffind_overlay) then
                !        self%real_values(i,j,1) = 0.0e0
                !        self%real_values(self%logical_dimensions(1)-i+1,j,1) = 0.0e0
                !    endif
                endif
                if (ctffind_overlay) then
                    if (current_spatial_frequency_squared .le. lowest_freq) then
                        self%real_values(i,j,1) = 0.0e0
                        self%real_values(self%logical_dimensions(1)-i+1,j,1) = 0.0e0
                    endif
                endif

            enddo
        enddo ! end of loop over half of the image

        !  Normalise the experimental power spectrum for display purposes
        if (ctffind_overlay) then
            !
            target_sigma = sqrt(values_in_rings%GetSampleVariance())
            do j=1,self%logical_dimensions(2)
                j_logi = real(j-self%physical_address_of_box_center(2))*inverse_logical_dimensions(2)
                j_logi_sq = j_logi**2
                do i=1,self%logical_dimensions(1)
                    i_logi = real(i-self%physical_address_of_box_center(1))*inverse_logical_dimensions(1)
                    i_logi_sq = i_logi**2
                    ! Where are we?
                    current_spatial_frequency_squared = j_logi_sq+i_logi_sq
                    if (i .gt. self%physical_address_of_box_center(1) .or. j .gt. self%physical_address_of_box_center(2)) then
                        ! In the experimental part of the diagnostic image
                        if (target_sigma .eq. 0.0) then
                            ! Blank spectrum, no rescaling to be done
                        else
                            self%real_values(i,j,1) = self%real_values(i,j,1) / target_sigma
                        endif
                    endif
                    !self%real_values(self%logical_dimensions(1)-i+1,self%logical_dimensions(2)-j+1,1) = &
                    !self%real_values(self%logical_dimensions(1)-i+1,self%logical_dimensions(2)-j+1,1) / target_sigma
                enddo
            enddo
            ! Loop over the outside of the theoretical CTF
            do j=1,self%physical_address_of_box_center(2)-1
                j_logi = real(j-self%physical_address_of_box_center(2))*inverse_logical_dimensions(2)
                j_logi_sq = j_logi**2
                do i=1,self%physical_address_of_box_center(1)
                    i_logi = real(i-self%physical_address_of_box_center(1))*inverse_logical_dimensions(1)
                    i_logi_sq = i_logi**2
                    ! Where are we?
                    current_spatial_frequency_squared = j_logi_sq+i_logi_sq
                    if (current_spatial_frequency_squared .gt. highest_freq) then
                        if (target_sigma .ne. 0.0) then
                            self%real_values(i,j,1) = self%real_values(i,j,1) / target_sigma
                            !self%real_values(self%logical_dimensions(1)-i+1,self%logical_dimensions(2)-j+1,1) = &
                            !self%real_values(self%logical_dimensions(1)-i+1,self%logical_dimensions(2)-j+1,1) / target_sigma
                        endif
                    endif
                enddo
            enddo
        endif

        ! The correlation terms
        if (present(correlation_terms)) then
            correlation_terms(1) = cross_product
            correlation_terms(2) = norm_image
            correlation_terms(3) = norm_ctf
        endif

        ! The astig penalty
        if (present(correlation) .or. present(correlation_terms)) then
            ! Restraint against crazy levels of astigmatism
            if (ctf%GetAstigmatismTolerance() .gt. 0.0) then
                astig_penalty = (ctf%GetAstigmatism())**2/2.0/(ctf%GetAstigmatismTolerance())**2/real(number_of_values)
            else
                astig_penalty = 0.0
            endif
            if (present(astigmatism_penalty)) astigmatism_penalty = astig_penalty
        endif

        ! The score
        if (present(correlation)) then
            correlation = cross_product / sqrt(norm_image*norm_ctf)
            correlation = correlation - astig_penalty
        endif
    end subroutine CTFOperation

    !>  \brief  Calculate an fsc with another image, and return it as a curve
    function GetFSCWith(self, other_image, number_of_bins_to_use, pixel_size_to_use,fsc_part_curve,fraction_of_voxels) &
                        result(fsc_curve)

        use Curves
        ! arguments

        class(Image),               intent(inout)   ::  self
        type(Image),                intent(inout)   ::  other_image

        integer,        optional,   intent(in)      ::  number_of_bins_to_use
        real,           optional,   intent(in)      ::  pixel_size_to_use
        type(Curve),    optional,   intent(inout)   ::  fsc_part_curve          !<  FSC curve with compensation for the fact that not all non-masked voxels are part of the object
        real,           optional,   intent(in)      ::  fraction_of_voxels      !<  Fraction of voxels which were not masked in the two volumes and which describe the object of interest (rather than solvent)

        ! private variables
        type(Curve)                             ::  fsc_curve

        integer                                 ::  used_number_of_bins

        integer                                 ::  bin_counter
        integer                                 ::  current_bin

        integer                                 ::  i
        integer                                 ::  j
        integer                                 ::  k

        real                                    ::  start_radius
        real                                    ::  end_radius
        real                                    ::  current_radius
        real                                    ::  bin_size

        real                                    ::  x
        real                                    ::  y
        real                                    ::  z

        integer                                 ::  j_logical
        integer                                 ::  k_logical

        logical                                 ::  must_fft
        logical                                 ::  have_pixel_size

        real, allocatable                       ::  fsc_data(:)
        real, allocatable                       ::  fsc_part_data(:)
        real, allocatable                       ::  first_denom(:)
        real, allocatable                       ::  second_denom(:)
        real, allocatable                       ::  resolution_data(:)

        real                                    ::  current_denom
        real                                    ::  ssnr


        ! do some checks

        if (.not. self%IsInSameSpaceAs(other_image)) then
            call this_program%TerminateWithFatalError('Image::GetFSCWith',&
                                                        'Images are in different spaces')
        endif
        ! check both images have same dimensions

        if (.not. self%HasSameDimensionsAs(other_image)) then
            call this_program%TerminateWithFatalError('Image::GetFSCWith',&
                                                        'Images have different dimensions')
        endif

        ! we may need to know the fraction of non-masked voxels in the object
        if (present(fsc_part_curve)) then
            if (.not. present(fraction_of_voxels)) then
                call this_program%TerminateWithFatalError('Image::GetFSCWith', &
                                                        'Need to know fraction of voxels in object')
            endif
        endif

        ! Work out the number of bins..

        if (present(number_of_bins_to_use)) then
            used_number_of_bins = number_of_bins_to_use
        else
            used_number_of_bins = self%GetMaximumRadius()
        endif

        if (present(pixel_size_to_use)) then
            have_pixel_size = .true.
        else
            have_pixel_size = .false.
        endif

        bin_size = 1. / used_number_of_bins

        ! allocate

        allocate(fsc_data(used_number_of_bins))
        allocate(fsc_part_data(used_number_of_bins))
        allocate(resolution_data(used_number_of_bins))
        allocate(first_denom(used_number_of_bins))
        allocate(second_denom(used_number_of_bins))

        must_fft = .false.

        if (self%IsInRealSpace()) then
            must_fft = .true.
            call self%ForwardFFT()
            call other_image%ForwardFFT()
        endif


        ! blank the arrays, and calculate resolution info

        do bin_counter = 1, used_number_of_bins

            fsc_data(bin_counter) = 0.0e0
            first_denom(bin_counter) = 0.0e0
            second_denom(bin_counter) = 0.0e0

            start_radius = real(bin_counter - 1) * real(bin_size)

            if (bin_counter .eq. used_number_of_bins) then
                end_radius = 1.0e0
            else
                end_radius = real(bin_counter) * real(bin_size)
            endif


            resolution_data(bin_counter) = ((start_radius + end_radius) / 2) * real(self%physical_upper_bound_complex(1) - 1) &
                                            * self%fourier_voxel_size(1)
            if (have_pixel_size) resolution_data(bin_counter) = 1. / (pixel_size_to_use / resolution_data(bin_counter))

        enddo

        ! Perform the actual calculation


        do k=1,self%logical_dimensions(3)
            k_logical = self%LogicalIndexGivenPhysicalIndexInFourierSpace(k, 3)
            z = real(k_logical) * self%fourier_voxel_size(3)
            z = z**2

            do j=1,self%logical_dimensions(2)
                j_logical = self%LogicalIndexGivenPhysicalIndexInFourierSpace(j, 2)
                y = real(j_logical)  * self%fourier_voxel_size(2)
                y = y**2

                do i=1,self%physical_upper_bound_complex(1)
                    x = real(i - 1) * self%fourier_voxel_size(1)
                    x = x**2

                    current_radius = sqrt(x + y + z)

                    !write(*,*) x, ',',y,',',z,' = ', current_radius

                    if (current_radius .lt. 0.5) then

                        ! work out which bin this radius is in..

                        current_bin = current_radius * real(used_number_of_bins) * 2.0e0
                        current_bin = current_bin + 1

                        fsc_data(current_bin) = fsc_data(current_bin) + real(self%complex_values(i, j, k) * &
                                                                        conjg(other_image%complex_values(i, j, k)))
                        first_denom(current_bin) = first_denom(current_bin) + abs(self%complex_values(i, j, k))**2
                        second_denom(current_bin) = second_denom(current_bin) + abs(other_image%complex_values(i, j, k))**2
                    endif

                enddo
            enddo
        enddo

        ! normalisation

        do bin_counter = 1, used_number_of_bins

            current_denom = sqrt(first_denom(bin_counter) * second_denom(bin_counter))

            if (current_denom .eq. 0.0e0) then
                fsc_data(bin_counter) = 0.0
            else
                fsc_data(bin_counter) = fsc_data(bin_counter) / current_denom
            endif

        enddo

        ! fsc part (following Frealign 9.09's implementation of Sindelar (2012))
        if (present(fsc_part_curve)) then
            do bin_counter = 1, used_number_of_bins
                ssnr = 2.0*fsc_data(bin_counter)/(1.0-fsc_data(bin_counter)+0.0001)
                if (fsc_data(bin_counter) .ge. 1.0) then
                    fsc_part_data(bin_counter) = 1.0
                else if (fraction_of_voxels .ne. 0.0) then
                    ssnr = ssnr/fraction_of_voxels
                    fsc_part_data(bin_counter) = ssnr/(abs(ssnr)+2.0)
                endif
            enddo
        endif

        if (must_fft) then
            call self%BackwardFFT()
            call other_image%BackwardFFT()
        endif

        ! copy into the output curve
        call fsc_curve%Init(size(resolution_data),'Spatial frequency','FSC')
        call fsc_curve%ClearData()
        call fsc_curve%AddPoints(resolution_data, fsc_data)

        !
        if (present(fsc_part_curve)) then
            call fsc_part_curve%Init(size(resolution_data),'Spatial frequency','FSC Part.')
            call fsc_part_curve%ClearData()
            call fsc_part_curve%AddPoints(resolution_data,fsc_part_data)
        endif

    end function GetFSCWith

    !>  \brief  In Fourier space, normalise voxel values such that the standard deviation of complex values within each resolution zone is 1.0
    subroutine ComputeSigmaOfResolutionZones(self,width_of_resolution_zones,sigma,number_of_terms)
        ! Arguments
        class(Image),                   intent(inout)   ::  self                        !<  Image or volume
        real,                           intent(in)      ::  width_of_resolution_zones   !<  In reciprocal pixels (1/pixels; 0.5 would mean there is only one resolution zone between the origin and Nyquist frequency)
        real(kind=8),                   intent(inout)   ::  sigma(:)                    !<  For each resolution zone, the standard deviations of the complex Fourier components (computed as the sqrt of the average squared amplitude)
        real(kind=8),                   intent(inout)   ::  number_of_terms(:)          !<  For each resolution zone, the number of terms which contributed to the statistics (Hermitian symmetry mates are not included)
        ! Private variables
        integer                         ::  address_x
        integer                         ::  address_y
        real                            ::  sq_spatial_freq_x
        real                            ::  sq_spatial_freq_y
        real                            ::  current_spatial_frequency
        integer                         ::  current_resolution_zone
        ! Start work

        ! Only meant for 2D images
        if (self%IsAVolume()) call this_program%TerminateWithFatalError('Image::ComputeSigmaOfResolutionZones',&
                                                                        'Not implemented for volumes')

        ! Make sure we're in Fourier space
        if (self%IsInRealSpace()) call self%ForwardFFT()


        ! Initialize
        number_of_terms = 0.0d0
        sigma           = 0.0d0


        ! First loop over image, to compute normalization factors
        do address_y=1,self%physical_upper_bound_complex(2)
            sq_spatial_freq_y = (self%LogicalIndexGivenPhysicalIndexInFourierSpace(address_y,2)*self%fourier_voxel_size(2))**2
            do address_x=1,self%physical_upper_bound_complex(1)
                sq_spatial_freq_x = ((address_x-1)*self%fourier_voxel_size(1))**2

                ! Some Fourier components have both Friedel mates explicitely described by the FFTW representation.
                ! We avoid counting them twice.
                if (self%FourierComponentIsRedundant(address_x,address_y,1)) cycle

                ! Current spatial frequency
                current_spatial_frequency = sqrt(sq_spatial_freq_x+sq_spatial_freq_y)

                ! Work out which resolution zone we are in
                current_resolution_zone = self%ResolutionZoneFromSpatialFrequency(current_spatial_frequency, &
                                                                                  width_of_resolution_zones)

                ! Accumulate stats
                number_of_terms(current_resolution_zone)    = number_of_terms(current_resolution_zone) + 1.0d0
                sigma(current_resolution_zone)              = sigma(current_resolution_zone)  &
                                                            + real(self%complex_values(address_x,address_y,1))**2 &
                                                            + aimag(self%complex_values(address_x,address_y,1))**2
            enddo
        enddo

        ! The standard deviation of a complex variable can be approximated as
        ! the sqrt of the average squared amplitude
        ! See derivation in notes 22-March-2015
        sigma = sqrt(sigma / number_of_terms)

    end subroutine ComputeSigmaOfResolutionZones

    pure integer function NumberOfResolutionZonesFromTheirWidth(self,width_of_resolution_zones) result(number_of_resolution_zones)
        class(Image),   intent(in)  ::  self
        real,           intent(in)  ::  width_of_resolution_zones   !<  In reciprocal pixels (1/pixels; 0.5 would mean there is only one resolution zone between the origin and Nyquist frequency)
        ! Start work
        number_of_resolution_zones = ceiling(self%GetMaximumDiagonalRadius() / width_of_resolution_zones)
    end function NumberOfResolutionZonesFromTheirWidth

    pure integer function ResolutionZoneFromSpatialFrequency(current_spatial_frequency,width_of_resolution_zones)
        real,           intent(in)  ::  current_spatial_frequency   !<  Spatial frequency for the Fourier component of interest (1/pixels)
        real,           intent(in)  ::  width_of_resolution_zones   !<  In reciprocal pixels (1/pixels; 0.5 would mean there is only one resolution zone between the origin and Nyquist frequency)
        ! Start work
        ResolutionZoneFromSpatialFrequency = int(current_spatial_frequency/width_of_resolution_zones) + 1
    end function ResolutionZoneFromSpatialFrequency


    !>  \brief  In Fourier space, normalise voxel values such that the standard deviation of complex values within each resolution zone is 1.0
    subroutine NormaliseSigmaOfResolutionZones(self,width_of_resolution_zones,number_of_terms)
        ! Arguments
        class(Image),                   intent(inout)   ::  self                        !<  Image or volume
        real,                           intent(in)      ::  width_of_resolution_zones   !<  In reciprocal pixels (1/pixels; 0.5 would mean there is only one resolution zone between the origin and Nyquist frequency)
        real(kind=8),                   intent(inout)   ::  number_of_terms(:)          !<  For each resolution zone, the number of terms which contributed to the statistics (Hermitian symmetry mates are not included)
        ! Private variables
        integer                         ::  address_x
        integer                         ::  address_y
        real                            ::  sq_spatial_freq_x
        real                            ::  sq_spatial_freq_y
        real                            ::  current_spatial_frequency
        integer                         ::  current_resolution_zone
        real(kind=8)                    ::  sigma(size(number_of_terms))
        ! Start work

        ! Only meant for 2D images
        if (self%IsAVolume()) call this_program%TerminateWithFatalError('Image::NormaliseSigmaOfResolutionZones',&
                                                                        'Not implemented for volumes')

        ! Make sure we're in Fourier space
        if (self%IsInRealSpace()) call self%ForwardFFT()

        ! Determine the sigma of complex values within each resolution zone
        call self%ComputeSigmaOfResolutionZones(width_of_resolution_zones,sigma,number_of_terms)

        ! Loop over image to apply the normalization
        do address_y=1,self%physical_upper_bound_complex(2)
            sq_spatial_freq_y = (self%LogicalIndexGivenPhysicalIndexInFourierSpace(address_y,2)*self%fourier_voxel_size(2))**2
            do address_x=1,self%physical_upper_bound_complex(1)
                sq_spatial_freq_x = ((address_x-1)*self%fourier_voxel_size(1))**2

                ! Current spatial frequency
                current_spatial_frequency = sqrt(sq_spatial_freq_x+sq_spatial_freq_y)

                ! Work out which resolution zone we are in
                current_resolution_zone = self%ResolutionZoneFromSpatialFrequency(current_spatial_frequency, &
                                                                                  width_of_resolution_zones)

                ! Apply normalization
                self%complex_values(address_x,address_y,1) = self%complex_values(address_x,address_y,1) &
                                                            / sigma(current_resolution_zone)
            enddo
        enddo
    end subroutine NormaliseSigmaOfResolutionZones

    !>  \brief  Compute the similarity between two images, following the Frealign 9 (2014/5) algorithm
    !!
    !!  The arrays which hold values per resolution zone should be allocated to the correct dimension (use the NumberOfResolutionZonesFromTheirWidth method)
    function ComputeFrealign9SimilarityMeasure( self,reference_image,low_spatial_freq_for_absolute,width_of_resolution_zones, &
                                                reference_sigma,number_of_terms, &
                                                numerator,denominator_self,denominator_reference) &
                                            result(similarity_measure)
        use Curves
        use ContrastTransferFunctions
        ! Arguments
        class(Image),                   intent(in)      ::  self                            !<  Experimental image. Must be in Fourier space and must be normalised by resolution zones, and multiplied by w=(SSNR)/(1+SSNR), where SSNR = CTF^2*PSSNR
        type(Image),                    intent(inout)   ::  reference_image                 !<  Reference image (e.g. projection from reference 3D). Must be in Fourier space, not yet normalized in resolution zones.
        real,                           intent(in)      ::  low_spatial_freq_for_absolute   !<  At and above this frequency (1/pixels), take absolute of cross-term
        real,                           intent(in)      ::  width_of_resolution_zones       !<  In reciprocal pixels (1/pixels; 0.5 would mean there is only one resolution zone between the origin and Nyquist frequency)
        real(kind=8),                   intent(inout)   ::  reference_sigma(:)              !<  For each resolution zone, the std dev of the Fourier components in the reference image
        real(kind=8),                   intent(inout)   ::  number_of_terms(:)              !<  For each resolution zone, the number of Fourier components (excluding Hermitian symmetry mates)
        real(kind=8),                   intent(inout)   ::  numerator(:)                    !<  For each resolution zone, the conjugate cross product of the Fourier components from the image and the reference
        real(kind=8),                   intent(inout)   ::  denominator_self(:)             !<  For each resolution zone, the correlation coefficient normalizing term corresponding to the experimental image
        real(kind=8),                   intent(inout)   ::  denominator_reference(:)        !<  For each resolution zone, the correlation coefficient normalizing term corresponding to the reference image
        ! Result
        real                        ::  similarity_measure
        ! Private variables
        integer                     ::  address_x
        integer                     ::  address_y
        real                        ::  sq_spatial_freq_x
        real                        ::  sq_spatial_freq_y
        real                        ::  current_spatial_frequency
        integer                     ::  current_resolution_zone
        real(kind=8)                ::  numerator_sum
        real(kind=8)                ::  denominator_self_sum
        real(kind=8)                ::  denominator_reference_sum

        ! Start work

        ! Initialise
        similarity_measure      =   0.0
        numerator               =   0.0d0
        denominator_self        =   0.0d0
        denominator_reference   =   0.0d0

        ! Compute normalisation factor for reference Fourier zones
        call reference_image%ComputeSigmaOfResolutionZones(width_of_resolution_zones,reference_sigma,number_of_terms)


        ! Loop over image
        do address_y=1,self%physical_upper_bound_complex(2)
            sq_spatial_freq_y = (self%LogicalIndexGivenPhysicalIndexInFourierSpace(address_y,2)*self%fourier_voxel_size(2))**2
            do address_x=1,self%physical_upper_bound_complex(1)
                sq_spatial_freq_x = ((address_x-1)*self%fourier_voxel_size(1))**2

                ! Some Fourier components have both Friedel mates explicitely described by the FFTW representation.
                ! We avoid counting them twice.
                if (self%FourierComponentIsRedundant(address_x,address_y,1)) cycle

                ! Current spatial frequency
                current_spatial_frequency = sqrt(sq_spatial_freq_x+sq_spatial_freq_y)

                ! Work out which resolution zone we are in
                current_resolution_zone = self%ResolutionZoneFromSpatialFrequency(current_spatial_frequency, &
                                                                                  width_of_resolution_zones)

                ! Accumulate terms of the correlation coefficient for this resolution shell
                numerator(current_resolution_zone) = numerator(current_resolution_zone) &
                                                    + self%complex_values(address_x,address_y,1) &
                                                        *conjg(reference_image%complex_values(address_x,address_y,1))

                denominator_self(current_resolution_zone) =     denominator_self(current_resolution_zone) &
                                                            +   real (self%complex_values(address_x,address_y,1))**2 &
                                                            +   aimag(self%complex_values(address_x,address_y,1))**2

                denominator_reference(current_resolution_zone) =    denominator_reference(current_resolution_zone) &
                                                            +   real (reference_image%complex_values(address_x,address_y,1))**2 &
                                                            +   aimag(reference_image%complex_values(address_x,address_y,1))**2


            enddo
        enddo

        ! Compute a score
        numerator_sum = 0.0d0
        denominator_self_sum = 0.0d0
        denominator_reference_sum = 0.0d0
        do current_resolution_zone=1,size(numerator)
            if (current_resolution_zone*width_of_resolution_zones .lt. low_spatial_freq_for_absolute) then
                ! No absolute
                numerator_sum = numerator_sum + numerator(current_resolution_zone)/reference_sigma(current_resolution_zone)&
                                                                                  /sqrt(number_of_terms(current_resolution_zone))
            else
                ! Take the absolute of the cross term
                numerator_sum = numerator_sum + abs(numerator(current_resolution_zone)/reference_sigma(current_resolution_zone)&
                                                                                 /sqrt(number_of_terms(current_resolution_zone)))
            endif
            denominator_self_sum = denominator_self_sum + denominator_self(current_resolution_zone)
            denominator_reference_sum = denominator_reference_sum + denominator_reference(current_resolution_zone)&
                                                                    /number_of_terms(current_resolution_zone)
        enddo
        similarity_measure = numerator_sum / (sqrt(denominator_self_sum)*sqrt(denominator_self_sum))


    end function ComputeFrealign9SimilarityMeasure

    !>  \brief  Compute the probability of an experimental image given a matching projection, following the formalism described in Lyumkis et al (2013)
    !!          and implemented in Frealign 9, particularly sigma.f
    !!
    !!          We are basically computing the log of the quantity defined in Equation 2 of Lyumkis et al JSB 183 (2013) 377-388
    subroutine ComputeFrealign9Probability(self,reference_image,mask_image)
        class(Image),       intent(in)  ::  self                !<  Experimental image
        type(Image),        intent(in)  ::  reference_image     !<  Projection image, assumed not filtered, but already phase shifted to match the experimental image
        type(Image),        intent(in)  ::  mask_image          !<  Mask image. Only non-zero pixels will be used to compute statistics
        ! Private variables
        real            ::  reference_image_scale_factor    ! alpha in Niko's notation
        real            ::  self_noise_variance             ! sigma^2 in Niko's notation
        real(kind=8)    ::  reference_image_average
        real(kind=8)    ::  reference_image_sum_of_squares  ! sum1 in Niko's notation
        real(kind=8)    ::  self_sum_of_squares             ! sum2 in Niko's notation
        real(kind=8)    ::  self_average
        real(kind=8)    ::  cross_product ! ccpart in Niko's notation
        integer         ::  address_x
        integer         ::  address_y
        integer         ::  number_of_pixels_in_mask
        real(kind=8)    ::  self_deviation
        real(kind=8)    ::  reference_image_deviation
        ! Start work

        ! Initialise values
        self_average = 0.0d0
        reference_image_average = 0.0d0
        number_of_pixels_in_mask = 0
        cross_product = 0.0d0
        reference_image_sum_of_squares = 0.0d0
        self_sum_of_squares = 0.0d0

        ! Compute averages within mask
        do address_y=1,self%logical_dimensions(2)
            do address_x=1,self%logical_dimensions(1)
                if (mask_image%real_values(address_x,address_y,1) .gt. 0.0) then
                    self_average = self_average + self%real_values(address_x,address_y,1)
                    reference_image_average = reference_image_average + reference_image%real_values(address_x,address_y,1)
                    number_of_pixels_in_mask = number_of_pixels_in_mask + 1
                endif
            enddo
        enddo
        self_average = self_average / number_of_pixels_in_mask
        reference_image_average = reference_image_average / number_of_pixels_in_mask

        ! Compute sums of square deviations and cross-product
        do address_y=1,self%logical_dimensions(2)
            do address_x=1,self%logical_dimensions(1)
                if (mask_image%real_values(address_x,address_y,1) .gt. 0.0) then
                    self_deviation = self%real_values(address_x,address_y,1) - self_average
                    reference_image_deviation = reference_image%real_values(address_x,address_y,1) - reference_image_average
                    cross_product = cross_product + self_deviation*reference_image_deviation
                    reference_image_sum_of_squares = reference_image_sum_of_squares + reference_image_deviation**2
                    self_sum_of_squares = self_sum_of_squares + self_deviation**2
                endif
            enddo
        enddo

        ! Compute scaling factors
        reference_image_scale_factor = abs(cross_product/reference_image_sum_of_squares)
        self_noise_variance = abs(self_sum_of_squares-reference_image_scale_factor**2*reference_image_sum_of_squares) &
                            / number_of_pixels_in_mask

    end subroutine ComputeFrealign9Probability

    !>  \brief  Calculate a cross correlation image with another image, this can then be searched for peaks
    subroutine CalculateCrossCorrelationImageWith(self, other_image)

        ! arguments

        class(Image),           intent(inout)   ::  self
        type(image),            intent(inout)   ::  other_image

        ! private variables
        logical ::  must_fft

        ! start work

        ! check both images are of same type

        if (.not. self%IsInSameSpaceAs(other_image)) then
            call this_program%TerminateWithFatalError('Image::CalculateCrossCorrelationImageWith',&
                                                        'Images are in different spaces')
        endif
        ! check both images have same dimensions

        if (.not. self%HasSameDimensionsAs(other_image)) then
            call this_program%TerminateWithFatalError('Image::CalculateCrossCorrelationImageWith',&
                                                        'Images have different dimensions')
        endif

        ! check whether we'll have to do fts internally

        must_fft = self%is_in_real_space

        ! fourier transform both inputs if necessary

        if (must_fft) then
            call self%ForwardFFT()
            call other_image%ForwardFFT()
        endif

        ! multiply by complex conjugate

        self%complex_values = self%complex_values * conjg(other_image%complex_values)

        if (self%object_is_centered_in_box) then

            ! we need object_is_centered_in_box will to be false, so that the correct phase shift is applied, we (Tim & Alexis) think
            self%object_is_centered_in_box = .false.

            ! swap the Quadrants so that the peak of a zero shift will be in the center

            call self%SwapRealSpaceQuadrants()


        endif


        ! reverse ft the image

        call self%BackwardFFT()

        ! if necessary, return the other image to real space
        if (must_fft) then
            call other_image%BackwardFFT()
        endif


    end subroutine CalculateCrossCorrelationImageWith

    !>  \brief  Write an image to disk, given a filename
    subroutine WriteToDiskGivenFilename(self,filename,wanted_location,pixel_size,print_file_info,delete_if_already_exists)
        use ImageFiles
        ! Arguments
        class(Image),               intent(inout)   ::  self
        character(len=*),           intent(in)      ::  filename                        !<  Name of file to write image to
        integer,        optional,   intent(in)      ::  wanted_location                 !<  Which location (section) in the file should we write to?
        real,           optional,   intent(in)      ::  pixel_size                      !<  Pixel size (in Angstroms) to be written to the header
        logical,        optional,   intent(in)      ::  print_file_info                 !<  Print to the terminal a summary of the file's properties
        logical,        optional,   intent(in)      ::  delete_if_already_exists        !<  Remove the file first if it already exists
        ! Private variables
        type(ImageFile)                             ::  my_imagefile
        ! Start work

        ! Initialise an ImageFile object
        call my_imagefile%Init(filename,pixel_size=pixel_size,delete_if_already_exists=delete_if_already_exists)

        ! Write the data
        call self%WriteToImageFile(my_imagefile,wanted_location)

        ! Print some info
        if (present(print_file_info)) then
            if (print_file_info) then
                call my_imagefile%PrintInfo()
            endif
        endif

        ! The destructor should take care of closing the file etc. (but with older versions of the compiler it does not?)
        call my_imagefile%Close()

    end subroutine WriteToDiskGivenFilename


    !>  \brief  Write an image to disk, given an initialised ImageFile object
    subroutine WriteToImageFile(self,my_imagefile,wanted_location,print_file_info)
        use ImageFiles
        ! Arguments
        class(Image),               intent(inout)   ::  self
        class(ImageFile),           intent(inout)   ::  my_imagefile                    !<  An initialised ImageFile object
        integer,        optional,   intent(in)      ::  wanted_location                 !<  Which location (section) in the file should we write to?
        logical,        optional,   intent(in)      ::  print_file_info                 !<  Print to the terminal a summary of the file's properties
        ! Private variables
        integer                                     ::  first_slice
        integer                                     ::  last_slice
        integer                                     ::  wwanted_location
        logical                                     ::  do_ft
        ! Start work

        ! Which location
        wwanted_location = 1
        if (present(wanted_location)) wwanted_location = wanted_location

        ! What slices/locations will we write to
        if (self%IsAVolume()) then
            !if (wwanted_location .gt. 1) then
            !    call this_program%TerminateWithFatalError('Image::WriteToImageFile','Stacks of volumes not supported yet')
            !endif
            first_slice = wwanted_location
            last_slice = self%logical_dimensions(3) + wwanted_location - 1
        else
            first_slice = wwanted_location
            last_slice = wwanted_location
        endif

        ! We always write real-space images to disk
        do_ft = .not. self%IsInRealSpace()
        if (do_ft) call self%BackwardFFT()

        ! Check for any NaNs
        if (self%HasNan()) then
            write(*,'(2a)') 'Warning(Image::WriteToImageFile): at least one NaN found in image to be written to ', &
                            my_imagefile%GetFilename()
            !call this_program%TerminateWithFatalError('Image::WriteToImageFile','NaN in image')
        endif

        ! Write to disk
        call my_imagefile%WriteSlicesToDisk(first_slice,last_slice,self%real_values,self%logical_dimensions(1))

        ! Print some info
        if (present(print_file_info)) then
            if (print_file_info) then
                call my_imagefile%PrintInfo()
            endif
        endif

        ! We always return the image in the space we found it. This is important when
        ! dumping files to disk for debugging purposes
        if (do_ft) call self%ForwardFFT()

    end subroutine WriteToImageFile

    !>  \brief  Read an image from disk, given a filename
    subroutine ReadFromDiskGivenFilename(self,filename,wanted_location,read_volume_not_slice)
        use ImageFiles
        ! Arguments
        class(Image),               intent(inout)   ::  self
        character(len=*),           intent(in)      ::  filename                        !<  Name of the file to read image from
        integer,        optional,   intent(in)      ::  wanted_location                 !<  Which location (section) in the file should we read into memory?
        logical,        optional,   intent(in)      ::  read_volume_not_slice           !<  Read a volume, not a slice
        ! Private variables
        type(ImageFile)                             ::  my_imagefile
        ! Start work

        ! Initialise an ImageFile object
        call my_imagefile%Init(filename)

        ! Read the data from disk
        call self%ReadFromImageFile(my_imagefile,wanted_location,read_volume_not_slice)
    end subroutine ReadFromDiskGivenFilename

    !>  \brief  Read an image from disk, given an initialised ImageFile object
    subroutine ReadFromImageFile(self,my_imagefile,wanted_location,read_volume_not_slice)
        use ImageFiles
        ! Arguments
        class(Image),               intent(inout)   ::  self
        class(ImageFile),           intent(inout)   ::  my_imagefile                    !<  An initialised ImageFile object
        integer,        optional,   intent(in)      ::  wanted_location                 !<  Which location (section) in the file should we read into memory?
        logical,        optional,   intent(in)      ::  read_volume_not_slice           !<  Read a volume, not a slice
        ! Private variables
        logical                                     ::  rread_volume_not_slice
        integer                                     ::  new_image_dimensions(3)
        integer                                     ::  first_slice
        integer                                     ::  last_slice
        integer                                     ::  wwanted_location
        ! Start work

        ! Reading 3D or 2D?
        rread_volume_not_slice = .false.
        if (present(read_volume_not_slice)) rread_volume_not_slice = read_volume_not_slice

        ! Which location?
        wwanted_location = 1
        if (present(wanted_location)) wwanted_location = wanted_location

        ! The image will have these dimensions
        new_image_dimensions = my_imagefile%GetDimensions()
        if (.not. rread_volume_not_slice) new_image_dimensions(3) = 1

        ! Allocate image
        call self%Allocate(dims=new_image_dimensions)

        ! The slices we will read
        if (rread_volume_not_slice) then
            if (wwanted_location .gt. 1) then
                call this_program%TerminateWithFatalError('Image::ReadFromImageFile','Stacks of volumes not supported yet')
            endif
            first_slice = 1
            last_slice = my_imagefile%GetStackSize()
        else
            first_slice = wwanted_location
            last_slice = wwanted_location
        endif

        ! Do the reading from disk
        call my_imagefile%ReadSlicesFromDisk(first_slice,last_slice,self%real_values)

        ! We assume that we read in real data
        self%is_in_real_space = .true.
    end subroutine ReadFromImageFile


    !>  \brief  Look for a near-meridonial cross-beta layer line
    subroutine LookForCrossBetaLayerLine(   self,expected_frequency,expected_frequency_uncertainty,expected_half_width,    &
                                            best_frequency,best_score)
        use UsefulFunctions, only : IsEven
        use EmpiricalDistributions
        ! arguments
        class(image),               intent(inout)   ::  self
        real,                       intent(in)      ::  expected_frequency              !<  Frequency in 1/pixel at which we expect to find the layer line
        real,                       intent(in)      ::  expected_frequency_uncertainty  !<  Half-range for the layer line search (in 1/pixel)
        real,                       intent(in)      ::  expected_half_width             !<  Half-width of the layer line in reciprocal pixels
        real,                       intent(out)     ::  best_frequency                  !<  Frequency at which we found the strongest layer line
        real,                       intent(out)     ::  best_score                      !<  The difference between the average value in the layer line and the average
                                                                                        !!  value in the corresponding ring, as a multiple of the standard deviation
                                                                                        !!  in the ring
        ! private variables
        integer                                     ::  i,j
        integer                                     ::  i_logi,j_logi
        real                                        ::  x_freq_sq,y_freq_sq
        real                                        ::  freq
        integer                                     ::  number_of_rings
        integer                                     ::  current_ring
        integer                                     ::  expected_ring
        logical                                     ::  in_layer_line_region
        real                                        ::  ring_width
        real                                        ::  expected_half_width_sq
        real                                        ::  expected_frequency_min_sq, expected_frequency_max_sq
        type(EmpiricalDistribution),    allocatable ::  amplitude_distribution(:), meridional_amplitude_distribution(:)
        real,                           allocatable ::  layer_line_score(:), ring_frequency(:)
        real                                        ::  current_pixel_amplitude
        integer,                        allocatable ::  best_ring(:)
        ! start work

        ! For later, we need squares of some of the arguments
        expected_half_width_sq = expected_half_width**2
        expected_frequency_min_sq = (expected_frequency-expected_frequency_uncertainty)**2
        expected_frequency_max_sq = (expected_frequency+expected_frequency_uncertainty)**2

        ! This is written with 2D images in mind
        if (self%IsAVolume()) call this_program%TerminateWithFatalError('Image::LookForCrossBetaLayerLine','Only for 2D images')

        ! Get Fourier transform
        if (self%is_in_real_space) call self%ForwardFFT()
        ! How wide will each ring be? (in reciprocal pixels)
        ring_width = 1.0 / minval(self%logical_dimensions(1:2))
        ! How many rings will concern ourselves with?
        number_of_rings = ceiling(2.0*expected_frequency_uncertainty/ring_width)
        if (IsEven(number_of_rings)) number_of_rings = number_of_rings + 1
        expected_ring = (number_of_rings-1)/2 + 1 ! we expect the layer line might be in this ring (the central one)
        ! Allocate memory for empirical distribution objects
        allocate(amplitude_distribution(number_of_rings),meridional_amplitude_distribution(number_of_rings))
        allocate(layer_line_score(number_of_rings))
        allocate(ring_frequency(number_of_rings))
        ! Initialise empirical distribution objects and remember what frequency each ring is at
#ifdef __INTEL_COMPILER
        do concurrent (current_ring=1:number_of_rings)
#else
        do current_ring=1,number_of_rings
#endif
            call amplitude_distribution(current_ring)%Init()
            call meridional_amplitude_distribution(current_ring)%Init()
            ring_frequency(current_ring) = expected_frequency + (current_ring-expected_ring)*ring_width
        enddo
        ! Loop over all elements of the amplitude spectrum. Compute the average amplitude average and std dev per ring.
        ! But when we are where a layer line could be, keep track of separate statistics.
        do j=1,self%physical_upper_bound_complex(2)
            j_logi = self%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2)
            y_freq_sq = (real(j_logi) / self%logical_dimensions(2))**2
            do i=1,self%physical_upper_bound_complex(1)
                i_logi = i-1
                x_freq_sq = (real(i_logi) / self%logical_dimensions(1))**2
                ! The spatial frequency in reciprocal pixels
                freq = sqrt(x_freq_sq+y_freq_sq)
                ! Which ring do we fall into?
                current_ring = expected_ring + nint((freq-expected_frequency)/ring_width)
                if (current_ring .ge. 1 .and. current_ring .le. number_of_rings) then
                    current_pixel_amplitude = cabs(self%complex_values(i,j,1))
                    ! Do we fall into where the layer line might be?
                    in_layer_line_region =  x_freq_sq .le. expected_half_width_sq .and. &
                                            y_freq_sq .ge. expected_frequency_min_sq .and. &
                                            y_freq_sq .le. expected_frequency_max_sq
                    if (in_layer_line_region) then
                        call meridional_amplitude_distribution(current_ring)%AddSampleValue(current_pixel_amplitude)
                    else
                        call amplitude_distribution(current_ring)%AddSampleValue(current_pixel_amplitude)
                    endif
                endif
            enddo
        enddo
        ! Compute layer line score for each ring
#ifdef __INTEL_COMPILER
        do concurrent (current_ring=1:number_of_rings)
#else
        do current_ring=1,number_of_rings
#endif
            !print *, 'ring ', current_ring, ' meridional mean ', meridional_amplitude_distribution(current_ring)%GetSampleMean()
            if (amplitude_distribution(current_ring)%GetSampleVariance() .gt. 0.0) then
                layer_line_score(current_ring) = ( meridional_amplitude_distribution(current_ring)%GetSampleMean() -&
                                                   amplitude_distribution(current_ring)%GetSampleMean() ) / &
                                                sqrt(amplitude_distribution(current_ring)%GetSampleVariance())
                !print *, 'ring ', current_ring, 'score ', layer_line_score(current_ring)
            else
                print *, 'sample variance of background is 0.0 for ring ', current_ring
            endif
        enddo
        ! Find the top candidate and return that
        best_ring      = maxloc(layer_line_score)
        best_frequency = ring_frequency(best_ring(1))
        best_score     = layer_line_score(best_ring(1))
    end subroutine LookForCrossBetaLayerLine

    !>  \brief  Compute the amplitude spectrum of the image. Can also compute phases, given the appropriate optional argument
    subroutine ComputeAmplitudeSpectrum(    self,amplitude_spectrum,    &
                                            compute_square,compute_log,compute_phases_instead,half_spectrum_only, &
                                            set_central_pixel_to_zero)
        use UsefulFunctions, only : RadiansToDegrees
        ! arguments
        class(image),               intent(inout)   ::  self
        type(image),                intent(inout)   ::  amplitude_spectrum
        logical,        optional,   intent(in)      ::  compute_square          !<  Compute square of amplitudes
        logical,        optional,   intent(in)      ::  compute_log             !<  Compute log of amplitudes
        logical,        optional,   intent(in)      ::  compute_phases_instead  !<  Compute spectrum of phases rather than amplitudes
        logical,        optional,   intent(in)      ::  half_spectrum_only      !<  Return a "half spectrum" with origin at first pixel
        logical,        optional,   intent(in)      ::  set_central_pixel_to_zero          ! private variables
        logical     ::  ccompute_phases_instead
        logical     ::  ccompute_square
        logical     ::  ccompute_log
        logical     ::  hhalf_spectrum_only
        logical     ::  sset_central_pixel_to_zero
        !type(image) ::  tmpimg
        integer     ::  half_spectrum_dimensions(3)
        ! debug
        integer :: i,j,k
        ! start work

        ! optional arguments
        ccompute_phases_instead = .false.
        if (present(compute_phases_instead)) ccompute_phases_instead = compute_phases_instead
        ccompute_log = .false.
        if (present(compute_log)) ccompute_log = compute_log
        ccompute_square = .false.
        if (present(compute_square)) ccompute_square = compute_square
        hhalf_spectrum_only = .false.
        if (present(half_spectrum_only)) hhalf_spectrum_only = half_spectrum_only
        sset_central_pixel_to_zero = .false.
        if (present(set_central_pixel_to_zero)) sset_central_pixel_to_zero = set_central_pixel_to_zero

        ! allocate output
        half_spectrum_dimensions = (/self%physical_upper_bound_complex(1),self%logical_dimensions(2),self%logical_dimensions(3)/)
        call amplitude_spectrum%Allocate(dims=half_spectrum_dimensions)

        ! copy to temporary image
        !tmpimg = self

        ! FFT
        if (self%is_in_real_space) call self%ForwardFFT()

        ! Take amplitudes
        if (ccompute_square) then
            amplitude_spectrum%real_values(1:amplitude_spectrum%logical_dimensions(1),:,:) = cabs(self%complex_values)**2
        else if (ccompute_log) then
            amplitude_spectrum%real_values(1:amplitude_spectrum%logical_dimensions(1),:,:) = log10(cabs(self%complex_values))
        else if (ccompute_phases_instead) then
            where(self%complex_values .ne. (0.0,0.0))
                amplitude_spectrum%real_values(1:amplitude_spectrum%logical_dimensions(1),:,:) &
                                            = RadiansToDegrees(atan2(aimag(self%complex_values),real(self%complex_values)))
            end where
        else
            do k=1,amplitude_spectrum%logical_dimensions(3)
                do j=1,amplitude_spectrum%logical_dimensions(2)
                    do i=1,amplitude_spectrum%logical_dimensions(1)
                        amplitude_spectrum%real_values(i,j,k) = cabs(self%complex_values(i,j,k))
                    enddo
                enddo
            enddo
            ! Original line below - for some reason (maybe compiler bug?) this segfaults with ifort 15.0.090
            !amplitude_spectrum%real_values(1:amplitude_spectrum%logical_dimensions(1),:,:) = cabs(self%complex_values)
        endif
        amplitude_spectrum%object_is_centered_in_box = .false.

        ! Convert to full image
        if (hhalf_spectrum_only) then
            if (sset_central_pixel_to_zero) amplitude_spectrum%real_values(1,1,1) = 0.0e0
        else
            call amplitude_spectrum%SpectrumFromHalfToFull()
            amplitude_spectrum%real_values( amplitude_spectrum%GetPhysicalIndexOfBoxCenter(1), &
                                            amplitude_spectrum%GetPhysicalIndexOfBoxCenter(2), &
                                            amplitude_spectrum%GetPhysicalIndexOfBoxCenter(3)) = 0.0e0
        endif

        !
        amplitude_spectrum%is_in_real_space = .true.
    end subroutine ComputeAmplitudeSpectrum

    !> \brief   Move the origin from (1,1,1) to the center of the image. Useful when the phase or amplitude
    !!          components of an image FT were used to create a real image, and one wants the full representation
    !!          of the FT, including negative frequencies of the first dimension. On output, the image is twice
    !!          as large in the first dimension.
    subroutine SpectrumFromHalfToFull(self)
        ! arguments
        class(image),    intent(inout)   ::  self
        ! private variables
        type(image) ::  imgout
        integer     ::  old_dim(3), new_dim(3)
        integer     ::  imgout_origin(3)
        integer     ::  i,j,k
        integer     ::  self_phys_addr(3)
        integer     ::  imgout_logi_addr(3)
        logical     ::  symmetry_mate

        ! start work

        ! work out old and new dimensions
        old_dim = self%logical_dimensions
        new_dim = [(old_dim(1)-1)*2,old_dim(2),old_dim(3)]
        if (new_dim(1) .eq. old_dim(1)-1) new_dim(1) = old_dim(1)

        ! prepare output image
        imgout%logical_dimensions = new_dim
        call imgout%Allocate()

        ! origin of output image
        imgout_origin = imgout%physical_address_of_box_center

        ! loop over output image
        do i=1,imgout%logical_dimensions(1)
            symmetry_mate = i .lt. imgout_origin(1) ! we are in the negative frequencies, so we'll have to go for the symmetry mate
            imgout_logi_addr(1) = i - imgout_origin(1) ! logical address
            if (symmetry_mate) imgout_logi_addr(1) = - imgout_logi_addr(1)
            do j=1,imgout%logical_dimensions(2)
                imgout_logi_addr(2) = j - imgout_origin(2) ! logical address
                if (symmetry_mate) imgout_logi_addr(2) = - imgout_logi_addr(2)
                do k=1,imgout%logical_dimensions(3)
                    imgout_logi_addr(3) = k - imgout_origin(3) ! logical address
                    if (symmetry_mate) imgout_logi_addr(3) = - imgout_logi_addr(3)
                    ! get the physical address in the input image corresponding to the output voxel
                    call imgout%PhysicalAddressGivenLogicalAddressInFourierSpace(imgout_logi_addr,self_phys_addr)
                    ! copy value over
                    imgout%real_values(i,j,k) = self%real_values(self_phys_addr(1),self_phys_addr(2),self_phys_addr(3))
                enddo
            enddo
        enddo

        ! replace input with output
        self = imgout

        self%object_is_centered_in_box = .true.

    end subroutine SpectrumFromHalfToFull

    !>  \brief  Randomise phases beyond a given spatial frequency
    subroutine RandomisePhases(self,threshold_spatial_frequency)
        class(Image),               intent(inout)   ::  self
        real,                       intent(in)      ::  threshold_spatial_frequency   !<  1/pixel (Nyquist is 0.5). Beyond this frequency, phases will be randomised.
        !
        integer ::  i,j,k, i_logi, j_logi, k_logi
        real    ::  x_freq_sq, y_freq_sq, z_freq_sq
        real    ::  threshold_spatial_frequency_sq
        real    ::  current_spatial_frequency_sq
        real    ::  current_phase, current_amplitude
        !
        threshold_spatial_frequency_sq = threshold_spatial_frequency**2
        ! We may need to be in Fourier space
        if (self%IsInRealSpace()) call self%ForwardFFT()
        ! Randomise the seed
        call random_seed()

        !
        do k=1,self%physical_upper_bound_complex(3)
            k_logi = self%LogicalIndexGivenPhysicalIndexInFourierSpace(k,3)
            z_freq_sq = (real(k_logi)/self%logical_dimensions(3))**2
            do j=1,self%physical_upper_bound_complex(2)
                j_logi = self%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2)
                y_freq_sq = (real(j_logi)/self%logical_dimensions(2))**2 + z_freq_sq
                do i=1,self%physical_upper_bound_complex(1)
                    i_logi=i-1
                    x_freq_sq = ((real(i_logi))/self%logical_dimensions(1))**2
                    !
                    current_spatial_frequency_sq = x_freq_sq + y_freq_sq
                    !
                    if (current_spatial_frequency_sq .gt. threshold_spatial_frequency_sq) then
                        call random_number(current_phase) ! between 0 and 1
                        current_phase = current_phase * 2.0 * pi
                        current_amplitude = cabs(self%complex_values(i,j,k))
                        self%complex_values(i,j,k) = cmplx(cos(current_phase)*current_amplitude, &
                                                           sin(current_phase)*current_amplitude)
                    endif
                enddo
            enddo
        enddo
    end subroutine RandomisePhases

    !>  \brief  Compute the phase spectrum of the image
    subroutine ComputePhaseSpectrum(self,phase_spectrum)
        ! arguments
        class(image),               intent(inout)   ::  self
        type(image),                intent(inout)   ::  phase_spectrum
        ! start work
        call self%ComputeAmplitudeSpectrum(phase_spectrum,compute_phases_instead=.true.)
    end subroutine ComputePhaseSpectrum

    !>  \brief  Return a real-space image with the Fourier amplitudes of the image object
    function GetAmplitudeSpectrum(self,half_spectrum_only) result(amplitude_spectrum)
        ! arguments
        class(image),               intent(inout)   ::  self
        logical,        optional,   intent(in)      ::  half_spectrum_only
        ! result
        type(image)                                 ::  amplitude_spectrum
        ! Start work
        call self%ComputeAmplitudeSpectrum(amplitude_spectrum,half_spectrum_only=half_spectrum_only)
    end function GetAmplitudeSpectrum

    !>  \brief  Return a real-space image with the Fourier phases of the image object
    function GetPhaseSpectrum(self) result(phase_spectrum)
        ! arguments
        class(image),               intent(inout)   ::  self
        ! result
        type(image)                                 ::  phase_spectrum
        ! Start work
        call self%ComputePhaseSpectrum(phase_spectrum)
    end function GetPhaseSpectrum

    !>  \brief  Return the logical index value (ie one of the 3 components of the address) given the physical index value. Applies only to dimensions 2 and 3. For dimension 1, just subtract 1.
    pure integer function LogicalIndexGivenPhysicalIndexInFourierSpace(self,physical_index,which_dimension) result(logical_index)
        ! arguments
        class(Image),   intent(in)  ::  self
        integer,        intent(in)  ::  physical_index      !<  The physical index along dimension 2 or 3 of the image object
        integer,        intent(in)  ::  which_dimension     !<  The dimension along which the conversion should be done
        ! start work
        if (physical_index .ge. self%physical_index_of_first_negative_frequency(which_dimension)) then
            logical_index = physical_index - (self%logical_dimensions(which_dimension)) - 1
        else
            logical_index = physical_index - 1
        endif
    end function LogicalIndexGivenPhysicalIndexInFourierSpace



    !>  \brief  Returns the physical address of voxel with given logical address.
    !!
    !!          The address returned is for a complex image with origin at address (1,1,1),
    !!          dimension 1 from 1 to dim(1)/2+1 (division by 2 rounded down), and other
    !!          2 dimensions from 1 to dim, folded halfway, the FFTW way.
    pure subroutine PhysicalAddressGivenLogicalAddressInFourierSpace(self,logi,phys)
        ! arguments
        class(Image),   intent(in)  ::  self
        integer,        intent(in)  ::  logi(3)
        ! result
        integer,        intent(out) ::  phys(3)
        ! private variables
        integer ::  i
        ! start work
        if (logi(1) .ge. 0) then
            ! first dimension
            phys(1) = logi(1) + 1
            ! second & third dimensions
            do i=2,3
                if (logi(i) .ge. 0) then
                    phys(i) = logi(i) + 1
                else
                    phys(i) = self%logical_dimensions(i) + 1 + logi(i)
                endif
            enddo
        else ! we need to grab the Friedel mate
            ! first dimension
            phys(1) = -logi(1) + 1
            ! second & third dimensions
            do i=2,3
                if (logi(i) .gt. 0) then
                    phys(i) = self%logical_dimensions(i) - logi(i) + 1
                else
                    phys(i) = - logi(i) + 1
                endif
            enddo
        endif
        if (self%logical_dimensions(3) .eq. 1) phys(3) = 1
    end subroutine PhysicalAddressGivenLogicalAddressInFourierSpace

    !>  \brief  Compute the physical address of the Hermitian symmetry mate. For most Fourier component,
    !!          the Hermitian symmetry mate is not explicitely addressable in the FFTW scheme and the physical
    !!          address is returned unchanged.
    pure subroutine PhysicalAddressOfHermitianMate(self,physical_address,physical_address_hermitian_mate)
        ! Arguments
        class(Image),   intent(in)  ::  self
        integer,        intent(in)  ::  physical_address(3)
        integer,        intent(out) ::  physical_address_hermitian_mate(3)
        ! temps
        integer ::  logical_address(3)
        integer ::  i
        !
        do i=1,3
            logical_address(i) = self%LogicalIndexGivenPhysicalIndexInFourierSpace(physical_address(i),i)
        enddo
        call self%PhysicalAddressGivenLogicalAddressInFourierSpace(-logical_address,physical_address_hermitian_mate)
    end subroutine PhysicalAddressOfHermitianMate

    !>  \brief  Compute the physical address of the Hermitian symmetry mate. For most Fourier component,
    !!          the Hermitian symmetry mate is not explicitely addressable in the FFTW scheme and the physical
    !!          address is returned unchanged.
    pure subroutine PhysicalIndicesOfHermitianMate(self,physical_index_i,physical_index_j,physical_index_k, &
                                                        physical_index_i_hermitian_mate, &
                                                        physical_index_j_hermitian_mate, &
                                                        physical_index_k_hermitian_mate)
        ! Arguments
        class(Image),   intent(in)  ::  self
        integer,        intent(in)  ::  physical_index_i, physical_index_j, physical_index_k
        integer,        intent(out) ::  physical_index_i_hermitian_mate
        integer,        intent(out) ::  physical_index_j_hermitian_mate,physical_index_k_hermitian_mate
        ! temps
        integer ::  logical_address(3)
        integer ::  physical_address_hermitian_mate(3)
        !
        logical_address(1) = self%LogicalIndexGivenPhysicalIndexInFourierSpace(physical_index_i,1)
        logical_address(2) = self%LogicalIndexGivenPhysicalIndexInFourierSpace(physical_index_j,2)
        logical_address(3) = self%LogicalIndexGivenPhysicalIndexInFourierSpace(physical_index_k,3)
        call self%PhysicalAddressGivenLogicalAddressInFourierSpace(-logical_address,physical_address_hermitian_mate)
        !
        physical_index_i_hermitian_mate = physical_address_hermitian_mate(1)
        physical_index_j_hermitian_mate = physical_address_hermitian_mate(2)
        physical_index_k_hermitian_mate = physical_address_hermitian_mate(3)
    end subroutine PhysicalIndicesOfHermitianMate


    !>  \brief  Update all properties related to looping & addressing in real & Fourier space, given the current logical dimensions.
    subroutine UpdateLoopingAndAddressing(self)
        use UsefulFunctions, only : IsEven
        class(image),               intent(inout)   ::  self
        ! private variables
        integer :: i
        ! start work

        ! The physical bounds of the complex address space
        self%physical_upper_bound_complex(2:3) = self%logical_dimensions(2:3)
        if (IsEven(self%logical_dimensions(1))) then
            self%physical_upper_bound_complex(1) = self%logical_dimensions(1)/2 + 1
        else
            self%physical_upper_bound_complex(1) = (self%logical_dimensions(1)-1)/2 + 1
        endif

        ! The address of the origin
        call self%UpdatePhysicalAddressOfBoxCenter()

        ! In each dimension, the physical index of the first pixel which stores negative frequencies
        ! Note that we never actually store negative frequencies for the first dimension. However, it is sometimes useful
        ! to pretend as though we did (for example when generating real-space images of CTFs).
        do i=1,3
            if (IsEven(self%logical_dimensions(i))) then
                self%physical_index_of_first_negative_frequency(i) =  self%logical_dimensions(i)/2 + 2
            else
                self%physical_index_of_first_negative_frequency(i) = (self%logical_dimensions(i)+3)/2
            endif
        enddo

        ! Update the Fourier voxel size
        self%fourier_voxel_size = 1.0d0 / dble(self%logical_dimensions)

        ! Update the bounds of the logical addresses
        do i=1,3
            if (IsEven(self%logical_dimensions(i))) then
                if (i .eq. 1) then
                    self%logical_lower_bound_complex(1) = -self%logical_dimensions(1)/2
                    self%logical_upper_bound_complex(1) =  self%logical_dimensions(1)/2
                else if (i .ge. 2) then
                    self%logical_lower_bound_complex(i) = -self%logical_dimensions(i)/2
                    self%logical_upper_bound_complex(i) =  self%logical_dimensions(i)/2 - 1
                endif
                self%logical_lower_bound_real(i)    = -self%logical_dimensions(i)/2
                self%logical_upper_bound_real(i)    =  self%logical_dimensions(i)/2-1
            else
                if (i .eq. 1) then
                    self%logical_lower_bound_complex(1) = -(self%logical_dimensions(1)-1)/2
                    self%logical_upper_bound_complex(1) =  (self%logical_dimensions(1)-1)/2
                else if (i .ge. 2) then
                    self%logical_lower_bound_complex(i) = -(self%logical_dimensions(i)-1)/2
                    self%logical_upper_bound_complex(i) =  (self%logical_dimensions(i)-1)/2
                endif
                self%logical_lower_bound_real(i)    = -(self%logical_dimensions(i)-1)/2
                self%logical_upper_bound_real(i)    =  (self%logical_dimensions(i)-1)/2
            endif
        enddo

    end subroutine UpdateLoopingAndAddressing

    !>  \brief  Returns the physical address of the image origin
    pure subroutine UpdatePhysicalAddressOfBoxCenter(self)
        ! arguments
        class(image),   intent(inout)  ::  self
        ! private variables
        integer ::  i
        ! start work
        do i=1,3
            self%physical_address_of_box_center(i) = IndexOfCentralPixelGivenLogicalDimension(self%logical_dimensions(i))
        enddo

    end subroutine UpdatePhysicalAddressOfBoxCenter

    !>  \brief Given a logical dimension of an image, return the index for the central pixel
    pure elemental function IndexOfCentralPixelGivenLogicalDimension(logical_dimension) result(i)
        use UsefulFunctions, only : IsEven
        ! arguments
        integer,    intent(in)  ::  logical_dimension   !<  logical dimension
        ! results
        integer                 ::  i                   !<  index of central voxel
        ! start work
        if (IsEven(logical_dimension)) then
            ! dimension is even
            i = logical_dimension/2+1
        else
            ! dimension is odd
            i = (logical_dimension-1)/2+1
        endif
    end function IndexOfCentralPixelGivenLogicalDimension

    !>  \brief  Work out whether a given Fourier component is the Hermitian mate of another component also described explicitely in the FFTW representation
    pure logical function FourierComponentIsRedundant(self,physical_index_x,physical_index_y,physical_index_z) result(is_redundant)
        use UsefulFunctions, only : IsEven
        class(Image),               intent(in)      ::  self
        integer,                    intent(in)      ::  physical_index_x
        integer,                    intent(in)      ::  physical_index_y
        integer,                    intent(in)      ::  physical_index_z
        ! Start work

        is_redundant = physical_index_x .eq. 1 .and. .not. ((physical_index_y .eq. 1) .and. (physical_index_z .eq. 1)) &
                                               .and. physical_index_y .ge. self%physical_index_of_first_negative_frequency(2)


        if (self%IsAVolume()) then
            is_redundant = is_redundant .and. physical_index_z .ge. self%physical_index_of_first_negative_frequency(3)
        endif

    end function FourierComponentIsRedundant

    !>  \brief  Work out whether a given Fourier component has a Hermitian mate which is also described explicitely by the FFTW
    pure logical function FourierComponentHasExplicitHermitianMate(self,physical_index_x,physical_index_y,physical_index_z) &
                            result(explicit_mate)
        use UsefulFunctions, only : IsEven
        class(Image),               intent(in)      ::  self
        integer,                    intent(in)      ::  physical_index_x
        integer,                    intent(in)      ::  physical_index_y
        integer,                    intent(in)      ::  physical_index_z
        ! Start work

        explicit_mate = physical_index_x .eq. 1 .and. .not. ((physical_index_y .eq. 1) .and. (physical_index_z .eq. 1))
        ! We assume that the Y dimension is the non-flat one
        if (IsEven(self%logical_dimensions(2))) then
            explicit_mate = explicit_mate .and. physical_index_y .ne. self%physical_index_of_first_negative_frequency(2)-1
        endif

        if (self%IsAVolume()) then
            if (IsEven(self%logical_dimensions(3))) then
                explicit_mate = explicit_mate .and. physical_index_z .ne. self%physical_index_of_first_negative_frequency(3)-1
            endif
        endif

    end function FourierComponentHasExplicitHermitianMate

    !>  \brief  Return the logical dimensions of the image
    pure function GetLogicalDimensions(self) result(logical_dimensions)
        class(image),               intent(in)      ::  self
        integer                                     ::  logical_dimensions(3)
        logical_dimensions = self%logical_dimensions
    end function GetLogicalDimensions

    !>  \brief  Return the logical dimension of the image along one of the axes
    pure integer function GetLogicalDimension(self,which_dimension)
        class(image),               intent(in)      ::  self
        integer,                    intent(in)      ::  which_dimension
        GetLogicalDimension = self%logical_dimensions(which_dimension)
    end function GetLogicalDimension

    pure integer function GetPhysicalUpperBoundComplex(self,which_dimension)
        class(image),               intent(in)      ::  self
        integer,                    intent(in)      ::  which_dimension
        GetPhysicalUpperBoundComplex = self%physical_upper_bound_complex(which_dimension)
    end function GetPhysicalUpperBoundComplex


    !>  \brief  Return the Fourier voxel size in the given dimension
    pure real function GetFourierVoxelSize(self,which_dimension)
        class(Image),               intent(in)      ::  self
        integer,                    intent(in)      ::  which_dimension
        GetFourierVoxelSize = self%fourier_voxel_size(which_dimension)
    end function GetFourierVoxelSize

    !>  \brief  Return the physical address of box center
    pure function GetPhysicalAddressOfBoxCenter(self) result(address)
        class(Image),               intent(in)      ::  self
        integer                                     ::  address(3)
        address = self%physical_address_of_box_center
    end function GetPhysicalAddressOfBoxCenter

    !>  \brief  Return one of the indices to the central pixel
    pure integer function GetPhysicalIndexOfBoxCenter(self,which_dimension)
        class(Image),               intent(in)      ::  self
        integer,                    intent(in)      ::  which_dimension
        GetPhysicalIndexOfBoxCenter = self%physical_address_of_box_center(which_dimension)
    end function GetPhysicalIndexOfBoxCenter

    !>  \brief  Change the logical dimensions of an image and update all related values
    pure subroutine SetLogicalDimensions_1(self,new_x,new_y,new_z)
        ! arguments
        class(image),               intent(inout)   ::  self
        integer,                    intent(in)      ::  new_x
        integer,                    intent(in)      ::  new_y
        integer,        optional,   intent(in)      ::  new_z
        !
        self%logical_dimensions(1) = new_x
        self%logical_dimensions(2) = new_y
        self%logical_dimensions(3) = 1
        if(present(new_z)) then
            self%logical_dimensions(3) = new_z
        endif
    end subroutine SetLogicalDimensions_1

    !>  \brief  Change the logical dimensions of an image
    pure subroutine SetLogicalDimensions_2(self,new_logical_dimensions)
        class(image),               intent(inout)   ::  self
        integer,                    intent(in)      ::  new_logical_dimensions(3)
        !
        call self%SetLogicalDimensions_1(new_logical_dimensions(1),new_logical_dimensions(2),new_logical_dimensions(3))
    end subroutine SetLogicalDimensions_2


    !>  \brief  Add an image to the end of an array of images
    subroutine AppendToArray(self,array_of_images)
        class(image),                   intent(in)      ::  self
        type(image),    allocatable,    intent(inout)   ::  array_of_images(:)  !<  Array of images. It will be grown to accomodate self.
        ! private variables
        type(image),    allocatable ::  tmp(:)
        ! start work
        allocate(tmp(size(array_of_images)+1))
        tmp(1:size(array_of_images)) = array_of_images(1:size(array_of_images))
        tmp(size(tmp)) = self
        call move_alloc(from=tmp,to=array_of_images) ! f2003 intrinsic
    end subroutine AppendToArray

    !>  \brief  Allocate memory for the Image object. Dimensions may be set in advance or supplied as arguments. Alternatively, a mould can be supplied.
    !!
    !!  If the object is already allocated with correct dimensions, nothing happens. Otherwise, object is deallocated first.
    subroutine Allocate(self,mould,dims, in_real_space)
        use iso_c_binding
        use omp_lib
        use fftw33, only : fftwf_plan_dft_r2c_3d, fftwf_plan_dft_r2c_2d
        use fftw33, only : fftwf_plan_dft_c2r_3d, fftwf_plan_dft_c2r_2d
        use fftw33, only : fftwf_alloc_complex, fftwf_malloc
        use fftw33, only : fftw_estimate, fftw_measure, fftw_patient, fftw_no_simd
        use UsefulFunctions, only : IsEven
        ! arguments
        class(image),               intent(inout)   ::  self
        type(image),    optional,   intent(in)      ::  mould       !<  The image will be allocated with this image's dimensions
        integer,        optional,   intent(in)      ::  dims(3)     !<  Dimensions of the image
        logical,        optional,   intent(in)      ::  in_real_space !< Should the image be in real space
        ! Private variables
        integer                 ::  array_shape(3)
        logical,    parameter   ::  debug   =   .false.
        logical                 ::  already_allocated_correctly

        ! Start work

        ! Do we need to do anything?
        already_allocated_correctly = self%is_in_memory
        if (present(dims)) already_allocated_correctly = already_allocated_correctly &
                                                        .and. all(self%logical_dimensions .eq. dims)
        if (present(mould)) already_allocated_correctly = already_allocated_correctly &
                                                        .and. (self .SameDimensions. mould)

        if (already_allocated_correctly) then
            if (present(mould)) then
                self%object_is_centered_in_box = mould%object_is_centered_in_box
                self%is_in_real_space = mould%is_in_real_space
            endif

            if (present(in_real_space)) self%is_in_real_space = in_real_space
        else
            ! Check whether memory is already allocated
            if (self%is_in_memory) then
                if (debug) write(*,'(2a)') '**debug(image_allocate): tyring to allocate memory to image already in memory, ',    &
                                            'will deallocate first'
                call self%Deallocate()
            endif

            if (present(mould)) then
                if (debug) write(*,'(a)') '**debug(image_allocate): using dimensions, type and object_is_centered_in_box from mould'
                call self%SetLogicalDimensions(mould%logical_dimensions(1),mould%logical_dimensions(2),mould%logical_dimensions(3))
                self%is_in_real_space = mould%is_in_real_space
                self%object_is_centered_in_box = mould%object_is_centered_in_box
            endif

            if (present(dims)) call self%SetLogicalDimensions(dims)


            if (present(in_real_space)) self%is_in_real_space = in_real_space

            ! Check the dimensions of the image have been set
            if (any(self%logical_dimensions .eq. 0)) then
                call this_program%TerminateWithFatalError('Image::Allocate','Image dimensions not defined')
            endif

            ! Update all looping- & addressing-related variables
            call self%UpdateLoopingAndAddressing()

            ! Workout the desired shape of the complex array
            if (IsEven(self%logical_dimensions(1))) then
                array_shape(1) = self%logical_dimensions(1)/2 + 1
            else
                array_shape(1) = (self%logical_dimensions(1)-1)/2 + 1
            endif
            array_shape(2) = self%logical_dimensions(2)
            array_shape(3) = self%logical_dimensions(3)

            ! Do the allocation, the first dimension gets a +2 because of the FFTW
            ! Letting FFTW do the allocation in C ensures that we will be using aligned memory
            !$omp   critical    (fftw_omp_crit)
            self%p = fftwf_malloc(int(product(array_shape)*8, c_size_t))
            ! Remember that we did the allocation
            self%is_in_memory = .true.
            !$omp   end critical    (fftw_omp_crit)

            ! Check whether the allocation went smoothly
            ! Don't know how to do this when using fftw_malloc

            ! Set up the complex_values pointer
            call c_f_pointer(self%p,self%complex_values,array_shape)

            ! Workout the shape of the real array
            array_shape(1) = 2 * array_shape(1)

            ! Set up the real_values pointer
            call c_f_pointer(self%p,self%real_values,array_shape)

            ! DO NOT Initialise the contents - bad for performance
            !self%complex_values = (0.0,0.0)

            ! Prepare the plans for FFTW
            if (.not. self%planned) then
                !$omp   critical    (fftw_omp_crit)
                if ( self%logical_dimensions(3) .gt. 1) then
                    self%plan_fwd = fftwf_plan_dft_r2c_3d(  self%logical_dimensions(3), &
                                                            self%logical_dimensions(2), &
                                                            self%logical_dimensions(1), &
                                                            self%real_values,           &
                                                            self%complex_values,        &
                                                            fftw_estimate)

                    self%plan_bwd = fftwf_plan_dft_c2r_3d(  self%logical_dimensions(3), &
                                                            self%logical_dimensions(2), &
                                                            self%logical_dimensions(1), &
                                                            self%complex_values,        &
                                                            self%real_values,           &
                                                            fftw_estimate)

                else
                    self%plan_fwd = fftwf_plan_dft_r2c_2d(  self%logical_dimensions(2), &
                                                            self%logical_dimensions(1), &
                                                            self%real_values,           &
                                                            self%complex_values,        &
                                                            fftw_estimate)

                    self%plan_bwd = fftwf_plan_dft_c2r_2d(  self%logical_dimensions(2), &
                                                            self%logical_dimensions(1), &
                                                            self%complex_values,        &
                                                            self%real_values,           &
                                                            fftw_estimate)

                endif
                self%planned = .true.
                !$omp   end critical    (fftw_omp_crit)
            else
                call this_program%TerminateWithFatalError('Image::Allocate','FFTW plans were already setup')
            endif
        endif

    end subroutine Allocate

    !>  \brief  Deallocate memory allocated to an image
    subroutine Deallocate(self)
        use fftw33
        ! arguments
        class(image),               intent(inout)   ::  self
        ! start work

        ! deallocate image data
        if (associated(self%complex_values)) nullify(self%complex_values)     ! dissociate the complex_values pointer from the memory
        if (associated(self%real_values)) nullify(self%real_values)   ! dissociate the actual memory
        ! free the memory
        !$omp   critical    (fftw_omp_crit)
        if (self%is_in_memory) then
            call fftwf_free(self%p)
            self%p = c_null_ptr
            self%is_in_memory   =   .false.
        endif
        if (self%planned) then
            call fftwf_destroy_plan(self%plan_fwd)
            self%plan_fwd = c_null_ptr
            call fftwf_destroy_plan(self%plan_bwd)
            self%plan_bwd = c_null_ptr
            self%planned     =   .false.
        endif
        !$omp   end critical    (fftw_omp_crit)
    end subroutine Deallocate

    !>  \brief  Check whether an image is allocated
    elemental pure logical function IsAllocated(self)
        class(image),               intent(in)      ::  self
        IsAllocated = self%is_in_memory
    end function IsAllocated

    !>  \brief  Check whether the object is centered in the box or has its center near the origin (first voxel)
    elemental pure logical function ObjectIsCenteredInBox(self)
        class(Image),   intent(in)  ::  self
        ObjectIsCenteredInBox = self%object_is_centered_in_box
    end function ObjectIsCenteredInBox

    !>  \brief  Reset all image properties
    !!  \todo   Make sure all properties are reset
    subroutine Reset(self)
        class(image),               intent(inout)   ::  self
        !
        self%logical_dimensions             =   0
        self%is_in_real_space               =   .true.
        self%object_is_centered_in_box      =   .true.
        self%physical_upper_bound_complex   =   0
        self%fourier_voxel_size             =   0.0
        self%p                              =   c_null_ptr
        self%real_values                    =>  null()
        self%complex_values                 =>  null()
        self%is_in_memory                   =   .false.
        self%plan_fwd                       =   c_null_ptr
        self%plan_bwd                       =   c_null_ptr
        self%planned                        =   .false.
    end subroutine Reset

    !>  \brief  Image destructor
    subroutine Destruct(self)
        type(image),                intent(inout)   ::  self
        call self%Deallocate()
    end subroutine Destruct

    !>  \brief  Image array destructor
    subroutine DestructArray(self)
        type(image),                intent(inout)   ::  self(:)
        ! private variables
        integer :: counter
        ! start work
        do counter=1,size(self)
            call self(counter)%Deallocate()
        enddo
    end subroutine DestructArray

    !>  \brief  print info about the image object. used mostly for debugging.
    subroutine PrintInfo(self)
        use fftw33, only : fftwf_print_plan
        ! arguments
        class(image),               intent(in)      ::  self
        ! private variables
        character(len=1),           parameter       ::  tab = '2'
        character(len=2),           parameter       ::  text_width = '50'
        logical,                    parameter       ::  print_fftw_plans = .true.
        logical,                    parameter       ::  print_rvalue_stats = .true.
        ! start work
        write(*,'(a)')                                      ' '
        write(*,'(a)')                                      'image info:'
        write(*,'('//tab//'x,a'//text_width//',3(i0,1x))')  'logical dimensions: ',     self%logical_dimensions
        write(*,'('//tab//'x,a'//text_width//',l1)')        'in real space? ',          self%is_in_real_space
        write(*,'('//tab//'x,a'//text_width//',l1)')        'Object centered in box? ',       self%object_is_centered_in_box
        write(*,'('//tab//'x,a'//text_width//',2l1)')       'real, complex arrays allocated, associated? ', &
                                                                                        associated(self%real_values), &
                                                                                        associated(self%complex_values)
        if (associated(self%real_values)) then
        write(*,'('//tab//'x,a'//text_width//',3(i0,1x))')  'real array dimensions: ',  size(self%real_values,dim=1), &
                                                                                        size(self%real_values,dim=2), &
                                                                                        size(self%real_values,dim=3)
        endif
        if (associated(self%complex_values)) then
        write(*,'('//tab//'x,a'//text_width//',3(i0,1x))')  'comp array dimensions: ',  size(self%complex_values,dim=1), &
                                                                                        size(self%complex_values,dim=2), &
                                                                                        size(self%complex_values,dim=3)
        endif
        write(*,'('//tab//'x,a'//text_width//',3(i0,1x))')  'physical upper bound comp.', self%physical_upper_bound_complex
        write(*,'('//tab//'x,a'//text_width//',3(i0,1x))')  'physical index of negative freq.', &
                                                                                self%physical_index_of_first_negative_frequency
        write(*,'('//tab//'x,a'//text_width//',l1)')        'image is in memory? ',     self%is_in_memory
        write(*,'('//tab//'x,a'//text_width//',l1)')        'fftw plans ready? ',       self%planned
#ifdef __INTEL_COMPILER
        write(*,'('//tab//'x,a'//text_width//',2l1)')       'real, complex contiguous? ',is_contiguous(self%real_values),&
                                                                                        is_contiguous(self%complex_values)
#endif
        if (print_fftw_plans) then
            flush(6)
            write(*,'(a)') 'forward plan:'
            call fftwf_print_plan(self%plan_fwd)
            flush(6)
            write(*,'(a)') 'backward plan:'
            call fftwf_print_plan(self%plan_bwd)
            flush(6)
        endif
        if (print_rvalue_stats .and. self%is_in_real_space) then
            write(*,'('//tab//'x,a'//text_width//',2(f0.3,1x))')  'real array min, max values: ',  self%GetMinimumValue(), &
                                                                                                   self%GetMaximumValue()
            write(*,'('//tab//'x,a'//text_width//',2(f0.3,1x))')  'real array mean, std values: ', self%GetAverageOfValues(), &
                                                                                                   self%GetSigmaOfValues()
            write(*,'('//tab//'x,a'//text_width//',l1)') 'any nan values? ', self%HasNan()
        endif
        write(*,'(a)')                                  ' '
    end subroutine PrintInfo

    !>  \brief Print out real or complex values
    subroutine PrintValues(self)
        class(Image),   intent(in)  ::  self
        ! Private variables
        integer ::  j,k

        ! Start work
        do k=1,self%logical_dimensions(3)
            do j=1,self%logical_dimensions(2)
                write(*,'(2(a,i0),a)',advance='no') 'slice ', k, ', line ', j, ': '
                if (self%IsInRealSpace()) then
                    write(*,*) self%real_values(1:self%logical_dimensions(1),j,k)
                else
                    write(*,*) self%complex_values(:,j,k)
                endif
            enddo
        enddo
    end subroutine PrintValues

    !>  \brief  Print out some statistics about the Fourier transform of the image
    subroutine PrintFourierStatistics(self)
        use EmpiricalDistributions
        ! Arguments
        class(Image),   intent(inout)   ::  self
        ! Private variables
        type(EmpiricalDistribution), allocatable ::  distributions_real(:), distributions_imaginary(:), distributions_amplitudes(:)
        integer ::  address_x,address_y,address_z
        real    ::  sq_spatial_freq_x, sq_spatial_freq_y, sq_spatial_freq_z
        real    ::  width_of_resolution_zones
        integer ::  number_of_resolution_zones
        integer ::  current_resolution_zone
        real    ::  current_spatial_frequency
        ! Start work

        ! How many resolution zones?
        if (self%IsAVolume()) then
            number_of_resolution_zones = nint(real(minval(self%logical_dimensions))*0.05) + 1
        else
            number_of_resolution_zones = nint(real(minval(self%logical_dimensions(1:2)))*0.05) + 1
        endif
        width_of_resolution_zones = 0.5 / (number_of_resolution_zones-1) ! in reciprocal pixels

        ! FT
        if (self%IsInRealSpace()) call self%ForwardFFT

        ! Allocate distributions
        allocate(distributions_real(number_of_resolution_zones), &
                 distributions_imaginary(number_of_resolution_zones), &
                 distributions_amplitudes(number_of_resolution_zones))

        ! Initialize distributions
        do current_resolution_zone=1,number_of_resolution_zones
            call distributions_real(current_resolution_zone)%Init()
            call distributions_imaginary(current_resolution_zone)%Init()
            call distributions_amplitudes(current_resolution_zone)%Init()
        enddo

        ! Loop over image
        do address_z=1,self%physical_upper_bound_complex(3)
            sq_spatial_freq_z = (self%LogicalIndexGivenPhysicalIndexInFourierSpace(address_z,3)*self%fourier_voxel_size(3))**2
            do address_y=1,self%physical_upper_bound_complex(2)
                sq_spatial_freq_y = (self%LogicalIndexGivenPhysicalIndexInFourierSpace(address_y,2)*self%fourier_voxel_size(2))**2
                do address_x=1,self%physical_upper_bound_complex(1)
                    sq_spatial_freq_x = ((address_x-1)*self%fourier_voxel_size(1))**2

                    ! Some Fourier components have both Friedel mates explicitely described by the FFTW representation.
                    ! We avoid counting them twice.
                    if (  self%FourierComponentIsRedundant(address_x,address_y,address_z)) cycle

                    ! Current spatial frequency
                    current_spatial_frequency = sqrt(sq_spatial_freq_x+sq_spatial_freq_y+sq_spatial_freq_z)

                    ! Work out which resolution zone we are in
                    current_resolution_zone = nint(current_spatial_frequency/width_of_resolution_zones) + 1

                    ! Add sample value
                    if (current_resolution_zone .le. number_of_resolution_zones) then
                        call distributions_real(current_resolution_zone)%AddSampleValue(&
                                                real(self%complex_values(address_x,address_y,address_z)))
                        call distributions_imaginary(current_resolution_zone)%AddSampleValue(&
                                                aimag(self%complex_values(address_x,address_y,address_z)))
                        call distributions_amplitudes(current_resolution_zone)%AddSampleValue(&
                                                cabs(self%complex_values(address_x,address_y,address_z)))
                    endif

                enddo
            enddo
        enddo

        ! Print out results
        write(*,'(15(a9,1x))')  'Zone #', 'Num pix', 'Freq', 'Ave(real)', 'Ave(r^2)','Ave(imag)', 'Ave(im^2)',&
                                'Ave(ampl)', 'Av^2(amp)','Var(ampl)', 'Std(ampl)','Av(amp^2)','S(av(a^2', &
                                'Var(F)', 'Std(F)'
        do current_resolution_zone=1,number_of_resolution_zones
            write(*,'(2(i9,1x),f9.2,1x,12(ES9.2,1x))')  &
                                                    current_resolution_zone, &
                                                    distributions_real(current_resolution_zone)%GetNumberOfSamples(), &
                                                    (current_resolution_zone-1)*width_of_resolution_zones, &
                                                    distributions_real(current_resolution_zone)%GetSampleMean(), &
                                                    distributions_real(current_resolution_zone)%GetSampleSumOfSquares()&
                                                    /distributions_real(current_resolution_zone)%GetNumberOfSamples(), &
                                                    distributions_imaginary(current_resolution_zone)%GetSampleMean(), &
                                                    distributions_imaginary(current_resolution_zone)%GetSampleSumOfSquares()&
                                                    /distributions_imaginary(current_resolution_zone)%GetNumberOfSamples(), &
                                                    distributions_amplitudes(current_resolution_zone)%GetSampleMean(), &
                                                    distributions_amplitudes(current_resolution_zone)%GetSampleMean()**2, &
                                                    distributions_amplitudes(current_resolution_zone)%GetSampleVariance(), &
                                                    sqrt(distributions_amplitudes(current_resolution_zone)%GetSampleVariance()), &
                                                    distributions_amplitudes(current_resolution_zone)%GetSampleSumOfSquares()&
                                                    /distributions_amplitudes(current_resolution_zone)%GetNumberOfSamples(), &
                                                    sqrt(distributions_amplitudes(current_resolution_zone)%GetSampleSumOfSquares()&
                                                    /distributions_amplitudes(current_resolution_zone)%GetNumberOfSamples()), &
                                                    distributions_real(current_resolution_zone)%GetSampleVariance() + &
                                                    distributions_imaginary(current_resolution_zone)%GetSampleVariance(), &
                                                    sqrt(distributions_real(current_resolution_zone)%GetSampleVariance() + &
                                                    distributions_imaginary(current_resolution_zone)%GetSampleVariance())

        enddo


    end subroutine PrintFourierStatistics


    pure logical function IsInRealSpace(self)
        class(image),   intent(in)  ::  self
        IsInRealSpace = self%is_in_real_space
    end function IsInRealSpace

    pure logical function IsInFourierSpace(self)
        class(image),   intent(in)  ::  self
        IsInFourierSpace = .not. self%is_in_real_space
    end function IsInFourierSpace

    !>  \brief  Compute the average value in the spectrum, ignoring the center and the cross
    subroutine ComputeAverageAndSigmaOfValuesInSpectrum(self,minimum_radius,maximum_radius,average,sigma,cross_half_width)
        use EmpiricalDistributions
        !
        class(Image),       intent(in)  ::  self
        real,               intent(in)  ::  minimum_radius
        real,               intent(in)  ::  maximum_radius
        real,               intent(out) ::  average
        real,               intent(out) ::  sigma
        real,   optional,   intent(in)  ::  cross_half_width    !<  Half width of the central cross to be ignored
        ! private variables
        integer                     ::  i,j
        real                        ::  x_sq, y_sq, rad_sq
        type(EmpiricalDistribution) ::  my_distribution
        real                        ::  min_rad_sq, max_rad_sq
        real                        ::  ccross_half_width_sq
        ! start work
        min_rad_sq = minimum_radius**2
        max_rad_sq = maximum_radius**2
        if (minimum_radius .ge. maximum_radius) then
            write(*,'(a,2(f0.3,1x))') '**error(Image::ComputeAverageAndSigmaOfValuesInSpectrum): bad min, max radii: ', &
                                            minimum_radius, maximum_radius
            call this_program%TerminateWithFatalError('Image::ComputeAverageAndSigmaOfValuesInSpectrum','Bad min,max radii')
        endif

        ! Half-width of the cross
        ccross_half_width_sq = 4.0
        if (present(cross_half_width)) ccross_half_width_sq = cross_half_width**2

        call my_distribution%Init()
        do j=1,self%logical_dimensions(2)
            y_sq = (j-self%GetPhysicalIndexOfBoxCenter(2))**2
            if (y_sq .le. ccross_half_width_sq) cycle
            do i=1,self%logical_dimensions(1)
                x_sq = (i-self%GetPhysicalIndexOfBoxCenter(1))**2
                rad_sq = x_sq + y_sq
                if ( x_sq .gt. ccross_half_width_sq .and. rad_sq .gt. min_rad_sq .and. rad_sq .lt. max_rad_sq) then
                    call my_distribution%AddSampleValue(self%real_values(i,j,1))
                endif
            enddo
        enddo
        average = my_distribution%GetSampleMean()
        sigma   = sqrt(my_distribution%GetSampleVariance())
    end subroutine ComputeAverageAndSigmaOfValuesInSpectrum

    !>  \brief  If values on the central cross are too high, threshold them
    subroutine  SetMaximumValueOnCentralCross(self,maximum_value)
        class(Image),   intent(inout)   ::  self
        real,           intent(in)      ::  maximum_value
        !
        !
        self%real_values(self%GetPhysicalIndexOfBoxCenter(1),:,:) = min(self%real_values(self%GetPhysicalIndexOfBoxCenter(1),:,:),&
                                                                        maximum_value)
        self%real_values(:,self%GetPhysicalIndexOfBoxCenter(2),:) = min(self%real_values(:,self%GetPhysicalIndexOfBoxCenter(2),:),&
                                                                        maximum_value)
    end subroutine SetMaximumValueOnCentralCross

    !>  \brief  return the sum of real values.
    real function GetSumOfValues(self)
        class(image),   intent(in)  ::  self
        GetSumOfValues = self%GetSumOfValuesDouble()
    end function GetSumOfValues

    !>  \brief Return the average of real values at the edges of the image
    pure real function GetAverageOfValuesOnEdges(self)
        class(image),    intent(in)  ::  self
        ! private variables
        integer         ::  i,j,k
        real(kind=8)    ::  s
        integer(kind=8) ::  c
        ! start work
        s = 0.0d0
        c = 0
        do k=1,self%logical_dimensions(3),max(self%logical_dimensions(3)-1,1)
            do j=1,self%logical_dimensions(2),max(self%logical_dimensions(2)-1,1)
                do i=1,self%logical_dimensions(1),max(self%logical_dimensions(1)-1,1)
                    s = s + self%real_values(i,j,k)
                    c = c + 1
                enddo
            enddo
        enddo
        GetAverageOfValuesOnEdges = s/c
    end function GetAverageOfValuesOnEdges

    !>  \brief  Return the avereage of real values within a ring specified by two radii
    pure real function GetAverageOfValuesInRing(self,radius_inner,radius_outer)
        class(Image),   intent(in)      ::  self
        real,           intent(in)      ::  radius_inner    !<  In pixels
        real,           intent(in)      ::  radius_outer    !<  In pixels
        ! private variables
        integer         ::  i,j,k
        real(kind=8)    ::  s
        integer(kind=8) ::  c
        real            ::  rad_x, rad_y, rad_z, rad
        real            ::  radius_inner_sq, radius_outer_sq
        ! start work
        radius_inner_sq = radius_inner**2
        radius_outer_sq = radius_outer**2
        s = 0.0d0
        c = 0
        do k=1,self%logical_dimensions(3)
            rad_z = (k-self%physical_address_of_box_center(3))**2
            do j=1,self%logical_dimensions(2)
                rad_y = (j-self%physical_address_of_box_center(2))**2 + rad_z
                do i=1,self%logical_dimensions(1)
                    rad_x = (i-self%physical_address_of_box_center(1))**2
                    !
                    rad = rad_x + rad_y
                    !
                    if ( rad .ge. radius_inner_sq .and. rad .le. radius_outer_sq ) then
                        !
                        s = s + self%real_values(i,j,k)
                        c = c + 1
                    endif
                enddo
            enddo
        enddo
        GetAverageOfValuesInRing = s/c
    end function GetAverageOfValuesInRing

    !>  \brief  return the average of real values of an image
    real function GetAverageOfValues(self)
        ! argument
        class(image),   intent(in)  ::  self
        ! start work
        GetAverageOfValues = self%GetSumOfValues() / real(product(self%logical_dimensions))
    end function GetAverageOfValues

    !>  \brief  Return the sum of real values as a DOUBLE PRECISION value.
    !!
    !!  The option ignore_friedel_mates is meant for use when summing CTF^2 volumes, which are real but represent reciprocal-space information arranged in the FFTW manner,
    !!  whereby in the first dimension, Friedel mates are described explicitely except for the Nyquist frequency element.
    double precision function GetSumOfValuesDouble(self,ignore_friedel_mates)
        ! arguments
        class(Image),               intent(in)      ::  self
        logical,        optional,   intent(in)      ::  ignore_friedel_mates    !<  if .true., friedel mates (which are descibed explicitely in the first dimension of a CTF^2 volume) will be ignored when computing the sum
        ! private variables
        integer ::  i,j,k
        logical ::  iignore_friedel_mates
        ! start work

        ! are we ignoring friedel mates?
        if (present(ignore_friedel_mates)) then
            iignore_friedel_mates = .true.
        else
            iignore_friedel_mates = .false.
        endif

        if (iignore_friedel_mates) then
            ! the logic for this is only written for even dimension images, because i'm too lazy to think about odd dimensions right now.
            if (mod(self%logical_dimensions(2),2) .ne. 0. .or. mod(self%logical_dimensions(2),2) .ne. 0.) then
                write(*,'(2a,3(i0,1x))')    '**error(SumDouble): images with odd 2nd or 3rd dimensions ', &
                                            'are not supported. Dimensions = ', self%logical_dimensions
                call this_program%TerminateWithFatalError(  'Image::GetSumOfValuesDouble',   &
                                                            'Odd dimensions not supported with ignore Friedel mates')
            endif
            ! this ignore_friedel_mates options is only really meant for ctf^2 volumes with specific dimensions. perhaps this restriction can be lifted, if needed.
            if (self%logical_dimensions(2) .ne. self%logical_dimensions(3) .or. &
                self%logical_dimensions(1) .ne. self%logical_dimensions(2)/2+1) then
                write(*,'(a,3(i5,1x))') '**error(SumDouble): unexpected input image dimensions: ', self%logical_dimensions
                call this_program%TerminateWithFatalError(  'Image::GetSumOfValuesDouble',   &
                                                            'Unexpected input image dimensions')
            endif
        endif

        GetSumOfValuesDouble = 0.0d0
        if (self%is_in_real_space) then
            do i=1,self%logical_dimensions(1)
                do j=1, self%logical_dimensions(2)
                    if (iignore_friedel_mates .and. i == 1 .and. j .gt. 1 .and. j .lt. self%logical_dimensions(2)/2+1) cycle
                    do k=1, self%logical_dimensions(3)
                        if (iignore_friedel_mates .and. i == 1 .and. j == 1 .and. k .gt. 1 &
                                                  .and. k .lt. self%logical_dimensions(3)/2+1) then
                            cycle
                        endif
                        GetSumOfValuesDouble = GetSumOfValuesDouble + dble(self%real_values(i,j,k))
                    enddo
                enddo
            enddo
        else
            call this_program%TerminateWithFatalError(  'Image::GetSumOfValuesDouble','Cannot operate in Fourier space')
        endif
    end function GetSumOfValuesDouble

    !>  \brief  return the minimum value in an image
    real function GetMinimumValue(self)
        ! arguments
        class(image),   intent(in)  ::  self
        ! private variables
        ! start work
        if (self%is_in_real_space) then
            GetMinimumValue = minval(self%real_values(   1:self%logical_dimensions(1),  &
                                                    1:self%logical_dimensions(2),   &
                                                    1:self%logical_dimensions(3)))
        else
            call this_program%TerminateWithFatalError(  'Image::GetMinimumValue','Complex images are not supported')
        endif
    end function GetMinimumValue

    !>  \brief  Set the maximum value in the image. Any values above the given maximum will be reset to the maximum.
    pure subroutine SetMaximumValue(self,new_maximum)
        class(Image),   intent(inout)   ::  self
        real,           intent(in)      ::  new_maximum
        if (self%IsInRealSpace()) then
            self%real_values = min(self%real_values,new_maximum)
        !else
        !    call this_program%TerminateWithFatalError('Image::SetMaximumValue','Only works in real space')
        endif
    end subroutine SetMaximumValue

    !>  \brief  Set the maximum and minimum values in the image. Any values above the given maximum will be reset to the maximum.
    pure subroutine SetMinimumAndMaximumValue(self,new_minimum,new_maximum)
        class(Image),   intent(inout)   ::  self
        real,           intent(in)      ::  new_minimum,new_maximum
        if (self%IsInRealSpace()) then
            self%real_values = min(self%real_values,new_maximum)
            self%real_values = max(self%real_values,new_minimum)
        !else
        !    call this_program%TerminateWithFatalError('Image::SetMaximumValue','Only works in real space')
        endif
    end subroutine SetMinimumAndMaximumValue

    !>  \brief  return the minimum value in an image
    real function GetMaximumValue(self,minimum_distance_from_center,minimum_distance_from_edge)
        ! arguments
        class(image),               intent(in)  ::  self
        integer,        optional,   intent(in)  ::  minimum_distance_from_center    !<  Only voxels at least this far away from the center in every dimension will be considered
        integer,        optional,   intent(in)  ::  minimum_distance_from_edge
        ! private variable
        integer ::  i,j,k
        integer ::  i_logi,j_logi,k_logi
        ! start work
        if (self%is_in_real_space) then
            GetMaximumValue = -huge(1.0e0)
            do k=1,self%logical_dimensions(3)
                k_logi = abs(k-self%physical_address_of_box_center(3))
                if (present(minimum_distance_from_center) .and. self%logical_dimensions(3) .gt. 1) then
                    if (k_logi .lt. minimum_distance_from_center) cycle
                endif
                if (present(minimum_distance_from_edge) .and. self%logical_dimensions(3) .gt. 1) then
                    if (k .lt. minimum_distance_from_edge) cycle
                    if (k .gt. self%logical_dimensions(3)-minimum_distance_from_edge) cycle
                endif
                do j=1,self%logical_dimensions(2)
                    j_logi = abs(j-self%physical_address_of_box_center(2))
                    if (present(minimum_distance_from_center)) then
                        if (j_logi .lt. minimum_distance_from_center) cycle
                    endif
                    if (present(minimum_distance_from_edge)) then
                        if (j .lt. minimum_distance_from_edge) cycle
                        if (j .gt. self%logical_dimensions(2)-minimum_distance_from_edge) cycle
                    endif
                    do i=1,self%logical_dimensions(1)
                        i_logi = abs(i-self%physical_address_of_box_center(1))
                        if (present(minimum_distance_from_center)) then
                            if (i_logi .lt. minimum_distance_from_center) cycle
                        endif
                        if (present(minimum_distance_from_edge)) then
                            if (i .lt. minimum_distance_from_edge) cycle
                            if (i .gt. self%logical_dimensions(1)-minimum_distance_from_edge) cycle
                        endif
                        !
                        if (self%real_values(i,j,k) .gt. GetMaximumValue) GetMaximumValue = self%real_values(i,j,k)
                    enddo
                enddo
            enddo
        else
            call this_program%TerminateWithFatalError(  'Image::GetMaximumValue','Complex images are not supported')
        endif
    end function GetMaximumValue

    !>  \brief  Return a histogram of real values in the image
    subroutine ComputeHistogramOfValues(self,my_histogram,number_of_bins)
        use Histograms
        class(Image),               intent(in)      ::  self
        type(Histogram),            intent(inout)   ::  my_histogram
        integer,                    intent(in)      ::  number_of_bins      !<  Number of bins for the histogram
        ! private variables
        integer ::  i,j,k
        ! start work

        call my_histogram%Init(self%GetMinimumValue(),self%GetMaximumValue(),number_of_bins,discard_extreme_values=.false.)
        do k=1,self%logical_dimensions(3)
            do j=1,self%logical_dimensions(2)
                do i=1,self%logical_dimensions(1)
                    call my_histogram%AddSampleValue(self%real_values(i,j,k))
                enddo
            enddo
        enddo
    end subroutine ComputeHistogramOfValues


    !>  \brief  return the maximum radius possible in an image in either X or Y.
    !!  The diagonal is NOT considered
    pure function GetMaximumRadius(self) result(maximum_radius)

        ! Arguments
        class(image),   intent(in)  ::  self

        ! Variables

        real                        ::  maximum_radius

        ! Start Work
        maximum_radius = maxval(self%physical_address_of_box_center) - 1


    end function GetMaximumRadius

    !>  \brief  return the maximum radius possible along a diagonal
    pure function GetMaximumDiagonalRadius(self)  result(maximum_radius)

        ! Arguments

        class(image),   intent(in)  ::  self

        ! Variables

        real                        ::  maximum_radius  !<  In pixels (real-space) or reciprocal pixels (Fourier space)

        ! Start work
        if (self%IsInRealSpace()) then
            maximum_radius = sqrt(real(sum((self%physical_address_of_box_center-1)**2)))
        else
            maximum_radius = sqrt(  (self%logical_upper_bound_complex(1)*self%fourier_voxel_size(1))**2 &
                                +   (self%logical_upper_bound_complex(2)*self%fourier_voxel_size(2))**2 &
                                +   (self%logical_upper_bound_complex(3)*self%fourier_voxel_size(3))**2)
        endif


    end function GetMaximumDiagonalRadius


    !>  \brief  return the maximum radius possible along a diagonal
    pure function GetRadiiGivenFraction(self, fraction_wanted)  result(returned_radii)

        ! Arguments

        class(image),   intent(in)  ::  self
        real,           intent(in)  ::  fraction_wanted

        ! Variables

        real                        ::  returned_radii(3)

        ! Start work

        returned_radii = (self%physical_address_of_box_center-1) * fraction_wanted

    end function GetRadiiGivenFraction

    !>  \brief  Return the variance of the real values of an image
    !!  Note: this algorithm is not precise and should not be used. See
    !!  http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    !!  \todo Replace with recommended robust algorithm
    pure real function GetVarianceOfValues(self)
        ! argument
        class(image),   intent(in)  ::  self
        ! private variables
        integer ::  i,j,k
        real(kind=8)    ::  sum, sum2
        ! start work
        sum = 0.0d0
        sum2 = 0.0d0
        do k=1,self%logical_dimensions(3)
            do j=1,self%logical_dimensions(2)
                do i=1,self%logical_dimensions(1)
                    sum = sum + self%real_values(i,j,k)
                    sum2 = sum2 + self%real_values(i,j,k)**2
                enddo
            enddo
        enddo
        GetVarianceOfValues = (sum2 - (sum**2/product(self%logical_dimensions)))/product(self%logical_dimensions)
    end function GetVarianceOfValues

    !>  \brief  Return the standard deviation of the real values of an image
    pure real function GetSigmaOfValues(self)
        implicit none
        ! argument
        class(image),   intent(in)  ::  self
        ! start work
        GetSigmaOfValues = sqrt(self%GetVarianceOfValues())
    end function GetSigmaOfValues

    !>  \brief  Check whether any of the volume voxels are NaN
    pure logical function HasNan(self)
#ifndef __GFORTRAN__
        use, intrinsic :: ieee_arithmetic
#endif
        ! arguments
        class(image),   intent(in)  ::  self
        ! private variables
        integer ::  i,j,k
        ! start work
        HasNan = .false.
        if (self%is_in_real_space) then
            do k=1,self%logical_dimensions(3)
                do j=1,self%logical_dimensions(2)
                    do i=1,self%logical_dimensions(1)
#ifndef __GFORTRAN__
                        if (ieee_is_nan(self%real_values(i,j,k))) then
#else
                        if (isnan(self%real_values(i,j,k))) then
#endif
                            HasNan = .true.
                            return
                        endif
                    enddo
                enddo
            enddo
        else
            do k=1,self%physical_upper_bound_complex(3)
                do j=1,self%physical_upper_bound_complex(2)
                    do i=1,self%physical_upper_bound_complex(1)
#ifndef __GFORTRAN__
                        if (ieee_is_nan(aimag(self%complex_values(i,j,k))) .or. ieee_is_nan(real(self%complex_values(i,j,k)))) then
#else
                        if (isnan(aimag(self%complex_values(i,j,k))) .or. isnan(real(self%complex_values(i,j,k)))) then
#endif
                            HasNan = .true.
                            return
                        endif
                    enddo
                enddo
            enddo
        endif
    end function HasNan

    !>  \brief  determine whether an image has been binarised (meaning it only contains 1's and 0's)
    logical function IsBinarised(self)
        class(image),    intent(in)  ::  self
        ! Private variables
        integer ::  i,j,k
        ! Start work
        ! Only makes sense for real images
        if (.not. self%is_in_real_space) then
            call this_program%TerminateWithFatalError('Image::IsBinarised','Image is in Fourier space')
        endif
        !
        IsBinarised = .true.
        ! if any value is other than 0.0 or 1.0, then not binarised
        do k=1,self%logical_dimensions(3)
            do j=1,self%logical_dimensions(2)
                do i=1,self%logical_dimensions(1)
                    if (    self%real_values(i,j,k) .ne. 0.0 .and. &
                            abs(self%real_values(i,j,k)-1.0) .gt. tiny(self%real_values(i,j,k))) then
                        IsBinarised = .false.
                        return
                    endif
                enddo
            enddo
        enddo
    end function IsBinarised

    !>  \brief  Determine whether an image is blank/zero
    pure logical function IsConstant(self)
        class(image),    intent(in)  ::  self
        ! start work
        if (self%is_in_real_space) then
            IsConstant = all(self%real_values(  1:self%logical_dimensions(1),   &
                                                1:self%logical_dimensions(2),   &
                                                1:self%logical_dimensions(3) )  &
                                                .eq. self%real_values(1,1,1) )
        else
            IsConstant = all(self%complex_values .eq. self%complex_values(1,1,1))
        endif
    end function IsConstant

    !>  \brief  Determine whether the image is square / cubic
    pure logical function IsSquare(self)
        class(image),   intent(in)  ::  self
        ! Start work
        IsSquare = self%logical_dimensions(1) .eq. self%logical_dimensions(2)
        if (self%IsAVolume()) IsSquare = IsSquare .and. (self%logical_dimensions(1) .eq. self%logical_dimensions(3))
    end function IsSquare

    !>  \brief  Determine whether the image is a volume (ie three-dimensional)
    pure logical function IsAVolume(self)
        class(image),   intent(in)  ::  self
        IsAVolume = all(self%logical_dimensions .gt. 1)
    end function IsAVolume

    !>  \brief  Are all logical dimensions even?
    pure logical function HasEvenDimensions(self)
        use UsefulFunctions, only : IsEven
        class(Image),   intent(in)  ::  self
        if (self%IsAVolume()) then
            HasEvenDimensions = all(IsEven(self%logical_dimensions))
        else
            HasEvenDimensions = all(IsEven(self%logical_dimensions(1:2)))
        endif
    end function HasEvenDimensions

    !>  \brief  Does the image have the same dimensions as another image
    pure logical function HasSameDimensionsAs(self,other)
        class(image),   intent(in)  ::  self
        class(image),   intent(in)  ::  other
        HasSameDimensionsAs = all(self%logical_dimensions .eq. other%logical_dimensions)
    end function HasSameDimensionsAs

    !>  \brief  Is the image in the same space (real or Fourier) as another image
    pure logical function IsInSameSpaceAs(self,other)
        class(image),   intent(in)  ::  self
        class(image),   intent(in)  ::  other
        IsInSameSpaceAs = (self%is_in_real_space .and. other%is_in_real_space) .or. &
                          (.not. self%is_in_real_space .and. .not. other%is_in_real_space)
    end function IsInSameSpaceAs

    !>  \brief  Get a 2d slice from a 3d volume along one of the axes
    subroutine GetSlice(self,slice,axis,address)
        ! Arguments
        class(image),   intent(in)      ::  self        !<  3d image
        type(image),    intent(inout)   ::  slice       !<  2d image
        integer,        intent(in)      ::  axis        !<  a slice perpendicular to this axis will be extracted
        integer,        intent(in)      ::  address     !<  the index of the slice along the chosen axis
        ! Start work

        ! Check the input image is 3d
        if (.not. self%IsAVolume()) call this_program%TerminateWithFatalError('Image::GetSlice','Input image is not a volume')

        ! Check axis number is 1, 2 or 3
        if (axis .lt. 1 .or. axis .gt. 3) then
            call this_program%TerminateWithFatalError('Image::GetSlice','Axis number must be between 1 and 3')
        endif

        ! Deallocate output image if necessary
        call slice%Deallocate()

        ! Set output image dimensions
        slice%logical_dimensions = 1
        select case (axis)
            case (1)
                slice%logical_dimensions(1) = self%logical_dimensions(2)
                slice%logical_dimensions(2) = self%logical_dimensions(3)
            case (2)
                slice%logical_dimensions(1) = self%logical_dimensions(1)
                slice%logical_dimensions(2) = self%logical_dimensions(3)
            case (3)
                slice%logical_dimensions(1) = self%logical_dimensions(1)
                slice%logical_dimensions(2) = self%logical_dimensions(2)
        end select

        ! Allocate output image
        call slice%Allocate()

        ! copy pixel data
        select case(axis)
            case (1)
                slice%real_values(1:slice%logical_dimensions(1),1:slice%logical_dimensions(2),1) &
                = self%real_values(address,1:self%logical_dimensions(2),1:self%logical_dimensions(3))
            case (2)
                slice%real_values(1:slice%logical_dimensions(1),1:slice%logical_dimensions(2),1) &
                = self%real_values(1:self%logical_dimensions(1),address,1:self%logical_dimensions(3))
            case (3)
                slice%real_values(1:slice%logical_dimensions(1),1:slice%logical_dimensions(2),1) &
                = self%real_values(1:self%logical_dimensions(1),1:self%logical_dimensions(2),address)
        end select
    end subroutine GetSlice

    !>  \brief  Extract a slice through the volume at an arbitrary orientation
    subroutine ExtractComplexSlice(self,euler_angles,slice,volume_z_shift,slice_x_shift,slice_y_shift,maximum_frequency)
        use OrientationsAndRotations, only : EulerAngles
        ! Arguments
        class(Image),       intent(in)      ::  self                !<  Must be in Fourier space
        type(EulerAngles),  intent(in)      ::  euler_angles
        type(Image),        intent(inout)   ::  slice
        real,               intent(in)      ::  volume_z_shift      !<  In pixels, the volume will be shifted by this amount on along the z axis before slice extraction (or this is what it'll look like)
        real,               intent(in)      ::  slice_x_shift       !<  In pixels, the slice will be shifted by this amount along the first dimension (after slice extraction)
        real,               intent(in)      ::  slice_y_shift       !<  In pixels, the slice will be shifted by this amount along the second dimension (after slice extraction)
        real,               intent(in)      ::  maximum_frequency   !<  In reciprocal pixels
        ! Private variables
        real    ::  spatial_frequency_in_slice(4)
        real    ::  spatial_frequency_in_volume(4)
        integer ::  i,j
        real    ::  y_sq,rad_sq
        real    ::  rotation_matrix(4,4)
        real    ::  total_phase_shift
        real    ::  volume_z_phase_shift
        real    ::  slice_x_phase_shift
        real    ::  slice_y_phase_shift
        real    ::  maximum_squared_frequency
        ! Start work

        if (.not. slice%IsAllocated()) then
            call this_program%TerminateWithFatalError('Image::ExtractComplexSlice','Slice not allocated')
        endif
        if (.not. self%IsInFourierSpace()) then
            call this_program%TerminateWithFatalError('Image::ExtractComplexSlice','Volume not in Fourier space')
        endif

        !
        slice%is_in_real_space = self%is_in_real_space
        slice%object_is_centered_in_box = self%object_is_centered_in_box

        ! Phase shifts
        volume_z_phase_shift = volume_z_shift * 2.0 * pi
        slice_x_phase_shift  = slice_x_shift  * 2.0 * pi
        slice_y_phase_shift  = slice_y_shift  * 2.0 * pi

        ! Maximum frequency
        maximum_squared_frequency = maximum_frequency**2

        ! Get the rotation matrix
        rotation_matrix = euler_angles%GetRotationMatrix(reverse=.true.)

        ! Loop over the slice
        spatial_frequency_in_slice = 0.0
        do j=1,slice%physical_upper_bound_complex(2)
            spatial_frequency_in_slice(2) = slice%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2) * slice%fourier_voxel_size(2)
            y_sq = spatial_frequency_in_slice(2)**2
            do i=1,slice%physical_upper_bound_complex(1)
                spatial_frequency_in_slice(1) = real(i-1)*slice%fourier_voxel_size(1)
                rad_sq = spatial_frequency_in_slice(1)**2 + y_sq
                if (rad_sq .lt. maximum_squared_frequency) then
                    spatial_frequency_in_volume = matmul(rotation_matrix,spatial_frequency_in_slice)
                    call self%GetComplexValueByLinearInterpolationNoBoundsCheckVolume(&
                                                    slice%complex_values(i,j,1), &
                                                    x=spatial_frequency_in_volume(1)*self%logical_dimensions(1), &
                                                    y=spatial_frequency_in_volume(2)*self%logical_dimensions(2), &
                                                    z=spatial_frequency_in_volume(3)*self%logical_dimensions(3) )
                    total_phase_shift = - volume_z_phase_shift* spatial_frequency_in_volume(3) &
                                        - slice_x_phase_shift * spatial_frequency_in_slice(1) &
                                        - slice_y_phase_shift * spatial_frequency_in_slice(2)
                    slice%complex_values(i,j,1) = slice%complex_values(i,j,1) &
                                                * cmplx(cos(total_phase_shift),sin(total_phase_shift))
                else
                    slice%complex_values(i,j,1) = (0.0,0.0)
                endif
            enddo
        enddo

    end subroutine ExtractComplexSlice

    !>  \brief  Insert a slice through the volume at an arbitrary orientation
    !!
    !!          The rotation matrix can be obtained from Euler angles using:
    !!          rotation_matrix = euler_angles%GetRotationMatrix(reverse=.true.)
    pure subroutine InsertComplexSlice(self,rotation_matrix,slice,volume_z_shift,slice_x_shift,slice_y_shift, &
                                        maximum_frequency,number_of_inserted_pixels)
        use OrientationsAndRotations, only : EulerAngles
        ! Arguments
        class(Image),       intent(inout)   ::  self                        !<  Must be in Fourier space
        real,               intent(in)      ::  rotation_matrix(4,4)
        type(Image),        intent(in)      ::  slice                       !<  Must be in Fourier space
        real,               intent(in)      ::  volume_z_shift              !<  In pixels, the volume will be shifted by this amount on along the z axis before slice insertion (or this is what it'll look like)
        real,               intent(in)      ::  slice_x_shift               !<  In pixels, the slice will be shifted by this amount along the first dimension (before slice insertion)
        real,               intent(in)      ::  slice_y_shift               !<  In pixels, the slice will be shifted by this amount along the second dimension (before slice insertion)
        real,               intent(in)      ::  maximum_frequency           !<  In reciprocal pixels
        integer,            intent(out)     ::  number_of_inserted_pixels   !<  Number of Fourier pixels from the slice inserted into the volume
        ! Private variables
        real    ::  spatial_frequency_in_slice(4)
        real    ::  spatial_frequency_in_volume(4)
        integer ::  i,j
        real    ::  y_sq,rad_sq
        real    ::  total_phase_shift
        real    ::  volume_z_phase_shift
        real    ::  slice_x_phase_shift
        real    ::  slice_y_phase_shift
        real    ::  maximum_squared_frequency
        ! Start work

        ! Phase shifts
        volume_z_phase_shift = volume_z_shift * 2.0 * pi
        slice_x_phase_shift  = slice_x_shift  * 2.0 * pi
        slice_y_phase_shift  = slice_y_shift  * 2.0 * pi

        ! Maximum frequency
        maximum_squared_frequency = maximum_frequency**2

        ! We will count how many slice pixels are inserted into the volume
        number_of_inserted_pixels = 0


        ! Loop over the slice
        spatial_frequency_in_slice(1:3) = 0.0
        spatial_frequency_in_slice(4)   = 1.0
        do j=1,slice%physical_upper_bound_complex(2)
            spatial_frequency_in_slice(2) = slice%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2) * slice%fourier_voxel_size(2)
            y_sq = spatial_frequency_in_slice(2)**2
            do i=1,slice%physical_upper_bound_complex(1)
                if (slice%FourierComponentIsRedundant(i,j,1)) cycle
                spatial_frequency_in_slice(1) = real(i-1)*slice%fourier_voxel_size(1)
                rad_sq = spatial_frequency_in_slice(1)**2 + y_sq
                if (rad_sq .lt. maximum_squared_frequency) then
                    number_of_inserted_pixels = number_of_inserted_pixels + 1
                    spatial_frequency_in_volume = matmul(rotation_matrix,spatial_frequency_in_slice)
                    total_phase_shift = - volume_z_phase_shift* spatial_frequency_in_volume(3) &
                                        - slice_x_phase_shift * spatial_frequency_in_slice(1) &
                                        - slice_y_phase_shift * spatial_frequency_in_slice(2)
                    call self%PutComplexValueByLinearInterpolationNoBoundsCheckVolume(&
                                                    slice%complex_values(i,j,1) &
                                                    * cmplx(cos(total_phase_shift),sin(total_phase_shift)), &
                                                    x=spatial_frequency_in_volume(1) * self%logical_dimensions(1), &
                                                    y=spatial_frequency_in_volume(2) * self%logical_dimensions(2), &
                                                    z=spatial_frequency_in_volume(3) * self%logical_dimensions(3))
                endif
            enddo
        enddo

    end subroutine InsertComplexSlice

    !>  \brief  Insert a slice through the volume at an arbitrary orientation
    !!
    !!          The rotation matrix can be obtained from Euler angles using:
    !!          rotation_matrix = euler_angles%GetRotationMatrix(reverse=.true.)
    !!
    !!  \todo   It may be more efficient to only insert slice within inner radius? That way we would save the matmuls for all the pixels beyond the inner radius
    pure subroutine InsertRealSlice(self,rotation_matrix,slice)
        use OrientationsAndRotations, only : EulerAngles
        ! Arguments
        class(Image),       intent(inout)   ::  self                !<  Must be in real space
        real,               intent(in)      ::  rotation_matrix(4,4)
        type(Image),        intent(in)      ::  slice               !<  Must be in real space
        ! Private variables
        real    ::  coordinates_in_slice(4)
        real    ::  coordinates_in_volume(4)
        integer ::  i,j
        real    ::  y_sq,rad_sq
        real    ::  physical_address_of_volume_center(3)
        real    ::  self_logical_dimensions_real(3)
        ! Start work

        ! We do the integer to real conversion just once
        physical_address_of_volume_center = self%physical_address_of_box_center
        self_logical_dimensions_real = real(self%logical_dimensions)

        ! Loop over the slice
        coordinates_in_slice = 0.0
        do j=1,slice%logical_dimensions(2)
            coordinates_in_slice(2) = real(j-slice%physical_address_of_box_center(2))
            y_sq = coordinates_in_slice(2)**2
            do i=1,slice%logical_dimensions(1)
                coordinates_in_slice(1) = real(i-slice%physical_address_of_box_center(1))
                rad_sq = coordinates_in_slice(1)**2 + y_sq

                coordinates_in_volume = matmul(rotation_matrix,coordinates_in_slice)
                coordinates_in_volume(1:3) = coordinates_in_volume(1:3) + physical_address_of_volume_center(1:3)

                if (    any(coordinates_in_volume(1:3) .lt. 1.01) &
                   .or. any(coordinates_in_volume(1:3) .gt. self_logical_dimensions_real - 0.01)) cycle


                call self%PutRealValueByLinearInterpolationNoBoundsCheckVolume(&
                                                slice%real_values(i,j,1),   &
                                                x=coordinates_in_volume(1), &
                                                y=coordinates_in_volume(2), &
                                                z=coordinates_in_volume(3))
            enddo
        enddo

    end subroutine InsertRealSlice

    !>  \brief  Extract a slice through the volume at an arbitrary orientation (attempt at vectorization)
    subroutine ExtractComplexSliceVector(self,euler_angles,slice,volume_z_shift,slice_x_shift,slice_y_shift,maximum_frequency)
        use OrientationsAndRotations, only : EulerAngles
        ! Arguments
        class(Image),       intent(in)      ::  self                !<  Must be in Fourier space
        type(EulerAngles),  intent(in)      ::  euler_angles
        type(Image),        intent(inout)   ::  slice
        real,               intent(in)      ::  volume_z_shift      !<  In pixels, the volume will be shifted by this amount on along the z axis before slice extraction (or this is what it'll look like)
        real,               intent(in)      ::  slice_x_shift       !<  In pixels, the slice will be shifted by this amount along the first dimension (after slice extraction)
        real,               intent(in)      ::  slice_y_shift       !<  In pixels, the slice will be shifted by this amount along the second dimension (after slice extraction)
        real,               intent(in)      ::  maximum_frequency   !<  In reciprocal pixels
        ! Private variables
        real    ::  spatial_frequency_in_slice(3,slice%physical_upper_bound_complex(1))
        real    ::  spatial_frequency_in_volume(3,slice%physical_upper_bound_complex(1))
        integer ::  i,j
        real    ::  y_sq
        real    ::  rotation_matrix(3,3),rotation_matrix_tmp(4,4)
        real    ::  total_phase_shift(slice%physical_upper_bound_complex(1))
        real    ::  volume_z_phase_shift
        real    ::  slice_x_phase_shift
        real    ::  slice_y_phase_shift
        real    ::  maximum_squared_frequency
        integer ::  physical_upper_bound_x
        real    ::  logical_coordinates_in_volume(3,slice%physical_upper_bound_complex(1))
        real    ::  logical_dimensions_r(3)
        ! Start work

        if (.not. slice%IsAllocated()) then
            call this_program%TerminateWithFatalError('Image::ExtractComplexSlice','Slice not allocated')
        endif
        if (.not. self%IsInFourierSpace()) then
            call this_program%TerminateWithFatalError('Image::ExtractComplexSlice','Volume not in Fourier space')
        endif


        !
        slice%is_in_real_space = self%is_in_real_space
        slice%object_is_centered_in_box = self%object_is_centered_in_box

        !
        logical_dimensions_r = real(self%logical_dimensions)

        ! Phase shifts
        volume_z_phase_shift = volume_z_shift * 2.0 * pi
        slice_x_phase_shift  = slice_x_shift  * 2.0 * pi
        slice_y_phase_shift  = slice_y_shift  * 2.0 * pi

        ! Maximum frequency
        maximum_squared_frequency = maximum_frequency**2

        ! Get the rotation matrix
        rotation_matrix_tmp = euler_angles%GetRotationMatrix(reverse=.true.)
        rotation_matrix = rotation_matrix_tmp(1:3,1:3)

        ! Loop over the slice
        spatial_frequency_in_slice = 0.0
        do j=1,slice%physical_upper_bound_complex(2)
            spatial_frequency_in_slice(2,:) = slice%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2) * slice%fourier_voxel_size(2)
            y_sq = spatial_frequency_in_slice(2,1)**2
            if (y_sq .gt. maximum_squared_frequency) cycle
            ! Work out how many pixels there are in the current row below the maximum frequency
            physical_upper_bound_x = min(int(sqrt(maximum_squared_frequency - y_sq)/slice%fourier_voxel_size(1)+1.0), &
                                         slice%physical_upper_bound_complex(1))
            ! Loop over the pixels in this row
            do i=1,physical_upper_bound_x
                spatial_frequency_in_slice(1,i)   = real(i-1)*slice%fourier_voxel_size(1)
                spatial_frequency_in_volume(:,i)= matmul(rotation_matrix,spatial_frequency_in_slice(:,i))
                logical_coordinates_in_volume(:,i)  = spatial_frequency_in_volume(:,i) * logical_dimensions_r
                total_phase_shift(i) = - volume_z_phase_shift* spatial_frequency_in_volume(3,i) &
                                    - slice_x_phase_shift * spatial_frequency_in_slice(1,i) &
                                    - slice_y_phase_shift * spatial_frequency_in_slice(2,i)


            enddo
            call self%GetComplexValueByLinearInterpolationNoBoundsCheckVolumeVector(&
                                                slice%complex_values, &
                                                physical_upper_bound_x, &
                                                j, &
                                                1, &
                                                logical_coordinates_in_volume)
            do i=1,physical_upper_bound_x
                slice%complex_values(i,j,1) = slice%complex_values(i,j,1) &
                                            * cmplx(cos(total_phase_shift(i)),sin(total_phase_shift(i)))
            enddo
            do i=physical_upper_bound_x+1,slice%physical_upper_bound_complex(1)
                slice%complex_values(i,j,1) = (0.0,0.0)
            enddo
        enddo

    end subroutine ExtractComplexSliceVector

    !>  \brief  Get a value from a non-integer physical "address" (i.e. where the first voxel is 1.0,1.0,1.0), using linear interpolation
    pure subroutine GetRealValueByLinearInterpolationNoBoundsCheckVolume(self,interpolated_value,x,y,z)
        ! Arguments
        class(Image),       intent(in)  ::  self
        real,               intent(out) ::  interpolated_value
        real,               intent(in)  ::  x
        real,               intent(in)  ::  y
        real,               intent(in)  ::  z
        ! Private variables
        integer ::  i_start,j_start,k_start
        real    ::  x_dist,y_dist,z_dist
        real    ::  x_dist_m,y_dist_m,z_dist_m
        ! Start work
        i_start = int(x)
        j_start = int(y)
        k_start = int(z)
        x_dist  = x-real(i_start,kind=4)
        y_dist  = y-real(j_start,kind=4)
        z_dist  = z-real(k_start,kind=4)
        x_dist_m = 1.0e0 - x_dist
        y_dist_m = 1.0e0 - y_dist
        z_dist_m = 1.0e0 - z_dist
        !
        ! TODO:
        !   - take care of edges of image
        !   - convert to taking (logical) coordinates as input
        !
        !
        interpolated_value =  self%real_values(i_start  ,j_start  ,k_start  ) * x_dist_m * y_dist_m * z_dist_m &
                            + self%real_values(i_start+1,j_start  ,k_start  ) * x_dist   * y_dist_m * z_dist_m &
                            + self%real_values(i_start  ,j_start+1,k_start  ) * x_dist_m * y_dist   * z_dist_m &
                            + self%real_values(i_start+1,j_start+1,k_start  ) * x_dist   * y_dist   * z_dist_m &
                            + self%real_values(i_start  ,j_start  ,k_start+1) * x_dist_m * y_dist_m * z_dist   &
                            + self%real_values(i_start+1,j_start  ,k_start+1) * x_dist   * y_dist_m * z_dist   &
                            + self%real_values(i_start  ,j_start+1,k_start+1) * x_dist_m * y_dist   * z_dist   &
                            + self%real_values(i_start+1,j_start+1,k_start+1) * x_dist   * y_dist   * z_dist

    end subroutine GetRealValueByLinearInterpolationNoBoundsCheckVolume

    !>  \brief  Get a value from a non-integer physical "address" (i.e. where the first voxel is 1.0,1.0,1.0), using linear interpolation
    pure subroutine GetRealValueByLinearInterpolationNoBoundsCheckImage(self,interpolated_value,x,y)
        ! Arguments
        class(Image),       intent(in)  ::  self
        real,               intent(out) ::  interpolated_value
        real,               intent(in)  ::  x
        real,               intent(in)  ::  y

        ! Private variables
        integer ::  i_start,j_start
        real    ::  x_dist,y_dist
        real    ::  x_dist_m,y_dist_m
        ! Start work
        i_start = int(x)
        j_start = int(y)

        x_dist  = x-real(i_start,kind=4)
        y_dist  = y-real(j_start,kind=4)

        x_dist_m = 1.0e0 - x_dist
        y_dist_m = 1.0e0 - y_dist


        !
        interpolated_value =  self%real_values(i_start  ,j_start  , 1 ) * x_dist_m * y_dist_m &
                            + self%real_values(i_start+1,j_start  , 1 ) * x_dist   * y_dist_m &
                            + self%real_values(i_start  ,j_start+1, 1 ) * x_dist_m * y_dist   &
                            + self%real_values(i_start+1,j_start+1, 1 ) * x_dist   * y_dist

    end subroutine GetRealValueByLinearInterpolationNoBoundsCheckImage

    !>  \brief  Put/spread a value from a non-integer physical "address" (i.e. where the first voxel is 1.0,1.0,1.0), using linear interpolation, into the nearest pixels
    pure subroutine PutRealValueByLinearInterpolationNoBoundsCheckImage(self,value_to_add,x,y)
        ! Arguments
        class(Image),       intent(inout)   ::  self
        real,               intent(in)      ::  value_to_add
        real,               intent(in)      ::  x
        real,               intent(in)      ::  y

        ! Private variables
        integer ::  i_start,j_start
        real    ::  x_dist,y_dist
        real    ::  x_dist_m,y_dist_m
        ! Start work
        i_start = int(x)
        j_start = int(y)

        x_dist  = x-real(i_start,kind=4)
        y_dist  = y-real(j_start,kind=4)

        x_dist_m = 1.0e0 - x_dist
        y_dist_m = 1.0e0 - y_dist

        !
        self%real_values(i_start  ,j_start  , 1 ) = self%real_values(i_start  ,j_start  , 1 ) &
                                                    + value_to_add  * x_dist_m * y_dist_m
        self%real_values(i_start+1,j_start  , 1 ) = self%real_values(i_start+1,j_start  , 1 ) &
                                                    + value_to_add  * x_dist   * y_dist_m
        self%real_values(i_start  ,j_start+1, 1 ) = self%real_values(i_start  ,j_start+1, 1 ) &
                                                    + value_to_add  * x_dist_m * y_dist
        self%real_values(i_start+1,j_start+1, 1 ) = self%real_values(i_start+1,j_start+1, 1 ) &
                                                    + value_to_add  * x_dist   * y_dist

    end subroutine PutRealValueByLinearInterpolationNoBoundsCheckImage

    !>  \brief  Put/spread a value from a non-integer physical "address" (i.e. where the first voxel is 1.0,1.0,1.0), using linear interpolation, into the nearest pixels
    pure subroutine PutRealValueByLinearInterpolationNoBoundsCheckVolume(self,value_to_add,x,y,z)
        ! Arguments
        class(Image),       intent(inout)   ::  self
        real,               intent(in)      ::  value_to_add
        real,               intent(in)      ::  x
        real,               intent(in)      ::  y
        real,               intent(in)      ::  z

        ! Private variables
        integer ::  i_start,j_start,k_start
        real    ::  x_dist,y_dist,z_dist
        real    ::  x_dist_m,y_dist_m,z_dist_m
        ! Start work
        i_start = int(x)
        j_start = int(y)
        k_start = int(z)

        x_dist  = x-real(i_start,kind=4)
        y_dist  = y-real(j_start,kind=4)
        z_dist  = z-real(k_start,kind=4)

        x_dist_m = 1.0e0 - x_dist
        y_dist_m = 1.0e0 - y_dist
        z_dist_m = 1.0e0 - z_dist

        !
        self%real_values(i_start  ,j_start  , k_start  )    = self%real_values(i_start  ,j_start  , k_start ) &
                                                            + value_to_add  * x_dist_m * y_dist_m * z_dist_m
        self%real_values(i_start+1,j_start  , k_start  )    = self%real_values(i_start+1,j_start  , k_start ) &
                                                            + value_to_add  * x_dist   * y_dist_m * z_dist_m
        self%real_values(i_start  ,j_start+1, k_start  )    = self%real_values(i_start  ,j_start+1, k_start ) &
                                                            + value_to_add  * x_dist_m * y_dist   * z_dist_m
        self%real_values(i_start+1,j_start+1, k_start  )    = self%real_values(i_start+1,j_start+1, k_start ) &
                                                            + value_to_add  * x_dist   * y_dist   * z_dist_m
        self%real_values(i_start  ,j_start  , k_start+1)    = self%real_values(i_start  ,j_start  , k_start+1 ) &
                                                            + value_to_add  * x_dist_m * y_dist_m * z_dist
        self%real_values(i_start+1,j_start  , k_start+1)    = self%real_values(i_start+1,j_start  , k_start+1 ) &
                                                            + value_to_add  * x_dist   * y_dist_m * z_dist
        self%real_values(i_start  ,j_start+1, k_start+1)    = self%real_values(i_start  ,j_start+1, k_start+1 ) &
                                                            + value_to_add  * x_dist_m * y_dist   * z_dist
        self%real_values(i_start+1,j_start+1, k_start+1)    = self%real_values(i_start+1,j_start+1, k_start+1 ) &
                                                            + value_to_add  * x_dist   * y_dist   * z_dist

    end subroutine PutRealValueByLinearInterpolationNoBoundsCheckVolume

    !>  \brief  Get a complex value from logical coordinates (i.e. origin at 0.0,0.0,0.0), using linear interpolation
    !!
    !!  \todo   This can be optimised further: do not do the Hermitian insertions. They can be done once all slices have been inserted.
    pure subroutine PutComplexValueByLinearInterpolationNoBoundsCheckVolume(self,value_to_add,x,y,z)
        ! Arguments
        class(Image),       intent(inout)   ::  self
        complex,            intent(in)      ::  value_to_add
        real,               intent(in)      ::  x                       !<  Logical coordinate, first dimension
        real,               intent(in)      ::  y                       !<  Logical coordinate, second dimension
        real,               intent(in)      ::  z                       !<  Logical coordinate, third dimension
        ! Private variables
        integer ::  i
        integer ::  i_start,j_start,k_start
        real    ::  x_dist,y_dist,z_dist
        real    ::  x_dist_m,y_dist_m,z_dist_m
        integer ::  logical_address(3,8), physical_address(3) ! The neighbours are numbered 1 to 8
        integer ::  hermitian_physical_address(3)
        real    ::  weight(8)
        ! Start work
        i_start  = floor(x)
        j_start  = floor(y)
        k_start  = floor(z)
        z_dist   = z-real(k_start)
        z_dist_m = 1.0e0 - z_dist
        y_dist   = y-real(j_start)
        y_dist_m = 1.0e0 - y_dist
        x_dist   = x-real(i_start)
        x_dist_m = 1.0e0 - x_dist

        ! Work out logical addresses of neighbors
        logical_address(1,1) = i_start
        logical_address(2,1) = j_start
        logical_address(3,1) = k_start
        logical_address(1,2) = i_start+1
        logical_address(2,2) = j_start
        logical_address(3,2) = k_start
        logical_address(1,3) = i_start
        logical_address(2,3) = j_start+1
        logical_address(3,3) = k_start
        logical_address(1,4) = i_start+1
        logical_address(2,4) = j_start+1
        logical_address(3,4) = k_start
        logical_address(1,5) = i_start
        logical_address(2,5) = j_start
        logical_address(3,5) = k_start+1
        logical_address(1,6) = i_start+1
        logical_address(2,6) = j_start
        logical_address(3,6) = k_start+1
        logical_address(1,7) = i_start
        logical_address(2,7) = j_start+1
        logical_address(3,7) = k_start+1
        logical_address(1,8) = i_start+1
        logical_address(2,8) = j_start+1
        logical_address(3,8) = k_start+1

        ! Work out weights of the neighbors
        weight(1) = x_dist_m * y_dist_m * z_dist_m
        weight(2) = x_dist   * y_dist_m * z_dist_m
        weight(3) = x_dist_m * y_dist   * z_dist_m
        weight(4) = x_dist   * y_dist   * z_dist_m
        weight(5) = x_dist_m * y_dist_m * z_dist
        weight(6) = x_dist   * y_dist_m * z_dist
        weight(7) = x_dist_m * y_dist   * z_dist
        weight(8) = x_dist   * y_dist   * z_dist

        ! Do the interpolation
        do i=1,8
            call self%PhysicalAddressGivenLogicalAddressInFourierSpace(logical_address(:,i),physical_address)
            if (logical_address(1,i) .ge. 0) then
                self%complex_values(physical_address(1),physical_address(2),physical_address(3)) = &
                self%complex_values(physical_address(1),physical_address(2),physical_address(3)) + weight(i)*value_to_add
            else
                self%complex_values(physical_address(1),physical_address(2),physical_address(3)) = &
                self%complex_values(physical_address(1),physical_address(2),physical_address(3)) + conjg(weight(i)*value_to_add)
            endif
            ! If the Hermitian mates are described explicitely, we need to insert this Fourier component a second time
            if (self%FourierComponentHasExplicitHermitianMate(physical_address(1),physical_address(2),physical_address(3))) then
                logical_address(2,i) = -1.0* logical_address(2,i)
                logical_address(3,i) = -1.0* logical_address(3,i)
                call self%PhysicalAddressGivenLogicalAddressInFourierSpace(logical_address(:,i),hermitian_physical_address)
                self%complex_values(hermitian_physical_address(1),hermitian_physical_address(2),hermitian_physical_address(3)) = &
                self%complex_values(hermitian_physical_address(1),hermitian_physical_address(2),hermitian_physical_address(3)) &
                + conjg(weight(i)*value_to_add)
            endif
        enddo

    end subroutine PutComplexValueByLinearInterpolationNoBoundsCheckVolume

    !>  \brief  Get a complex value from logical coordinates (i.e. origin at 0.0,0.0,0.0), using linear interpolation
    subroutine GetComplexValueByNearestNeighborNoBoundsCheckVolume(self,interpolated_value,x,y,z)
        ! Arguments
        class(Image),       intent(in)  ::  self
        complex,            intent(out) ::  interpolated_value
        real,               intent(in)  ::  x                       !<  Logical coordinate, first dimension
        real,               intent(in)  ::  y                       !<  Logical coordinate, second dimension
        real,               intent(in)  ::  z                       !<  Logical coordinate, third dimension
        ! Private variables
        integer ::  logical_address(3), physical_address(3) ! The neighbours are numbered 1 to 8
        ! Start work

        logical_address(1) = anint(x)
        logical_address(2) = anint(y)
        logical_address(3) = anint(z)

        call self%PhysicalAddressGivenLogicalAddressInFourierSpace(logical_address,physical_address)


        interpolated_value = self%complex_values(physical_address(1),physical_address(2),physical_address(3))

    end subroutine GetComplexValueByNearestNeighborNoBoundsCheckVolume


    !>  \brief  Get a complex value from logical coordinates (i.e. origin at 0.0,0.0,0.0), using linear interpolation
    subroutine GetComplexValueByLinearInterpolationNoBoundsCheckVolume(self,interpolated_value,x,y,z)
        ! Arguments
        class(Image),       intent(in)  ::  self
        complex,            intent(out) ::  interpolated_value
        real,               intent(in)  ::  x                       !<  Logical coordinate, first dimension
        real,               intent(in)  ::  y                       !<  Logical coordinate, second dimension
        real,               intent(in)  ::  z                       !<  Logical coordinate, third dimension
        ! Private variables
        integer ::  i
        integer ::  i_start,j_start,k_start
        real    ::  x_dist,y_dist,z_dist
        real    ::  x_dist_m,y_dist_m,z_dist_m
        integer ::  logical_address(3,8), physical_address(3) ! The neighbours are numbered 1 to 8
        real    ::  weight(8)
        ! Start work
        i_start  = int(x)
        j_start  = int(y)
        k_start  = int(z)
        z_dist   = z-real(k_start)
        z_dist_m = 1.0e0 - z_dist
        y_dist   = y-real(j_start)
        y_dist_m = 1.0e0 - y_dist
        x_dist   = x-real(i_start)
        x_dist_m = 1.0e0 - x_dist

        ! Work out logical addresses of neighbors
        logical_address(1,1) = i_start
        logical_address(2,1) = j_start
        logical_address(3,1) = k_start
        logical_address(1,2) = i_start+1
        logical_address(2,2) = j_start
        logical_address(3,2) = k_start
        logical_address(1,3) = i_start
        logical_address(2,3) = j_start+1
        logical_address(3,3) = k_start
        logical_address(1,4) = i_start+1
        logical_address(2,4) = j_start+1
        logical_address(3,4) = k_start
        logical_address(1,5) = i_start
        logical_address(2,5) = j_start
        logical_address(3,5) = k_start+1
        logical_address(1,6) = i_start+1
        logical_address(2,6) = j_start
        logical_address(3,6) = k_start+1
        logical_address(1,7) = i_start
        logical_address(2,7) = j_start+1
        logical_address(3,7) = k_start+1
        logical_address(1,8) = i_start+1
        logical_address(2,8) = j_start+1
        logical_address(3,8) = k_start+1

        ! Work out weights of the neighbors
        weight(1) = x_dist_m * y_dist_m * z_dist_m
        weight(2) = x_dist   * y_dist_m * z_dist_m
        weight(3) = x_dist_m * y_dist   * z_dist_m
        weight(4) = x_dist   * y_dist   * z_dist_m
        weight(5) = x_dist_m * y_dist_m * z_dist
        weight(6) = x_dist   * y_dist_m * z_dist
        weight(7) = x_dist_m * y_dist   * z_dist
        weight(8) = x_dist   * y_dist   * z_dist

        ! Do the interpolation
        interpolated_value = 0.0e0
        do i=1,8
            call self%PhysicalAddressGivenLogicalAddressInFourierSpace(logical_address(:,i),physical_address)
            if (logical_address(1,i) .ge. 0) then
                interpolated_value = interpolated_value + weight(i) &
                                    *       self%complex_values(physical_address(1),physical_address(2),physical_address(3))
            else
                interpolated_value = interpolated_value + weight(i) &
                                    * conjg(self%complex_values(physical_address(1),physical_address(2),physical_address(3)))
            endif
        enddo

    end subroutine GetComplexValueByLinearInterpolationNoBoundsCheckVolume

        !>  \brief  Get a complex value from logical coordinates (i.e. origin at 0.0,0.0,0.0), using linear interpolation
    pure subroutine GetComplexValueByLinearInterpolationNoBoundsCheckVolumeVector(self,interpolated_values, &
                                                                                number_of_values_to_interpolate, &
                                                    row_to_interpolate_into,slice_to_interpolate_into,logical_coordinates)
        ! Arguments
        class(Image),       intent(in)  ::  self
        complex,            intent(inout)   ::  interpolated_values(:,:,:)
        integer,            intent(in)  ::  number_of_values_to_interpolate
        integer,            intent(in)  ::  row_to_interpolate_into
        integer,            intent(in)  ::  slice_to_interpolate_into
        real,               intent(in)  ::  logical_coordinates(3,number_of_values_to_interpolate) !< Logical coordinates at which to interpolate
        ! Private variables
        integer ::  i
        integer ::  i_start,j_start,k_start
        real    ::  x_dist,y_dist,z_dist
        real    ::  x_dist_m,y_dist_m,z_dist_m
        integer ::  logical_address(3,8,number_of_values_to_interpolate)
        integer ::  physical_address(3,8,number_of_values_to_interpolate) ! The neighbours are numbered 1 to 8
        real    ::  weight(8,number_of_values_to_interpolate)
        integer ::  counter
        ! Start work


        ! Loop over the values to interpolate
        do counter = 1, number_of_values_to_interpolate
            i_start  = int(logical_coordinates(1,counter))
            j_start  = int(logical_coordinates(2,counter))
            k_start  = int(logical_coordinates(3,counter))
            z_dist   = logical_coordinates(3,counter)-real(k_start)
            z_dist_m = 1.0e0 - z_dist
            y_dist   = logical_coordinates(2,counter)-real(j_start)
            y_dist_m = 1.0e0 - y_dist
            x_dist   = logical_coordinates(1,counter)-real(i_start)
            x_dist_m = 1.0e0 - x_dist

            ! Work out logical addresses of neighbors
            logical_address(1,1,counter) = i_start
            logical_address(2,1,counter) = j_start
            logical_address(3,1,counter) = k_start
            logical_address(1,2,counter) = i_start+1
            logical_address(2,2,counter) = j_start
            logical_address(3,2,counter) = k_start
            logical_address(1,3,counter) = i_start
            logical_address(2,3,counter) = j_start+1
            logical_address(3,3,counter) = k_start
            logical_address(1,4,counter) = i_start+1
            logical_address(2,4,counter) = j_start+1
            logical_address(3,4,counter) = k_start
            logical_address(1,5,counter) = i_start
            logical_address(2,5,counter) = j_start
            logical_address(3,5,counter) = k_start+1
            logical_address(1,6,counter) = i_start+1
            logical_address(2,6,counter) = j_start
            logical_address(3,6,counter) = k_start+1
            logical_address(1,7,counter) = i_start
            logical_address(2,7,counter) = j_start+1
            logical_address(3,7,counter) = k_start+1
            logical_address(1,8,counter) = i_start+1
            logical_address(2,8,counter) = j_start+1
            logical_address(3,8,counter) = k_start+1

            ! Work out weights of the neighbors
            weight(1,counter) = x_dist_m * y_dist_m * z_dist_m
            weight(2,counter) = x_dist   * y_dist_m * z_dist_m
            weight(3,counter) = x_dist_m * y_dist   * z_dist_m
            weight(4,counter) = x_dist   * y_dist   * z_dist_m
            weight(5,counter) = x_dist_m * y_dist_m * z_dist
            weight(6,counter) = x_dist   * y_dist_m * z_dist
            weight(7,counter) = x_dist_m * y_dist   * z_dist
            weight(8,counter) = x_dist   * y_dist   * z_dist

            ! Do the interpolation
            interpolated_values(counter,row_to_interpolate_into,slice_to_interpolate_into) = 0.0e0
            do i=1,8
                call self%PhysicalAddressGivenLogicalAddressInFourierSpace(logical_address(:,i,counter),&
                                                                           physical_address(:,i,counter))
                if (logical_address(1,i,counter) .ge. 0) then
                    interpolated_values(counter,row_to_interpolate_into,slice_to_interpolate_into) = &
                    interpolated_values(counter,row_to_interpolate_into,slice_to_interpolate_into) + weight(i,counter) &
                                        *       self%complex_values(physical_address(1,i,counter), &
                                                                    physical_address(2,i,counter), &
                                                                    physical_address(3,i,counter))
                else
                    interpolated_values(counter,row_to_interpolate_into,slice_to_interpolate_into) = &
                    interpolated_values(counter,row_to_interpolate_into,slice_to_interpolate_into) + weight(i,counter) &
                                        * conjg(self%complex_values(physical_address(1,i,counter), &
                                                                    physical_address(2,i,counter), &
                                                                    physical_address(3,i,counter)))
                endif
            enddo
        enddo

    end subroutine GetComplexValueByLinearInterpolationNoBoundsCheckVolumeVector


    !> \brief   Apply a forward FT to the Image object. The FT is scaled.
    !!          The DC component is at (self%DIM(1)/2+1,self%DIM(2)/2+1,self%DIM(3)/2+1) (assuming even dimensions) or at (1,1,1) by default.
    !!
    !!
    !! For details on FFTW, see http://www.fftw.org/
    !! A helpful page for understanding the output format: http://www.dsprelated.com/showmessage/102675/1.php
    !! A helpful page to learn about vectorization and FFTW benchmarking: http://www.dsprelated.com/showmessage/76779/1.php
    !!  \todo   Check this: http://objectmix.com/fortran/371439-ifort-openmp-fftw-problem.html
#ifdef SKIP_RUNTIME_CHECKS
    pure &
#endif
    subroutine ForwardFFT(self,scale)
        use fftw33
        ! Arguments
        class(image),           intent(inout)   ::  self        !<  Image object to FT
        logical,    optional,   intent(in)      ::  scale       !<  Whether to scale the values so that a BwdFT of the output would yield the original image back, rather than a scaled version.
        ! Private variables
        logical                 ::  sscale
        ! Start work
#ifndef SKIP_RUNTIME_CHECKS
        ! Check the image is in memory
        if (.not. self%is_in_memory) then
            call this_program%TerminateWithFatalError('Image::ForwardFFT','Image not allocated')
        endif

        ! check that the image is currently in real space
        if (.not. self%is_in_real_space) then
            call this_program%TerminateWithFatalError('Image::ForwardFFT','Forward FFT called when already in Fourier space')
        endif
#endif

        ! Scale?
        sscale = .true.
        if (present(scale)) sscale = scale

        ! Do the actual FT
        call fftwf_execute_dft_r2c(self%plan_fwd,self%real_values,self%complex_values)

        ! Scale the Fourier transform
        if (sscale) then
            !self%complex_values = self%complex_values / sqrt(real(self%dim(1) * self%dim(2)* self%dim(3)))
            self%complex_values = self%complex_values / real(product(self%logical_dimensions))
        endif

        ! The image is now in Fourier space
        self%is_in_real_space = .false.

    end subroutine ForwardFFT

    !>  \brief Apply a backward FT to the Image object.
    !!
#ifdef SKIP_RUNTIME_CHECKS
    pure &
#endif
    subroutine BackwardFFT(self,scale)
        use fftw33
        ! Arguments
        class(image),           intent(inout)   ::  self
        logical,    optional,   intent(in)      ::  scale  !<  Whether to scale the values so that a BwdFT of the output would yield the original image back, rather than a scaled version.
        ! Private variables
        logical :: sscale
        ! Start work
#ifndef SKIP_RUNTIME_CHECKS
        ! Check the image is in memory
        if (.not. self%is_in_memory) then
            call this_program%TerminateWithFatalError('Image::BackwardFFT','Image not allocated')
        endif

        ! Check that the image is currently complex
        if (self%is_in_real_space) then
            call this_program%TerminateWithFatalError('Image::BackwardFFT','Backward FFT called when already in real space')
        endif
#endif

        ! Scale?
        sscale = .true.
        if (present(scale)) sscale = scale


        ! Do the actual FT
        call fftwf_execute_dft_c2r(self%plan_bwd,self%complex_values,self%real_values)


        ! Scale the output
        if (sscale) then
            !self%complex_values = self%complex_values / sqrt(real(self%dim(1) * self%dim(2)* self%dim(3)))
        endif

        ! Set the image type
        self%is_in_real_space = .true.

    end subroutine BackwardFFT

    !>  \brief  Set all pixels to given complex value
    subroutine AssignComplexToImage(self,complex_value)
        class(image),   intent(inout)   ::  self
        complex,        intent(in)      ::  complex_value
        if (self%is_in_real_space) then
            call this_program%TerminateWithFatalError('Image::AssignComplexToImage','Image is in real space')
        endif
        self%complex_values = complex_value
    end subroutine AssignComplexToImage

    !>  \brief  Set all pixels to given real value
    subroutine AssignRealToImage(self,real_value)
        class(image),   intent(inout)   ::  self
        real,           intent(in)      ::  real_value
        if (.not. self%is_in_real_space) then
            call this_program%TerminateWithFatalError('Image::AssignRealToImage','Image is in Fourier space')
        endif
        self%real_values = real_value
    end subroutine AssignRealToImage

    !>  \brief  Assignment (=) operation. self_rhs is copied to self_lhs
    subroutine AssignImageToImage(self_lhs,self_rhs)
        class(image),   intent(inout)   ::  self_lhs
        type(image),    intent(in)      ::  self_rhs
        integer :: i,j,k
        call self_lhs%Allocate(mould=self_rhs)
        if (.not. associated(self_rhs%real_values) .or. .not. associated(self_lhs%real_values)) then
            call this_program%TerminateWithFatalError('AssignImageToImage','real_values array is not associated')
        endif
        ! With Intel compiler (version 14), the explicit looping below is significantly
        ! faster than the more modern array assignment syntax.
        do k=1,size(self_lhs%real_values,3)
            do j=1,size(self_lhs%real_values,2)
                do i=1,size(self_lhs%real_values,1)
                    self_lhs%real_values(i,j,k) = self_rhs%real_values(i,j,k)
                enddo
            enddo
        enddo
    end subroutine AssignImageToImage

    !
    ! Unit tests
    !
    subroutine image_base_unit_tests()
        implicit none
        !
        call image_base_unit_test_1()
        call image_base_unit_test_2()
        call image_base_unit_test_3()
        call image_base_unit_test_4()
    end subroutine image_base_unit_tests

    !>  \brief  check image allocation and finalisation
    subroutine image_base_unit_test_1()
        implicit none
        ! private variables
        type(image) ::  self,self_copy
        ! start work
        write(*,'(a)') '**info(image_base_unit_test_1): entering'
        self%logical_dimensions = [32,32,32]
        call self%Allocate()
        self_copy = self
        write(*,'(a)') '**info(image_base_unit_test_1): exiting'
    end subroutine image_base_unit_test_1

    !>  \brief  check image allocate, deallocation and assignements
    subroutine image_base_unit_test_2()
        implicit none
        ! private variables
        type(image) ::  self, self_copy
        !
        write(*,'(a)') '**info(image_base_unit_test_2): checking allocation of even-dimension 3d image'
        self%logical_dimensions = [32,32,32]
        call self%Allocate()
        call self%Deallocate()
        !
        write(*,'(a)') '**info(image_base_unit_test_2): checking allocation of odd-dimension 3d image'
        self%logical_dimensions = [55,55,55]
        call self%Allocate()
        call self%Deallocate()
        !
        write(*,'(a)') '**info(image_base_unit_test_2): checking allocation of arbirtrary dimension 3d image'
        self%logical_dimensions = [32,55,120]
        call self%Allocate()
        call self%Deallocate()
        !
        write(*,'(a)') '**info(image_base_unit_test_2): checking image assignments (1)'
        call self%Allocate()
        self_copy = self
        call self%Deallocate()
        call random_number(self_copy%real_values) ! should still be able to use self_copy
        call self_copy%Deallocate()
        !
        write(*,'(a)') '**info(image_base_unit_test_2): checking image assignments (2)'
        call self%Allocate()
        self_copy = self
        call self_copy%Deallocate()
        call random_number(self%real_values) ! should still be able to use self
        call self%Deallocate()
        !
        write(*,'(a)') '**info(image_base_unit_test_2): checking image assignments (3)'
        self%logical_dimensions = [32,32,1]
        call self%Allocate()
        self_copy = self
        call self_copy%Deallocate()
        call random_number(self%real_values) ! should still be able to use self
        call self%Deallocate()
    end subroutine image_base_unit_test_2

    !>  \brief  check that an image has fsc of 1.0 with itself
    subroutine image_base_unit_test_3()
        implicit none
        ! private variables
        type(image) ::  self, self_copy
        !real, allocatable :: fsc_values(:)
        ! start work
        write(*,'(a)') '**info(image_base_unit_test_3): entering'
        self%logical_dimensions = [13,13,13]
        call self%Allocate()
        call random_number(self%real_values)
        self_copy = self
        !
        !call compute_fsc(self,self_copy,fsc_values)
        !allocate(fsc_values(7))
        !fsc_values = 1.0
        !
        !if (any(fsc_values .lt. 0.999)) then
        !    write(*,'(a)') '**error(fsc_unit_test_1): not all fsc values are 1.0'
        !    write(*,*) fsc_values
        !    call terminate('failed unit test')
        !else
        !    write(*,'(a)') '**info(fsc_unit_test_1): passed unit test'
        !endif
        write(*,'(a)') '**info(image_base_unit_test_3): exiting'
    end subroutine image_base_unit_test_3

    !>  \brief  test functionalistyt of the image_fft module
    subroutine image_fft_unit_tests()
        implicit none
        call image_fft_unit_test_1()
        call image_fft_unit_test_2()
    end subroutine image_fft_unit_tests

    !>  \brief  check that we can go back and forth between real and complex image and that values are kept
    subroutine image_fft_unit_test_1()
        implicit none
        ! private variables
        type(image) ::  self, self_bak
        real, allocatable :: diff(:,:,:)
        ! start work
        write(*,'(a)') '**info(image_fft_unit_test_1): entering'
        self%logical_dimensions = [ 32,32,1 ]
        call self%Allocate()
        call self_bak%Allocate(mould=self)
        ! use fortran's intrinsic random number generator
        call random_number(self%real_values)
        ! back up initial state
        self_bak = self
        call self%WriteToDisk('image_fft_unit_1_in.mrc')
        ! go to fourier space and back
        call self%ForwardFFT()
        call self%BackwardFFT()
        !
        call self%WriteToDisk('image_fft_unit_1_out.mrc')
        ! check we still have the same thing on return
        allocate(diff(self%logical_dimensions(1),self%logical_dimensions(2),self%logical_dimensions(3)))
        diff =  abs(self    %real_values(1:self%logical_dimensions(1),1:self%logical_dimensions(2),1:self%logical_dimensions(3)) &
                -   self_bak%real_values(1:self%logical_dimensions(1),1:self%logical_dimensions(2),1:self%logical_dimensions(3)) )
        if (any(diff .gt. 0.001)) then
            write(*,'(a)') '**error(image_fft_unit_test_1): test failed'
            call this_program%TerminateWithFatalError('Image::fft_unit_test_1','Test failed')
        else
            write(*,'(a)') '**info(image_fft_unit_test_1): test passed'
        endif
        !
        write(*,'(a)') '**info(image_fft_unit_test_1): exiting'
    end subroutine image_fft_unit_test_1

    !> \brief   check we can do an ft on an image that was allocated by assignment
    subroutine image_fft_unit_test_2()
        implicit none
        ! private variables
        type(image) ::  self,self_copy
        ! start work
        write(*,'(a)') '**info(image_fft_unit_test_2): entering'
        self%logical_dimensions = [32,32,1]
        call self%Allocate()
        self_copy = self
        call self%Deallocate()
        !
        write(*,'(a)') '**info(image_fft_unit_test_2): checking ft works on image allocated via assignment after source deallocated'
        call self_copy%ForwardFFT()
        call self_copy%BackwardFFT()
        !
        write(*,'(a)') '**info(image_fft_unit_test_2): exiting'
    end subroutine image_fft_unit_test_2

    !>  \brief  Checking addressing scheme, looping convience functions, definitions of Hermitian mates etc.
    subroutine image_base_unit_test_4()
        type(Image) ::  self
        !
        write(*,'(a)') '**info(image_base_unit_test_4): entering'
        ! Even-dimensions
        call self%Allocate(dims=[4,4,1],in_real_space=.false.)
        call CheckHermitianMateAddresses(self)
        call self%Allocate(dims=[4,4,4],in_real_space=.false.)
        call CheckHermitianMateAddresses(self)
        ! Odd dimension
        call self%Allocate(dims=[5,5,1],in_real_space=.false.)
        call CheckHermitianMateAddresses(self)
        call self%Allocate(dims=[5,5,5],in_real_space=.false.)
        call CheckHermitianMateAddresses(self)
        ! Mixed dimensions
        call self%Allocate(dims=[4,5,1],in_real_space=.false.)
        call CheckHermitianMateAddresses(self)
        call self%Allocate(dims=[4,5,6],in_real_space=.false.)
        call CheckHermitianMateAddresses(self)


        write(*,'(a)') '**info(image_base_unit_test_4): exiting'
    end subroutine image_base_unit_test_4

    subroutine CheckHermitianMateAddresses(self)
        type(Image),    intent(in)   ::  self
        integer ::  i,j,k
        integer ::  address(3)
        integer ::  hermitian_address(3)
        !
        do k=1,self%physical_upper_bound_complex(3)
            do j=1,self%physical_upper_bound_complex(2)
                do i=1,self%physical_upper_bound_complex(1)
                    if (self%FourierComponentHasExplicitHermitianMate(i,j,k)) then
                        address = [i,j,k]
                        call self%PhysicalAddressOfHermitianMate(address,hermitian_address)
                        ! Check the hermitian address is valid
                        if (any(hermitian_address .lt. 1) .or. any(hermitian_address .gt. self%physical_upper_bound_complex)) then
                            print *, self%logical_dimensions
                            print *, address
                            print *, hermitian_address
                            call this_program%TerminateWithFatalError('Image::CheckHermitianMatesAddresses', &
                                                                      'Failed: invalid addres')
                        endif
                        ! If we landed somewhere else in memory, check this is known to be redundant
                        if (any(address .ne. hermitian_address)) then
                            if (.not. self%FourierComponentHasExplicitHermitianMate(hermitian_address(1), &
                                                                                    hermitian_address(2), &
                                                                                    hermitian_address(3))) then
                                print *, self%logical_dimensions
                                print *, address
                                print *, hermitian_address
                                call this_program%TerminateWithFatalError('Image::CheckHermitianMatesAddresses', &
                                                'Failed: hermitian address not known to have explicit Hermitian mate')
                            endif
                        endif
                    endif
                    if (self%FourierComponentIsRedundant(i,j,k)) then
                        address = [i,j,k]
                        call self%PhysicalAddressOfHermitianMate(address,hermitian_address)
                        if (.not. self%FourierComponentHasExplicitHermitianMate(hermitian_address(1),&
                                                                                hermitian_address(2),&
                                                                                hermitian_address(3))) then
                            print *, self%logical_dimensions
                            print *, address
                            print *, hermitian_address
                            call this_program%TerminateWithFatalError('Image::CheckHermitianMatesAddresses', &
                                                            'Failed: hermitian address not known to be have redundant mate')
                        endif
                        if (any(hermitian_address .lt. 1) .or. any(hermitian_address .gt. self%physical_upper_bound_complex)) then
                            print *, self%logical_dimensions
                            print *, address
                            print *, hermitian_address
                            call this_program%TerminateWithFatalError('Image::CheckHermitianMatesAddresses', &
                                                                      'Failed: invalid address for mate of redundant')
                        endif
                    endif
                enddo
            enddo
        enddo
    end subroutine
end module Images
