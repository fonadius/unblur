!>  \brief  Class to deal with TIFF image files
module TiffImageFiles
    use libtiff
    use Globals
    use iso_c_binding, only : c_ptr, c_null_ptr
    implicit none
    private

    type, public :: TiffImageFile
        private
        type(c_ptr) ::  tiff = c_null_ptr
        contains
            procedure,  public  ::  Init
            procedure,  public  ::  Open
            procedure,  public  ::  GetDimensions
            procedure,  public  ::  ReadImageFromDisk
            procedure,  public  ::  Close
            procedure,  public  ::  PrintInfo
            procedure           ::  IsTiled
            procedure           ::  NumberOfStrips
            procedure           ::  GetStackSize
    end type TiffImageFile

    contains

    subroutine PrintInfo(self)
        class(TiffImageFile),                   intent(in)      ::  self
        call TIFFPrintInfo(self%tiff)
    end subroutine

    subroutine Close(self)
        use iso_c_binding, only : c_null_ptr, c_associated
        class(TiffImageFile),                   intent(inout)   ::  self
        if (c_associated(self%tiff)) then
            call TIFFClose(self%tiff)
            self%tiff = c_null_ptr
        endif
    end subroutine Close

    subroutine ReadImageFromDisk(self,image_number,output_array)
        use StringManipulations, only : IntegerToString
        use iso_c_binding, only : c_int32_t, c_f_pointer, c_int16_t
        ! Arguments
        class(TiffImageFile),                   intent(in)      ::  self
        integer,                                intent(in)      ::  image_number        !<  Number of image to be read in. First image is numbered 1.
        real,                                   intent(inout)   ::  output_array(:,:)
        ! Private variables
        integer         ::  strip_counter
        type(c_ptr)     ::  strip_buffer_cptr
        integer         ::  success
        integer         ::  number_of_bytes_placed_in_buffer
        integer         ::  tiff_dimensions(3)
        integer(kind=1),    pointer ::  temp_byte_array(:)
        integer         ::  rows_per_strip
        integer         ::  current_first_row, current_last_row
        integer         ::  current_first_byte, current_last_byte
        integer         ::  row_counter
        integer         ::  previous_number_of_bytes_placed_in_buffer
        ! Start work

        ! What are the dimensions of the image stack held by the the TIFF file
        tiff_dimensions = self%GetDimensions()

        ! Make sure we're on the correct TIFF dictionary
        if (image_number .gt. tiff_dimensions(3)) then
            call this_program%TerminateWithFatalError('TiffImageFiles::ReadImageFromDisk','Invalid image number')
        endif
        success = TIFFSetDirectory(self%tiff,int(image_number-1,kind=c_int16_t))
        if (success .ne. 1) then
            call this_program%TerminateWithFatalError('TiffImageFiles::ReadImageFromDisk',&
                                                      'Failed to set directory to : '//&
                                                      IntegerToString(image_number-1))
        endif

        ! We don't support tiled organization, only strips
        if (self%IsTiled()) then
            call this_program%TerminateWithFatalError('TiffImageFiles::ReadImageFromDisk', &
                                                      'Tile-based TIFF files not supported')
        endif

        ! We only support 8-bit grayscale TIFFs
        if (TIFFGetBitsPerSample(self%tiff) .ne. 8) then
            call this_program%TerminateWithFatalError('TiffImageFiles::ReadImageFromDisk',&
                                                      'Unsupported bit depth: '//&
                                                      IntegerToString(TIFFGetBitsPerSample(self%tiff)))
        endif
        if (TIFFGetSamplesPerPixel(self%tiff) .ne. 1) then
            call this_program%TerminateWithFatalError('TiffImageFiles::ReadImageFromDisk', &
                                                     'Unsupported number of samples per pixel: '//&
                                                     IntegerToString(TIFFGetSamplesPerPixel(self%tiff)))
        endif
        ! The defines below copied from tiff.h on 150317
#define	    SAMPLEFORMAT_UINT		1	/* !unsigned integer data */
#define	    SAMPLEFORMAT_INT		2	/* !signed integer data */
#define	    SAMPLEFORMAT_IEEEFP		3	/* !IEEE floating point data */
#define	    SAMPLEFORMAT_VOID		4	/* !untyped data */
#define	    SAMPLEFORMAT_COMPLEXINT	5	/* !complex signed int */
#define	    SAMPLEFORMAT_COMPLEXIEEEFP	6	/* !complex ieee floating */
        if (TIFFGetSampleFormat(self%tiff) .ne. SAMPLEFORMAT_UINT) then
            call this_program%TerminateWithFatalError('TiffImageFiles::ReadImageFromDisk',&
                                                      'Unsupported sample format: '//&
                                                      IntegerToString(TIFFGetSampleFormat(self%tiff)))
        endif

        ! Check the data array has correct dimensions
        if (size(output_array,1) .ne. tiff_dimensions(1) .or. size(output_array,2) .ne. tiff_dimensions(2)) then
            print *, size(output_array,1), size(output_array,2)
            print *, tiff_dimensions(1:2)
            call this_program%TerminateWithFatalError('TiffImageFiles::ReadImageFromDisk',&
                                                      'Data array has wrong dimensions')
        endif

        ! How many rows of pixels per strip?
        rows_per_strip = TIFFGetRowsPerStrip(self%tiff)


        ! Allocate memory to hold data read from strip
        strip_buffer_cptr = TIFFAllocateStripBuffer(self%tiff)

        !
        previous_number_of_bytes_placed_in_buffer = 0

        ! Loop over strips
        do strip_counter=1,self%NumberOfStrips()
            ! Read a strip from the file
            number_of_bytes_placed_in_buffer = TIFFReadEncodedStrip(self%tiff, &
                                                                    int(strip_counter-1,kind=c_int32_t), &
                                                                    strip_buffer_cptr, &
                                                                    int(-1,kind=c_int32_t))
            if (number_of_bytes_placed_in_buffer .le. 0) then
                call this_program%TerminateWithFatalError('TiffImageFiles::ReadImageFromDisk', &
                                                          'Failed to read from strip '//&
                                                          IntegerToString(strip_counter))
            endif
            ! Setup F pointer to the strip data
            if (number_of_bytes_placed_in_buffer .ne. previous_number_of_bytes_placed_in_buffer) then
                call c_f_pointer(strip_buffer_cptr,temp_byte_array,[number_of_bytes_placed_in_buffer])
            endif
            ! Sanity check
            if (number_of_bytes_placed_in_buffer .ne. rows_per_strip * size(output_array,1)) then
                call this_program%TerminateWithFatalError('TiffImageFiles::ReadImageFromDisk', &
                                                          'Unexpected number of bytes in buffer ')
            endif
            ! Which rows of the image did we just get
            current_first_row = rows_per_strip*(strip_counter-1) + 1
            current_last_row  = rows_per_strip*strip_counter
            ! Copy (and cast) the data to the output array
            do row_counter=current_first_row,current_last_row
                current_first_byte = (row_counter-current_first_row)*size(output_array,1)+1
                current_last_byte  = current_first_byte + size(output_array,1) - 1
                output_array(1:size(output_array,1),size(output_array,2)-row_counter+1) = &
                                                temp_byte_array(current_first_byte:current_last_byte)
            enddo
        enddo

        ! Deallocate strip buffer
        success = TIFFDeallocateStripBuffer(strip_buffer_cptr)

    end subroutine ReadImageFromDisk

    logical function IsTiled(self)
        class(TiffImageFile),   intent(in)      ::  self
        IsTiled = TIFFIsTiled(self%tiff) .ne. 0
    end function IsTiled

    integer function NumberOfStrips(self)
        class(TiffImageFile),   intent(in)      ::  self
        NumberOfStrips = TIFFNumberOfStrips(self%tiff)
    end function NumberOfStrips



    !>  \brief  Return the logical dimensions of the image data
    function GetDimensions(self)
        class(TiffImageFile),   intent(in)      ::  self
        integer                                 ::  GetDimensions(3)
        ! Start work
        GetDimensions(1) = TIFFGetWidth(self%tiff)
        GetDimensions(2) = TIFFGetLength(self%tiff)
        GetDimensions(3) = self%GetStackSize()
    end function GetDimensions

    !>  \brief  Work out how many images are in the tiff file
    integer function GetStackSize(self)
        class(TiffImageFile),   intent(in)  ::  self
        !
        integer ::  success
        ! Start work
        do
            if (TIFFLastDirectory(self%tiff) .ne. 0) exit
            success = TIFFReadDirectory(self%tiff)
            if (success .ne. 1) stop 'Error setting tiff directory, or already at last directory'
        enddo
        GetStackSize = TIFFCurrentDirectory(self%tiff) + 1
    end function GetStackSize


    !>  \brief  Initialise
    subroutine Init(self)
        class(TiffImageFile),   intent(inout)   ::  self
        !
    end subroutine Init

    !>  \brief  Open a file on disk
    subroutine Open(self,filename,delete_if_already_exists)
        use iso_c_binding
        use StringManipulations, only : GetCStringFromFString
        use libtiff, only : TIFFOpen
        ! Arguments
        class(TiffImageFile),               intent(inout)   ::  self
        character(len=*),                   intent(in)      ::  filename
        logical,                optional,   intent(in)      ::  delete_if_already_exists
        ! Private variables
        character(len=:),   allocatable ::  open_mode
        !type(c_ptr)::  my_tiff
        character(kind=c_char),   allocatable ::  filename_c(:), open_mode_c(:)
        ! Start work
        open_mode = 'rc' ! 'r+' for read-write, I think
        if (present(delete_if_already_exists)) then
            if (delete_if_already_exists) open_mode = 'w'
        endif
        filename_c =  GetCStringFromFString(filename)
        open_mode_c = GetCStringFromFString(open_mode)
        self%tiff = TIFFOpen(filename_c,open_mode_c)
    end subroutine Open

end module TiffImageFiles
