!>  \brief  Class to deal with image files on disks
!!
!!  The following file formats are (will be) supported:
!!  - Imagic: http://imagescience.de/formats/index.htm
!!  - Spider: http://www.wadsworth.org/spider_doc/spider/docs/image_doc.html
!!  - MRC: http://www2.mrc-lmb.cam.ac.uk/image2000.html
!!
!!
module ImageFiles ! When ready to take over from frealix, do not call this core any more
    use ImageHeaders
    use Globals
    use TiffImageFiles
    implicit none
    private
    public :: ImageFile

    type ImageFile
        private
        character(len=filename_max_len) ::  filename                        !<  Filename
        logical                         ::  initialised     =   .false.     !<  Set to true when the object has been initialised, false when it's closed.
        logical                         ::  was_written_to  =   .false.     !<  Whether data was written to the file since it was opened
        character(len=1)                ::  header_format                   !<  'I' (imagic),'M' (mrc),'S' (spider) or 'T' (tiff) file format.
        integer                         ::  lun_img         =   0           !<  Logical unit number for the image file. Set to 0 unless the file is currently open.
        integer                         ::  lun_hed         =   0           !<  Logical unit number for the header file (in the case of MRC or Spider file, this is the same as the image file, so that LUN_HED == LUN_IMG). Set to 0 unless currently open.
        class(ImageHeader), allocatable ::  header                          !<  Image header object
        type(TiffImageFile)             ::  tiff                            !<  TIFF image file object
        contains
            procedure           ::  ImageFileInit
            procedure, public   ::  Init => ImageFileInit
            procedure           ::  FormatFromFilename
            procedure           ::  Exists
            procedure           ::  Open
            procedure           ::  Close
            procedure, public   ::  ReadSlicesFromDisk
            procedure, public   ::  ReadSliceFromDisk
            procedure, public   ::  WriteSlicesToDisk
            procedure, public   ::  WriteSliceToDisk
            procedure, public   ::  PrintInfo
!            procedure, public   ::  ReadImage
!            procedure, public   ::  WriteImage
            procedure, public   ::  GetStackSize
            procedure, public   ::  GetDimensions
            procedure, public   ::  GetDimension
!            procedure, public   ::  GetDim
            procedure, public   ::  GetPixelSize
            procedure, public   ::  GetFilename
            final               ::  Destructor
    end type ImageFile

    contains

        !>  \brief  Print out basic information about the file
        subroutine PrintInfo(self)
            ! Arguments
            class(ImageFile),           intent(in)      ::  self
            ! Start work
            write(*,'(/2a)') 'Summary information for file ', trim(adjustl(self%filename))
            if (self%header_format .eq. 'T') then
                call self%tiff%PrintInfo()
            else
                call self%header%PrintInfo()
            endif
            write(*,'(a)') ' '
        end subroutine PrintInfo

        !>  \brief  Return the dimension of the image stack
        function GetDimensions(self)
            ! Arguments
            class(ImageFile),           intent(in)      ::  self
            ! Result
            integer                                     ::  GetDimensions(3)
            ! Start work
            if (self%header_format .eq. 'T') then
                GetDimensions = self%tiff%GetDimensions()
            else
                GetDimensions = self%header%GetDimensions()
            endif
        end function GetDimensions

        !>  \brief  Return one of the dimensio of the image stack
        function GetDimension(self,which_dimension)
            ! Arguments
            class(ImageFile),           intent(in)      ::  self
            integer,                    intent(in)      ::  which_dimension
            ! Result
            integer                                     ::  GetDimension
            ! Private variables
            integer ::  dimensions(3)
            ! Start work
            if (self%header_format .eq. 'T') then
                dimensions = self%tiff%GetDimensions()
                GetDimension = dimensions(which_dimension)
            else
                GetDimension = self%header%GetDimension(which_dimension)
            endif
        end function GetDimension


        !>  \brief Return the filename
        function GetFilename(self)
            ! Arguments
            class(ImageFile),           intent(in)      ::  self
            ! Result
            character(len=filename_max_len)             ::  GetFilename
            GetFilename = self%filename
        end function GetFilename

        !>  \brief  Return the pixel size of the image data (in Angstroms)
        real function GetPixelSize(self)
            ! Arguments
            class(ImageFile),           intent(in)      ::  self
            ! Start work
            GetPixelSize = self%header%GetPixelSize()
        end function GetPixelSize

        !>  \brief  Return the number of 2D images in the stack
        integer function GetStackSize(self)
            ! Arguments
            class(ImageFile),           intent(in)      ::  self
            ! Private variables
            integer ::  dimensions(3)
            ! Start work
            if (self%header_format .eq. 'T') then
                dimensions = self%tiff%GetDimensions()
                GetStackSize = dimensions(3)
            else
                GetStackSize = self%header%GetStackSize()
            endif
        end function GetStackSize

        !>  \brief  Read a slice of the image file from disk into memory
        subroutine WriteSliceToDisk(self,slice_number,real_array,logical_dimension_1)
            ! Arguments
            class(ImageFile),               intent(inout)   ::  self
            integer,                        intent(in)      ::  slice_number        !<  Number of the slice to read in (the first slice in the file is numbered 1)
            real,                           intent(in)      ::  real_array(:,:,:)   !<  Array of reals. Will be (re)allocated if needed
            integer,                        intent(in)      ::  logical_dimension_1 !<  Logical size of the array in the first dimension. This will be written to disk: real_array(1:logical_dimension_1,:,:).
            ! Start work
            call self%WriteSlicesToDisk(slice_number,slice_number,real_array,logical_dimension_1)
        end subroutine WriteSliceToDisk

        !>  \brief  Write a set of contiguous slices to an image file on disk. The array of reals should have +2 elements in the first dimension if logical_dimension_1 is even.
        subroutine WriteSlicesToDisk(self,first_slice,last_slice,real_array,logical_dimension_1)
            use StringManipulations, only : IntegerToString
#ifndef __GFORTRAN__
            use ieee_arithmetic, only : ieee_is_nan
#endif
            ! Arguments
            class(ImageFile),               intent(inout)   ::  self
            integer,                        intent(in)      ::  first_slice         !<  First slice to write in (the first slice in the file is numbered 1)
            integer,                        intent(in)      ::  last_slice          !<  Last slice to write in
            real,                           intent(in)      ::  real_array(:,:,:)   !<  Array of reals.
            integer,                        intent(in)      ::  logical_dimension_1 !<  Logical size of the array in the first dimension. This will be written to disk: real_array(1:logical_dimension_1,:,:).
            ! Private variables
            integer                     ::  io_status
            integer                     ::  dimensions(3)
            integer(kind=8)             ::  first_byte
            real(kind=4),   allocatable ::  temp_array(:,:,:)
            real                        ::  min_value, max_value
            integer                     ::  i,j,k
            logical,    parameter       ::  survey_pixel_values = .true.
            ! Start work

            ! Get the dimensions of the image file
            dimensions = self%header%GetDimensions()

            ! Check that the first and last slice numbers given make sense
            if (first_slice .gt. last_slice) then
                call this_program%TerminateWithFatalError('ImageFile_WriteSlicesToDisk','Last < first slice')
            endif

            ! Check that the array dimensions and the file dimensions are compatible.
            ! When writing to the first location, we redefine the file dimensions.
            if (first_slice .eq. 1) then
                dimensions(1) = logical_dimension_1 !size(real_array,1)-2
                dimensions(2) = size(real_array,2)
                dimensions(3) = size(real_array,3)
                call self%header%SetDimensions(dimensions)
            else
                if (logical_dimension_1  .ne. dimensions(1) .or. &
                    size(real_array,2)   .ne. dimensions(2) .or. &
                    size(real_array,3)   .gt. dimensions(3) ) then
                    call this_program%TerminateWithFatalError('ImageFile_WriteSlicesToDisk',    &
                                                              'Data array and image file have incompatible dimensions')
                endif
            endif

            ! Work out the position of the first byte we want to write to
            first_byte =    int(self%header%FirstDataByte(),kind=8) +   &
                            int((first_slice-1),kind=8) *               &
                            int(product(dimensions(1:2)),kind=8) *      &
                            int(self%header%BytesPerPixel(),kind=8)

            ! Writing to TIFF not yet supported
            if (self%header_format .eq. 'T') then
                call this_program%TerminateWithFatalError('ImageFile::WriteSlicesToDisk','Writing to TIFF files not yet supported')
            endif

            ! Write data to disk
            select case (self%header%BytesPerPixel())
                case(4)
                    temp_array = real_array(1:dimensions(1),:,:)
                    write(unit=self%lun_img,pos=first_byte,iostat=io_status) temp_array
                case default
                    call this_program%TerminateWithFatalError('ImageFile_WriteSlicesToDisk', &
                                                      'Unsupported bit-depth: '//IntegerToString(self%header%BytesPerPixel()))
            end select

            ! Check the write was successful
            if (io_status .ne. 0) then
                write(*,'(a,i0,2a)') '**ERROR(WriteSlicesFromDisk): I/O error ', io_status, ' when writing to: ', self%filename
                call this_program%TerminateWithFatalError('ImageFiles_WriteSlicesToDisk','I/O error')
            endif

            ! May need to update file dimensions
            dimensions(3) = max(dimensions(3),last_slice)
            call self%header%SetDimensions(dimensions)

            ! May need to update min and max in the header
            if (survey_pixel_values) then
                min_value =  huge(1.0e0)
                max_value = -huge(1.0e0)
                do k=1,size(real_array,3)
                    do j=1,size(real_array,2)
                        do i=1,logical_dimension_1
#ifndef __GFORTRAN__
                            if (.not. ieee_is_nan(real_array(i,j,k))) then
#else
                            if (.not. isnan(real_array(i,j,k))) then
#endif
                                if (real_array(i,j,k) .gt. max_value) max_value = real_array(i,j,k)
                                if (real_array(i,j,k) .lt. min_value) min_value = real_array(i,j,k)
                            endif
                        enddo
                    enddo
                enddo
                if (min_value .lt. self%header%GetMinimumPixelValue()) call self%header%SetMinimumPixelValue(min_value)
                if (max_value .gt. self%header%GetMaximumPixelValue()) call self%header%SetMaximumPixelValue(max_value)
            endif

            ! Remember that we wrote to the file
            self%was_written_to = .true.

        end subroutine WriteSlicesToDisk

        !>  \brief  Read a slice of the image file from disk into memory
        subroutine ReadSliceFromDisk(self,slice_number,real_array)
            ! Arguments
            class(ImageFile),           intent(in)      ::  self
            integer,                        intent(in)      ::  slice_number        !<  Number of the slice to read in (the first slice in the file is numbered 1)
            real,       allocatable,        intent(inout)   ::  real_array(:,:,:)   !<  Array of reals. Will be (re)allocated if needed
            ! Start work
            call self%ReadSlicesFromDisk(slice_number,slice_number,real_array)
        end subroutine ReadSliceFromDisk

        !>  \brief  Read a set of contiguous slices of the image file from disk into memory. The array of reals should have +2 elements in the first dimension.
        subroutine ReadSlicesFromDisk(self,first_slice,last_slice,real_array)
            use UsefulFunctions, only : IsOdd
            use iso_c_binding
            use ImageHeaders, only : data_are_bytes, data_are_integer, data_are_float
            ! Arguments
            class(ImageFile),               intent(in)      ::  self
            integer,                        intent(in)      ::  first_slice         !<  First slice to read in (the first slice in the file is numbered 1)
            integer,                        intent(in)      ::  last_slice          !<  Last slice to read in
            real,                           intent(inout)   ::  real_array(:,:,:)   !<  Array of reals. Will be (re)allocated if needed
            ! Private variables
            integer ::  io_status
            integer ::  current_slice
            integer ::  dimensions(3)
            integer(kind=8) ::  first_byte
            logical ::  array_is_ready
            real(kind=4), allocatable       ::  temp_32bit_float_array(:,:,:)
            integer(kind=1), allocatable    ::  temp_byte_array(:,:,:)
            integer(kind=2), allocatable    ::  temp_16bit_int_array(:,:,:)
            character(len=100)              ::  io_message
            ! Start work

            ! Get the dimensions of the image file
            dimensions = self%GetDimensions()

            ! Check that the first and last slice numbers given make sense
            if (first_slice .gt. last_slice) then
                call this_program%TerminateWithFatalError('ImageFile::ReadSlicesFromDisk','Last < first slice')
            endif
            if (last_slice .gt. dimensions(3)) then
                call this_program%TerminateWithFatalError('ImageFile::ReadSlicesFromDisk','Last slice > dim_z')
            endif

            ! Work out the dimensions of the array we're about to read in
            dimensions(3) = last_slice - first_slice + 1

            ! Check the array is properly allocated
            if (IsOdd(dimensions(1))) then
                array_is_ready =    size(real_array,1) .eq. dimensions(1) + 1
            else
                array_is_ready =    size(real_array,1) .eq. dimensions(1) + 2
            endif
            array_is_ready      =   array_is_ready                              &
                             .and. (size(real_array,2) .eq. dimensions(2))      &
                             .and. (size(real_array,3) .eq. dimensions(3))

            ! We need the memory to have been properly allocated
            if (.not. array_is_ready) then
                call this_program%TerminateWithFatalError('ImageFile::ReadSlicesFromDisk','Array is not properly allocated')
            endif

            ! If we're dealing with a TIFF file, we rely on the libtiff library to do the reading, decompressing etc for us
            if (self%header_format .eq. 'T') then
                allocate(temp_32bit_float_array(dimensions(1),dimensions(2),dimensions(3)))
                do current_slice=first_slice,last_slice
                    call self%tiff%ReadImageFromDisk(current_slice,temp_32bit_float_array(:,:,current_slice-first_slice+1))
                enddo
                real_array(1:dimensions(1),:,:) = temp_32bit_float_array
                deallocate(temp_32bit_float_array)
            else
                ! Not TIFF

                ! Work out the position of the first byte we want to read in
                first_byte =    int(self%header%FirstDataByte(), kind=8) +                              &
                                int((first_slice-1), kind=8)*int(product(dimensions(1:2)), kind=8) *    &
                                int(self%header%BytesPerPixel(), kind=8)

                ! Read data from disk
                ! Todo, make sure we first query from the header whether data is in real/float or integers and whether it is signed or unsigned
                select case (self%header%BytesPerPixel())
                    case(1) ! Byte data
                        allocate(temp_byte_array(dimensions(1),dimensions(2),dimensions(3)))
                        read(unit=self%lun_img,pos=first_byte,iostat=io_status,iomsg=io_message) temp_byte_array
                        ! Conversion from unsigned byte integer (which MRC appears to be) is tricky because Fortran doesn't do unsigned integer natively.
                        ! The following IAND trick is courtesy of Jim Dempsey at http://software.intel.com/en-us/forums/showthread.php?t=64400
                        ! Confusingly, the MRC format documentation implies that one should expect signed integers, which seems to be incorrect: http://www2.mrc-lmb.cam.ac.uk/image2000.html
                        ! IMOD documentation indicates that prior to IMOD 4.2.23, unsigned bytes were used and that one needs to inspect the imodStamp header to check
                        if (self%header%PixelDataAreSigned()) then
                            real_array(1:dimensions(1),:,:) = temp_byte_array(:,:,:)
                        else
                            real_array(1:dimensions(1),:,:) = real(iand(int(temp_byte_array(:,:,:),kind=4),int(255,kind=4)))
                        endif
                        deallocate(temp_byte_array)
                    case(2) ! 16-bit data
                        select case (self%header%GetPixelDataType())
                            case(data_are_integer)
                                allocate(temp_16bit_int_array(dimensions(1),dimensions(2),dimensions(3)))
                                read(unit=self%lun_img,pos=first_byte,iostat=io_status,iomsg=io_message) temp_16bit_int_array
                                if (self%header%PixelDataAreSigned()) then
                                    real_array(1:dimensions(1),:,:) = real(temp_16bit_int_array(:,:,:))
                                else
                                    real_array(1:dimensions(1),:,:) = real(iand(int(temp_16bit_int_array(:,:,:),kind=4), &
                                                                                int(huge(int(1,kind=2)),kind=4)))
                                endif
                                deallocate(temp_16bit_int_array)
                            case default
                                call this_program%TerminateWithFatalError('ImageFile::ReadSlicesFromDisk', &
                                                                            'Non-integer 16-bit data not supported')
                        end select
                    case(4) ! 32-bit data
                        allocate(temp_32bit_float_array(dimensions(1),dimensions(2),dimensions(3)))
                        read(unit=self%lun_img,pos=first_byte,iostat=io_status,iomsg=io_message) temp_32bit_float_array
                        real_array(1:dimensions(1),:,:) = temp_32bit_float_array

                        deallocate(temp_32bit_float_array)
                    case default
                        write(*,'(2a)') 'filename: ', self%filename
                        write(*,'(a,i0,a)') 'bit depth: ', self%header%BytesPerPixel(), ' bytes'
                        call this_program%TerminateWithFatalError('ImageFile::ReadSlicesFromDisk','Unsupported bit-depth')
                end select

                ! Make sure we set the non-used part of the array to 0.0
                if (.not. self%header%PixelDataAreComplex()) then
                    real_array(dimensions(1)+1:,:,:) = 0.0e0
                endif

                ! Check the read was successful
                if (io_status .ne. 0) then
                    write(*,'(a,i0,2a)') '**ERROR(ReadSlicesFromDisk): I/O error ', io_status, ' when reading from: ', self%filename
                    write(*,'(2a)') 'IO error message was: ', io_message
                    call this_program%TerminateWithFatalError('ImageFiles::ReadSlicesFromDisk','I/O error')
                endif

                ! Is this file from a machine with opposite endinaness?
                if (.not. self%header%HasLocalEndianess()) then
                    call this_program%TerminateWithFatalError('ImageFile::ReadSlicesFromDisk', &
                                'Files created by machines with the opposite endianess are not supported')
                endif
            endif ! end of test for tiff

        end subroutine ReadSlicesFromDisk

        !>  \brief  Initialise an ImageFile object
        subroutine ImageFileInit(self,filename,mould,dim_x,dim_y,dim_z,pixel_size,delete_if_already_exists)
            ! Arguments
            class(ImageFile),               intent(inout)   ::  self        !<  Imagefile object to be initialised
            character(len=*),               intent(in)      ::  filename    !<  Filename
            type(ImageFile),    optional,   intent(in)      ::  mould       !<  An imagefile object, which the new object will be initialised to resemble.
            integer,            optional,   intent(in)      ::  dim_x       !<  Dimension in X
            integer,            optional,   intent(in)      ::  dim_y       !<  Dimension in Y
            integer,            optional,   intent(in)      ::  dim_z       !<  Dimension in Z
            real,               optional,   intent(in)      ::  pixel_size  !<  Pixel size of image data (in Angstroms)
            logical,            optional,   intent(in)      ::  delete_if_already_exists    !<  If the file already exists on disk, replace it
            ! Private variables
            !logical ::  file_already_exists
            ! Start work
            if (self%initialised) then
                call this_program%TerminateWithFatalError('ImageFileCore_Init', &
                                                          'Attempt to initialise an ImageFile which is already initialised')
            endif

            ! Remove leading blanks in filename
            self%filename = adjustl(filename)

            !! Does the file already exist on disk?
            !file_already_exists = self%Exists()


            !
            ! Work out the file format the user wants to use.
            ! First priority is given to any extension present in the filename.
            ! Next, we check if a mould was given
            ! Otherwise, we revert to the default from the runtime parameters.
            !
            !> \todo check the header format conforms, or whether the header format tells us what kind of file
            self%header_format = default_file_format ! this comes from Globals

            if (self%FormatFromFilename() .ne. 'n') then
                self%header_format = self%FormatFromFilename()
            else
                if (present(mould)) then
                    self%header_format = mould%header_format
                endif
            endif

            ! Allocate & init the header object
            select case(self%header_format)
                case ('M')
                    allocate(MrcImageHeader :: self%header)
                    call self%header%Init()
                case ('T')
                    ! No need to allocate a header for TIFF files, instead we initialise the tif object
                    call self%tiff%Init()
                case default
                    call this_program%TerminateWithFatalError('ImageFile_Init','File format not supported yet')
            end select


            ! Open the file
            call self%Open(delete_if_already_exists)

            if (self%header_format .ne. 'T') then
                ! If the file exists, read the header, if not, setup the header
                ! with default values and the given file dimensions
                if (self%Exists()) then
                    call self%header%ReadFromDisk(self%lun_hed)
                else
                    call self%header%ResetToDefaults()
                    if (present(mould)) then
                        self%header = mould%header
                    else if (present(dim_x) .and. present(dim_y) .and. present(dim_z)) then
                        call self%header%SetDimensions([dim_x,dim_y,dim_z])
                    endif
                endif

                ! Pixel size
                if (present(pixel_size)) call self%header%SetPixelSize(pixel_size)
            endif

            ! The file is now initialised
            self%initialised = .true.

        end subroutine ImageFileInit

        !>  \brief  Clean up when we're done
        subroutine Destructor(self)
            ! Arguments
            type(ImageFile),    intent(inout)   ::  self
            ! Start work
            call self%Close()
            self%initialised = .false.
        end subroutine Destructor


        !>  \brief Open the file(s) for the ImageFile
        subroutine Open(self,delete_if_already_exists)
            use StringManipulations, only : FilenameReplaceExtension
            ! Arguments
            class(ImageFile),               intent(inout)   ::  self
            logical,            optional,   intent(in)      ::  delete_if_already_exists
            ! Private variables
            character(len=9)    ::  rw_string
            character(len=7)    ::  status_string
            ! Start work


            if (self%header_format .eq. 'T') then
                !
                call self%tiff%Open(self%filename,delete_if_already_exists)
            else

                ! We need to prepare a string for the Open statement
                rw_string = 'READWRITE'

                ! What is the status of the file?
                status_string = 'UNKNOWN'
                if (present(delete_if_already_exists)) then
                    if (delete_if_already_exists) status_string  = 'REPLACE'
                endif

                ! Get an IO unit number for the header file
                self%lun_hed = this_program%GetAvailableUnit()

                !
                if (self%header_format .eq. 'I') call FilenameReplaceExtension(self%filename,'hed')
                open(unit=self%lun_hed,access='STREAM',file=self%filename,action=rw_string,status=status_string)
                if (self%header_format .eq. 'I') then
                    call FilenameReplaceExtension(self%filename,'img')
                    self%lun_img = this_program%GetAvailableUnit()
                    open(unit=self%lun_img,access='STREAM',file=self%filename,action=rw_string,status=status_string)
                else
                   self%lun_img = self%lun_hed
                endif
            endif

            !
            self%was_written_to = .false.
        end subroutine Open

        !>  \brief  Close the file(s) and "de-initialise" the ImageFile object
        subroutine Close(self)
            use UsefulFunctions, only : UnitIsOpen
            ! Arguments
            class(ImageFile),   intent(inout)   ::  self
            ! Start work
            if (self%initialised) then
                if (UnitIsOpen(self%lun_hed)) then
                    if (self%was_written_to) call self%header%WriteToDisk(self%lun_hed)
                    close(self%lun_hed)
                    call this_program%ReleaseUnit(self%lun_hed)
                endif
                if (self%header_format .eq. 'I') then
                    if (UnitIsOpen(self%lun_img)) then
                        close(self%lun_img)
                        call this_program%ReleaseUnit(self%lun_img)
                    endif
                endif
            endif
            if (allocated(self%header)) call self%header%Destroy()
            if (allocated(self%header)) deallocate(self%header)
            if (self%header_format .eq. 'T') call self%tiff%Close()
            self%initialised = .false.
        end subroutine Close

        !>  \brief  Check whether the file exists on disk
        logical function Exists(self)
            use UsefulFunctions, only : FileExists, FileSizeFromFilename
            use StringManipulations, only : FilenameReplaceExtension
            class(ImageFile),   intent(inout)  ::  self
            if (self%header_format .eq. 'I') then
                call FilenameReplaceExtension(self%filename,'img')
                Exists = FileExists(self%filename) .and. FileSizeFromFilename(self%filename) .gt. 0
                call FilenameReplaceExtension(self%filename,'hed')
                Exists = FileExists(self%filename) .and. FileSizeFromFilename(self%filename) .gt. 0.and. Exists
            else
                Exists = FileExists(self%filename)
				if (Exists) Exists = Exists .and. FileSizeFromFilename(self%filename) .gt. 0
            endif
        end function Exists

        !>  \brief Return a one letter code for the file format designated by the extension in the filename
        !!
        !!  If .mrc: M
        !!  If .spi: S
        !!  If .img: I
        !!  If .hed: I
        !!  If .tif: T
        !!  Else: N
        pure function FormatFromFilename(self)
            use StringManipulations, only : ExtensionFromFilename
            ! Arguments
            class(ImageFile),   intent(in)  ::  self
            ! Return value
            character(len=1)    ::  FormatFromFilename
            ! Private variable
            character(len=3)    ::  extension
            ! Start work
            extension = ExtensionFromFilename(self%filename)
            select case(extension)
                case ('img','hed')
                    FormatFromFilename = 'I'
                case ('mrc','map','st','ctf','mrcs')
                    FormatFromFilename = 'M'
                case ('spi')
                    FormatFromFilename = 'S'
                case ('tif','tiff')
                    FormatFromFilename = 'T'
                case default
                    FormatFromFilename = 'N'
            end select
        end function FormatFromFilename
end module ImageFiles
