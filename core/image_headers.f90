!>  \brief  Module with type and routine definitions to deal with image file headers
!!
!!  The following file formats are supported:
!!  - imagic: http://imagescience.de/formats/index.htm
!!  - spider: http://www.wadsworth.org/spider_doc/spider/docs/image_doc.html
!!  - mrc: http://www2.mrc-lmb.cam.ac.uk/image2000.html
!!
!!  \todo Add destructor(s)
module ImageHeaders

    use Globals
    use ImageHeaderRecords
    use iso_c_binding

    implicit none

    private
    public  ::  ImageHeader, MrcImageHeader

    integer,    parameter,  public  ::  data_are_bytes      =   1
    integer,    parameter,  public  ::  data_are_integer    =   2
    integer,    parameter,  public  ::  data_are_float      =   3

    type ImageHeader
        private
        integer                                     ::  length                  !<  length (in bytes) of the header
        integer(kind=1),            allocatable     ::  byte_array(:)           !<  array of bytes (why is this public??)
        contains
            procedure, public   ::  Init
            procedure, public   ::  Destroy
            procedure, public   ::  PrintInfo
            procedure, public   ::  ReadFromDisk
            procedure, public   ::  WriteToDisk
            procedure, public   ::  ResetToDefaults
            procedure, public   ::  HasLocalEndianess
            procedure           ::  SetMachineStamp
            procedure, public   ::  GetPixelSize
            procedure, public   ::  SetPixelSize
            procedure, public   ::  GetMinimumPixelValue
            procedure, public   ::  GetMaximumPixelValue
            procedure, public   ::  SetMinimumPixelValue
            procedure, public   ::  SetMaximumPixelValue
            procedure, public   ::  GetStackSize
            procedure, public   ::  GetDimensions
            procedure, public   ::  GetDimension
            procedure, public   ::  SetDimensions
            procedure, public   ::  FirstDataByte
            procedure, public   ::  Assign
            procedure, public   ::  IsAllocated
            procedure, public   ::  BytesPerPixel
            procedure, public   ::  PixelDataAreSigned
            procedure, public   ::  GetPixelDataType
            procedure, public   ::  PixelDataAreComplex
            generic,   public   ::  assignment(=)   =>  Assign
    end type ImageHeader

    type, extends(ImageHeader) :: MrcImageHeader
        type(IntgImageHeaderRecord)                 ::  nx      !<  number of columns (fastest changing in map)
        type(IntgImageHeaderRecord)                 ::  ny      !<  number of rows
        type(IntgImageHeaderRecord)                 ::  nz      !<  number of sections (slowest changing in map)
        type(IntgImageHeaderRecord)                 ::  mode    !<  data type:  0 image: signed 8-bit bytes rante -128 to 127
                                                                !!              1 image: 16-bit halfwords
                                                                !!              2 image: 32-bit reals
                                                                !!              3 transform: complex 16-bit integers
                                                                !!              4 transform: complex 32-bit reals
        type(IntgImageHeaderRecord)                 ::  nxstart !<  number of first column in map (default = 0)
        type(IntgImageHeaderRecord)                 ::  nystart !<  number of first row in map
        type(IntgImageHeaderRecord)                 ::  nzstart !<  number of first section in map
        type(IntgImageHeaderRecord)                 ::  mx      !<  number of intervals along x
        type(IntgImageHeaderRecord)                 ::  my      !<  number of intervals along y
        type(IntgImageHeaderRecord)                 ::  mz      !<  number of intervals along z
        type(RealImageHeaderRecord)                 ::  cella1  !<  cell dimensions in angstroms
        type(RealImageHeaderRecord)                 ::  cella2
        type(RealImageHeaderRecord)                 ::  cella3
        type(RealImageHeaderRecord)                 ::  cellb1  !<  cell angles in degrees
        type(RealImageHeaderRecord)                 ::  cellb2
        type(RealImageHeaderRecord)                 ::  cellb3
        type(IntgImageHeaderRecord)                 ::  mapc    !<  axis corresponding to columns  (1,2,3 for x,y,z)
        type(IntgImageHeaderRecord)                 ::  mapr    !<  axis corresponding to rows     (1,2,3 for x,y,z)
        type(IntgImageHeaderRecord)                 ::  maps    !<  axis corresponding to sections (1,2,3 for x,y,z)
        type(RealImageHeaderRecord)                 ::  dmin    !<  minimum density value
        type(RealImageHeaderRecord)                 ::  dmax    !<  maximum density value
        type(RealImageHeaderRecord)                 ::  dmean   !<  mean density value
        type(IntgImageHeaderRecord)                 ::  ispg    !<  space group number 0 or 1 (default=0)
        type(IntgImageHeaderRecord)                 ::  nsymbt  !<  number of bytes used for symmetry data (0 or 80)
        type(IntgImageHeaderRecord)                 ::  extra   !<  extra space used for anything - 0 by default
        type(IntgImageHeaderRecord)                 ::  originx !<  origin in x used for transforms
        type(IntgImageHeaderRecord)                 ::  originy !<  origin in y used for transforms
        type(IntgImageHeaderRecord)                 ::  originz !<  origin in z used for transforms
        type(CharImageHeaderRecord)                 ::  map     !<  character string 'map' to identify file type
        type(IntgImageHeaderRecord)                 ::  machst  !<  machine stamp
        type(RealImageHeaderRecord)                 ::  rms     !<  rms deviation of map from mean density
        type(IntgImageHeaderRecord)                 ::  nlabl   !<  number of labels being used
        type(CharImageHeaderRecord)                 ::  label(10)!<  10 80-character text labels
        contains
            final                                   ::  DestructorMRC
    end type MrcImageHeader


    contains

        subroutine Destroy(self)
            class(ImageHeader), intent(inout)   ::  self
            select type(self)
                type is (MrcImageHeader)
                    call DestructorMRC(self)
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::Destroy','Unsupported header type')
            end select
        end subroutine Destroy



        subroutine DestructorMRC(self)
            type(MrcImageHeader),  intent(inout)   ::  self
            !
            integer ::  current_label
            !
            call self%nx     %Destroy()
            call self%ny     %Destroy()
            call self%nz     %Destroy()
            call self%mode   %Destroy()
            call self%nxstart%Destroy()
            call self%nystart%Destroy()
            call self%nzstart%Destroy()
            call self%mx     %Destroy()
            call self%my     %Destroy()
            call self%mz     %Destroy()
            call self%cella1 %Destroy()
            call self%cella2 %Destroy()
            call self%cella3 %Destroy()
            call self%cellb1 %Destroy()
            call self%cellb2 %Destroy()
            call self%cellb3 %Destroy()
            call self%mapc   %Destroy()
            call self%mapr   %Destroy()
            call self%maps   %Destroy()
            call self%dmin   %Destroy()
            call self%dmax   %Destroy()
            call self%dmean  %Destroy()
            call self%ispg   %Destroy()
            call self%nsymbt %Destroy()
            call self%extra  %Destroy()
            call self%originx%Destroy()
            call self%originy%Destroy()
            call self%originz%Destroy()
            call self%map    %Destroy()
            call self%machst %Destroy()
            call self%rms    %Destroy()
            call self%nlabl  %Destroy()
            do current_label=1,10
                call self%label(current_label)  %Destroy()
            enddo
            if (allocated(self%byte_array)) deallocate(self%byte_array)
        end subroutine DestructorMRC



        !>  \brief  Return the local machine's "machinestamp", which is endian-specific
        function GetLocalMachineStamp()
            !use iso_c_binding
            integer(kind=4) ::  GetLocalMachineStamp
            ! Private variables
            integer(kind=1)                 ::  machst(4)
            integer(kind=4),  parameter     ::  a0  =   48
            integer(kind=4),  parameter     ::  a1  =   49
            integer(kind=4),  parameter     ::  a2  =   50
            integer(kind=4),  parameter     ::  a3  =   51
            integer(kind=4)                 ::  i
            character(len=4)                ::  ich
            !type(c_ptr)                     ::  i_cptr
            logical,            parameter   ::  debug = .false.
            ! Start work
            i=a0+a1*256+a2*(256**2)+a3*(256**3) !  = 858927408 (decimal)
                                                !  = 0011 0011 0011 0010 0011 0001 0011 0000 (binary, little endian)
                                                !  when this is converted to ASCII characters (1 byte per character, with the most significant bit always 0)
                                                ! this will give different results on little- and big-endian machines
                                                ! For example, '0' in ASCII has decimal value 48 and bit value 011 0000
                                                ! '3' in ASCII has decimal value 51 and bit value 011 0011
                                                ! Therefore the value computed above, when converted to bytes will have the first byte
                                                ! read off as ASCII character '0' on little-endian and '3' on big-endian machines

            ! Take the bit pattern over from the 4byte integer to an array of 4 characters (each 1 byte)
            ich = transfer(i,ich)
            !! Create a c pointer
            !i_cptr = c_loc(i)
            !! convert to a f pointer
            !call c_f_pointer(i_cptr,ich)

            if (ich.eq.'0123') then
                if (debug) write(*,'(a)') '**debug(GetLocalMachineStamp): machine is little-endian (dec/osf, intel, amd ...)'
                !0100 0100
                machst(1)=68
                !0100 0001
                machst(2)=65
                machst(3)=0
                machst(4)=0
            elseif (ich.eq.'3210') then
                if (debug) write(*,'(a)') '**debug(GetLocalMachineStamp): machine is big-endian (sgi, sun, hp, ibm)'
                !0001 0001
                machst(1)=17
                !0001 0001
                machst(2)=17
                machst(3)=0
                machst(4)=0
            else
                if (debug) write(*,'(a)') '**debug(GetLocalMachineStamp): mixed endianity machine (vax)'
                !0010 0010
                machst(1)=34
                !0010 0001
                machst(2)=33
                machst(3)=0
                machst(4)=0
            endif
            ! Convert machst (4 bytes) to a 4-byte integer
            GetLocalMachineStamp = transfer(machst,GetLocalMachineStamp)
        end function GetLocalMachineStamp



        !>  \brief  Check allocation status of byte array
        pure logical function IsAllocated(self)
            ! Arguments
            class(ImageHeader),     intent(in)      ::  self
            ! Start work
            IsAllocated = allocated(self%byte_array)
        end function

        !>  \brief  Return the number of bytes per pixel
        integer function BytesPerPixel(self)
            ! Argument
            class(ImageHeader),     intent(in)      ::  self
            ! Start work
            select type(self)
                type is (MrcImageHeader)
                    select case(self%mode%Get())
                        case(0)
                            BytesPerPixel = 1
                        case(1,3,6)
                            BytesPerPixel = 2
                        case(2,4)
                            BytesPerPixel = 4
                        case default
                            call this_program%TerminateWithFatalError('ImageHeader::BytesPerPixel','Nonsensical MRC mode number')
                    end select
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::BytesPerPixel','Format not supported')
            end select
        end function BytesPerPixel

        !>  \brief  Does the header indicate that pixel density values are stored in unsigned format?
        logical function PixelDataAreSigned(self)
            ! Argument
            class(ImageHeader),     intent(in)      ::  self
            ! Start work
            select type(self)
                type is (MrcImageHeader)
                    select case (self%mode%Get())
                        case (0)
                            ! Note that MRC mode 0 is sometimes signed, sometimes unsigned. TODO: sort this out by checking the imodStamp header
                            PixelDataAreSigned = .true.
                        case (1,2,3,4)
                            PixelDataAreSigned = .true.
                        case (6,16)
                            PixelDataAreSigned = .false.
                        case default
                            call this_program%TerminateWithFatalError('ImageHeader::PixelDataAreSigned', &
                                                                        'Nonsensical MRC mode number')
                    end select
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::PixelDataAreSigned','Format not supported')
            end select
        end function PixelDataAreSigned

        !>  \brief  Work out whether pixel data are byte, integer or float
        integer function GetPixelDataType(self)
            ! Argument
            class(ImageHeader),     intent(in)      ::  self
            ! Start work
            select type(self)
                type is (MrcImageHeader)
                    select case (self%mode%Get())
                        case (0)
                            GetPixelDataType = data_are_bytes
                        case (1,3,6)
                            GetPixelDataType = data_are_integer
                        case (2,4)
                            GetPixelDataType = data_are_float
                        case default
                            call this_program%TerminateWithFatalError('ImageHeader::GetPixelDataType','Nonsensical MRC mode number')
                    end select
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::GetPixelDataType','Format not supported')
            end select
        end function GetPixelDataType

        !>  \brief  Does the header indicate that pixel density values are complex numbers
        logical function PixelDataAreComplex(self)
            ! Argument
            class(ImageHeader),     intent(in)      ::  self
            ! Start work
            select type(self)
                type is (MrcImageHeader)
                    select case (self%mode%Get())
                        case (0,1,2,6)
                            PixelDataAreComplex = .false.
                        case (3,4)
                            PixelDataAreComplex = .true.
                        case default
                            call this_program%TerminateWithFatalError('ImageHeader::PixelDataAreComplex', &
                                                                      'Nonsensical MRC mode number')
                    end select
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::PixelDataAreComplex','Format not supported')
            end select
        end function PixelDataAreComplex


        !>  \brief  Return the index of the first byte containing image data
        integer function FirstDataByte(self)
            ! Argument
            class(ImageHeader),     intent(in)      ::  self
            ! Start work
            select type(self)
                type is (MrcImageHeader)
                    FirstDataByte = 1025 + self%nsymbt%Get()
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::FirstDataByte','Format not supported')
            end select
        end function FirstDataByte

        !>  \brief  Assignment overloader
        subroutine Assign(lhs,rhs)
            ! Arguments
            class(ImageHeader),     intent(inout)   ::  lhs
            class(ImageHeader),     intent(in)      ::  rhs
            ! Start work
            select type(lhs)
                type is (MrcImageHeader)
                    select type(rhs)
                        type is (MrcImageHeader)
                            call lhs%Init()
                            lhs%byte_array = rhs%byte_array
                        class default
                            call this_program%TerminateWithFatalError('ImageHeader::Assign','Format not supported (RHS)')
                    end select
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::Assign','Format not supported (LHS)')
            end select
        end subroutine Assign


        !>  \brief  Read the header data from disk
        subroutine ReadFromDisk(self,lun)
            ! arguments
            class(ImageHeader),     intent(inout)   ::  self
            integer,                intent(in)      ::  lun
            ! private variables
            integer ::  io_status
            character(len=512) :: io_message
            ! start work
            select type(self)
                type is (MrcImageHeader)
                    read(unit=lun,pos=1,iostat=io_status,iomsg=io_message) self%byte_array
                    if (io_status .ne. 0) then
                        write(*,'(a,i0,2a)') '**error(ImageHeader::ReadFromDisk): error ', io_status, &
                                             ' when reading header bytes from disk: ', trim(io_message)
                        call this_program%TerminateWithFatalError('ImageHeader::ReadFromDisk','I/O error')
                    endif
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::ReadFromDisk','Format not supported')
            end select
        end subroutine ReadFromDisk

        !>  \brief  Read the header data from disk
        subroutine WriteToDisk(self,lun)
            ! arguments
            class(ImageHeader),     intent(inout)   ::  self
            integer,                intent(in)      ::  lun
            ! private variables
            integer ::  io_status
            ! start work
            call self%SetMachineStamp()
            write(unit=lun,pos=1,iostat=io_status) self%byte_array
            if (io_status .ne. 0) then
                write(*,'(a,i0,a)') '**error(ImageHeader::WriteToDisk): error ', io_status, ' when writing header bytes to disk'
                print *, allocated(self%byte_array)
                print *, lun
                call this_program%TerminateWithFatalError('ImageHeader::WriteToDisk','I/O error')
            endif
        end subroutine WriteToDisk

        !>  \brief  Print out information contained in the header to the terminal
        subroutine PrintInfo(self)
            ! arguments
            class(ImageHeader),     intent(in)  ::  self
            ! private variables
            integer ::  current_label
            real    ::  pixel_size(3)
            integer ::  bytes_per_pixel
            ! Start work
            select type(self)
                type is (MrcImageHeader)
                    write(*,'(a,3(i0,1x))')     'Number of columns, rows, sections: ', self%nx%Get(), self%ny%Get(), self%nz%Get()
                    write(*,'(a,i0)')           'MRC data mode: ',  self%mode%Get()
                    bytes_per_pixel = self%BytesPerPixel()
                    write(*,'(a,i0)')           'Bit depth: ',      bytes_per_pixel*8
                    pixel_size = 0.0
                    if (self%mx%Get() .ne. 0) pixel_size(1) = self%cella1%Get()/self%mx%Get()
                    if (self%my%Get() .ne. 0) pixel_size(2) = self%cella2%Get()/self%my%Get()
                    if (self%mz%Get() .ne. 0) pixel_size(3) = self%cella3%Get()/self%mz%Get()
                    write(*,'(a,3(f0.3,1x))')   'Pixel size: ', pixel_size
                    do current_label=1,self%nlabl%Get()
                        write(*,'(a,i2,2a)') 'Label ', current_label, ': ', self%label(current_label)%Get()
                    enddo
                class default
                    call this_program%TerminateWithFatalError('ImageHeaders::PrintInfo','Unsupported file format')
            end select
        end subroutine PrintInfo

        !>  \brief  Reset all the values in the header to default values
        subroutine ResetToDefaults(self)
            ! Arguments
            class(ImageHeader), intent(inout)   ::  self
            ! Private variables
            character(len=1)    ::  blank_string
            integer             ::  i
            ! Start work
            blank_string = ' '
            select type(self)
                type is (MrcImageHeader)
                    self%nx          = 0
                    self%ny          = 0
                    self%nz          = 0
                    self%mode        = 2
                    self%nxstart     = 0
                    self%nystart     = 0
                    self%nzstart     = 0
                    self%mx          = 1
                    self%my          = 1
                    self%mz          = 1
                    self%cella1      = 1.0
                    self%cella2      = 1.0
                    self%cella3      = 1.0
                    self%cellb1      = 90.0
                    self%cellb2      = 90.0
                    self%cellb3      = 90.0
                    self%mapc        = 1
                    self%mapr        = 2
                    self%maps        = 3
                    self%dmin        = 0.0
                    self%dmax        = 0.0
                    self%dmean       = 0.0
                    self%ispg        = 0
                    self%nsymbt      = 0
                    self%extra       = 0
                    self%originx     = 0
                    self%originy     = 0
                    self%originz     = 0
                    self%map         = 'MAP '
                    self%machst      = 0
                    self%rms         = 0.0
                    self%nlabl       = 0
                    do i=1,10
                    self%label(i) = ' '
                    enddo
                class default
                    call this_program%TerminateWithFatalError('ResetToDefaults','Format not supported')
            end select
        end subroutine ResetToDefaults

        !>  \brief  Return the minimum pixel value
        real function GetMinimumPixelValue(self)
            class(ImageHeader), intent(in)  ::  self
            !
            select type(self)
                type is (MrcImageHeader)
                    GetMinimumPixelValue = self%dmin%Get()
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::GetMinimumPixelValue','Format not supported')
            end select
        end function GetMinimumPixelValue

        !>  \brief  Return the maximum pixel value
        real function GetMaximumPixelValue(self)
            class(ImageHeader), intent(in)  ::  self
            !
            select type(self)
                type is (MrcImageHeader)
                    GetMaximumPixelValue = self%dmax%Get()
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::GetMinimumPixelValue','Format not supported')
            end select
        end function GetMaximumPixelValue

        !>  \brief  Set the minimum pixel value
        subroutine SetMinimumPixelValue(self,new_value)
            class(ImageHeader), intent(inout)   ::  self
            real,               intent(in)      ::  new_value
            select type(self)
                type is (MrcImageHeader)
                    self%dmin = new_value
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::SetMinimumPixelValue','Format not supported')
            end select
        end subroutine SetMinimumPixelValue

        !>  \brief  Set the maximum pixel value
        subroutine SetMaximumPixelValue(self,new_value)
            class(ImageHeader), intent(inout)   ::  self
            real,               intent(in)      ::  new_value
            select type(self)
                type is (MrcImageHeader)
                    self%dmax = new_value
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::SetMinimumPixelValue','Format not supported')
            end select
        end subroutine SetMaximumPixelValue

        !>  \brief  Return the pixel size
        real function GetPixelSize(self)
            ! Arguments
            class(ImageHeader), intent(in)  ::  self
            ! Start work
            select type(self)
                type is (MrcImageHeader)
                    if (self%mx%Get() .ne. 0) then
                        GetPixelSize = self%cella1%Get()/self%mx%Get()
                    else
                        GetPixelSize = 0.0
                    endif
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::GetPixelSize','Format not supported')
            end select
        end function GetPixelSize

        !>  \brief  Set the pixel size
        subroutine SetPixelSize(self,pixel_size)
            ! Arguments
            class(ImageHeader), intent(inout)   ::  self
            real,               intent(in)      ::  pixel_size
            ! Start work
            select type(self)
                type is (MrcImageHeader)
                    self%cella1 = pixel_size * self%mx%Get()
                    self%cella2 = pixel_size * self%my%Get()
                    self%cella3 = pixel_size * self%mz%Get()
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::SetPixelSize','Format not supported')
            end select
        end subroutine SetPixelSize

        !> \brief   Set the machine stamp
        subroutine SetMachineStamp(self)
            ! Arguments
            class(ImageHeader), intent(inout)   ::  self
            ! Start work
            select type(self)
                type is (MrcImageHeader)
                    self%machst = GetLocalMachineStamp()
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::SetMachineStamp','Format not supported')
            end select
        end subroutine SetMachineStamp

        !>  \brief  Check whether the header was created by a machine with the same endianess as the machine we're on now
        !!          If the header doesn't have a machine stamp (or it is 0), we assume it's in the local endianess.
        logical function HasLocalEndianess(self)
            class(ImageHeader), intent(in)      ::  self

            !select type(self)
            !    type is (MrcImageHeader)
            !        HasLocalEndianess = (GetLocalMachineStamp() .eq. self%machst%Get()) .or. self%machst%Get() .eq. 0
            !    class default
            !        call this_program%TerminateWithFatalError('ImageHeader::SetMachineStamp','Format not supported')
            !end select

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SETTING THIS TO ALWAYS RETURN TRUE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! It seems most programs ignore the endianess setting in the header, which should basically all be the !
            ! same these days. As such I am setting this function to always return true, the code is being
            ! left in place so that we can revert this back at a later date if required.
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SETTING THIS TO ALWAYS RETURN TRUE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            HasLocalEndianess = .true.
        end function HasLocalEndianess

        !>  \brief  Return the number of 2D images in the stack
        integer function GetStackSize(self) result(stack_size)
            ! Arguments
            class(ImageHeader), intent(in)      ::  self
            ! Start work
            select type(self)
                type is (MrcImageHeader)
                    stack_size = self%nz%Get()
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::GetStackSize','Format not supported')
            end select
        end function GetStackSize


        !>  \brief  Return the dimensions stored in the header
        function GetDimensions(self) result(dimensions)
            ! Arguments
            class(ImageHeader), intent(in)      ::  self
            ! Result
            integer ::  dimensions(3)
            ! Start work
            select type(self)
                type is (MrcImageHeader)
                    dimensions = [self%nx%Get(),self%ny%Get(),self%nz%Get()]
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::GetDimensions','Format not supported')
            end select
        end function GetDimensions

        !>  \brief  Return one of the dimensions stored in the header
        function GetDimension(self,which_dimension) result(dimension)
            ! Arguments
            class(ImageHeader), intent(in)      ::  self
            integer,            intent(in)      ::  which_dimension
            ! Result
            integer ::  dimension
            ! Start work
            select type(self)
                type is (MrcImageHeader)
                    select case (which_dimension)
                        case (1)
                            dimension = self%nx%Get()
                        case (2)
                            dimension = self%ny%Get()
                        case (3)
                            dimension = self%nz%Get()
                        case default
                            call this_program%TerminateWithFatalError('ImageHeader::GetDimension','Dimension should be 1, 2 or 3')
                    end select
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::GetDimension','Format not supported')
            end select
        end function GetDimension

        !>  \brief  Set the dimensions stored in the header
        subroutine SetDimensions(self,dimensions)
            ! Arguments
            class(ImageHeader), intent(inout)   ::  self
            integer,            intent(in)      ::  dimensions(3)
            ! Start work
            select type(self)
                type is (MrcImageHeader)
                    self%nx = dimensions(1)
                    self%ny = dimensions(2)
                    self%nz = dimensions(3)
                    self%mx = dimensions(1)
                    self%my = dimensions(2)
                    self%mz = dimensions(3)
                class default
                    call this_program%TerminateWithFatalError('ImageHeader::SetDimensions','Format not supported')
            end select
        end subroutine SetDimensions

        !>  \brief  initialise a header
        subroutine Init(self,length)
            ! arguments
            class(ImageHeader),     target,     intent(inout)   ::  self                !<  image header object
            integer,                optional,   intent(in)      ::  length              !<  length of the header record (number of bytes). defaults to 1024 if not given.
            ! private variables
            integer             ::  llength
            integer             ::  ierr
            integer             ::  current_label
            character(len=100)  ::  error_message
            ! start work

            ! length
            if (present(length)) then
                llength = length
            else
                llength = 1024
            endif

            ! deallocate if necessary
            if (allocated(self%byte_array)) then
                if (size(self%byte_array) .ne. llength) deallocate(self%byte_array)
            endif

            ! allocate memory for the byte array
            if (.not. allocated(self%byte_array)) then
                allocate(self%byte_array(llength),stat=ierr,errmsg=error_message)
                if (ierr .ne. 0) then
                    write(*,'(a,i0,2a)') '**error(ImageHeader::Init): memory allocation failed with error ', ierr, ': ', &
                                                                                            trim(adjustl(error_message))
                    call this_program%TerminateWithFatalError('ImageHeader::Init', 'Failed memory allocation')
                endif
            endif

            ! Zero the byte array
            self%byte_array = 0

            !
            ! do work specific to the different types of headers
            !
            select type(self)
                type is (MrcImageHeader)
                    ! initialise
                    call self%nx     %init  (1,     1,    self%byte_array)
                    call self%ny     %init  (2,     5,    self%byte_array)
                    call self%nz     %init  (3,     9,    self%byte_array)
                    call self%mode   %init  (4,     13,   self%byte_array)
                    call self%nxstart%init  (5,     17,   self%byte_array)
                    call self%nystart%init  (6,     21,   self%byte_array)
                    call self%nzstart%init  (7,     25,   self%byte_array)
                    call self%mx     %init  (8,     29,   self%byte_array)
                    call self%my     %init  (9,     33,   self%byte_array)
                    call self%mz     %init  (10,    37,   self%byte_array)
                    call self%cella1 %init  (11,    41,   self%byte_array)
                    call self%cella2 %init  (12,    45,   self%byte_array)
                    call self%cella3 %init  (13,    49,   self%byte_array)
                    call self%cellb1 %init  (14,    53,   self%byte_array)
                    call self%cellb2 %init  (15,    57,   self%byte_array)
                    call self%cellb3 %init  (16,    61,   self%byte_array)
                    call self%mapc   %init  (17,    65,   self%byte_array)
                    call self%mapr   %init  (18,    69,   self%byte_array)
                    call self%maps   %init  (19,    73,   self%byte_array)
                    call self%dmin   %init  (20,    77,   self%byte_array)
                    call self%dmax   %init  (21,    81,   self%byte_array)
                    call self%dmean  %init  (22,    85,   self%byte_array)
                    call self%ispg   %init  (23,    89,   self%byte_array)
                    call self%nsymbt %init  (24,    93,   self%byte_array)
                    call self%extra  %init  (25,    97,   self%byte_array)
                    call self%originx%init  (50,    197,  self%byte_array)
                    call self%originy%init  (51,    201,  self%byte_array)
                    call self%originz%init  (52,    205,  self%byte_array)
                    call self%map    %init  (53,    209,  self%byte_array, length=4)
                    call self%machst %init  (54,    213,  self%byte_array)
                    call self%rms    %init  (55,    217,  self%byte_array)
                    call self%nlabl  %init  (56,    221,  self%byte_array)
                    do current_label=1,10
                        call self%label(current_label)  %init  (57,    225,  self%byte_array, length=80)
                    enddo

                    ! Initialise the area of memory to all 0
                    self%byte_array = 0
            class default
                call this_program%TerminateWithFatalError('ImageHeader::Init', 'Unsupported image file format')
            end select
        end subroutine Init



end module ImageHeaders
