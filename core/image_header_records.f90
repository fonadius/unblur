!>  \brief  Image header records are contained in image headers. Each contains a piece of information regarding the image file.
module ImageHeaderRecords
    use Globals
    implicit none
    private
    public :: IntgImageHeaderRecord, RealImageHeaderRecord, CharImageHeaderRecord

    !>  A header record is a labelled value which describes a characteristic of an imagefile
    type :: ImageHeaderRecord
        private
        integer                                     ::  index_position                      !<  the position of the record within the file header. starting at 1 and incrementing
        integer                                     ::  byte_position                       !<  the position of the first byte of the record within the header.
        integer(kind=1),    pointer                 ::  byte_array(:)   =>  null()          !<  pointer to the array of bytes containing the actual header values

        contains
            private
            procedure,  public                      ::  Init
            procedure,  public                      ::  Destroy
    end type ImageHeaderRecord

    type, extends(ImageHeaderRecord) :: IntgImageHeaderRecord
        contains
            private
            procedure                               ::  GetIntg
            procedure,  public                      ::  Get => GetIntg
            procedure                               ::  SetIntg
            generic,    public                      ::  assignment(=)   =>  SetIntg
            final                                   ::  DestructorIntg
    end type

    type, extends(ImageHeaderRecord) :: RealImageHeaderRecord
        contains
            private
            procedure                               ::  GetReal
            procedure,  public                      ::  Get => GetReal
            procedure                               ::  SetReal
            generic,    public                      ::  assignment(=)   =>  SetReal
            final                                   ::  DestructorReal
    end type

    type, extends(ImageHeaderRecord) :: CharImageHeaderRecord
        integer                                     ::  length = 0                          !<  length of the string of characters
        contains
            private
            procedure                               ::  GetChar
            procedure,  public                      ::  Get => GetChar
            procedure                               ::  SetChar
            generic,    public                      ::  assignment(=)   =>  SetChar
            final                                   ::  DestructorChar
    end type


    interface assignment(=)
        module procedure  GetIntgAssign
        module procedure  GetRealAssign
        module procedure  GetCharAssign
    end interface

    contains

        !>  \brief initialise a header record
        subroutine Init(self,index_position,byte_position,header_byte_array,length)
            ! arguments
            class(ImageHeaderRecord),       intent(inout)   ::  self                    !<  header record
            integer,                        intent(in)      ::  index_position
            integer,                        intent(in)      ::  byte_position           !<
            integer(kind=1),    target,     intent(in)      ::  header_byte_array(:)    !<  byte array to point to
            integer,            optional,   intent(in)      ::  length                  !<  length of character string
            ! private variables

            ! start work
            self%index_position =   index_position
            self%byte_position  =   byte_position
            self%byte_array     =>  header_byte_array

            ! if a character
            select type(self)
                type is (CharImageHeaderRecord)
                    self%length = 1
                    if (present(length)) self%length = length
            end select
        end subroutine Init

        !>  \brief Generic destructor
        subroutine Destroy(self)
            class(ImageHeaderRecord),   intent(inout)   ::  self
            select type(self)
                type is (RealImageHeaderRecord)
                    call DestructorReal(self)
                type is (IntgImageHeaderRecord)
                    call DestructorIntg(self)
                type is (CharImageHeaderRecord)
                    call DestructorChar(self)
                class default
                    call this_program%TerminateWithFatalError('ImageHeaderRecord::Destroy','Unknown ImageHeaderRecord type')
            end select
        end subroutine Destroy

        subroutine DestructorReal(self)
            type(RealImageHeaderRecord),    intent(inout)   ::  self
            if (associated(self%byte_array)) nullify(self%byte_array)
        end subroutine DestructorReal
        subroutine DestructorIntg(self)
            type(IntgImageHeaderRecord),    intent(inout)   ::  self
            if (associated(self%byte_array)) nullify(self%byte_array)
        end subroutine DestructorIntg
        subroutine DestructorChar(self)
            type(CharImageHeaderRecord),    intent(inout)   ::  self
            if (associated(self%byte_array)) nullify(self%byte_array)
        end subroutine DestructorChar

        !>  \brief  Clean-up when we're done
        subroutine DestructorSlave(self)
            ! Arguments
            class(ImageHeaderRecord),   intent(inout)   ::  self
            ! Start work
            if (associated(self%byte_array)) nullify(self%byte_array)
        end subroutine DestructorSlave


        !>  \brief  Write value to header record
        subroutine SetIntg(self,value)
            use iso_c_binding
            ! arguments
            class(IntgImageHeaderRecord),       intent(inout)   ::  self    !<  header record
            integer(kind=4),        target,     intent(in)      ::  value   !<  integer value to be written to header
            ! private variables
            integer(kind=1),        pointer     ::  byte_ptr(:)
            type(c_ptr)                         ::  cptr
            integer                             ::  ptr_shape(1)
            ! start work

            ! prepare a byte pointer to the value
            cptr = c_loc(value)
            ptr_shape = [4]
            call c_f_pointer(cptr,byte_ptr,ptr_shape)

            ! copy the bytes over
            self%byte_array(self%byte_position:self%byte_position+3) = byte_ptr(1:4)

            ! nullify pointers
            nullify(byte_ptr)
        end subroutine SetIntg

        !>  \brief  Write value to header record
        subroutine SetReal(self,value)
            use iso_c_binding
            ! arguments
            class(RealImageHeaderRecord),       intent(inout)   ::  self    !<  header record
            real(kind=4),           target,     intent(in)      ::  value   !<  integer value to be written to header
            ! private variables
            integer(kind=1),        pointer     ::  byte_ptr(:)
            type(c_ptr)                         ::  cptr
            integer                             ::  ptr_shape(1)
            ! start work

            ! prepare a byte pointer to the value
            cptr = c_loc(value)
            ptr_shape = [4]
            call c_f_pointer(cptr,byte_ptr,ptr_shape)

            ! copy the bytes over
            self%byte_array(self%byte_position:self%byte_position+3) = byte_ptr(1:4)

            ! nullify pointers
            nullify(byte_ptr)
        end subroutine SetReal

        !>  \brief  Write value to header record
        subroutine SetChar(self,value)
            use StringManipulations, only : IntegerToString
            ! arguments
            class(CharImageHeaderRecord),       intent(inout)   ::  self    !<  header record
            character(len=*),                   intent(in)      ::  value   !<  integer value to be written to header
            ! private variables
            character(len=:), allocatable   ::  tmp_string
            ! start work
            allocate(character(len=self%length) :: tmp_string)
            write(tmp_string,'(a'//IntegerToString(self%length)//')') value

            ! copy the bytes over
            self%byte_array(self%byte_position:self%byte_position+self%length-1) = transfer(tmp_string,self%byte_array)

        end subroutine SetChar

        !>  \brief  Read value from header record (to be used on RHS of assignments)
        subroutine GetIntgAssign(value,self)
            class(IntgImageHeaderRecord),   intent(in)      ::  self
            integer(kind=4),                intent(out)     ::  value
            ! start work
            value =  self%GetIntg()
        end subroutine GetIntgAssign

        !>  \brief  Read value from header record (to be used on RHS of assignments)
        subroutine GetRealAssign(value,self)
            class(RealImageHeaderRecord),   intent(in)      ::  self
            real(kind=4),                   intent(out)     ::  value
            ! start work
            value =  self%GetReal()
        end subroutine GetRealAssign

        !>  \brief  Read value from header record (to be used on RHS of assignments)
        subroutine GetCharAssign(value,self)
            class(CharImageHeaderRecord),   intent(in)      ::  self
            character(len=:),   allocatable,intent(inout)   ::  value
            ! start work
            if (allocated(value)) deallocate(value)
            value =  self%GetChar()
        end subroutine GetCharAssign

        !>  \brief  Read value from header record
        integer function GetIntg(self)
            use iso_c_binding
            ! arguments
            class(IntgImageHeaderRecord),   intent(in)      ::  self    !<  header record
            ! private variables
            integer(kind=1),    pointer     ::  byte_ptr(:)
            type(c_ptr)                     ::  cptr
            integer(kind=4),    target      ::  tmpval
            integer                         ::  ptr_shape(1)
            ! start work

            ! prepare a byte pointer to the temp value
            cptr = c_loc(tmpval)
            ptr_shape = [4]
            call c_f_pointer(cptr,byte_ptr,ptr_shape)

            ! copy bytes over
            byte_ptr(1:4) = self%byte_array(self%byte_position:self%byte_position+3)

            ! nullify pointers
            nullify(byte_ptr)

            ! the integer is now copied over to the result
            GetIntg = tmpval
        end function GetIntg

        !>  \brief  Read value from header record
        real function GetReal(self)
            use iso_c_binding
            ! arguments
            class(RealImageHeaderRecord),   intent(in)      ::  self    !<  header record
            ! private variables
            integer(kind=1),    pointer     ::  byte_ptr(:)
            type(c_ptr)                     ::  cptr
            real(kind=4),       target      ::  tmpval
            integer                         ::  ptr_shape(1)
            ! start work

            ! prepare a byte pointer to the temp value
            cptr = c_loc(tmpval)
            ptr_shape = [4]
            call c_f_pointer(cptr,byte_ptr,ptr_shape)

            ! copy bytes over
            byte_ptr(1:4) = self%byte_array(self%byte_position:self%byte_position+3)

            ! nullify pointers
            nullify(byte_ptr)

            ! the integer is now copied over to the result
            GetReal = tmpval
        end function GetReal

        !>  \brief  Read value from header record
        function GetChar(self)
            use iso_c_binding
            ! arguments
            class(CharImageHeaderRecord),   intent(in)      ::  self    !<  header record
            ! result
            character(len=:), allocatable   ::  GetChar
            ! start work

            ! allocate the output
            allocate(character(len=self%length) :: GetChar)

            GetChar = transfer(self%byte_array(self%byte_position:self%byte_position+self%length-1),GetChar)

        end function GetChar



end module ImageHeaderRecords
