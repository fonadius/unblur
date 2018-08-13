!>  \brief  class & methods defining run-time parameters and their behaviour
module UserSuppliedParameters

    implicit none

    private

    type,   abstract,   public  ::  UserSuppliedParameter
        character(len=:), allocatable       ::  keyword                 !<  keyword used in command files to identify the parameter
        character(len=:), allocatable       ::  description             !<  description
        logical                             ::  set_by_user =   .false. !<  flag to indicate whether the value has been set by the user, or is default

    contains

        procedure                           ::  PrintInfo
    end type

    type, extends(UserSuppliedParameter), public :: UserSuppliedReal
        real                           ::  value
    end type

    type, extends(UserSuppliedParameter), public :: UserSuppliedInteger
        integer                        ::  value
     end type

    type, extends(UserSuppliedParameter), public :: UserSuppliedLogical
        logical                        ::  value
    end type

    type, extends(UserSuppliedParameter), public :: UserSuppliedFilename
        character(len=:), allocatable  ::  value
    end type


    contains

    !>  \brief  print info about the parameter, over multiple lines if necessary
    subroutine PrintInfo(self)

        ! arguments
        class(UserSuppliedParameter),     intent(in)  ::  self
        ! private variables
        integer,            parameter ::  keyword_column_width = 16
        character(len=4),   parameter ::  keyword_column_width_c   =   '16'
        integer,            parameter ::  value_column_width = 45
        character(len=4),   parameter ::  value_column_width_c   =   '45'
        integer,            parameter ::  setbyuser_column_width = 7
        character(len=4),   parameter ::  setbyuser_column_width_c   =   '7'
        integer,            parameter ::  description_column_width = 64
        character(len=4),   parameter ::  description_column_width_c   =   '64'

        ! start work

        ! First off just write out the keyword..

        write(*,'(a'//keyword_column_width_c//' ,1x)',advance='no') trim(adjustl(self%keyword))

        select type (self)
!            type is (UserSuppliedParameter)
!            ! Do Nothing - this should never be called
            class is (UserSuppliedReal)
                write(*,'(f'//value_column_width_c//'.4  ,1x)',advance='no') self%value
            class is (UserSuppliedInteger)
                write(*,'(i'//value_column_width_c//'  ,1x)',advance='no') self%value
            class is (UserSuppliedLogical)
                write(*,'(l'//value_column_width_c//'  ,1x)',advance='no') self%value
            class is (UserSuppliedFilename)
                write(*,'(a'//value_column_width_c//'  ,1x)',advance='no') self%value
        end select

        write(*,'(l'//setbyuser_column_width_c//'  ,1x)',advance='no') self%set_by_user
        write(*,'(a)') ' '
    end subroutine PrintInfo
end module UserSuppliedParameters
