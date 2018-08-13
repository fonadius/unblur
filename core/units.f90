module Units
    use Globals
    private

    public  ::  unit_to_string, unit_conversion, convert

    ! A few units we might use
    integer,            parameter,  public  ::  millimeters                 =   0
    integer,            parameter,  public  ::  microns                     =   1
    integer,            parameter,  public  ::  angstroms                   =   2
    integer,            parameter,  public  ::  pixels                      =   3
    integer,            parameter,  public  ::  degrees                     =   4
    integer,            parameter,  public  ::  radians                     =   5
    integer,            parameter,  public  ::  reciprocal_angstroms        =   6
    integer,            parameter,  public  ::  reciprocal_pixels           =   7
    character(len=2),   parameter           ::  millimeters_text            =   'mm'
    character(len=2),   parameter           ::  microns_text                =   'um'
    character(len=1),   parameter           ::  angstroms_text              =   'A'
    character(len=6),   parameter           ::  pixels_text                 =   'pixels'
    character(len=7),   parameter           ::  degrees_text                =   'degrees'
    character(len=7),   parameter           ::  radians_text                =   'radians'
    character(len=3),   parameter           ::  reciprocal_angstroms_text   =   '1/A'
    character(len=8),   parameter           ::  reciprocal_pixels_text      =   '1/pixels'

    contains

    !>  \brief  Return converted value
    real elemental function convert(current_value,current_unit,new_unit,pixel_size) result(output_value)
        real,               intent(in)      ::  current_value
        integer,            intent(in)      ::  current_unit
        integer,            intent(in)      ::  new_unit
        real,   optional,   intent(in)      ::  pixel_size
        ! private variable
        integer ::  temp_unit
        ! start work
        output_value = current_value
        temp_unit = current_unit
        call unit_conversion(output_value,temp_unit,new_unit,pixel_size)
    end function convert

    !>  \brief  Return a converted value
    elemental subroutine unit_conversion(current_value,current_unit,new_unit,pixel_size)
        ! arguments
        real,               intent(inout)   ::  current_value
        integer,            intent(inout)   ::  current_unit    !<  Will be changed by the routine to the desired new unit
        integer,            intent(in)      ::  new_unit
        real,   optional,   intent(in)      ::  pixel_size      !<  Size of a pixel in Angstroms. May be required for the conversion.
        ! private variables
        logical ::  conversion_was_successful
        logical ::  pixel_size_missing
        ! start work
        conversion_was_successful = .true.
        pixel_size_missing = .false.

        ! If we need the pixel size, check we have it
        if (current_unit .ne. new_unit) then
            if (     current_unit   .eq. pixels             &
                .or. new_unit       .eq. pixels             &
                .or. current_unit   .eq. reciprocal_pixels  &
                .or. new_unit       .eq. reciprocal_pixels) then
                ! we need to know the pixel size
                pixel_size_missing = .not. present(pixel_size)
                if (pixel_size_missing) then
                    conversion_was_successful = .false.
                    !write(*,'(a)') '**error(unit_conversion): pixel size was not supplied'
                endif
            endif
        endif

        if (current_unit .ne. new_unit .and. .not. pixel_size_missing) then
            select case (current_unit)
            case (microns)
                select case (new_unit)
                case (pixels)
                    ! microns to pixels
                    current_value = current_value * 1.0e4 / pixel_size
                case default
                    conversion_was_successful = .false.
                end select
            case (millimeters)
                select case (new_unit)
                case (pixels)
                    ! mm to pixels
                    current_value = current_value * 1.0e7 / pixel_size
                case default
                    conversion_was_successful = .false.
                end select
            case (angstroms)
                select case (new_unit)
                case (pixels)
                    ! A to pixels
                    current_value = current_value / pixel_size
                case default
                    conversion_was_successful = .false.
                end select
            case (reciprocal_angstroms)
                select case (new_unit)
                case (reciprocal_pixels)
                    ! 1/A to 1/pixels
                    current_value = current_value * pixel_size
                case default
                    conversion_was_successful = .false.
                end select
            case (pixels)
                select case (new_unit)
                case (angstroms)
                    ! pix to A
                    current_value = current_value * pixel_size
                case default
                    conversion_was_successful = .false.
                end select
            case (degrees)
                select case (new_unit)
                case (radians)
                    ! degrees to radians
                    current_value = current_value / 180.0e0 * pi
                case default
                    conversion_was_successful = .false.
                end select
            case (radians)
                select case (new_unit)
                case (degrees)
                    ! rad to deg
                    current_value = current_value / pi * 180.0e0
                case default
                    conversion_was_successful = .false.
                end select
            case default
                conversion_was_successful = .false.
            end select
        endif

        if (conversion_was_successful) then
            current_unit = new_unit
        else
            !write(*,'(4a)') '**error(unit_conversion): don''t know how to convert ', unit_to_string(current_unit), ' to ', unit_to_string(new_unit)
            !call this_program%TerminateWithFatalError('unit_conversion','failed conversion')
        endif
    end subroutine unit_conversion

    !>  \brief  Return a string useful describing the units
    function unit_to_string(unit_identifier) result(string)
        integer,            intent(in)  ::  unit_identifier
        character(len=:),   allocatable ::  string
        !
        select case (unit_identifier)
            case (millimeters)
                string = millimeters_text
            case (microns)
                string = microns_text
            case (angstroms)
                string = angstroms_text
            case (pixels)
                string = pixels_text
            case (degrees)
                string = degrees_text
            case (radians)
                string = radians_text
            case (reciprocal_angstroms)
                string = reciprocal_angstroms_text
            case (reciprocal_pixels)
                string = reciprocal_pixels_text
            case default
                call this_program%TerminateWithFatalError('unit_string','unknown unit')
        end select
    end function unit_to_string


end module Units
