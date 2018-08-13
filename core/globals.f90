module Globals



    use ProgramInstances

    type(ProgramInstance)                   ::  this_program

    integer,                parameter       ::  filename_max_len        =   200             !< Maximum number of characters for a file name
    integer,                parameter       ::  line_max_len            =   8192            !< Maximum number of characters for a line
    integer,                parameter       ::  user_input_max_len      =   50              !< Length used for user input questions
    character(len=1),       parameter       ::  default_file_format     =   'M'             !< The default file format:
                                                                                            !! 'I', 'M' or 'S' for imagic, mrc, spider

    ! Constants
    real,                   parameter       ::  pi                      =   acos(-1.0e0)
    real(kind=8),           parameter       ::  dpi                     =   acos(-1.0d0)
    real,                   parameter       ::  protein_density         =   0.86e0          !<    Protein density in Dalton / cubic Angstrom. Converted from 1.43 g/cm^3, see
                                                                                            !!  Quillin, M. L. and Matthews, B. W., "Accurate calculation of the density of proteins.", Acta Crystallogr D Biol Crystallogr, vol. 56, pp. 791--794, 2000

end module
