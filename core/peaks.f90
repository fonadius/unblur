module Peaks

    implicit none

    private

    type, public :: Peak

        real       ::     CoOrdinates(3)  =   [0.0e0,0.0e0,0.0e0] !<  X, Y, Z coordinates
        real       ::     PeakHeight                              !<  Value of the peak / peak height
        real       ::     PeakHeightInSigmas                      !<  Number of times above the standard deviation
        contains
            procedure   ::  PrintInfo
    end type Peak

    contains

    subroutine PrintInfo(self)
        class(Peak),    intent(in)  ::  self
        write(*,'(5(f0.3,1x))') self%CoOrdinates, self%PeakHeight, self%PeakHeightInSigmas
    end subroutine PrintInfo

end module
