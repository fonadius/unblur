!>  \brief  Classes & methods to deal with electron dose, optimal filters etc
module ElectronDoses

    use Globals

    implicit none

    real,   public, parameter       ::  critical_dose_at_dc = huge(1.0e0) * 0.001

    type,   public ::   ElectronDose
        private
        real    ::  acceleration_voltage    =   300.0                               !<  300 or 200 keV

        ! Critical dose function parameters
        !real    ::  critical_dose_scaling =  0.0
        !real    ::  critical_dose_power   =  0.0

        real    ::  critical_dose_a  =  0.0
        real    ::  critical_dose_b  =  0.0
        real    ::  critical_dose_c  =  0.0

        real    ::  voltage_scaling_factor = 0.0


        contains
            procedure,  public  ::  Init
            procedure,  public  ::  ApplyDoseFilterToImage
            procedure,  public  ::  GetRealDoseFilterImageAsIfComplexSize
            procedure           ::  CriticalDose
            procedure           ::  OptimalDoseGivenCriticalDose
            procedure           ::  SignalToNoiseFromDoseGivenCriticalDose
            procedure           ::  electron_dose_assign
            procedure,  public, nopass  ::  DoseFilter
            generic,    public  ::  assignment(=) => electron_dose_assign
    end type

    contains

        !>  \brief  Initialise an ElectronDose object
        subroutine Init(self,acceleration_voltage)
            use StringManipulations, only : RealToString
            ! arguments
            class(ElectronDose),    intent(inout)   ::  self
            real,                   intent(in)      ::  acceleration_voltage    !<  300.0 or 200.0 keV
            ! private variables

            ! start work

            ! Acceleration voltage
            if (acceleration_voltage .lt. 301.0 .and. acceleration_voltage .gt. 299.0) then
                self%acceleration_voltage = 300.0
                self%voltage_scaling_factor = 1.0e0
            else if (acceleration_voltage .lt. 201.0 .and. acceleration_voltage .gt. 199.0) then
                self%acceleration_voltage = 200.0
                self%voltage_scaling_factor = 0.8e0
            else
                call this_program%TerminateWithFatalError('ElectronDose::Init',&
                            'Bad acceleration voltage: '//RealToString(acceleration_voltage,2))
            endif

            ! Set up the critical curve function

            self%critical_dose_a =  0.24499
            self%critical_dose_b =  -1.6649
            self%critical_dose_c =  2.8141

            !self%critical_dose_power   = -1.22776
            !self%critical_dose_scaling =  0.7078384

            !if (self%acceleration_voltage .gt. 250.0) then
                ! Accel voltage is 300 kV
           !     self%critical_dose_scaling = self%critical_dose_scaling * 1.25
            !endif

        end subroutine Init


        !>  \brief  Apply a dose filter to the image
        subroutine ApplyDoseFilterToImage(self,img,dose_start,dose_finish,pixel_size)
            use Images
            ! arguments
            class(ElectronDose),    intent(in)      ::  self
            type(Image),            intent(inout)   ::  img
            real,                   intent(in)      ::  dose_start  !<  Number of electrons per square angstroms incident on the specimen at beginning of exposure
            real,                   intent(in)      ::  dose_finish !<  Number of electrons per square angstroms incident on the specimen at end of exposure
            real,                   intent(in)      ::  pixel_size  !<  In Angstroms, the size of the image pixel
            ! private variables
            integer ::  i,j
            real    ::  x,y
            real    ::  current_critical_dose
            real    ::  current_optimal_dose
            real    ::  filter_value
            real,   parameter   ::  critical_dose_at_dc = huge(1.0e0) * 0.001
            ! start work

            !write(*,'(a,2(f0.2,1x))') '**debug(ApplyDoseFilterToImage): dose start, finish: ', dose_start, dose_finish


            ! Loop over image
            if  (img%IsInRealSpace()) then
                do j=1,img%GetLogicalDimension(2)
                    y = ((j-img%GetPhysicalIndexOfBoxCenter(2)) * img%fourier_voxel_size(2))**2
                    do i=1,img%GetLogicalDimension(1)
                        x = ((i-img%GetPhysicalIndexOfBoxCenter(1)) * img%fourier_voxel_size(1))**2
                        !
                        if (i .eq. img%GetPhysicalIndexOfBoxCenter(1) .and. j .eq. img%GetPhysicalIndexOfBoxCenter(2)) then
                            current_critical_dose = critical_dose_at_dc
                        else
                            current_critical_dose = self%CriticalDose(sqrt(x+y)/pixel_size)
                        endif
                        !
                        current_optimal_dose = self%OptimalDoseGivenCriticalDose(current_critical_dose)
                        !
                        if (abs(dose_finish-current_optimal_dose) .lt. abs(dose_start-current_optimal_dose)) then
                            ! We scale the current Fourier component by a value between 0.0 and 1.0,
                            ! where 1.0 can only be obtained if the current frame was the first in the
                            ! movie (i.e. there was no pre-exposure of the specimen).
                            ! The scaling factor is related to the SNR of the image, which in the case of a frame
                            ! from a dose-fractionated exposure (movie) can be computed as the SNR of the movie at the end
                            ! of the frame minus the SNR of the movie at the beginning of the frame
                            !filter_value = (self%SignalToNoiseFromDoseGivenCriticalDose(dose_finish,current_critical_dose) &
                            !               -self%SignalToNoiseFromDoseGivenCriticalDose(dose_start,current_critical_dose)) &
                            !            /   self%SignalToNoiseFromDoseGivenCriticalDose(dose_finish-dose_start,&
                            !                                                            current_critical_dose)
                            !if (filter_value .lt. 0.0e0) then
                            !    call this_program%TerminateWithFatalError('ApplyDoseFilterToImage','Negative filter value')
                            !endif
                            img%real_values(i,j,1) = img%real_values(i,j,1) * DoseFilter(dose_finish,current_critical_dose)
                        else
                            ! If we are beyond the optimal dose, the frame is just adding noise to the current
                            ! Fourier component, therefore we don't let it contribute since otherwise the SNR
                            ! would be made worse
                            img%real_values(i,j,1) = 0.0e0
                        endif
                    enddo
                enddo
            else
                do j=1,img%GetLogicalDimension(2)
                    y = (img%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2) * img%fourier_voxel_size(2))**2
                    do i=1,img%physical_upper_bound_complex(1)
                        x = ((i-1) * img%fourier_voxel_size(1))**2
                        !
                        if (i .eq. 1 .and. j .eq. 1) then
                            current_critical_dose = critical_dose_at_dc
                        else
                            current_critical_dose = self%CriticalDose(sqrt(x+y)/pixel_size)
                        endif
                        !
                        current_optimal_dose = self%OptimalDoseGivenCriticalDose(current_critical_dose)
                        !
                        if (abs(dose_finish-current_optimal_dose) .lt. abs(dose_start-current_optimal_dose)) then
                            ! We scale the current Fourier component by a value between 0.0 and 1.0,
                            ! where 1.0 can only be obtained if the current frame was the first in the
                            ! movie (i.e. there was no pre-exposure of the specimen).
                            ! The scaling factor is related to the SNR of the image, which in the case of a frame
                            ! from a dose-fractionated exposure (movie) can be computed as the SNR of the movie at the end
                            ! of the frame minus the SNR of the movie at the beginning of the frame
                            !filter_value = (self%SignalToNoiseFromDoseGivenCriticalDose(dose_finish,current_critical_dose) &
                            !               -self%SignalToNoiseFromDoseGivenCriticalDose(dose_start,current_critical_dose)) &
                            !            /   self%SignalToNoiseFromDoseGivenCriticalDose(dose_finish-dose_start, &
                            !                                                            current_critical_dose)
                            !if (filter_value .lt. 0.0e0) then
                            !    call this_program%TerminateWithFatalError('ApplyDoseFilterToImage','Negative filter value')
                            !endif
                            img%complex_values(i,j,1) = img%complex_values(i,j,1) * DoseFilter(dose_finish,current_critical_dose)
                        else
                            ! If we are beyond the optimal dose, the frame is just adding noise to the current
                            ! Fourier component, therefore we don't let it contribute since otherwise the SNR
                            ! would be made worse
                            img%complex_values(i,j,1) = 0.0e0
                        endif
                    enddo
                enddo
            endif
        end subroutine ApplyDoseFilterToImage

        !>  \brief  Return the dose filter, into a real image for later multiplication

        subroutine GetRealDoseFilterImageAsIfComplexSize(self,img,dose_start,dose_finish,pixel_size)
            use Images
            ! arguments
            class(ElectronDose),    intent(in)      ::  self
            type(Image),            intent(inout)   ::  img
            real,                   intent(in)      ::  dose_start  !<  Number of electrons per square angstroms incident on the specimen at beginning of exposure
            real,                   intent(in)      ::  dose_finish !<  Number of electrons per square angstroms incident on the specimen at end of exposure
            real,                   intent(in)      ::  pixel_size  !<  In Angstroms, the size of the image pixel
            ! private variables
            integer ::  i,j
            real    ::  x,y
            real    ::  current_critical_dose
            real    ::  current_optimal_dose
            real    ::  filter_value
            real,   parameter   ::  critical_dose_at_dc = huge(1.0e0) * 0.001
            ! start work

            !write(*,'(a,2(f0.2,1x))') '**debug(ApplyDoseFilterToImage): dose start, finish: ', dose_start, dose_finish


            ! Loop over image
            if  (img%IsInRealSpace()) then
                do j=1,img%logical_dimensions(2)
                    y = (img%LogicalIndexGivenPhysicalIndexInFourierSpace(j,2) * img%fourier_voxel_size(2))**2
                    do i=1,img%logical_dimensions(1)
                        x = ((i-1) * (img%fourier_voxel_size(1)*0.5))**2
                        !
                        if (i .eq. 1 .and. j .eq. 1) then
                            current_critical_dose = critical_dose_at_dc
                        else
                            current_critical_dose = self%CriticalDose(sqrt(x+y)/pixel_size)
                        endif
                        !
                        current_optimal_dose = self%OptimalDoseGivenCriticalDose(current_critical_dose)
                        !
                        if (abs(dose_finish-current_optimal_dose) .lt. abs(dose_start-current_optimal_dose)) then
                            ! We scale the current Fourier component by a value between 0.0 and 1.0,
                            ! where 1.0 can only be obtained if the current frame was the first in the
                            ! movie (i.e. there was no pre-exposure of the specimen).
                            ! The scaling factor is related to the SNR of the image, which in the case of a frame
                            ! from a dose-fractionated exposure (movie) can be computed as the SNR of the movie at the end
                            ! of the frame minus the SNR of the movie at the beginning of the frame
                            !filter_value = (self%SignalToNoiseFromDoseGivenCriticalDose(dose_finish,current_critical_dose) &
                            !               -self%SignalToNoiseFromDoseGivenCriticalDose(dose_start,current_critical_dose)) &
                            !            /   self%SignalToNoiseFromDoseGivenCriticalDose(dose_finish-dose_start, &
                            !                                                            current_critical_dose)
                            !if (filter_value .lt. 0.0e0) then
                            !    call this_program%TerminateWithFatalError('ApplyDoseFilterToImage','Negative filter value')
                            !endif
                            img%real_values(i,j,1) = DoseFilter(dose_finish,current_critical_dose)
                        else
                            ! If we are beyond the optimal dose, the frame is just adding noise to the current
                            ! Fourier component, therefore we don't let it contribute since otherwise the SNR
                            ! would be made worse
                            img%real_values(i,j,1) = 0.0e0
                        endif
                    enddo
                enddo

            else

            call this_program%TerminateWithFatalError('ElectronDoses:GetRealDoseFilterImageAsIfComplexSize','Image not real')

            endif

        end subroutine GetRealDoseFilterImageAsIfComplexSize


        !>  \brief  Compute the dose filter, which is the signal attenuation factor due to radiation damage
        !!  See notes of 23-Oct-2014 for almost-complete derivation. This is the same function as Niko (and others, presumably) have used and
        !!  is based on assumption of exponential amplitude decay, as first described by Hayward & Glaeser (1975), I believe.
        pure real function DoseFilter(dose_at_end_of_frame,critical_dose)
            ! arguments
            real,       intent(in)  ::  dose_at_end_of_frame
            real,       intent(in)  ::  critical_dose
            !
            DoseFilter = exp((-0.5 * dose_at_end_of_frame)/critical_dose)
        end function DoseFilter

        !>  \brief  Given a spatial frequency, return the critical dose in electrons per square Angstroms
        pure real function CriticalDose(self,spatial_frequency)
            ! arguments
            class(ElectronDose),    intent(in)  ::  self
            real,                   intent(in)  ::  spatial_frequency   !<  In reciprocal Angstroms
            !
            ! Obtained this by fitting data obtained by Tim Grant from his DLP project
            !CriticalDose = self%critical_dose_scaling * spatial_frequency**self%critical_dose_power
             CriticalDose = (self%critical_dose_a * spatial_frequency**self%critical_dose_b + self%critical_dose_c) &
                            * self%voltage_scaling_factor
            !CriticalDose = ((self%critical_dose_a / (spatial_frequency - self%critical_dose_b)) - self%critical_dose_c) * self%voltage_scaling_factor
        end function CriticalDose

        !>  \brief  Given a number of electrons and a critical dose, return the signal-to-noise ratio
        pure real function SignalToNoiseFromDoseGivenCriticalDose(self,dose,critical_dose) result(snr)
            ! arguments
            class(ElectronDose),    intent(in)  ::  self
            real,                   intent(in)  ::  dose            !<  Number of electrons per unit area
            real,                   intent(in)  ::  critical_dose   !<  In electrons per unit area
            !
            ! From Hayward & Glaeser (1979) via Baker et al (JSB 2010; Eq 6)
            if (dose .eq. 0.0e0) then
                snr = 0.0e0
            else
                snr = (1-exp(-dose * 0.5 / critical_dose))**2/dose
            endif
        end function SignalToNoiseFromDoseGivenCriticalDose

        !>  \brief  Given the critical dose, return an estimate of the optimal dose (at which the SNR is maximised)
        pure real function OptimalDoseGivenCriticalDose(self,critical_dose) result(optimal_dose)
            !
            class(ElectronDose),    intent(in)  ::  self
            real,                   intent(in)  ::  critical_dose
            ! There is actually an analytical solution, found by Wolfram Alpha:
            ! optimal_dose = -critical_dose - 2*critical_dose*W(-1)(-1/(2*sqrt(e)))
            ! where W(k) is the analytic continuation of the product log function
            ! http://mathworld.wolfram.com/LambertW-Function.html
            ! However, there is an acceptable numerical approximation, which I
            ! checked using a spreadsheet and the above formula.
            ! Here, we use the numerical approximation:
            optimal_dose = 2.51284 * critical_dose
        end function OptimalDoseGivenCriticalDose


        !>  \brief  Assignment overload
        subroutine electron_dose_assign(lhs,rhs)
            class(ElectronDose),    intent(inout)   ::  lhs
            class(ElectronDose),    intent(in)      ::  rhs
            !
            lhs%acceleration_voltage  = rhs%acceleration_voltage

            !lhs%critical_dose_scaling = rhs%critical_dose_scaling
            !lhs%critical_dose_power   = rhs%critical_dose_power

            lhs%critical_dose_a = rhs%critical_dose_a
            lhs%critical_dose_b = rhs%critical_dose_b
            lhs%critical_dose_c = rhs%critical_dose_c
            lhs%voltage_scaling_factor = rhs%voltage_scaling_factor


        end subroutine electron_dose_assign


end module ElectronDoses
