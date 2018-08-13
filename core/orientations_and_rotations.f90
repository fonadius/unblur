!>  \brief  Set of classes and methods to deal with orientations and rotations in three dimensions
module OrientationsAndRotations
    use Units
    implicit none
    private

    !>  \brief  Euler angles, following the Frealign / Spider convention
    type, public ::  EulerAngles
        private
        integer ::  angle_units =   degrees
        real    ::  phi         =   0.0     !<  Phi angle (Spider/Frealign convention) - twist  (equiv. to Imagic's Gamma)
        real    ::  theta       =   0.0     !<  Theta angle (Spider/Frealign convention) - out-of-plane (equiv. to Imagic's Beta)
        real    ::  psi         =   0.0     !<  Psi angle, (Spider/Frealign convention) - in-plane (equiv. to Imagic's Alpha)
        contains
            procedure,  public  ::  Init
            procedure,  public  ::  SetAnglesInDegrees
            procedure           ::  ConvertToRadians
            procedure,  public  ::  GetRotationMatrix
    end type

    contains

    !>  \brief  Initialise a set of Euler angles
    pure subroutine Init(self)
        class(EulerAngles), intent(inout)   ::  self
        call self%SetAnglesInDegrees(0.0,0.0,0.0)
    end subroutine Init

    !>  \brief  Set Euler angles to given values in degrees
    pure subroutine SetAnglesInDegrees(self,phi,theta,psi)
        class(EulerAngles), intent(inout)   ::  self
        real,               intent(in)      ::  phi,theta,psi
        self%phi            =   phi
        self%theta          =   theta
        self%psi            =   psi
        self%angle_units    =   degrees
    end subroutine SetAnglesInDegrees

    !>  \brief  Convert Euler angles to radians
    pure subroutine ConvertToRadians(self)
        class(EulerAngles), intent(inout)   ::  self
        if (self%angle_units .ne. radians) then
            self%phi            =   convert(self%phi,self%angle_units,radians)
            self%theta          =   convert(self%theta,self%angle_units,radians)
            self%psi            =   convert(self%psi,self%angle_units,radians)
            self%angle_units    =   radians
        endif
    end subroutine ConvertToRadians

    !>  \brief  Return a 4x4 rotation matrix, given Euler angles following the Spider / Frealign convention,
    !!          as decribed in http://www.wadsworth.org/spider_doc/spider/docs/euler.html
    !!
    !!
    !!      By Niko, on the lab forum:
    !!      "Here are the Frealign conventions. The directions i use here assume that the particle image origin is in the lower left corner (as is the case in the program Ximdisp).
    !!      1) In-plane rotation (psi angle in Frealign): a positive value describes a counterclockwise rotation of the particle. To align it with the reference, the negative value of this angle needs to be applied, i.e. a clockwise rotation.
    !!      2) X,y shifts: Positive values describe shifts to the right and up. The negative values have to be applied to align the particle with a reference.
    !!      3) Frealign applies the shifts first and then the rotations. There are no mirrors used in Frealign."
    !!
    !!      The reverse flag can be used to obtain the reverse transformation described by the Euler angles. In other words,
    !!      multiplying the matrix returned when using REVERSE=.FALSE. by the matrix returned when using REVERSE=.TRUE. should
    !!      give the identity matrix.
    !!
    !!
    pure function GetRotationMatrix(self,shift,reverse) result(rotmat)
        ! Arguments
        class(EulerAngles),     intent(in)      ::  self
        real,       optional,   intent(in)      ::  shift(3)    !<  Shift (applied after the rotation in forward direction)
        logical,    optional,   intent(in)      ::  reverse     !<  Whether the inverse matrix should be returned. See subroutine documentation for more detail.
        ! Result
        real(kind=8)    ::  rotmat(4,4)
        ! Private variables
        real    ::  pphi,ttheta,ppsi
        real    ::  rotmat_phi(4,4), rotmat_theta(4,4), rotmat_psi(4,4)
        logical ::  rreverse
        integer ::  i
        ! Start work

        ! Initialise resut matrix (to the identity matrix)
        rotmat = 0.0
        do i=1,4
            rotmat(i,i) = 1.0
        enddo

        ! Check whether we'll be reversing (by default, we don't)
        rreverse = .false.
        if (present(reverse)) rreverse = reverse

        ! Let's make sure we're in radians
        pphi    = convert(self%phi,  self%angle_units,radians)
        ttheta  = convert(self%theta,self%angle_units,radians)
        ppsi    = convert(self%psi,  self%angle_units,radians)


        ! Reverse the angles  if this was asked.
        !
        ! Note the way frealign thinks about Euler angles:
        ! To it, you apply the Euler angles to the 3d, then shift the projection to obtain the original image,
        ! whereas to Imagic (and Spider i think), the Euler angles are applied after the image is shifted to
        ! back-project it into the 3D.
        if (rreverse) then
            pphi = -pphi
            ttheta = -ttheta
            ppsi = -ppsi
        endif

        ! Set the matrix elements.
        rotmat_psi(1,:)   = (/ cos(ppsi),     sin(ppsi), 0.0,         0.0/)
        rotmat_psi(2,:)   = (/-sin(ppsi),     cos(ppsi), 0.0,         0.0/)
        rotmat_psi(3,:)   = (/ 0.0,           0.0,       1.0,         0.0/)
        rotmat_psi(4,:)   = (/ 0.0,           0.0,       0.0,         1.0/)

        rotmat_theta(1,:) = (/ cos(ttheta),   0.0,      -sin(ttheta), 0.0/)
        rotmat_theta(2,:) = (/ 0.0,           1.0,       0.0,         0.0/)
        rotmat_theta(3,:) = (/ sin(ttheta),   0.0,       cos(ttheta), 0.0/)
        rotmat_theta(4,:) = (/ 0.0,           0.0,       0.0,         1.0/)

        rotmat_phi(1,:)   = (/ cos(pphi),     sin(pphi), 0.0,         0.0/)
        rotmat_phi(2,:)   = (/-sin(pphi),     cos(pphi), 0.0,         0.0/)
        rotmat_phi(3,:)   = (/ 0.0,           0.0,       1.0,         0.0/)
        rotmat_phi(4,:)   = (/ 0.0,           0.0,       0.0,         1.0/)

        if (rreverse) then
            if (present(shift)) then
                rotmat_psi(1:3,4) = -shift(1:3)
            endif
            rotmat = matmul(rotmat_theta,rotmat_psi)
            rotmat = matmul(rotmat_phi,rotmat)
        else
            rotmat = matmul(rotmat_theta,rotmat_phi)
            rotmat = matmul(rotmat_psi,rotmat)
            if (present(shift)) then
                rotmat(1:3,4) = shift(1:3)
            endif
        endif
    end function GetRotationMatrix



end module OrientationsAndRotations
