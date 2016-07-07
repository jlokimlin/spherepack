module type_Workspace

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: Workspace

    ! Declare derived data type
    type, abstract, public :: Workspace
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        logical,                public :: initialized = .false.
        real (wp), allocatable, public :: legendre_workspace(:)
        real (wp), allocatable, public :: forward_scalar(:)
        real (wp), allocatable, public :: forward_vector(:)
        real (wp), allocatable, public :: backward_scalar(:)
        real (wp), allocatable, public :: backward_vector(:)
        real (wp), allocatable, public :: real_harmonic_coefficients(:,:)
        real (wp), allocatable, public :: imaginary_harmonic_coefficients(:,:)
        real (wp), allocatable, public :: real_polar_harmonic_coefficients(:,:)
        real (wp), allocatable, public :: imaginary_polar_harmonic_coefficients(:,:)
        real (wp), allocatable, public :: real_azimuthal_harmonic_coefficients(:,:)
        real (wp), allocatable, public :: imaginary_azimuthal_harmonic_coefficients(:,:)
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure, public :: destroy_workspace
        procedure, public :: copy_workspace
        !----------------------------------------------------------------------
    end type Workspace


contains



    subroutine copy_workspace(this, object_to_be_copied)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (Workspace), intent (out) :: this
        class (Workspace), intent (in)  :: object_to_be_copied
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if (.not.object_to_be_copied%initialized) then
            error stop 'Uninitialized object of class (Workspace): '&
                //'in assignment (=) '
        end if

        !
        !==> Make copies
        !
        this%initialized = object_to_be_copied%initialized
        this%legendre_workspace = object_to_be_copied%legendre_workspace
        this%forward_scalar = object_to_be_copied%forward_scalar
        this%forward_vector = object_to_be_copied%forward_vector
        this%backward_scalar = object_to_be_copied%backward_scalar
        this%backward_vector = object_to_be_copied%backward_vector
        this%real_harmonic_coefficients = object_to_be_copied%real_harmonic_coefficients
        this%imaginary_harmonic_coefficients = object_to_be_copied%imaginary_harmonic_coefficients
        this%real_polar_harmonic_coefficients = object_to_be_copied%real_polar_harmonic_coefficients
        this%imaginary_polar_harmonic_coefficients = object_to_be_copied%imaginary_polar_harmonic_coefficients
        this%real_azimuthal_harmonic_coefficients = object_to_be_copied%real_azimuthal_harmonic_coefficients
        this%imaginary_azimuthal_harmonic_coefficients = object_to_be_copied%imaginary_azimuthal_harmonic_coefficients

    end subroutine copy_workspace



    subroutine destroy_workspace(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Workspace), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if (.not.this%initialized) return

        !
        !==> Release memory
        !
        if (allocated(this%legendre_workspace)) then
            deallocate(this%legendre_workspace)
        end if

        if (allocated(this%forward_scalar)) then
            deallocate(this%forward_scalar)
        end if

        if (allocated(this%forward_vector)) then
            deallocate(this%forward_vector)
        end if

        if (allocated(this%backward_scalar)) then
            deallocate(this%backward_scalar)
        end if

        if (allocated(this%backward_vector)) then
            deallocate(this%backward_vector)
        end if

        if (allocated(this%real_harmonic_coefficients)) then
            deallocate(this%real_harmonic_coefficients)
        end if

        if (allocated(this%imaginary_harmonic_coefficients)) then
            deallocate(this%imaginary_harmonic_coefficients)
        end if

        if (allocated(this%real_polar_harmonic_coefficients)) then
            deallocate(this%real_polar_harmonic_coefficients)
        end if

        if (allocated(this%imaginary_polar_harmonic_coefficients)) then
            deallocate(this%imaginary_polar_harmonic_coefficients)
        end if

        if (allocated(this%real_azimuthal_harmonic_coefficients)) then
            deallocate(this%real_azimuthal_harmonic_coefficients)
        end if

        if (allocated(this%imaginary_azimuthal_harmonic_coefficients)) then
            deallocate(this%imaginary_azimuthal_harmonic_coefficients)
        end if

        ! Reset flag
        this%initialized = .false.

    end subroutine destroy_workspace


end module type_Workspace
