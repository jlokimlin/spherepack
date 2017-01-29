module type_Workspace

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: Workspace

    type, public, abstract :: Workspace
        ! Type components
        logical,               public :: initialized = .false.
        real(wp), allocatable, public :: legendre_workspace(:)
        real(wp), allocatable, public :: forward_scalar(:)
        real(wp), allocatable, public :: forward_vector(:)
        real(wp), allocatable, public :: backward_scalar(:)
        real(wp), allocatable, public :: backward_vector(:)
        real(wp), allocatable, public :: real_harmonic_coefficients(:,:)
        real(wp), allocatable, public :: imaginary_harmonic_coefficients(:,:)
        real(wp), allocatable, public :: real_polar_harmonic_coefficients(:,:)
        real(wp), allocatable, public :: imaginary_polar_harmonic_coefficients(:,:)
        real(wp), allocatable, public :: real_azimuthal_harmonic_coefficients(:,:)
        real(wp), allocatable, public :: imaginary_azimuthal_harmonic_coefficients(:,:)
    contains
        ! Type-bound procedures
        procedure, public :: destroy_workspace
        procedure, public :: copy_workspace
    end type Workspace

contains

    subroutine copy_workspace(self, other)

        ! Dummy arguments
        class(Workspace), intent(out) :: self
        class(Workspace), intent(in)  :: other

        ! Check if object is usable
        if (.not.other%initialized) then
            error stop 'Uninitialized object of class(Workspace): '&
                //'in assignment (=) '
        end if

        !  Make copies
        self%initialized = other%initialized
        self%legendre_workspace = other%legendre_workspace
        self%forward_scalar = other%forward_scalar
        self%forward_vector = other%forward_vector
        self%backward_scalar = other%backward_scalar
        self%backward_vector = other%backward_vector
        self%real_harmonic_coefficients = other%real_harmonic_coefficients
        self%imaginary_harmonic_coefficients = other%imaginary_harmonic_coefficients
        self%real_polar_harmonic_coefficients = other%real_polar_harmonic_coefficients
        self%imaginary_polar_harmonic_coefficients = other%imaginary_polar_harmonic_coefficients
        self%real_azimuthal_harmonic_coefficients = other%real_azimuthal_harmonic_coefficients
        self%imaginary_azimuthal_harmonic_coefficients = other%imaginary_azimuthal_harmonic_coefficients

    end subroutine copy_workspace

    subroutine destroy_workspace(self)

        ! Dummy arguments
        class(Workspace), intent(inout) :: self

        ! Check flag
        if (.not.self%initialized) return

        !  Release memory
        if (allocated(self%legendre_workspace)) then
            deallocate(self%legendre_workspace)
        end if

        if (allocated(self%forward_scalar)) then
            deallocate(self%forward_scalar)
        end if

        if (allocated(self%forward_vector)) then
            deallocate(self%forward_vector)
        end if

        if (allocated(self%backward_scalar)) then
            deallocate(self%backward_scalar)
        end if

        if (allocated(self%backward_vector)) then
            deallocate(self%backward_vector)
        end if

        if (allocated(self%real_harmonic_coefficients)) then
            deallocate(self%real_harmonic_coefficients)
        end if

        if (allocated(self%imaginary_harmonic_coefficients)) then
            deallocate(self%imaginary_harmonic_coefficients)
        end if

        if (allocated(self%real_polar_harmonic_coefficients)) then
            deallocate(self%real_polar_harmonic_coefficients)
        end if

        if (allocated(self%imaginary_polar_harmonic_coefficients)) then
            deallocate(self%imaginary_polar_harmonic_coefficients)
        end if

        if (allocated(self%real_azimuthal_harmonic_coefficients)) then
            deallocate(self%real_azimuthal_harmonic_coefficients)
        end if

        if (allocated(self%imaginary_azimuthal_harmonic_coefficients)) then
            deallocate(self%imaginary_azimuthal_harmonic_coefficients)
        end if

        ! Reset flag
        self%initialized = .false.

    end subroutine destroy_workspace

end module type_Workspace
