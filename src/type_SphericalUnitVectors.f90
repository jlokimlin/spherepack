module type_SphericalUnitVectors

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_Grid, only: &
        SphericalGrid

    use type_TrigonometricFunctions, only: &
        TrigFunctions => TrigonometricFunctions

    use type_ThreeDimensionalVector, only: &
        Vector => ThreeDimensionalVector, &
        assignment(=), &
        operator(*)

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: SphericalUnitVectors

    ! Declare derived data type
    type, public :: SphericalUnitVectors
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        logical,                    public :: initialized = .false.
        integer (ip),               public :: NUMBER_OF_LONGITUDES = 0
        integer (ip),               public :: NUMBER_OF_LATITUDES = 0
        type (Vector), allocatable, public :: radial(:, :)
        type (Vector), allocatable, public :: polar(:, :)
        type (Vector), allocatable, public :: azimuthal(:, :)
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure, public :: create => create_spherical_unit_vectors
        procedure, public :: destroy => destroy_spherical_unit_vectors
        procedure, public :: get_spherical_angle_components
        final             :: finalize_spherical_unit_vectors
        !----------------------------------------------------------------------
    end type SphericalUnitVectors


contains


    subroutine create_spherical_unit_vectors( this, grid_type, trig_functions )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalUnitVectors),            intent (in out) :: this
        class (SphericalGrid),                            intent (in out) :: grid_type
        class (TrigFunctions), optional, target, intent (in out) :: trig_functions
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)                   :: k,  l !! Counters
        type (TrigFunctions), pointer :: ptr => null()
        !----------------------------------------------------------------------

        ! Check if object is usable
        call this%destroy()

        ! Check if polymorphic argument is usable
        if ( grid_type%initialized .eqv. .false. ) then
            error stop 'TYPE(SphericalUnitVectors): '&
                //'initialized polymorphic argument CLASS(Grid)'
        end if

        ! Associate array extends
        associate( &
            nlat => grid_type%NUMBER_OF_LATITUDES, &
            nlon => grid_type%NUMBER_OF_LONGITUDES &
            )

            this%NUMBER_OF_LATITUDES = nlat
            this%NUMBER_OF_LONGITUDES = nlon

            ! Allocate memory
            allocate( this%radial(nlat, nlon) )
            allocate( this%polar(nlat, nlon) )
            allocate( this%azimuthal(nlat, nlon) )

            ! Address optional argument
            if (present(trig_functions)) then

                ! Check if polymorphic argument is usable
                if ( trig_functions%initialized .eqv. .false. ) then
                    error stop 'TYPE(SphericalUnitVectors): '&
                        //'initialized polymorphic argument CLASS(Grid)'
                else
                    ! Assign pointer
                    ptr => trig_functions
                end if
            else
                ! Compute trigonometric functions from scratch
                allocate( TrigFunctions :: ptr )
                call ptr%create( grid_type )
            end if

            ! Compute spherical unit vectors
            associate( &
                r => this%radial, &
                theta => this%polar, &
                phi => this%azimuthal, &
                sint => ptr%sint, &
                cost => ptr%cost, &
                sinp => ptr%sinp, &
                cosp => ptr%cosp &
                )
                do l = 1, nlon
                    do k = 1, nlat

                        ! set radial unit vector
                        r(k, l) = &
                            Vector( &
                            x = sint(k) * cosp(l), &
                            y = sint(k) * sinp(l), &
                            z = cost(k) &
                            )

                        ! set polar unit vector
                        theta(k, l) = &
                            Vector( &
                            x = cost(k) * cosp(l), &
                            y = cost(k) * sinp(l), &
                            z = -sint(k) &
                            )

                        ! set azimuthal unit vector
                        phi(k, l) = &
                            Vector( &
                            x = -sinp(l), &
                            y =  cosp(l), &
                            z =  0.0_wp &
                            )
                    end do
                end do
            end associate
        end associate

        ! Garbage collection
        if (present(trig_functions)) then
            nullify( ptr )
        else
            deallocate( ptr )
        end if

        ! Set flag
        this%initialized = .true.

    end subroutine create_spherical_unit_vectors
    


    subroutine destroy_spherical_unit_vectors( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalUnitVectors), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if ( this%initialized .eqv. .false. ) return

        ! Release memory
        if ( allocated( this%radial ) ) deallocate ( this%radial )
        if (allocated(this%polar)) deallocate(  this%polar )
        if (allocated(this%azimuthal)) deallocate( this%azimuthal )

        ! Reset constants
        this%NUMBER_OF_LONGITUDES = 0
        this%NUMBER_OF_LATITUDES = 0

        ! Reset flag
        this%initialized = .false.

    end subroutine destroy_spherical_unit_vectors



    subroutine get_spherical_angle_components( this, &
        vector_function, polar_component, azimuthal_component )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalUnitVectors), intent (in out) :: this
        real (wp),                    intent (in)     :: vector_function(:, :, :)
        real (wp),                    intent (out)    :: polar_component(:, :)
        real (wp),                    intent (out)    :: azimuthal_component(:, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)  :: k, l !! Counters
        type (Vector) :: vector_field !! To cast array to vector
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(SphericalUnitVectors): '&
                //'uninitialized object in GET_SPHERICAL_ANGLE_COMPONENTS'
        end if

        ! Calculate the spherical angle components
        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            nlon => this%NUMBER_OF_LONGITUDES &
            )
            do l = 1, nlon
                do k = 1, nlat
                    ! Cast array to vector
                    vector_field = vector_function(:, k, l)
                    associate( &
                        theta => this%polar(k, l), &
                        phi => this%azimuthal(k, l) &
                        )
                        ! set the theta component
                        polar_component( k, l ) = theta.dot.vector_field
                        ! set the azimuthal_component
                        azimuthal_component( k, l ) = phi.dot.vector_field
                    end associate
                end do
            end do
        end associate

    end subroutine get_spherical_angle_components


    subroutine finalize_spherical_unit_vectors( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (SphericalUnitVectors), intent (in out)    :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_spherical_unit_vectors


end module type_SphericalUnitVectors
