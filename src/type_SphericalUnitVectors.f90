module type_SphericalUnitVectors

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_SphericalGrid, only: &
        SphericalGrid

    use type_TrigonometricFunctions, only: &
        TrigonometricFunctions

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
        type (Vector), allocatable, public :: radial(:,:)
        type (Vector), allocatable, public :: polar(:,:)
        type (Vector), allocatable, public :: azimuthal(:,:)
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



    ! Declare constructor
    interface SphericalUnitVectors
        module procedure spherical_unit_vectors_constructor
    end interface



contains



    function spherical_unit_vectors_constructor(grid) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalGrid), intent (in out) :: grid
        type (SphericalUnitVectors)            :: return_value
        !----------------------------------------------------------------------

        call return_value%create(grid)

    end function spherical_unit_vectors_constructor



    subroutine create_spherical_unit_vectors(this, grid)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalUnitVectors),  intent (in out) :: this
        class (SphericalGrid),         intent (in out) :: grid
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)                  :: k,  l !! Counters
        type (TrigonometricFunctions) :: trig_func
        !----------------------------------------------------------------------

        ! Check if object is usable
        call this%destroy()

        ! Check if polymorphic argument is usable
        if ( grid%initialized .eqv. .false.) then
            error stop 'Object of class (SphericalUnitVectors): '&
                //'uninitialized polymorphic argument of class (SphericalGrid) '&
                //'in create_spherical_unit_vectors'
        end if


        associate( &
            nlat => grid%NUMBER_OF_LATITUDES, &
            nlon => grid%NUMBER_OF_LONGITUDES &
            )

            ! Set constants
            this%NUMBER_OF_LATITUDES = nlat
            this%NUMBER_OF_LONGITUDES = nlon

            !
            !==> Allocate memory
            !
            allocate(this%radial(nlat, nlon) )
            allocate(this%polar(nlat, nlon) )
            allocate(this%azimuthal(nlat, nlon) )

            ! Compute required trigonometric functions
            trig_func = TrigonometricFunctions(grid)

            ! Compute spherical unit vectors
            associate( &
                r => this%radial, &
                theta => this%polar, &
                phi => this%azimuthal, &
                sint => trig_func%sint, &
                cost => trig_func%cost, &
                sinp => trig_func%sinp, &
                cosp => trig_func%cosp &
                )

                do l = 1, nlon
                    do k = 1, nlat

                        ! set radial unit vector
                        r(k, l) = &
                            Vector(&
                            sint(k) * cosp(l), &
                            sint(k) * sinp(l), &
                            cost(k) &
                            )

                        ! set polar unit vector
                        theta(k, l) = &
                            Vector( &
                            cost(k) * cosp(l), &
                            cost(k) * sinp(l), &
                            -sint(k) &
                            )

                        ! set azimuthal unit vector
                        phi(k, l) = &
                            Vector( &
                            -sinp(l), &
                            cosp(l), &
                            0.0_wp &
                            )
                    end do
                end do
            end associate
        end associate

        ! Set flag
        this%initialized = .true.

    end subroutine create_spherical_unit_vectors
    


    subroutine destroy_spherical_unit_vectors(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalUnitVectors), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if (.not.this%initialized) return

        ! Release memory
        if (allocated(this%radial)) deallocate( this%radial )
        if (allocated(this%polar)) deallocate(  this%polar )
        if (allocated(this%azimuthal)) deallocate( this%azimuthal )

        ! Reset constants
        this%NUMBER_OF_LONGITUDES = 0
        this%NUMBER_OF_LATITUDES = 0

        ! Reset flag
        this%initialized = .false.

    end subroutine destroy_spherical_unit_vectors



    subroutine get_spherical_angle_components(this, &
        vector_function, polar_component, azimuthal_component )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalUnitVectors), intent (in out) :: this
        real (wp),                    intent (in)     :: vector_function(:,:,:)
        real (wp),                    intent (out)    :: polar_component(:,:)
        real (wp),                    intent (out)    :: azimuthal_component(:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)  :: k, l !! Counters
        type (Vector) :: vector_field !! To cast array to vector
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.this%initialized) then
            error stop 'Object of class (SphericalUnitVectors): '&
                //'uninitialized object in get_spherical_angle_components'
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
                        polar_component(k, l) = theta.dot.vector_field
                        ! set the azimuthal_component
                        azimuthal_component(k, l) = phi.dot.vector_field
                    end associate
                end do
            end do
        end associate

    end subroutine get_spherical_angle_components




    subroutine get_vector_function(this, &
        radial_component, polar_component, azimuthal_component, vector_function)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalUnitVectors), intent (in out) :: this
        real (wp),                    intent (in)     :: radial_component(:,:)
        real (wp),                    intent (in)     :: polar_component(:,:)
        real (wp),                    intent (in)     :: azimuthal_component(:,:)
        real (wp),                    intent (out)    :: vector_function(:,:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)  :: k, l !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.this%initialized) then
            error stop 'TYPE(SphericalUnitVectors): '&
                //'uninitialized object in GET_VECTOR_FUNCTION'
        end if


        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            nlon => this%NUMBER_OF_LONGITUDES &
            )
            do l = 1, nlon
                do k = 1, nlat
                    associate( &
                        r => this%radial(k,l), &
                        theta => this%polar(k, l), &
                        phi => this%azimuthal(k, l) &
                        )
                        !
                        !==> Calculate the spherical angle components
                        !
                        vector_function(:, k, l ) = &
                            r * radial_component(k, l) &
                            + theta * polar_component(k, l) &
                            + phi * azimuthal_component(k,l)

                    end associate
                end do
            end do
        end associate

    end subroutine get_vector_function



    subroutine finalize_spherical_unit_vectors(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (SphericalUnitVectors), intent (in out)    :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_spherical_unit_vectors



end module type_SphericalUnitVectors
