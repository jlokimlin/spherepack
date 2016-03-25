module type_GaussianSphere

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stderr => ERROR_UNIT

    use type_Sphere, only: &
        Sphere

    use type_GaussianWorkspace, only: &
        GaussianWorkspace

    use type_GaussianGrid, only: &
        GaussianGrid

    use type_SphericalUnitVectors, only: &
        SphericalUnitVectors

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
    public :: GaussianSphere

    !----------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !----------------------------------------------------------------------
    character (len=250) :: error_message     !! Probably long enough
    integer (ip)        :: allocate_status   !! To check allocation status
    integer (ip)        :: deallocate_status !! To check deallocation status
    !----------------------------------------------------------------------

    ! Declare derived data type
    type, extends (Sphere), public :: GaussianSphere
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure, public  :: create => create_gaussian_sphere
        procedure, public  :: destroy => destroy_gaussian_sphere
        procedure, public  :: get_colatitude_derivative !! Vtsgs
        procedure, public  :: get_vector_laplacian !! Vlapgs
        procedure, public  :: invert_vector_laplacian !! Ivlapgs
        procedure, public  :: get_stream_function_and_velocity_potential
        procedure, public  :: invert_stream_function_and_velocity_potential
        procedure, public  :: perform_grid_transfers
        procedure, public  :: perform_geo_math_coordinate_transfers
        procedure, public  :: perform_scalar_analysis => gaussian_scalar_analysis
        procedure, public  :: perform_scalar_synthesis => gaussian_scalar_synthesis
        procedure, public  :: perform_scalar_projection !! Shpg
        procedure, public  :: perform_vector_analysis => gaussian_vector_analysis
        procedure, public  :: perform_vector_synthesis => gaussian_vector_synthesis
        procedure, public  :: get_Legendre_functions
        procedure, public  :: compute_surface_integral
        final              :: finalize_gaussian_sphere
        !----------------------------------------------------------------------
    end type GaussianSphere


contains


    subroutine create_gaussian_sphere( this, nlat, nlon, isym, itype, isynt, rsphere )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out)        :: this
        integer (ip),           intent (in)            :: nlat
        integer (ip),           intent (in)            :: nlon
        integer (ip),           intent (in), optional  :: isym      !! Either 0, 1, or 2
        integer (ip),           intent (in), optional  :: itype     !! Either 0, 1, 2, 3, ..., 8
        integer (ip),           intent (in), optional  :: isynt
        real (wp),              intent (in), optional  :: rsphere
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: scalar_sym
        integer (ip) :: vector_sym
        integer (ip) :: num_synt
        real (wp)    :: radius
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Allocate polymorphic components
        allocate( GaussianGrid :: this%grid )
        allocate( GaussianWorkspace :: this%workspace )

        ! Initialize polymorphic types
        associate( &
            grid => this%grid, &
            workspace => this%workspace &
            )
            ! Initialize gaussian grid
            select type (grid)
                class is (GaussianGrid)
                call grid%create( nlat, nlon )
            end select
            ! Initialize gaussian workspace
            select type (workspace)
                class is (GaussianWorkspace)
                call workspace%create( nlat, nlon )
            end select
        end associate

        ! Initialize constants
        scalar_sym = 0
        vector_sym = 0
        num_synt = 1
        radius = 1.0_wp

        ! Address optional arguments
        if (present(isym)) scalar_sym = isym
        if (present(itype)) vector_sym = itype
        if (present(isynt)) num_synt = isynt
        if (present(rsphere)) radius = rsphere

        ! Create parent type
        call this%create_sphere( nlat, nlon, scalar_sym, vector_sym, num_synt, radius )

        ! Set initialization flag
        this%initialized = .true.
        
    end subroutine create_gaussian_sphere


    subroutine destroy_gaussian_sphere( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if ( this%initialized .eqv. .false. ) return

        ! Release memory from parent type
        call this%destroy_sphere()

        ! Reset initialization flag
        this%initialized = .false.

    end subroutine destroy_gaussian_sphere



    function compute_surface_integral( this, scalar_function ) result( return_value )
        !
        !< Purpose:
        !
        ! computes the (scalar) surface integral on the sphere (S^2):
        !
        ! * Trapezoidal rule    in phi:   0 <=  phi  <= 2*pi
        ! * Gaussian quadrature in theta: 0 <= theta <= pi
        !
        !   \int_{S^2} f( theta, phi ) dS
        !
        !   where
        !
        !   dS = sin(theta) dtheta dphi
        !
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        real (wp),              intent (in)     :: scalar_function(:,:)
        real (wp)                               :: return_value
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)           :: k  !! counter
        real (wp), allocatable :: summation(:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(GaussianSphere): '&
                //'uninitialized object in GET_SURFACE_INTEGRAL'
        end if

        ! Allocate memory
        associate( nlat => this%NUMBER_OF_LATITUDES )
            allocate(summation(nlat))
        end associate

        ! compute the integrant
        associate( grid => this%grid )
            select type(grid)
                class is (GaussianGrid)
                associate( &
                    nlat => grid%NUMBER_OF_LATITUDES, &
                    dphi => grid%LONGITUDINAL_MESH, &
                    wts => grid%gaussian_weights, &
                    f => scalar_function &
                    )
                    ! Apply trapezoidal rule
                    do k = 1, nlat
                        summation(k) = sum(f(k,:)) * dphi
                    end do
                    ! Apply gaussian quadrature
                    summation = summation * wts
                end associate
            end select
        end associate

        ! Set integral \int_{S^2} f( theta, phi ) dS
        return_value = sum( summation )

        ! Release memory
        deallocate( summation )

    end function compute_surface_integral
    


    subroutine compute_first_moment( this, scalar_function, first_moment )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere),      intent (in out)  :: this
        real (wp),                      intent (in)      :: scalar_function(:,:)
        type (Vector),  intent (out)     :: first_moment
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)            :: k, l !! Counters
        real (wp), allocatable :: integrant(:,:,:)
        !----------------------------------------------------------------------

        !----------------------------------------------------------------------
        ! Check if object is usable
        !----------------------------------------------------------------------

        call this%assert_initialized()

        !----------------------------------------------------------------------
        ! Allocate array
        !----------------------------------------------------------------------

        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            nlon => this%NUMBER_OF_LONGITUDES &
            )

            ! Allocate arrays
            allocate( &
                integrant( nlat, nlon, 3 ), &
                stat=allocate_status &
                )

            ! Check allocation status
            if ( allocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (SpherepackWrapper)'
                write( stderr, '(A)' ) 'Allocation failed in COMPUTE_FIRST_MOMENT'

            end if

        end associate

        !----------------------------------------------------------------------
        ! compute integrant
        !----------------------------------------------------------------------

        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            nlon => this%NUMBER_OF_LONGITUDES &
            )

            ! compute integrant
            do l = 1, nlon
                do k = 1, nlat

                    associate( &
                        u => this%unit_vectors%radial( k, l ), &
                        f => scalar_function(k, l) &
                        )

                        integrant(k, l, 1) = u%x * f
                        integrant(k, l, 2) = u%y * f
                        integrant(k, l, 3) = u%z * f

                    end associate
                end do
            end do
        end associate

        !----------------------------------------------------------------------
        ! compute first moment
        !----------------------------------------------------------------------


        associate( &
            M  => first_moment, &
            f1 => integrant(:,:, 1), &
            f2 => integrant(:,:, 2), &
            f3 => integrant(:,:, 3) &
            )

            m%x = this%compute_surface_integral( f1 )

            m%y = this%compute_surface_integral( f2 )

            m%z = this%compute_surface_integral( f3 )

        end associate

        !----------------------------------------------------------------------
        ! Deallocate array
        !----------------------------------------------------------------------

        deallocate( &
            integrant, &
            stat=deallocate_status &
            )

        ! Check deallocate status
        if ( deallocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (SpherepackWrapper)'
            write( stderr, '(A)' ) 'Deallocation failed in COMPUTE_FIRST_MOMENT'


        end if

    end subroutine compute_first_moment
    !
    
    !
    ! Public SPHEREPACK 3.2 methods
    !
    
    !
    subroutine get_colatitude_derivative( this, polar_component, azimuthal_component )
        !
        !< Purpose:
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#vtsgs.html
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        real (wp),        intent (out)      :: polar_component(:)     !! vt
        real (wp),        intent (out)      :: azimuthal_component(:) !! wt
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: ierror
        !----------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        !        ! TODO incorporate Vtsgsi into type(workspace)
        !        subroutine vtsgsi(nlat, nlon, wvts, lwvts, work, lwork, dwork, ldwork, ierror)

        !        call Vtsgs( &
        !            size(this%grid%latitudes ), size(this%grid%longitudes ),  this%ityp, 1, &
        !            polar_component, azimuthal_component, &
        !            this%NUMBER_OF_LATITUDES, this%NUMBER_OF_LONGITUDES, &
        !            this%workspace%real_polar_harmonic_coefficients, this%workspace%imaginary_polar_harmonic_coefficients, &
        !            this%workspace%real_azimuthal_harmonic_coefficients, this%workspace%imaginary_azimuthal_harmonic_coefficients, &
        !            size(this%workspace%real_polar_harmonic_coefficients, dim=1), size(this%workspace%real_polar_harmonic_coefficients, dim=2), &
        !            this%workspace%wvts, this%workspace%size(wvts), &
        !            this%workspace%legendre_workspace, size(this%workspace%legendre_workspace), ierror)

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine get_colatitude_derivative

!
!    subroutine invert_helmholtz( this, &
!        helmholtz_constant, source_term, solution, perturbation_optional )
!        !----------------------------------------------------------------------
!        ! Dictionary: calling arguments
!        !----------------------------------------------------------------------
!        class (GaussianSphere), intent (in out)        :: this
!        real (wp),        intent (in)            :: helmholtz_constant
!        real (wp),        intent (in)            :: source_term(:,:)
!        real (wp),        intent (out)           :: solution(:,:)
!        real (wp),        intent (out), optional :: perturbation_optional
!        !----------------------------------------------------------------------
!        ! Dictionary: local variables
!        !----------------------------------------------------------------------
!        real (wp):: perturbation
!        integer (ip):: error_flag
!        !----------------------------------------------------------------------
!
!        !----------------------------------------------------------------------
!        ! Check if object is usable
!        !----------------------------------------------------------------------
!
!        call this%assert_initialized()
!
!        !----------------------------------------------------------------------
!        ! Set (real) scalar spherica harmonic coefficients
!        !----------------------------------------------------------------------
!
!        call this%perform_scalar_analysis( source_term )
!
!        !----------------------------------------------------------------------
!        ! Invoke SPHEREPACK 3.2
!        !----------------------------------------------------------------------
!
!        associate( &
!            nlat   => this%NUMBER_OF_LATITUDES, &
!            nlon   => this%NUMBER_OF_LONGITUDES, &
!            isym   => this%SCALAR_SYMMETRIES, &
!            nt     => this%NUMBER_OF_SYNTHESES, &
!            xlmbda => helmholtz_constant, &
!            sf     => solution, &
!            ids    => size(solution, dim=1), &
!            jds    => size(solution, dim=2), &
!            a      => this%workspace%real_harmonic_coefficients, &
!            b      => this%workspace%imaginary_harmonic_coefficients, &
!            mdab   => size(this%workspace%real_harmonic_coefficients, dim=1 ), &
!            ndab   => size(this%workspace%real_harmonic_coefficients, dim=2 ), &
!            wshsgs => this%workspace%backward_scalar, &
!            lshsgs => size(this%workspace%backward_scalar ), &
!            work   => this%workspace%legendre_workspace, &
!            lwork  => size(this%workspace%legendre_workspace ), &
!            pertrb => perturbation, &
!            ierror => error_flag &
!            )
!
!            call Islapgs( nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
!                mdab, ndab, wshsgs, lshsgs, work, lwork, pertrb, ierror )
!
!        end associate
!
!        !----------------------------------------------------------------------
!        ! Address the error flag
!        !----------------------------------------------------------------------
!
!        if ( error_flag == 0 ) then
!
!            return
!
!        else
!
!            write( stderr, '(A)') 'SPHEREPACK 3.2 error: invert_HELMHOLTZ'
!
!            if ( error_flag == 1 ) then
!
!                write( stderr, '(A)') 'Error in the specification of NLAT'
!
!            else if ( error_flag == 2 ) then
!
!                write( stderr, '(A)') 'Error in the specification of NLON'
!
!            else if ( error_flag == 3 ) then
!
!                write( stderr, '(A)') 'Error in the specification of VECTOR_SYMMETRIES'
!
!            else if (error_flag == 4) then
!
!                write( stderr, '(A)') 'Error in the specification of NUMBER_OF_synthESES'
!
!            else if ( error_flag == 5 ) then
!
!                write( stderr, '(A)') 'Invalid extent for SOLUTION'
!                write( stderr, '(A)') 'size(SOLUTION, dim=1 )'
!
!            else if ( error_flag == 6 ) then
!
!                write( stderr, '(A)') 'Invalid extent for SOLUTION'
!                write( stderr, '(A)') 'size(SOLUTION, dim=2 )'
!
!            else if ( error_flag == 7 ) then
!
!                write( stderr, '(A)') 'Invalid extent for '
!                write( stderr, '(A)') 'REAL_HARMONIC_COEFFICIENTS (A) '
!                write( stderr, '(A)') 'or'
!                write( stderr, '(A)') 'IMAGINARY_HARMONIC_COEFFICIENTS (B)'
!                write( stderr, '(A)') 'size(SOLUTION, dim=1 )'
!
!            else if ( error_flag == 8 ) then
!
!                write( stderr, '(A)') 'Invalid extent for '
!                write( stderr, '(A)') 'REAL_HARMONIC_COEFFICIENTS (A) '
!                write( stderr, '(A)') 'or'
!                write( stderr, '(A)') 'IMAGINARY_HARMONIC_COEFFICIENTS (B)'
!                write( stderr, '(A)') 'size(SOLUTION, dim=2 )'
!
!            else if ( error_flag == 9 ) then
!
!                write( stderr, '(A)') 'Invalid extent for WSHSGS'
!                write( stderr, '(A)') 'size(WSHSGS )'
!
!            else if ( error_flag == 10 ) then
!
!                write( stderr, '(A)') 'Invalid extent for WORK'
!                write( stderr, '(A)') 'size(WORK )'
!
!            else
!
!                write( stderr, '(A)') 'Undetermined error flag'
!
!            end if
!        end if
!
!        !----------------------------------------------------------------------
!        ! Address optional arguments
!        !----------------------------------------------------------------------
!
!        if ( present( perturbation_optional ) ) then
!
!            perturbation_optional = perturbation
!
!        end if
!
!    end subroutine invert_helmholtz
    !
    
    ! TODO
    subroutine get_vector_laplacian( this )
        !
        !< Purpose:
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: ierror
        !----------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine get_vector_laplacian
    !
    
    ! TODO
    subroutine invert_vector_laplacian( this )
        !
        !< Purpose:
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: ierror
        !----------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine invert_vector_laplacian
    !
    
    ! TODO
    subroutine get_stream_function_and_velocity_potential( this )
        !
        !< Purpose:
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: ierror
        !----------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine get_stream_function_and_velocity_potential
    !
    
    ! TODO
    subroutine invert_stream_function_and_velocity_potential( this )
        !
        !< Purpose:
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: ierror
        !----------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine invert_stream_function_and_velocity_potential
    !
    
    ! TODO
    subroutine perform_grid_transfers( this )
        !
        !< Purpose:
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: ierror
        !----------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine perform_grid_transfers
    !
    
    ! TODO
    subroutine perform_geo_math_coordinate_transfers( this )
        !
        !< Purpose:
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: ierror
        !----------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine perform_geo_math_coordinate_transfers
    


    subroutine gaussian_scalar_analysis( this, scalar_function )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        real (wp),              intent (in)     :: scalar_function(:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: error_flag
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(GaussianSphere): '&
                //'uninitialized object in GAUSSIAN_SCALAR_ANALYSIS'
        end if

        select type (this)
            class is (GaussianSphere)
            associate( workspace => this%workspace )
                select type (workspace)
                    class is (GaussianWorkspace)
                    ! perform the (real) spherical harmonic analysis
                    associate( &
                        nlat => this%NUMBER_OF_LATITUDES, &
                        nlon => this%NUMBER_OF_LONGITUDES, &
                        isym => this%SCALAR_SYMMETRIES, &
                        nt => this%NUMBER_OF_SYNTHESES, &
                        g => scalar_function, &
                        idg => this%NUMBER_OF_LATITUDES, &
                        jdg => this%NUMBER_OF_LONGITUDES, &
                        a => workspace%real_harmonic_coefficients, &
                        b => workspace%imaginary_harmonic_coefficients, &
                        mdab => this%NUMBER_OF_LATITUDES, &
                        ndab => this%NUMBER_OF_LATITUDES, &
                        wshags => workspace%forward_scalar, &
                        lshags => size(workspace%forward_scalar), &
                        work => workspace%legendre_workspace, &
                        lwork => size(workspace%legendre_workspace), &
                        ierror => error_flag &
                        )
                        call shags( nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                            wshags, lshags, work, lwork, ierror )
                    end associate
                end select
            end associate
        end select

        !----------------------------------------------------------------------
        ! Address the error flag
        !----------------------------------------------------------------------

        if ( error_flag == 0 ) then

            return

        else
            write( stderr, '(A)') 'SPHEREPACK 3.2 error: perform_SCALAR_ANALYSIS'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'


            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Invalid extent for WSHAGS'
                write( stderr, '(A)') 'size(WSHAGS )'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Invalid extent for WORK'
                write( stderr, '(A)') 'size(WORK )'

            else if ( error_flag == 5 ) then

                write( stderr, '(A)') 'Invalid extent for DWORK'
                write( stderr, '(A)') 'size(DWORK )'

            else if ( error_flag == 6 ) then

                write( stderr, '(A)') 'Failure in GAQD to compute gaussian points '
                write( stderr, '(A)') '(due to failure in eigenvalue routine)'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine gaussian_scalar_analysis
    


    subroutine gaussian_scalar_synthesis( this, scalar_function )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        real (wp),              intent (out)    :: scalar_function(:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: error_flag
        !----------------------------------------------------------------------

                ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(GaussianSphere): '&
                //'uninitialized object in GAUSSIAN_SCALAR_ANALYSIS'
        end if

        select type (this)
            class is (GaussianSphere)
            associate( workspace => this%workspace )
                select type (workspace)
                    class is (GaussianWorkspace)
                    ! perform (real) spherical harmonic synthesis
                    associate( &
                        nlat => this%NUMBER_OF_LATITUDES, &
                        nlon => this%NUMBER_OF_LONGITUDES, &
                        isym => this%SCALAR_SYMMETRIES, &
                        nt => this%NUMBER_OF_SYNTHESES, &
                        g => scalar_function, &
                        idg => size(scalar_function, dim=1), &
                        jdg => size(scalar_function, dim=2), &
                        a => workspace%real_harmonic_coefficients, &
                        b => workspace%imaginary_harmonic_coefficients, &
                        mdab => size(workspace%real_harmonic_coefficients, dim=1), &
                        ndab => size(workspace%real_harmonic_coefficients, dim=2), &
                        wshsgs => workspace%backward_scalar, &
                        lshsgs => size(workspace%backward_scalar ), &
                        work => workspace%legendre_workspace, &
                        lwork => size(workspace%legendre_workspace ), &
                        ierror => error_flag &
                        )
                        call shsgs(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                            wshsgs, lshsgs, work, lwork, ierror)
                    end associate
                end select
            end associate
        end select

        !----------------------------------------------------------------------
        ! Address the error flag
        !----------------------------------------------------------------------

        if ( error_flag == 0 ) then

            return

        else

            write( stderr, '(A)') 'SPHEREPACK 3.2 error: perform_SCALAR_synthESIS'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'


            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Invalid extent for WSHSGS'
                write( stderr, '(A)') 'size(WSHSGS )'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Invalid extent for WORK'
                write( stderr, '(A)') 'size(WORK )'

            else if ( error_flag == 5 ) then

                write( stderr, '(A)') 'Invalid extent for DWORK'
                write( stderr, '(A)') 'size(DWORK )'

            else if ( error_flag == 6 ) then

                write( stderr, '(A)') 'Failure in GAQD to compute gaussian points '
                write( stderr, '(A)') '(due to failure in eigenvalue routine)'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine gaussian_scalar_synthesis
    !
    
    ! TODO
    subroutine perform_scalar_projection( this, scalar_function, scalar_projection )
        !
        !< Purpose:
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#shpg.html
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out)        :: this
        real (wp), dimension (:,:), intent (in)  :: scalar_function
        real (wp), dimension (:,:), intent (out) :: scalar_projection
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: ierror
        !----------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ! TODO: Include driver program into type(workspace)
        !        call Shpgi( &
        !            this%NUMBER_OF_LATITUDES, this%NUMBER_OF_LONGITUDES, this%SCALAR_SYMMETRIES, this%TRIANGULAR_TRUNCATION_LIMIT, &
        !            this%workspace%wshp, size(this%workspace%wshp ), &
        !            this%workspace%iwshp, size(this%workspace%iwshp ), &
        !            this%workspace%legendre_workspace, size(this%workspace%legendre_workspace ), ierror )
        !
        !        call Shpg( &
        !            this%NUMBER_OF_LATITUDES, this%NUMBER_OF_LONGITUDES, this%SCALAR_SYMMETRIES, this%TRIANGULAR_TRUNCATION_LIMIT, &
        !            scalar_function, scalar_projection, this%NUMBER_OF_LATITUDES, &
        !            this%workspace%wshp, size(this%workspace%wshp ), &
        !            this%workspace%iwshp, size(this%workspace%iwshp ), &
        !            this%workspace%legendre_workspace, size(this%workspace%legendre_workspace ), ierror )

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine perform_scalar_projection
    !
    
    !
    subroutine gaussian_vector_analysis( this, vector_field )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        real (wp),              intent (in)     :: vector_field(:,:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        real (wp), allocatable :: polar_component(:,:)
        real (wp), allocatable :: azimuthal_component(:,:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(GaussianSphere): '&
                //'uninitialized object in GAUSSIAN_VECTOR_ANALYSIS'
        end if

        ! Allocate memory
        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            nlon => this%NUMBER_OF_LONGITUDES &
            )
            allocate( polar_component(nlat, nlon) )
            allocate( azimuthal_component(nlat, nlon) )
        end associate

        ! compute the spherical angle components
        associate( &
            F => vector_field, &
            v => polar_component, &
            w => azimuthal_component &
            )
            call this%unit_vectors%get_spherical_angle_components( F, v, w )
        end associate

        ! Perform vector analysis
        select type (this)
            class is (GaussianSphere)
            associate( workspace => this%workspace )
                select type (workspace)
                    class is (GaussianWorkspace)
                    associate( &
                        nlat => this%NUMBER_OF_LATITUDES, &
                        nlon => this%NUMBER_OF_LONGITUDES, &
                        ityp => this%VECTOR_SYMMETRIES, &
                        nt => this%NUMBER_OF_SYNTHESES, &
                        v => polar_component, &
                        w => azimuthal_component, &
                        idvw => size(polar_component, dim=1 ), &
                        jdvw => size(polar_component, dim=2 ), &
                        br => workspace%real_polar_harmonic_coefficients, &
                        bi => workspace%imaginary_polar_harmonic_coefficients, &
                        cr => workspace%real_azimuthal_harmonic_coefficients, &
                        ci => workspace%imaginary_azimuthal_harmonic_coefficients, &
                        mdab => size(workspace%real_polar_harmonic_coefficients, dim=1 ), &
                        ndab => size(workspace%real_polar_harmonic_coefficients, dim=2 ), &
                        wvhags => workspace%forward_vector, &
                        lvhags => size(workspace%forward_vector ), &
                        work => workspace%legendre_workspace, &
                        lwork => size(workspace%legendre_workspace ), &
                        ierror => error_flag &
                        )
                        call vhags( nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                            mdab, ndab, wvhags, lvhags, work, lwork, ierror )
                    end associate
                end select
            end associate
        end select

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Error in the specification of NLAT'
            case(2)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Error in the specification of NLON'
            case(3)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Error in the specification of VECTOR_SYMMETRIES'
            case(4)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Error in the specification of NUMBER_OF_synthESES'
            case(5)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Invalid DIM=1 extent for '&
                    //'POLAR_COMPONENT (THETA) or AZIMUTHAL_COMPONENT (PHI)'
            case(6)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Invalid DIM=2 extent '&
                    //'POLAR_COMPONENT (THETA) or AZIMUTHAL_COMPONENT (PHI)'
            case(7)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Invalid DIM=1 extent for BR or CR'
            case(8)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Invalid DIM=1 extent for BI or CI'
            case(9)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Invalid extent for FORWARD_VECTOR'
            case(10)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Invalid extent for LEGENDRE_WORKSPACE'
            case default
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Undetermined error flag'
        end select

        ! Release memory
        deallocate( polar_component)
        deallocate( azimuthal_component)

    end subroutine gaussian_vector_analysis



    subroutine gaussian_vector_synthesis( this, polar_component, azimuthal_component )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        real (wp),              intent (out)    :: polar_component(:,:)
        real (wp),              intent (out)    :: azimuthal_component(:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: error_flag
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(GaussianSphere): '&
                //'uninitialized object in GAUSSIAN_VECTOR_SYNTHESIS'
        end if

        ! Perform vector analysis
        select type (this)
            class is (GaussianSphere)
            associate( workspace => this%workspace )
                select type (workspace)
                    class is (GaussianWorkspace)
                    associate( &
                        nlat => this%NUMBER_OF_LATITUDES, &
                        nlon => this%NUMBER_OF_LONGITUDES, &
                        ityp => this%VECTOR_SYMMETRIES, &
                        nt => this%NUMBER_OF_SYNTHESES, &
                        v => polar_component, &
                        w => azimuthal_component, &
                        idvw => size(polar_component, dim=1 ),  &
                        jdvw => size(polar_component, dim=2 ),  &
                        br => workspace%real_polar_harmonic_coefficients, &
                        bi => workspace%imaginary_polar_harmonic_coefficients, &
                        cr => workspace%real_azimuthal_harmonic_coefficients, &
                        ci => workspace%imaginary_azimuthal_harmonic_coefficients, &
                        mdab => size(workspace%real_polar_harmonic_coefficients, dim=1 ), &
                        ndab => size(workspace%real_polar_harmonic_coefficients, dim=2 ), &
                        wvhsgs => workspace%backward_vector, &
                        lvhsgs => size(workspace%backward_vector ), &
                        work => workspace%legendre_workspace, &
                        lwork => size(workspace%legendre_workspace ), &
                        ierror => error_flag &
                        )
                        call vhsgs( nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                            mdab, ndab, wvhsgs, lvhsgs, work, lwork, ierror )
                    end associate
                end select
            end associate
        end select

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_SYNTHESIS'&
                    //'Error in the specification of NLAT'
            case(2)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_SYNTHESIS'&
                    //'Error in the specification of NLON'
            case(3)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_SYNTHESIS'&
                    //'Error in the specification of VECTOR_SYMMETRIES'
            case(4)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_SYNTHESIS'&
                    //'Error in the specification of NUMBER_OF_synthESES'
            case(5)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_SYNTHESIS'&
                    //'Invalid DIM=1 extent for '&
                    //'POLAR_COMPONENT (THETA) or AZIMUTHAL_COMPONENT (PHI)'
            case(6)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_SYNTHESIS'&
                    //'Invalid DIM=2 extent '&
                    //'POLAR_COMPONENT (THETA) or AZIMUTHAL_COMPONENT (PHI)'
            case(7)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_SYNTHESIS'&
                    //'Invalid DIM=1 extent for BR or CR'
            case(8)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_SYNTHESIS'&
                    //'Invalid DIM=1 extent for BI or CI'
            case(9)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_SYNTHESIS'&
                    //'Invalid extent for BACKWARD_VECTOR'
            case(10)
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_SYNTHESIS'&
                    //'Invalid extent for LEGENDRE_WORKSPACE'
            case default
                error stop 'TYPE(GaussianSphere) in GAUSSIAN_VECTOR_SYNTHESIS'&
                    //'Undetermined error flag'
        end select

    end subroutine gaussian_vector_synthesis
    !
    
    ! TODO
    subroutine get_legendre_functions( this )
        !
        !< Purpose:
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: ierror
        !----------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine get_legendre_functions
    !
    
    ! TODO
    subroutine Icosahedral_geodesic( this )
        !
        !< Purpose:
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: ierror
        !----------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Icosahedral_geodesic
    !
    
    ! TODO
    subroutine perform_multiple_ffts( this )
        !
        !< Purpose:
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: ierror
        !----------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine perform_multiple_ffts


    subroutine finalize_gaussian_sphere( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (GaussianSphere), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_gaussian_sphere


end module type_GaussianSphere
