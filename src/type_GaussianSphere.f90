module type_GaussianSphere

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

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

    use type_Vector3D, only: &
        Vector => Vector3D, &
        assignment(=), &
        operator(*)
    
    use scalar_analysis_routines, only: &
        ShagsAux

    use scalar_synthesis_routines, only: &
        ShsgsAux

    use module_vhags, only: &
        VhagsAux

    use module_vhsgs, only: &
        VhsgsAux

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: GaussianSphere


    ! Declare derived data type
    type, extends (Sphere), public :: GaussianSphere
        !----------------------------------------------------------------------
        ! Type components
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Type-bound procedures
        !----------------------------------------------------------------------
        procedure, public  :: create => create_gaussian_sphere
        procedure, public  :: destroy => destroy_gaussian_sphere
        procedure, public  :: perform_scalar_analysis => gaussian_scalar_analysis
        procedure, public  :: perform_scalar_synthesis => gaussian_scalar_synthesis
        procedure, public  :: vector_analysis_from_spherical_components => gaussian_vector_analysis
        procedure, public  :: perform_vector_synthesis => gaussian_vector_synthesis
        procedure, public  :: compute_surface_integral
        procedure, public  :: compute_first_moment
        final              :: finalize_gaussian_sphere
        !----------------------------------------------------------------------
    end type GaussianSphere


    ! Declare constructor
    interface GaussianSphere
        module procedure gaussian_sphere_constructor
    end interface



contains



    function gaussian_sphere_constructor(nlat, nlon) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in) :: nlat !! number of latitudinal points 0 <= theta <= pi
        integer (ip), intent (in) :: nlon !! number of longitudinal points 0 <= phi <= 2*pi
        type (GaussianSphere)     :: return_value
        !----------------------------------------------------------------------

        call return_value%create(nlat, nlon)

    end function gaussian_sphere_constructor



    subroutine create_gaussian_sphere(this, nlat, nlon, ntrunc, isym, itype, nt, rsphere)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out)        :: this
        integer (ip),           intent (in)            :: nlat
        integer (ip),           intent (in)            :: nlon
        integer (ip),           intent (in), optional  :: ntrunc
        integer (ip),           intent (in), optional  :: isym      !! Either 0, 1, or 2
        integer (ip),           intent (in), optional  :: itype     !! Either 0, 1, 2, 3, ..., 8
        integer (ip),           intent (in), optional  :: nt
        real (wp),              intent (in), optional  :: rsphere
        !--------------------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: ntrunc_op
        integer (ip) :: isym_op
        integer (ip) :: ityp_op
        integer (ip) :: nt_op
        real (wp)    :: rsphere_op
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        !
        !==> Allocate polymorphic components
        !
        allocate( this%grid, source=GaussianGrid(nlat,nlon) )
        allocate( this%workspace, source=GaussianWorkspace(nlat, nlon) )

        !
        !==> Address optional arguments
        !
        if (present(ntrunc)) then
            ntrunc_op = ntrunc
        else
            ntrunc_op = nlat - 1
        end if

        if (present(isym)) then
            isym_op = isym
        else
            isym_op = 0
        end if

        if (present(itype)) then
            ityp_op = itype
        else
            ityp_op = 0
        end if

        if (present(nt)) then
            nt_op = nt
        else
            nt_op = 1
        end if

        if (present(rsphere)) then
            rsphere_op = rsphere
        else
            rsphere_op = 1.0_wp
        end if

        !
        !==> Create parent type
        !
        call this%create_sphere(nlat, nlon, ntrunc_op, isym_op, ityp_op, nt_op, rsphere_op)

        ! Set flag
        this%initialized = .true.
        
    end subroutine create_gaussian_sphere


    subroutine destroy_gaussian_sphere(this)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if (.not.this%initialized) return

        ! Release memory from parent type
        call this%destroy_sphere()

        ! Reset initialization flag
        this%initialized = .false.

    end subroutine destroy_gaussian_sphere



    subroutine gaussian_scalar_analysis(this, scalar_function)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        real (wp),              intent (in)     :: scalar_function(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)    :: error_flag
        type (ShagsAux) :: aux
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.this%initialized) then
            error stop 'Uninitialized object of class(GaussianSphere): '&
                //'in gaussian_scalar_analysis'
        end if


        select type (this)
            class is (GaussianSphere)
            associate( workspace => this%workspace )
                select type (workspace)
                    class is (GaussianWorkspace)
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
                        !
                        !==> perform the (real) spherical harmonic scalar analysis
                        !
                        call aux%shags(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                            wshags, lshags, work, lwork, ierror)

                    end associate
                end select
            end associate
        end select

        ! Address error flag
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_analysis'&
                    //'invalid specification of NUMBER_OF_LATITUDES'
            case (2)
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_analysis'&
                    //'invalid specification of NUMBER_OF_LONGITUDES'
            case (3)
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_analysis'&
                    //'invalid specification for scalar_forward'
            case (4)
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_analysis'&
                    //'invalid specification for legendre_workspace'
            case (5)
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_analysis'&
                    //'invalid specification for dwork'
            case (6)
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_analysis'&
                    //' failure in compute_gaussian_latitudes_and_weights due to failure in eigenvalue routine'
            case default
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_analysis'&
                    //' undetermined error'
        end select

    end subroutine gaussian_scalar_analysis
    


    subroutine gaussian_scalar_synthesis(this, scalar_function)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        real (wp),              intent (out)    :: scalar_function(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)    :: error_flag
        type (ShsgsAux) :: aux
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.this%initialized) then
            error stop 'Uninitialized object of class(GaussianSphere): '&
                //'in gaussian_scalar_synthesis'
        end if

        select type (this)
            class is (GaussianSphere)
            associate( workspace => this%workspace )
                select type (workspace)
                    class is (GaussianWorkspace)
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
                        lshsgs => size(workspace%backward_scalar), &
                        work => workspace%legendre_workspace, &
                        lwork => size(workspace%legendre_workspace), &
                        ierror => error_flag &
                        )
                        !
                        !==> Perform gaussian scalar synthesis
                        !
                        call aux%shsgs(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                            wshsgs, lshsgs, work, lwork, ierror)

                    end associate
                end select
            end associate
        end select

        !
        !==> Address error flag
        !
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_synthesis'&
                    //'invalid specification of NUMBER_OF_LATITUDES'
            case (2)
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_synthesis'&
                    //'invalid specification of NUMBER_OF_LONGITUDES'
            case (3)
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_synthesis'&
                    //'invalid specification for scalar_backward'
            case (4)
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_synthesis'&
                    //'invalid specification for legendre_workspace'
            case (5)
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_synthesis'&
                    //'invalid specification for dwork'
            case (6)
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_synthesis'&
                    //'Failure in compute_gaussian_latitudes_and_weights due to failure in eigenvalue routine'
            case default
                error stop 'Object of class(GaussianSphere) in gaussian_scalar_synthesis'&
                    //'Undetermined error flag'
        end select


    end subroutine gaussian_scalar_synthesis



    subroutine gaussian_vector_analysis(this, polar_component, azimuthal_component)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        real (wp),              intent (in)     :: polar_component(:,:)
        real (wp),              intent (in)     :: azimuthal_component(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)    :: error_flag
        type (VhagsAux) :: aux
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.this%initialized) then
            error stop 'Uninitialized object of class(GaussianSphere): '&
                //'in gaussian_vector_analysis'
        end if

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
                        idvw => size(polar_component, dim=1), &
                        jdvw => size(polar_component, dim=2), &
                        br => workspace%real_polar_harmonic_coefficients, &
                        bi => workspace%imaginary_polar_harmonic_coefficients, &
                        cr => workspace%real_azimuthal_harmonic_coefficients, &
                        ci => workspace%imaginary_azimuthal_harmonic_coefficients, &
                        mdab => size(workspace%real_polar_harmonic_coefficients, dim=1), &
                        ndab => size(workspace%real_polar_harmonic_coefficients, dim=2), &
                        wvhags => workspace%forward_vector, &
                        lvhags => size(workspace%forward_vector), &
                        work => workspace%legendre_workspace, &
                        lwork => size(workspace%legendre_workspace), &
                        ierror => error_flag &
                        )

                        !
                        !==> Perform gaussian vector analysis
                        !
                        call aux%vhags(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                            mdab, ndab, wvhags, lvhags, work, lwork, ierror)

                    end associate
                end select
            end associate
        end select

        !
        !==> Address error flag
        !
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop 'Object of class(GaussianSphere) in '&
                    //'gaussian_vector_analysis'&
                    //'invalid specification of NUMBER_OF_LATITUDES'
            case (2)
                error stop 'Object of class(GaussianSphere) in '&
                    //'gaussian_vector_analysis'&
                    //'invalid specification of NUMBER_OF_LONGITUDES'
            case (3)
                error stop 'Object of class(GaussianSphere) in '&
                    //'gaussian_vector_analysis'&
                    //'invalid specification of VECTOR_SYMMETRIES'
            case (4)
                error stop 'Object of class(GaussianSphere) in '&
                    //'gaussian_vector_analysis'&
                    //'invalid specification of NUMBER_OF_SYNTHESES'
            case (5)
                error stop 'Object of class(GaussianSphere) in '&
                    //'gaussian_vector_analysis'&
                    //'invalid dim=1 extent for '&
                    //'polar_component (theta) or azimuthal_component (phi)'
            case (6)
                error stop 'Object of class(GaussianSphere) in '&
                    //'gaussian_vector_analysis'&
                    //'invalid dim=2 extent '&
                    //'polar_component (theta) or azimuthal_component (phi)'
            case (7)
                error stop 'Object of class(GaussianSphere) in '&
                    //'gaussian_vector_analysis'&
                    //'invalid dim=1 extent for br or cr'
            case (8)
                error stop 'Object of class(GaussianSphere) in '&
                    //'gaussian_vector_analysis'&
                    //'invalid dim=1 extent for bi or ci'
            case (9)
                error stop 'Object of class(GaussianSphere) in '&
                    //'gaussian_vector_analysis'&
                    //'invalid specification for forward_vector'
            case (10)
                error stop 'Object of class(GaussianSphere) in '&
                    //'gaussian_vector_analysis'&
                    //'invalid specification for legendre_workspace'
            case default
                error stop 'Object of class(GaussianSphere) in '&
                    //'gaussian_vector_analysis'&
                    //'Undetermined error flag'
        end select

    end subroutine gaussian_vector_analysis



    subroutine gaussian_vector_synthesis(this, polar_component, azimuthal_component)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        real (wp),              intent (out)    :: polar_component(:,:)
        real (wp),              intent (out)    :: azimuthal_component(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)    :: error_flag
        type (VhsgsAux) :: aux
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.this%initialized) then
            error stop 'Uninitialized object of class(GaussianSphere): '&
                //'in gaussian_vector_synthesis'
        end if

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
                        idvw => size(polar_component, dim=1),  &
                        jdvw => size(polar_component, dim=2),  &
                        br => workspace%real_polar_harmonic_coefficients, &
                        bi => workspace%imaginary_polar_harmonic_coefficients, &
                        cr => workspace%real_azimuthal_harmonic_coefficients, &
                        ci => workspace%imaginary_azimuthal_harmonic_coefficients, &
                        mdab => size(workspace%real_polar_harmonic_coefficients, dim=1), &
                        ndab => size(workspace%real_polar_harmonic_coefficients, dim=2), &
                        wvhsgs => workspace%backward_vector, &
                        lvhsgs => size(workspace%backward_vector), &
                        work => workspace%legendre_workspace, &
                        lwork => size(workspace%legendre_workspace), &
                        ierror => error_flag &
                        )
                        !
                        !==> Perform gaussian vector synthesis
                        !
                        call aux%vhsgs(nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
                            br, bi, cr, ci, mdab, ndab, wvhsgs, lvhsgs, work, lwork, ierror)

                    end associate
                end select
            end associate
        end select

        !
        !==> Address error flag
        !
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop 'Object of class(GaussianSphere) in gaussian_vector_synthesis'&
                    //'invalid specification of NUMBER_OF_LATITUDES'
            case (2)
                error stop 'Object of class(GaussianSphere) in gaussian_vector_synthesis'&
                    //'invalid specification of NUMBER_OF_LONGITUDES'
            case (3)
                error stop 'Object of class(GaussianSphere) in gaussian_vector_synthesis'&
                    //'invalid specification of VECTOR_SYMMETRIES'
            case (4)
                error stop 'Object of class(GaussianSphere) in gaussian_vector_synthesis'&
                    //'invalid specification of NUMBER_OF_SYNTHESES'
            case (5)
                error stop 'Object of class(GaussianSphere) in gaussian_vector_synthesis'&
                    //'invalid dim=1 extent for '&
                    //'polar_component (theta) or azimuthal_component (phi)'
            case (6)
                error stop 'Object of class(GaussianSphere) in gaussian_vector_synthesis'&
                    //'invalid dim=2 extent '&
                    //'polar_component (theta) or azimuthal_component (phi)'
            case (7)
                error stop 'Object of class(GaussianSphere) in gaussian_vector_synthesis'&
                    //'invalid dim=1 extent for br or cr'
            case (8)
                error stop 'Object of class(GaussianSphere) in gaussian_vector_synthesis'&
                    //'invalid dim=1 extent for bi or ci'
            case (9)
                error stop 'Object of class(GaussianSphere) in gaussian_vector_synthesis'&
                    //'invalid specification for backward_vector'
            case (10)
                error stop 'Object of class(GaussianSphere) in gaussian_vector_synthesis'&
                    //'invalid specification for legendre_workspace'
            case default
                error stop 'Object of class(GaussianSphere) in gaussian_vector_synthesis'&
                    //'Undetermined error flag'
        end select

    end subroutine gaussian_vector_synthesis
    


    function compute_surface_integral(this, scalar_function) result (return_value)
        !
        ! Purpose:
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
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianSphere), intent (in out) :: this
        real (wp),              intent (in)     :: scalar_function(:,:)
        real (wp)                               :: return_value
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)           :: k  !! counter
        real (wp), allocatable :: summation(:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.this%initialized) then
            error stop 'Uninitialized object of class(GaussianSphere): '&
                //'in get_surface_integral'
        end if

        !
        !==> Allocate memory
        !
        associate( nlat => this%NUMBER_OF_LATITUDES )

            allocate(summation(nlat))

        end associate

        !
        !==> compute the integrant
        !
        associate( grid => this%grid )
            select type(grid)
                class is (GaussianGrid)
                associate( &
                    nlat => grid%NUMBER_OF_LATITUDES, &
                    dphi => grid%LONGITUDINAL_MESH, &
                    wts => grid%gaussian_weights, &
                    f => scalar_function &
                    )
                    !
                    !==> Apply trapezoidal rule
                    !
                    do k = 1, nlat
                        summation(k) = sum(f(k, :)) * dphi
                    end do
                    !
                    !==> Apply gaussian quadrature
                    !
                    summation = summation * wts
                end associate
            end select
        end associate

        !
        !==> Set integral \int_{S^2} f( theta, phi ) dS
        !
        return_value = sum(summation)

        !
        !==> Release memory
        !
        deallocate(summation)

    end function compute_surface_integral



    subroutine compute_first_moment(this, scalar_function, first_moment)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianSphere),  intent (in out) :: this
        real (wp),               intent (in)     :: scalar_function(:,:)
        class (Vector),          intent (out)    :: first_moment
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)           :: k, l !! Counters
        real (wp), allocatable :: integrant(:,:,:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.this%initialized) then
            error stop 'Uninitialized object of class(GaussianSphere): '&
                //'in compute_first_moment'
        end if


        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            nlon => this%NUMBER_OF_LONGITUDES &
            )
            !
            !==> Allocate memory
            !
            allocate( integrant(nlat, nlon, 3) )


            do l = 1, nlon
                do k = 1, nlat
                    associate( &
                        u => this%unit_vectors%radial(k, l), &
                        f => scalar_function(k, l) &
                        )
                        !
                        !==> Compute integrant
                        !
                        integrant(k, l, 1) = u%x * f
                        integrant(k, l, 2) = u%y * f
                        integrant(k, l, 3) = u%z * f

                    end associate
                end do
            end do
        end associate


        associate( &
            m => first_moment, &
            f1 => integrant(:,:,1), &
            f2 => integrant(:,:,2), &
            f3 => integrant(:,:,3) &
            )
            !
            !==> Compute first moment
            !
            m%x = this%compute_surface_integral(f1)
            m%y = this%compute_surface_integral(f2)
            m%z = this%compute_surface_integral(f3)

        end associate

        ! Release memory
        deallocate( integrant )

    end subroutine compute_first_moment

    

    subroutine finalize_gaussian_sphere(this)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        type (GaussianSphere), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_gaussian_sphere


end module type_GaussianSphere
