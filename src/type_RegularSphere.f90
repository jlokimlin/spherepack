module type_RegularSphere

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_Sphere, only: &
        Sphere

    use type_RegularWorkspace, only: &
        RegularWorkspace

    use type_RegularGrid, only: &
        RegularGrid

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
    public :: RegularSphere


    ! Declare derived data type
    type, extends (Sphere), public :: RegularSphere
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure, public  :: create => create_regular_sphere
        procedure, public  :: destroy => destroy_regular_sphere
        procedure, public  :: perform_scalar_analysis => regular_scalar_analysis
        procedure, public  :: perform_scalar_synthesis => regular_scalar_synthesis
        procedure, public  :: vector_analysis_from_spherical_components => regular_vector_analysis
        procedure, public  :: perform_vector_synthesis => regular_vector_synthesis
        final              :: finalize_regular_sphere
        !----------------------------------------------------------------------
    end type RegularSphere


contains


    subroutine create_regular_sphere( this, nlat, nlon, ntrunc, isym, itype, isynt, rsphere )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularSphere), intent (in out)        :: this
        integer (ip),          intent (in)            :: nlat
        integer (ip),          intent (in)            :: nlon
        integer (ip),          intent (in), optional  :: ntrunc
        integer (ip),          intent (in), optional  :: isym      !! Either 0, 1, or 2
        integer (ip),          intent (in), optional  :: itype     !! Either 0, 1, 2, 3, ..., 8
        integer (ip),          intent (in), optional  :: isynt
        real (wp),             intent (in), optional  :: rsphere
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: ntrunc_op
        integer (ip) :: isym_op
        integer (ip) :: ityp_op
        integer (ip) :: isynt_op
        real (wp)    :: rsphere_op
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Allocate polymorphic components
        allocate( RegularGrid :: this%grid )
        allocate( RegularWorkspace :: this%workspace )

        ! Initialize polymorphic types
        associate( &
            grid => this%grid, &
            workspace => this%workspace &
            )
            ! Initialize Regular grid
            select type (grid)
                class is (RegularGrid)
                call grid%create( nlat, nlon )
            end select
            ! Initialize Regular workspace
            select type (workspace)
                class is (RegularWorkspace)
                call workspace%create( nlat, nlon )
            end select
        end associate

        ! Initialize constants
        ntrunc_op = nlat - 1
        isym_op = 0
        ityp_op = 0
        isynt_op = 1
        rsphere_op = 1.0_wp

        ! Address optional arguments
        if (present(ntrunc)) ntrunc_op = ntrunc
        if (present(isym)) isym_op = isym
        if (present(itype)) ityp_op = itype
        if (present(isynt)) isynt_op = isynt
        if (present(rsphere)) rsphere_op = rsphere

        ! Create parent type
        call this%create_sphere( &
        nlat=nlat, nlon=nlon, ntrunc=ntrunc_op, &
        isym=isym_op, itype=ityp_op, isynt=isynt_op, rsphere=rsphere_op )

        ! Set initialization flag
        this%initialized = .true.
        
    end subroutine create_regular_sphere



    subroutine destroy_regular_sphere( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularSphere), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if ( this%initialized .eqv. .false. ) return

        ! Release memory from parent type
        call this%destroy_sphere()

        ! Reset initialization flag
        this%initialized = .false.

    end subroutine destroy_regular_sphere



    subroutine regular_scalar_analysis( this, scalar_function )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularSphere), intent (in out) :: this
        real (wp),             intent (in)     :: scalar_function(:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: error_flag
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(RegularSphere): '&
                //'uninitialized object in regular_SCALAR_ANALYSIS'
        end if

        select type (this)
            class is (RegularSphere)
            associate( workspace => this%workspace )
                select type (workspace)
                    class is (RegularWorkspace)
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
                        wshaes => workspace%forward_scalar, &
                        lshaes => size(workspace%forward_scalar), &
                        work => workspace%legendre_workspace, &
                        lwork => size(workspace%legendre_workspace), &
                        ierror => error_flag &
                        )
                        call shaes( nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                            wshaes, lshaes, work, lwork, ierror )
                    end associate
                end select
            end associate
        end select

        ! Address the error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'TYPE(RegularSphere) in regular_SCALAR_ANALYSIS'&
                    // 'Error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'TYPE(RegularSphere) in regular_SCALAR_ANALYSIS'&
                    //'Error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'TYPE(RegularSphere) in regular_SCALAR_ANALYSIS'&
                    //'Invalid extent for SCALAR_FORWARD'
            case(4)
                error stop 'TYPE(RegularSphere) in regular_SCALAR_ANALYSIS'&
                    //'Invalid extent for LEGENDRE_WORKSPACE'
            case(5)
                error stop 'TYPE(RegularSphere) in regular_SCALAR_ANALYSIS'&
                    //'Invalid extent for DWORK'
            case default
                error stop 'TYPE(RegularSphere) in regular_SCALAR_ANALYSIS'&
                    // 'Undetermined error flag'
        end select

    end subroutine regular_scalar_analysis
    


    subroutine regular_scalar_synthesis( this, scalar_function )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularSphere), intent (in out) :: this
        real (wp),             intent (out)    :: scalar_function(:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: error_flag
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(RegularSphere): '&
                //'uninitialized object in regular_SCALAR_SYNTHESIS'
        end if

        ! Perform (real) spherical harmonic synthesis
        select type (this)
            class is (RegularSphere)
            associate( workspace => this%workspace )
                select type (workspace)
                    class is (RegularWorkspace)
                    ! Associate various quantities
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
                        wshses => workspace%backward_scalar, &
                        lshses => size(workspace%backward_scalar ), &
                        work => workspace%legendre_workspace, &
                        lwork => size(workspace%legendre_workspace ), &
                        ierror => error_flag &
                        )
                        call shses(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                            wshses, lshses, work, lwork, ierror)
                    end associate
                end select
            end associate
        end select

        ! Address the error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'TYPE(RegularSphere) in regular_SCALAR_SYNTHESIS'&
                    // 'Error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'TYPE(RegularSphere) in regular_SCALAR_SYNTHESIS'&
                    //'Error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'TYPE(RegularSphere) in regular_SCALAR_SYNTHESIS'&
                    //'Invalid extent for SCALAR_FORWARD'
            case(4)
                error stop 'TYPE(RegularSphere) in regular_SCALAR_SYNTHESIS'&
                    //'Invalid extent for LEGENDRE_WORKSPACE'
            case(5)
                error stop 'TYPE(RegularSphere) in regular_SCALAR_SYNTHESIS'&
                    //'Invalid extent for DWORK'
            case default
                error stop 'TYPE(RegularSphere) in regular_SCALAR_SYNTHESIS'&
                    // 'Undetermined error flag'
        end select

    end subroutine regular_scalar_synthesis



    subroutine regular_vector_analysis( this, &
        polar_component, azimuthal_component )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularSphere), intent (in out) :: this
        real (wp),             intent (in)     :: polar_component(:,:)
        real (wp),             intent (in)     :: azimuthal_component(:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: error_flag
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(RegularSphere): '&
                //'uninitialized object in REGULAR_VECTOR_ANALYSIS'
        end if

        ! Perform vector analysis
        select type (this)
            class is (RegularSphere)
            associate( workspace => this%workspace )
                select type (workspace)
                    class is (RegularWorkspace)
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
                        wvhaes => workspace%forward_vector, &
                        lvhaes => size(workspace%forward_vector ), &
                        work => workspace%legendre_workspace, &
                        lwork => size(workspace%legendre_workspace ), &
                        ierror => error_flag &
                        )
                        call vhaes(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                            mdab, ndab, wvhaes, lvhaes, work, lwork, ierror )
                    end associate
                end select
            end associate
        end select

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'TYPE(RegularSphere) in '&
                    //'REGULAR_VECTOR_ANALYSIS'&
                    //'Error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'TYPE(RegularSphere) in '&
                    //'REGULAR_VECTOR_ANALYSIS'&
                    //'Error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'TYPE(RegularSphere) in '&
                    //'REGULAR_VECTOR_ANALYSIS'&
                    //'Error in the specification of VECTOR_SYMMETRIES'
            case(4)
                error stop 'TYPE(RegularSphere) in '&
                    //'REGULAR_VECTOR_ANALYSIS'&
                    //'Error in the specification of NUMBER_OF_SYNTHESES'
            case(5)
                error stop 'TYPE(RegularSphere) in '&
                    //'REGULAR_VECTOR_ANALYSIS'&
                    //'Invalid DIM=1 extent for '&
                    //'POLAR_COMPONENT (THETA) or AZIMUTHAL_COMPONENT (PHI)'
            case(6)
                error stop 'TYPE(RegularSphere) in '&
                    //'REGULAR_VECTOR_ANALYSIS'&
                    //'Invalid DIM=2 extent '&
                    //'POLAR_COMPONENT (THETA) or AZIMUTHAL_COMPONENT (PHI)'
            case(7)
                error stop 'TYPE(RegularSphere) in '&
                    //'REGULAR_VECTOR_ANALYSIS'&
                    //'Invalid DIM=1 extent for BR or CR'
            case(8)
                error stop 'TYPE(RegularSphere) in '&
                    //'REGULAR_VECTOR_ANALYSIS'&
                    //'Invalid DIM=1 extent for BI or CI'
            case(9)
                error stop 'TYPE(RegularSphere) in '&
                    //'REGULAR_VECTOR_ANALYSIS'&
                    //'Invalid extent for FORWARD_VECTOR'
            case(10)
                error stop 'TYPE(RegularSphere) in '&
                    //'REGULAR_VECTOR_ANALYSIS'&
                    //'Invalid extent for LEGENDRE_WORKSPACE'
            case default
                error stop 'TYPE(RegularSphere) in '&
                    //'REGULAR_VECTOR_ANALYSIS'&
                    //'Undetermined error flag'
        end select

    end subroutine regular_vector_analysis



    subroutine regular_vector_synthesis( this, polar_component, azimuthal_component )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularSphere), intent (in out) :: this
        real (wp),             intent (out)    :: polar_component(:,:)
        real (wp),             intent (out)    :: azimuthal_component(:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: error_flag
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(RegularSphere): '&
                //'uninitialized object in regular_VECTOR_SYNTHESIS'
        end if

        ! Perform vector analysis
        select type (this)
            class is (RegularSphere)
            associate( workspace => this%workspace )
                select type (workspace)
                    class is (RegularWorkspace)
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
                        wvhses => workspace%backward_vector, &
                        lvhses => size(workspace%backward_vector ), &
                        work => workspace%legendre_workspace, &
                        lwork => size(workspace%legendre_workspace ), &
                        ierror => error_flag &
                        )
                        call vhses(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                            mdab, ndab, wvhses, lvhses, work, lwork, ierror )
                    end associate
                end select
            end associate
        end select

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'TYPE(RegularSphere) in REGULAR_VECTOR_SYNTHESIS'&
                    //'Error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'TYPE(RegularSphere) in REGULAR_VECTOR_SYNTHESIS'&
                    //'Error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'TYPE(RegularSphere) in REGULAR_VECTOR_SYNTHESIS'&
                    //'Error in the specification of VECTOR_SYMMETRIES'
            case(4)
                error stop 'TYPE(RegularSphere) in REGULAR_VECTOR_SYNTHESIS'&
                    //'Error in the specification of NUMBER_OF_synthESES'
            case(5)
                error stop 'TYPE(RegularSphere) in REGULAR_VECTOR_SYNTHESIS'&
                    //'Invalid DIM=1 extent for '&
                    //'POLAR_COMPONENT (THETA) or AZIMUTHAL_COMPONENT (PHI)'
            case(6)
                error stop 'TYPE(RegularSphere) in REGULAR_VECTOR_SYNTHESIS'&
                    //'Invalid DIM=2 extent '&
                    //'POLAR_COMPONENT (THETA) or AZIMUTHAL_COMPONENT (PHI)'
            case(7)
                error stop 'TYPE(RegularSphere) in REGULAR_VECTOR_SYNTHESIS'&
                    //'Invalid DIM=1 extent for BR or CR'
            case(8)
                error stop 'TYPE(RegularSphere) in REGULAR_VECTOR_SYNTHESIS'&
                    //'Invalid DIM=1 extent for BI or CI'
            case(9)
                error stop 'TYPE(RegularSphere) in REGULAR_VECTOR_SYNTHESIS'&
                    //'Invalid extent for BACKWARD_VECTOR'
            case(10)
                error stop 'TYPE(RegularSphere) in REGULAR_VECTOR_SYNTHESIS'&
                    //'Invalid extent for LEGENDRE_WORKSPACE'
            case default
                error stop 'TYPE(RegularSphere) in REGULAR_VECTOR_SYNTHESIS'&
                    //'Undetermined error flag'
        end select

    end subroutine regular_vector_synthesis
    

    subroutine finalize_regular_sphere( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (RegularSphere), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_regular_sphere


end module type_RegularSphere
