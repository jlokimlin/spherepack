module spherepack

    use, intrinsic :: ISO_Fortran_env, only: &
        stderr => ERROR_UNIT, &
        stdout => OUTPUT_UNIT

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        HALF_PI, &
        PI, &
        TWO_PI

    use divergence_routines, only: &
        divec, dives, divgc, divgs, &
        idivec, idives, idivgc, idivgs

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    use coordinate_transfer_routines, only: &
        geo2maths, math2geos, geo2mathv, math2geov

    use gradient_routines, only: &
        gradec, grades, gradgc, gradgs, &
        igradec, igrades, igradgc, igradgs

    use module_idvtec, only: &
        idvtec

    use module_idvtes, only: &
        idvtes

    use module_idvtgc, only: &
        idvtgc

    use module_idvtgs, only: &
        idvtgs

    use icosahedral_geodesic_routines, only: &
        ihgeod

    use module_isfvpec, only: &
        isfvpec

    use module_isfvpes, only: &
        isfvpes

    use module_isfvpgc, only: &
        isfvpgc

    use module_isfvpgs, only: &
        isfvpgs

    use module_sfvpec, only: &
        sfvpec

    use module_sfvpes, only: &
        sfvpes

    use module_sfvpgc, only: &
        sfvpgc

    use module_sfvpgs, only: &
        sfvpgs

    use scalar_analysis_routines, only: &
        shaes, shaesi, &
        shaec, shaeci, &
        shags, shagsi, &
        shagc, shagci

    use scalar_projection_routines, only: &
        shpe, shpei, &
        shpg, shpgi

    use scalar_synthesis_routines, only: &
        shsec, shseci, &
        shses, shsesi, &
        shsgc, shsgci, &
        shsgs, shsgsi

    use scalar_laplacian_routines, only: &
        slapec, slapes, slapgc, slapgs, &
        islapec, islapes, islapgc, islapgs

    use module_sshifte, only: &
        sshifte, sshifti

    use module_trssph, only: &
        trssph

    use module_trvsph, only: &
        trvsph

    use vector_analysis_routines, only: &
        vhaec, vhaeci, &
        vhaes, vhaesi, &
        vhagc, vhagci, &
        vhags, vhagsi

    use vector_synthesis_routines, only: &
        vhsec, vhseci, &
        vhses, vhsesi, &
        vhsgc, vhsgci, &
        vhsgs, vhsgsi

    !    use module_visequ, only: &
    !        visequ
    !
    !    use module_visgau, only: &
    !        visgau
    !
    !    use module_visgeo, only: &
    !        visgeo

    use vector_laplacian_routines, only: &
        vlapec, vlapes, vlapgc, vlapgs, &
        ivlapec, ivlapes, ivlapgc, ivlapgs

    use vorticity_routines, only: &
        vrtec, vrtes, vrtgc, vrtgs, &
        ivrtec, ivrtes, ivrtgc, ivrtgs

    use module_vshifte, only: &
        vshifte, vshifti

    use colatitudinal_derivative_routines, only: &
        vtsec, vtses, vtsgc, vtsgs, &
        vtseci, vtsesi, vtsgci, vtsgsi
    
    use type_FFTpack, only: &
        FFTpack

    use type_RealPeriodicFastFourierTransform, only: &
        RealPeriodicFastFourierTransform, &
        hrffti, hrfftf, hrfftb

    use type_LegendreAux, only: &
        LegendreAux, &
        alfk, lfp, lfpt, lfim, lfin

    use type_Vector3D, only: &
        Vector3D, &
        Vector3DPointer, &
        assignment(=), &
        operator(*)
    
    use type_Sphere, only: &
        Sphere

    use type_GaussianGrid, only:&
        GaussianGrid
        
    use type_GaussianSphere, only: &
        GaussianSphere
    
    use type_RegularGrid, only: &
        RegularGrid

    use type_RegularSphere, only: &
        RegularSphere

    use type_SpherepackUtility, only: &
        SpherepackUtility

    ! Explicit typing only
    implicit none
    
    ! Everything is private unless stated otherwise
    private

    ! Constants
    public :: wp, ip
    public :: HALF_PI, PI, TWO_PI

    ! Derived data types
    public :: FFTpack
    public :: RealPeriodicFastFourierTransform
    public :: LegendreAux
    public :: GaussianGrid
    public :: GaussianSphere
    public :: RegularGrid
    public :: RegularSphere
    public :: Sphere
    public :: Vector3D
    public :: Vector3DPointer
    public :: assignment(=)
    public :: operator(*)

    ! Colatitude derivative
    public :: vtsec, vtses, vtsgc, vtsgs
    public :: vtseci, vtsesi, vtsgci, vtsgsi

    ! Gradient
    public :: gradec, grades, gradgc, gradgs

    ! Inverse gradient
    public :: igradec, igrades, igradgc, igradgs

    ! Divergence
    public :: divec, dives, divgc, divgs

    ! Inverse divergence
    public :: idivec, idives, idivgc, idivgs

    ! Vorticity
    public :: vrtec, vrtes, vrtgc, vrtgs

    ! Inverse vorticity
    public :: ivrtec, ivrtes, ivrtgc, ivrtgs

    ! Gaussian wts & pts
    public :: compute_gaussian_latitudes_and_weights

    ! Geo/math coordinate transfers
    public :: geo2maths, math2geos, geo2mathv, math2geov

    ! Multiple ffts
    public :: hrffti, hrfftf, hrfftb

    public :: idvtec, idvtes, idvtgc, idvtgs

    public :: ihgeod
    public :: isfvpec, isfvpes, isfvpgc, isfvpgs
    public :: islapec, islapes, islapgc, islapgs
    public :: ivlapec, ivlapes, ivlapgc, ivlapgs

    public :: sfvpec, sfvpes, sfvpgc, sfvpgs
    public :: shaec, shaes, shagc, shags
    public :: shaeci, shaesi, shagci, shagsi
    public :: shpe, shpei
    public :: shpg, shpgi
    public :: shsec, shses, shsgc, shsgs
    public :: shseci, shsesi, shsgci, shsgsi
    public :: slapec, slapes, slapgc, slapgs
    public :: sshifte, sshifti
    public :: trssph, trvsph
    public :: vhaec, vhaes, vhagc, vhags
    public :: vhaeci, vhaesi, vhagci, vhagsi
    public :: vhsec, vhses, vhsgc, vhsgs
    public :: vhseci, vhsesi, vhsgci, vhsgsi
    !    public :: visequ
    !    public :: visgau
    !    public :: visgeo
    public :: vlapec, vlapes, vlapgc, vlapgs

    public :: vshifte, vshifti

    public :: alfk, lfp, lfpt, lfim, lfin

    ! Temporary solution for testing
    public :: vecout, iout, name, vout, check_error

    interface vecout
        module procedure vecout_array
        module procedure vecout_scalar
    end interface vecout

    interface vout
        module procedure vout_array
        module procedure vout_scalar
    end interface vout

    type, public :: SpherepackWaveTable
        ! Type components
        real(wp), allocatable :: scalar_forward_regular(:)
        real(wp), allocatable :: scalar_backward_regular(:)
        real(wp), allocatable :: scalar_forward_regular_saved(:)
        real(wp), allocatable :: scalar_backward_regular_saved(:)
        real(wp), allocatable :: scalar_forward_gaussian(:)
        real(wp), allocatable :: scalar_backward_gaussian(:)
        real(wp), allocatable :: scalar_forward_gaussian_saved(:)
        real(wp), allocatable :: scalar_backward_gaussian_saved(:)
        real(wp), allocatable :: vector_forward_regular(:)
        real(wp), allocatable :: vector_backward_regular(:)
        real(wp), allocatable :: vector_forward_regular_saved(:)
        real(wp), allocatable :: vector_backward_regular_saved(:)
        real(wp), allocatable :: vector_forward_gaussian(:)
        real(wp), allocatable :: vector_backward_gaussian(:)
        real(wp), allocatable :: vector_forward_gaussian_saved(:)
        real(wp), allocatable :: vector_backward_gaussian_saved(:)
        logical :: first_call_scalar_forward_regular = .true.
        logical :: first_call_scalar_backward_regular = .true.
        logical :: first_call_scalar_forward_regular_saved = .true.
        logical :: first_call_scalar_backward_regular_saved = .true.
        logical :: first_call_scalar_forward_gaussian = .true.
        logical :: first_call_scalar_backward_gaussian = .true.
        logical :: first_call_scalar_forward_gaussian_saved = .true.
        logical :: first_call_scalar_backward_gaussian_saved = .true.
        logical :: first_call_vector_forward_regular = .true.
        logical :: first_call_vector_backward_regular = .true.
        logical :: first_call_vector_forward_regular_saved = .true.
        logical :: first_call_vector_backward_regular_saved = .true.
        logical :: first_call_vector_forward_gaussian = .true.
        logical :: first_call_vector_backward_gaussian = .true.
        logical :: first_call_vector_forward_gaussian_saved = .true.
        logical :: first_call_vector_backward_gaussian_saved = .true.
    contains
        ! Type-bound procedures
        procedure :: initialize_scalar_forward_gaussian_saved
    end type SpherepackWaveTable

    type, public :: SpherepackScalarHarmonic
        ! Type components
        logical,               public :: initialized = .false.
        integer(ip),           public :: NUMBER_OF_LONGITUDES = 0
        integer(ip),           public :: NUMBER_OF_LATITUDES = 0
        integer(ip),           public :: NUMBER_OF_SYNTHESES = 0
        real(wp), allocatable, public :: real_component(:, :, :)
        real(wp), allocatable, public :: imaginary_component(:, :, :)
    contains
        ! Type-bound procedures
        procedure, public  :: create => create_spherepack_scalar_harmonic
        procedure, public  :: destroy => destroy_spherepack_scalar_harmonic
        procedure, private :: copy_spherepack_scalar_harmonic
        ! Generic type-bound procedures
        generic, public :: assignment(=) => copy_spherepack_scalar_harmonic
    end type SpherepackScalarHarmonic

    ! Declare user-defined constructor
    interface SpherepackScalarHarmonic
        module procedure spherepack_scalar_harmonic_constructor
    end interface

contains

    function spherepack_scalar_harmonic_constructor(nlat, nlon, nt) &
        result (return_value)

        ! Dummy arguments
        integer(ip),           intent(in) :: nlat ! Number of latitudinal points 0 <= theta <= pi
        integer(ip),           intent(in) :: nlon ! Number of longitudinal points 0 <= phi <= 2*pi
        integer(ip), optional, intent(in) :: nt ! Number of syntheses
        type(SpherepackScalarHarmonic)    :: return_value

        ! Local variables
        integer(ip) :: number_of_syntheses

        ! Address optional argument
        if (present(nt)) then
            number_of_syntheses = nt
        else
            number_of_syntheses = 1
        end if

        call return_value%create(nlat, nlon, number_of_syntheses)

    end function spherepack_scalar_harmonic_constructor

    subroutine create_spherepack_scalar_harmonic(self, nlat, nlon, nt)

        ! Dummy arguments
        class(SpherepackScalarHarmonic), intent(inout) :: self
        integer(ip),                     intent(in)    :: nlat
        integer(ip),                     intent(in)    :: nlon
        integer(ip),                     intent(in)    :: nt

        ! Local variables
        integer(ip) :: mdab, ndab

        ! Ensure that object is usable
        call self%destroy()

        !  Set constants
        self%NUMBER_OF_LATITUDES = nlat
        self%NUMBER_OF_LONGITUDES = nlon
        self%NUMBER_OF_SYNTHESES = nt

        !  Set upper limit for vector m subscript
        select case(mod(nlon, 2))
            case(0)
                mdab = min(nlat, (nlon + 2)/2)
            case default
                mdab = min(nlat, (nlon + 1)/2)
        end select

        !mdab = nlat
        !ndab = nlat

        !  Allocate memory
        allocate(self%real_component(mdab, ndab, nt))
        allocate(self%imaginary_component(mdab, ndab, nt))

        ! Set flag
        self%initialized = .true.

    end subroutine create_spherepack_scalar_harmonic

    subroutine copy_spherepack_scalar_harmonic(self, other)

        ! Dummy arguments
        class(SpherepackScalarHarmonic), intent(out) :: self
        class(SpherepackScalarHarmonic), intent(in)  :: other

        ! Check if object is usable
        if (.not.other%initialized) then
            error stop 'Uninitialized object of class(SpherepackScalarHarmonic): '&
                //'in assignment(=) '
        end if

        !  Make copies
        self%initialized = other%initialized
        self%NUMBER_OF_LATITUDES = other%NUMBER_OF_LATITUDES
        self%NUMBER_OF_LONGITUDES = other%NUMBER_OF_LONGITUDES
        self%NUMBER_OF_SYNTHESES = other%NUMBER_OF_SYNTHESES
        self%real_component = other%real_component
        self%imaginary_component = other%imaginary_component

    end subroutine copy_spherepack_scalar_harmonic

    subroutine destroy_spherepack_scalar_harmonic(self)

        ! Dummy arguments
        class(SpherepackScalarHarmonic), intent(inout) :: self

        ! Check flag
        if (.not.self%initialized) return

        if (allocated(self%real_component)) then
            deallocate (self%real_component)
        end if

        if (allocated(self%imaginary_component)) then
            deallocate (self%imaginary_component)
        end if

        !  Reset constants
        self%NUMBER_OF_LONGITUDES = 0
        self%NUMBER_OF_LATITUDES = 0
        self%NUMBER_OF_SYNTHESES = 0

        ! Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_spherepack_scalar_harmonic

    subroutine initialize_scalar_forward_gaussian_saved(self, nlat, nlon)
        class(SpherepackWaveTable), intent(inout) :: self
        integer(ip),                intent(in)    :: nlat
        integer(ip),                intent(in)    :: nlon

        ! Local variables
        type(SpherepackUtility) :: util
        integer(ip)             :: lshags, lwork, ldwork, error_flag

        ! Release memory if necessary
        if (allocated(self%scalar_forward_gaussian_saved)) then
            deallocate (self%scalar_forward_gaussian_saved)
        end if

        ! Compute required wavetable size
        lshags = util%get_lshags(nlat, nlon)

        ! Allocate memory
        allocate (self%scalar_forward_gaussian_saved(lshags))

        ! Compute required workspace sizes
        lwork = util%get_lwork_for_shagsi(nlat)
        ldwork = util%get_ldwork_for_shagsi(nlat)

        ! Initialize forward wavetable
        block
            real(wp) :: work(lwork), dwork(ldwork)
            associate (wshags => self%scalar_forward_gaussian_saved)
                call shagsi(nlat, nlon, wshags, error_flag)
            end associate
        end block

    end subroutine initialize_scalar_forward_gaussian_saved

    subroutine spherepack_shags(symmetries, scalar_function, wavetable, scalar_harmonic)
        integer(ip),                     intent(in)    :: symmetries
        real(wp),                        intent(in)    :: scalar_function(:,:,:)
        class(SpherepackWaveTable),      intent(inout) :: wavetable
        class(SpherepackScalarHarmonic), intent(inout) :: scalar_harmonic

        ! Local variables
        type(SpherepackUtility) :: util
        integer(ip)             :: lshags, lwork, ldwork, error_flag

        if (.not. scalar_harmonic%initialized) then
            error stop 'Unusable scalar_harmonic in spherepack_shags'
        end if

        associate (&
            nlat => scalar_harmonic%NUMBER_OF_LATITUDES, &
            nlon => scalar_harmonic%NUMBER_OF_LONGITUDES, &
            nt => scalar_harmonic%NUMBER_OF_SYNTHESES, &
            a => scalar_harmonic%real_component, &
            b => scalar_harmonic%imaginary_component, &
            mdab => size(scalar_harmonic%real_component, dim=1), &
            ndab => size(scalar_harmonic%real_component, dim=2), &
            g => scalar_function, &
            idg => size(scalar_function, dim=1), &
            jdg => size(scalar_function, dim=2), &
            isym => symmetries, &
            ierror => error_flag &
           )

            ! Ensure that forward wavetable is usable
            if (wavetable%first_call_scalar_forward_gaussian_saved) then
                call wavetable%initialize_scalar_forward_gaussian_saved(nlat, nlon)
            end if

            ! Compute required workspace size
            lwork = util%get_lwork_for_gaussian_saved(isym, nt, nlat, nlon)

            block
                real(wp) :: work(lwork)
                associate (wshags => wavetable%scalar_forward_gaussian_saved)
                    call shags(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                        wshags, lshags, work, lwork, ierror)
                end associate
            end block
        end associate

        ! Set call flag
        wavetable%first_call_scalar_forward_gaussian_saved = .false.

    end subroutine spherepack_shags

    subroutine vecout_array(vec, nam, vec_size)

        ! Dummy arguments
        real(wp),         intent(in) ::vec(vec_size)
        character(len=*), intent(in) :: nam
        integer(ip),      intent(in) :: vec_size

        ! Local variables
        integer(ip) :: i

        write( stdout, 109) nam, (vec(i), i=1, vec_size)
109     format(1h a4, /(1h 8e11.4))

    end subroutine vecout_array

    subroutine vecout_scalar(vec, nam, vec_size)

        ! Dummy arguments
        real(wp),         intent(in) :: vec
        integer(ip),      intent(in) :: vec_size
        character(len=*), intent(in) ::  nam

        write( stdout, '(a, 8e11.4)') nam, vec

    end subroutine vecout_scalar

    subroutine vout_scalar(var, nam)

        ! Dummy arguments
        real(wp),         intent(in) :: var
        character(len=*), intent(in) :: nam

        write( stdout, '(a, e12.5)') nam, var

    end subroutine vout_scalar

    subroutine vout_array(var, nam)

        ! Dummy arguments
        real(wp),         intent(in) :: var(:)
        character(len=*), intent(in) :: nam

        write( stdout, *) nam, var

    end subroutine vout_array

    subroutine iout(ivar, nam)

        ! Dummy arguments
        integer(ip),      intent(in) :: ivar
        character(len=*), intent(in) :: nam

        write( stdout, '(a, i5)') nam, ivar

    end subroutine iout

    subroutine name(routine_name)

        ! Dummy arguments
        character(len=*), intent(in) :: routine_name

        write(stdout, '(a)') routine_name

    end subroutine name

    subroutine check_error(ierror)

        ! Dummy arguments
        integer(ip), intent(in) :: ierror

        if (ierror /= 0) write(stderr, '(a, i5)') '   ierror', ierror

    end subroutine check_error

end module spherepack
