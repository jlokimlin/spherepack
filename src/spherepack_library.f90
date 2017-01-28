module spherepack_library

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

    use geo2math_coordinate_transfer_routines, only: &
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

    use module_shpe, only: &
        shpe, shpei

    use module_shpg, only: &
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

    use type_RealPeriodicTransform, only: &
        RealPeriodicTransform, &
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

    ! Explicit typing only
    implicit none
    
    ! Everything is private unless stated otherwise
    private

    ! Constants
    public :: wp, ip
    public :: HALF_PI, PI, TWO_PI

    ! Derived data types
    public :: FFTpack
    public :: RealPeriodicTransform
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

    ! Procedural methods
    public :: divec, dives, divgc, divgs
    public :: compute_gaussian_latitudes_and_weights
    public :: geo2maths, math2geos, geo2mathv, math2geov
    public :: gradec, grades, gradgc, gradgs
    public :: hrffti, hrfftf, hrfftb
    public :: idivec, idives, idivgc, idivgs
    public :: idvtec, idvtes, idvtgc, idvtgs
    public :: igradec, igrades, igradgc, igradgs
    public :: ihgeod
    public :: isfvpec, isfvpes, isfvpgc, isfvpgs
    public :: islapec, islapes, islapgc, islapgs
    public :: ivlapec, ivlapes, ivlapgc, ivlapgs
    public :: ivrtec, ivrtes, ivrtgc, ivrtgs
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
    public :: vrtec, vrtes, vrtgc, vrtgs
    public :: vshifte, vshifti
    public :: vtsec, vtses, vtsgc, vtsgs
    public :: vtseci, vtsesi, vtsgci, vtsgsi
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

contains

    subroutine vecout_array(vec,nam,vec_size)

        ! Dummy arguments
        real(wp),         intent(in) ::vec(vec_size)
        character(len=*), intent(in) :: nam
        integer(ip),      intent(in) :: vec_size

        ! Local variables
        integer(ip) :: i

        write( stdout, 109) nam, (vec(i),i=1,vec_size)
109     format(1h a4,/(1h 8e11.4))

    end subroutine vecout_array

    subroutine vecout_scalar(vec,nam,vec_size)

        ! Dummy arguments
        real(wp),         intent(in) :: vec
        integer(ip),      intent(in) :: vec_size
        character(len=*), intent(in) ::  nam

        write( stdout, '(a, 8e11.4)') nam, vec

    end subroutine vecout_scalar

    subroutine vout_scalar(var,nam)

        ! Dummy arguments
        real(wp),         intent(in) :: var
        character(len=*), intent(in) :: nam

        write( stdout,'(a, e12.5)') nam, var

    end subroutine vout_scalar

    subroutine vout_array(var,nam)

        ! Dummy arguments
        real(wp),         intent(in) :: var(:)
        character(len=*), intent(in) :: nam

        write( stdout, *) nam, var

    end subroutine vout_array

    subroutine iout(ivar,nam)

        ! Dummy arguments
        integer(ip),      intent(in) :: ivar
        character(len=*), intent(in) :: nam

        write( stdout, '(a,i5)') nam, ivar

    end subroutine iout

    subroutine name(routine_name)

        ! Dummy arguments
        character(len=*), intent(in) :: routine_name

        write(stdout,'(a)') routine_name

    end subroutine name

    subroutine check_error(ierror)

        ! Dummy arguments
        integer(ip), intent(in) :: ierror

        if (ierror /= 0) write(stderr, '(a,i5)') '   ierror', ierror

    end subroutine check_error

end module spherepack_library
