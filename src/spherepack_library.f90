module spherepack_library

    use, intrinsic :: ISO_Fortran_env, only: &
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

    use module_geo2math, only: &
        geo2maths, math2geos, geo2mathv, math2geov

    use module_gradec, only: &
        gradec

    use module_grades, only: &
        grades

    use module_gradgc, only: &
        gradgc

    use module_gradgs, only: &
        gradgs

    use module_idvtec, only: &
        idvtec

    use module_idvtes, only: &
        idvtes

    use module_idvtgc, only: &
        idvtgc

    use module_idvtgs, only: &
        idvtgs

    use module_igradec, only: &
        igradec

    use module_igrades, only: &
        igrades

    use module_igradgc, only: &
        igradgc

    use module_igradgs, only: &
        igradgs

    use module_ihgeod, only: &
        ihgeod

    use module_isfvpec, only: &
        isfvpec

    use module_isfvpes, only: &
        isfvpes

    use module_isfvpgc, only: &
        isfvpgc

    use module_isfvpgs, only: &
        isfvpgs

    use module_islapec, only: &
        islapec

    use module_islapes, only: &
        islapes

    use module_islapgc, only: &
        islapgc

    use module_islapgs, only: &
        islapgs

    use module_ivlapec, only: &
        ivlapec

    use module_ivlapes, only: &
        ivlapes

    use module_ivlapgc, only: &
        ivlapgc

    use module_ivlapgs, only: &
        ivlapgs

    use module_ivrtec, only: &
        ivrtec

    use module_ivrtes, only: &
        ivrtes

    use module_ivrtgc, only: &
        ivrtgc

    use module_ivrtgs, only: &
        ivrtgs

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

    use module_slapec, only: &
        slapec

    use module_slapes, only: &
        slapes

    use module_slapgc, only: &
        slapgc

    use module_slapgs, only: &
        slapgs

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

    use module_vlapec, only: &
        vlapec

    use module_vlapes, only: &
        vlapes

    use module_vlapgc, only: &
        vlapgc

    use module_vlapgs, only: &
        vlapgs

    use module_vrtec, only: &
        vrtec

    use module_vrtes, only: &
        vrtes

    use module_vrtgc, only: &
        vrtgc

    use module_vrtgs, only: &
        vrtgs

    use module_vshifte, only: &
        vshifte, vshifti

    use module_vtsec, only: &
        vtsec, vtseci

    use module_vtses, only: &
        vtses, vtsesi

    use module_vtsgc, only: &
        vtsgc, vtsgci

    use module_vtsgs, only: &
        vtsgs, vtsgsi
    
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
    public :: vecout

    interface vecout
        module procedure vecout_array
        module procedure vecout_scalar
    end interface vecout

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

end module spherepack_library
