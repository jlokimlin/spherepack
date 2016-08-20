module spherepack_library

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, &
        TWO_PI

    use module_divec, only: &
        divec

    use module_dives, only: &
        dives

    use module_divgc, only: &
        divgc

    use module_divgs, only: &
        divgs

    use module_gaqd, only: &
        gaqd

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

    use module_idivec, only: &
        idivec

    use module_idives, only: &
        idives

    use module_idivgc, only: &
        idivgc

    use module_idivgs, only: &
        idivgs

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

    use module_shaec, only: &
        shaec, shaeci

    use scalar_analysis_routines, only: &
        shaes, shaesi, &
        shags, shagsi

    use module_shagc, only: &
        shagc, shagci

    use module_shpe, only: &
        shpe, shpei

    use module_shpg, only: &
        shpg, shpgi

    use module_shsec, only: &
        shsec, shseci

    use module_shses, only: &
        shses, shsesi

    use module_shsgc, only: &
        shsgc, shsgci

    use module_shsgs, only: &
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

    use module_vhaec, only: &
        vhaec, vhaeci

    use module_vhaes, only: &
        vhaes, vhaesi

    use module_vhagc, only: &
        vhagc, vhagci

    use module_vhags, only: &
        vhags, vhagsi

    use module_vhsec, only: &
        vhsec, vhseci

    use module_vhses, only: &
        vhses, vhsesi

    use module_vhsgc, only: &
        vhsgc, vhsgci

    use module_vhsgs, only: &
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

    use type_HFFTpack, only: &
        HFFTpack, &
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
    public :: PI, TWO_PI
    ! Classes
    public :: FFTpack
    public :: HFFTpack
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
    public :: gaqd
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

contains





end module spherepack_library
