!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                      SPHEREPACK version 3.2                   *
!     *                                                               *
!     *       A Package of Fortran Subroutines and Programs           *
!     *                                                               *
!     *              for Modeling Geophysical Processes               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *                  John Adams and Paul Swarztrauber             *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
!
! ... file vhags.f90
!
!     this file contains code and documentation for subroutines
!     vhags and vhagsi
!
! ... files which must be loaded with vhags.f90
!
!     type_SpherepackAux.f90, type_HFFTpack.f90, compute_gaussian_latitudes_and_weights.f90
!
!     subroutine vhags(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci,
!    +                 mdab, ndab, wvhags, lvhags, work, lwork, ierror)
!
!     subroutine vhags performs the vector spherical harmonic analysis
!     on the vector field (v, w) and stores the result in the arrays
!     br, bi, cr, and ci. v(i, j) and w(i, j) are the colatitudinal
!     (measured from the north pole) and east longitudinal components
!     respectively, located at the gaussian colatitude point theta(i)
!     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
!     representation of (v, w) is given at output parameters v, w in
!     subroutine vhses.  
!
!     input parameters
!
!     nlat   the number of points in the gaussian colatitude grid on the
!            full sphere. these lie in the interval (0, pi) and are computed
!            in radians in theta(1) <...< theta(nlat) by subroutine compute_gaussian_latitudes_and_weights.
!            if nlat is odd the equator will be included as the grid point
!            theta((nlat+1)/2).  if nlat is even the equator will be
!            excluded as a grid point and will lie half way between
!            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
!            note: on the half sphere, the number of grid points in the
!            colatitudinal direction is nlat/2 if nlat is even or
!            (nlat+1)/2 if nlat is odd.
!
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     ityp   = 0  no symmetries exist about the equator. the analysis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
!                 j=1, ..., nlon.
!
!            = 1  no symmetries exist about the equator. the analysis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
!                 j=1, ..., nlon. the curl of (v, w) is zero. that is,
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 2  no symmetries exist about the equator. the analysis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
!                 j=1, ..., nlon. the divergence of (v, w) is zero. i.e.,
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!            = 3  v is symmetric and w is antisymmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 4  v is symmetric and w is antisymmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the curl of (v, w) is zero. that is,
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 5  v is symmetric and w is antisymmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the divergence of (v, w) is zero. i.e.,
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!            = 6  v is antisymmetric and w is symmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 7  v is antisymmetric and w is symmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the curl of (v, w) is zero. that is,
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 8  v is antisymmetric and w is symmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the divergence of (v, w) is zero. i.e.,
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!
!     nt     the number of analyses.  in the program that calls vhags,
!            the arrays v, w, br, bi, cr, and ci can be three dimensional
!            in which case multiple analyses will be performed.
!            the third index is the analysis index which assumes the 
!            values k=1, ..., nt.  for a single analysis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that all the arrays are two
!            dimensional.
!
!     v, w    two or three dimensional arrays (see input parameter nt)
!            that contain the vector function to be analyzed.
!            v is the colatitudnal component and w is the east 
!            longitudinal component. v(i, j), w(i, j) contain the
!            components at the gaussian colatitude point theta(i)
!            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges
!            are defined above at the input parameter ityp.
!
!     idvw   the first dimension of the arrays v, w as it appears in
!            the program that calls vhags. if ityp .le. 2 then idvw
!            must be at least nlat.  if ityp .gt. 2 and nlat is
!            even then idvw must be at least nlat/2. if ityp .gt. 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls vhags. jdvw must be at least nlon.
!
!     mdab   the first dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhags. mdab must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhags. ndab must be at
!            least nlat.
!
!     wvhags an array which must be initialized by subroutine vhgsi.
!            once initialized, wvhags can be used repeatedly by vhags
!            as long as nlon and nlat remain unchanged.  wvhags must
!            not be altered between calls of vhags.
!
!     lvhags the dimension of the array wvhags as it appears in the
!            program that calls vhags. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhags must be at least
!
!            l1*l2(nlat+nlat-l1+1)+nlon+15
!
!        ??? (nlat+1)*(nlat+1)*nlat/2 + nlon + 15
!
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhags. define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if ityp .le. 2 then lwork must be at least
!            the larger of the two quantities
!
!               3*nlat*(nlat+1)+2  (required by vhagsi)
!
!            and
!
!               (2*nt+1)*nlat*nlon
!
!            if ityp .gt. 2 then lwork must be at least
!            the larger of the two quantities
!
!               3*nlat*(nlat+1)+2  (required by vhagsi)
!
!            and
!
!              (2*nt+1)*l2*nlon
!
!
!     **************************************************************
!
!     output parameters
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!     cr, ci  that contain the vector spherical harmonic coefficients
!            in the spectral representation of v(i, j) and w(i, j) given
!            in the discription of subroutine vhses. br(mp1, np1),
!            bi(mp1, np1), cr(mp1, np1), and ci(mp1, np1) are computed
!            for mp1=1, ..., mmax and np1=mp1, ..., nlat except for np1=nlat
!            and odd mp1. mmax=min(nlat, nlon/2) if nlon is even or
!            mmax=min(nlat, (nlon+1)/2) if nlon is odd.
!      
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of ityp
!            = 4  error in the specification of nt
!            = 5  error in the specification of idvw
!            = 6  error in the specification of jdvw
!            = 7  error in the specification of mdab
!            = 8  error in the specification of ndab
!            = 9  error in the specification of lvhags
!            = 10 error in the specification of lwork
!
!
!     subroutine vhagsi(nlat, nlon, wvhags, lvhags, dwork, ldwork, ierror)
!
!     subroutine vhagsi initializes the array wvhags which can then be
!     used repeatedly by subroutine vhags until nlat or nlon is changed.
!
!     input parameters
!
!     nlat   the number of points in the gaussian colatitude grid on the
!            full sphere. these lie in the interval (0, pi) and are computed
!            in radians in theta(1) <...< theta(nlat) by subroutine compute_gaussian_latitudes_and_weights.
!            if nlat is odd the equator will be included as the grid point
!            theta((nlat+1)/2).  if nlat is even the equator will be
!            excluded as a grid point and will lie half way between
!            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
!            note: on the half sphere, the number of grid points in the
!            colatitudinal direction is nlat/2 if nlat is even or
!            (nlat+1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     lvhags the dimension of the array wvhags as it appears in the
!            program that calls vhagsi. lvhags must be at least
!
!               3*nlat*(nlat+1)+2  (required by vhagsi)
!
!     dwork  a real work space that does not need to be saved
!
!     ldwork the dimension of the array dwork as it appears in the
!            program that calls vhagsi. ldwork must be at least
!
!                   (3*nlat*(nlat+3)+2)/2
!
!     **************************************************************
!
!     output parameters
!
!     wvhags an array which is initialized for use by subroutine vhags.
!            once initialized, wvhags can be used repeatedly by vhags
!            as long as nlat and nlon remain unchanged.  wvhags must not
!            be altered between calls of vhags.
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lvhags
!            = 4  error in the specification of ldwork
!
module module_vhags

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI

    use type_HFFTpack, only: &
        HFFTpack

    use type_SpherepackAux, only: &
        SpherepackAux

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: vhags
    public :: vhagsi
    public :: VhagsAux

    ! Declare derived data type
    type, public :: VhagsAux
        !-----------------------------------------
        ! Type components
        !-----------------------------------------
    contains
        !-----------------------------------------
        ! Type-bound procedures
        !-----------------------------------------
        procedure, nopass :: vhags
        procedure, nopass :: vhagsi
        procedure, nopass :: get_lvhags
        procedure, nopass :: get_ldwork
        procedure, nopass :: get_legendre_workspace_size
        !-----------------------------------------
    end type VhagsAux


contains


    pure function get_lvhags(nlat, nlon) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip)         :: l1, l2
        type(SpherepackAux) :: sphere_aux
        !----------------------------------------------------------------------

        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = ((nlat+1)**2)*nlat/2+nlon+15
        !return_value = max(3*nlat*(nlat+1)+2, l1*l2*(2*nlat-l1+1)+nlon+15+2*nlat)

    end function get_lvhags




    pure function get_ldwork(nlat) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)  :: nlat
        integer(ip)               :: return_value
        !----------------------------------------------------------------------

        return_value = (3*nlat*(nlat+3)+2)/2

    end function get_ldwork




    pure function get_legendre_workspace_size(nlat, nlon, nt, ityp) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip),           intent(in) :: nlat
        integer(ip),           intent(in) :: nlon
        integer(ip), optional, intent(in) :: nt
        integer(ip), optional, intent(in) :: ityp
        integer(ip)                        :: return_value
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: nt_op, ityp_op, l2
        !----------------------------------------------------------------------

        !
        !==> Address optional arguments
        !
        if (present(nt)) then
            nt_op = nt
        else
            nt_op = 1
        end if

        if (present(ityp)) then
            ityp_op = ityp
        else
            ityp_op = 0
        end if

        !
        !==> Compute workspace size
        !
        if (ityp <= 2) then
            ! Set workspace size
            return_value = max(3*nlat*(nlat+1)+2, (2*nt_op+1)*nlat*nlon)
        else
            ! Compute parity
            select case (mod(nlat, 2))
                case (0)
                    l2 = nlat/2
                case default
                    l2 = (nlat + 1)/2
            end select
            ! Set workspace size
            return_value = max(3*nlat*(nlat+1)+2, (2*nt_op+1)*l2*nlon)
        end if


    end function get_legendre_workspace_size




    subroutine vhags(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
        mdab, ndab, wvhags, lvhags, work, lwork, ierror)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: ityp
        integer(ip), intent(in)  :: nt
        real(wp),    intent(in)  :: v(idvw, jdvw, nt)
        real(wp),    intent(in)  :: w(idvw, jdvw, nt)
        integer(ip), intent(in)  :: idvw
        integer(ip), intent(in)  :: jdvw
        real(wp),    intent(out) :: br(mdab,ndab,nt)
        real(wp),    intent(out) :: bi(mdab, ndab,nt)
        real(wp),    intent(out) :: cr(mdab,ndab,nt)
        real(wp),    intent(out) :: ci(mdab, ndab,nt)
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        real(wp),    intent(in)  :: wvhags(lvhags)
        integer(ip), intent(in)  :: lvhags
        real(wp),    intent(out) :: work(lwork)
        integer(ip), intent(in)  :: lwork
        integer(ip), intent(out) :: ierror
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: idv, idz, imid, ist
        integer(ip) :: lnl, lzimn, mmax
        integer(ip) :: workspace_indices(7)
        !----------------------------------------------------------------------

        mmax = min(nlat, (nlon+1)/2)
        idz = (mmax*(2*nlat-mmax+1))/2
        imid = (nlat+1)/2
        lzimn = idz*imid

        if (ityp <= 2) then
            ist = imid
            idv = nlat
        else
            ist = 0
            idv = imid
        end if

        lnl = nt*idv*nlon

        !
        !==> Check validity of input arguments
        !
        if (nlat < 3) then
            ierror = 1
            return
        else if (nlon < 1) then
            ierror = 2
            return
        else if (ityp < 0 .or. ityp > 8) then
            ierror = 3
            return
        else if (nt < 0) then
            ierror = 4
            return
        else if ( &
            (ityp <= 2 .and. idvw < nlat) &
            .or. &
            (ityp > 2 .and. idvw < imid) &
            ) then
            ierror = 5
            return
        else if (jdvw < nlon) then
            ierror = 6
            return
        else if (mdab < mmax) then
            ierror = 7
            return
        else if (ndab < nlat) then
            ierror = 8
            return
        else if (lvhags < 2*lzimn+nlon+15) then
            ierror = 9
            return
        else if (lwork < 2*lnl+idv*nlon) then
            ierror = 10
            return
        else
            ierror = 0
        end if

        !
        !==> Compute workspace pointers
        !
        workspace_indices = get_workspace_indices(nlat, imid, ist, lnl)

        associate( &
            jw1 => workspace_indices(1), &
            jw2 => workspace_indices(2), &
            jw3 => workspace_indices(3), &
            iw1 => workspace_indices(4), &
            iw2 => workspace_indices(5), &
            iw3 => workspace_indices(6), &
            iw4 => workspace_indices(7) &
            )

            call vhags1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
                br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
                work(iw4), idz, wvhags(jw1), wvhags(jw2), wvhags(jw3))

        end associate


    contains


        pure function get_workspace_indices(nlat, imid, ist, lnl) result (return_value)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: imid
            integer(ip), intent(in)  :: ist
            integer(ip), intent(in)  :: lnl
            integer(ip)               :: return_value(7)
            !----------------------------------------------------------------------
            ! Local variables
            !----------------------------------------------------------------------
            integer(ip) :: lmn
            !----------------------------------------------------------------------

            associate( i => return_value )
                !
                !==> set wvhags pointers
                !
                lmn = nlat*(nlat+1)/2
                i(1) = 1
                i(2) = i(1)+imid*lmn
                i(3) = i(2)+imid*lmn
                !
                !==> set work pointers
                !
                i(4) = ist+1
                i(5) = lnl+1
                i(6) = i(5)+ist
                i(7) = i(5)+lnl

            end associate

        end function get_workspace_indices



        subroutine vhags1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
            ndab, br, bi, cr, ci, idv, ve, vo, we, wo, work, idz, vb, wb, wrfft)

            real(wp) :: bi
            real(wp) :: br
            real(wp) :: ci
            real(wp) :: cr
            real(wp) :: fsn
            integer(ip) :: i
            integer(ip) :: idv
            integer(ip) :: idvw
            integer(ip) :: idz
            integer(ip) :: imid
            integer(ip) :: imm1
            integer(ip) :: ityp
            integer(ip) :: j
            integer(ip) :: jdvw
            integer(ip) :: k
            integer(ip) :: m
            integer(ip) :: mb
            integer(ip) :: mdab
            integer(ip) :: mlat
            integer(ip) :: mlon
            integer(ip) :: mmax
            integer(ip) :: mp1
            integer(ip) :: mp2
            integer(ip) :: ndab
            integer(ip) :: ndo1
            integer(ip) :: ndo2
            integer(ip) :: nlat
            integer(ip) :: nlon
            integer(ip) :: nlp1
            integer(ip) :: np1
            integer(ip) :: nt
            real(wp) :: tsn
            real(wp) :: v
            real(wp) :: vb
            real(wp) :: ve
            real(wp) :: vo
            real(wp) :: w
            real(wp) :: wb
            real(wp) :: we
            real(wp) :: wo
            real(wp) :: work
            real(wp) :: wrfft
            dimension v(idvw, jdvw, *), w(idvw, jdvw, *), br(mdab, ndab,*), &
                bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *), &
                ve(idv, nlon, *), vo(idv, nlon, *), we(idv, nlon, *), &
                wo(idv, nlon, *), work(*), &
                vb(imid,*), wb(imid,*), wrfft(*)

            type(HFFTpack)      :: hfft

            nlp1 = nlat+1
            tsn = 2.0_wp/nlon
            fsn = 4.0_wp/nlon
            mlat = mod(nlat, 2)
            mlon = mod(nlon, 2)
            mmax = min(nlat, (nlon+1)/2)

            select case (mlat)
                case (0)
                    imm1 = imid
                    ndo1 = nlat
                    ndo2 = nlat-1
                case default
                    imm1 = imid-1
                    ndo1 = nlat-1
                    ndo2 = nlat
            end select

            if (ityp <= 2) then
                do k=1, nt
                    do i=1, imm1
                        do j=1, nlon
                            ve(i, j, k) = tsn*(v(i, j, k)+v(nlp1-i, j, k))
                            vo(i, j, k) = tsn*(v(i, j, k)-v(nlp1-i, j, k))
                            we(i, j, k) = tsn*(w(i, j, k)+w(nlp1-i, j, k))
                            wo(i, j, k) = tsn*(w(i, j, k)-w(nlp1-i, j, k))
                        end do
                    end do
                end do
            else
                do k=1, nt
                    do i=1, imm1
                        do j=1, nlon
                            ve(i, j, k) = fsn*v(i, j, k)
                            vo(i, j, k) = fsn*v(i, j, k)
                            we(i, j, k) = fsn*w(i, j, k)
                            wo(i, j, k) = fsn*w(i, j, k)
                        end do
                    end do
                end do
            end if

            if (mlat /= 0) then
                do k=1, nt
                    do j=1, nlon
                        ve(imid, j, k) = tsn*v(imid, j, k)
                        we(imid, j, k) = tsn*w(imid, j, k)
                    end do
                end do
            end if

            do k=1, nt
                call hfft%forward(idv, nlon, ve(1, 1, k), idv, wrfft, work)
                call hfft%forward(idv, nlon, we(1, 1, k), idv, wrfft, work)
            end do

            !
            !==> Set polar coefficients to zero
            !
            select case (ityp)
                case (0:1,3:4,6:7)
                    do k=1, nt
                        do mp1=1, mmax
                            do np1=mp1, nlat
                                br(mp1, np1, k) = 0.0_wp
                                bi(mp1, np1, k) = 0.0_wp
                            end do
                        end do
                    end do
            end select

            !
            !==> Set azimuthal coefficients to zero
            !
            select case (ityp)
                case (0,2:3,5:6,8)
                    do k=1, nt
                        do mp1=1, mmax
                            do np1=mp1, nlat
                                cr(mp1, np1, k) = 0.0_wp
                                ci(mp1, np1, k) = 0.0_wp
                            end do
                        end do
                    end do
            end select

            select case (ityp)
                case (0)
                    !
                    !==> case ityp=0 ,  no symmetries
                    !
                    !    case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*ve(i, 1, k)
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*we(i, 1, k)
                            end do
                        end do
                    end do

                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*vo(i, 1, k)
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*wo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 <= ndo1) then

                            do k=1, nt
                                do i=1, imm1
                                    do np1=mp1, ndo1, 2
                                        br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-2, k) &
                                            +wb(i, np1+mb)*we(i, 2*mp1-1, k)
                                        bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-1, k) &
                                            -wb(i, np1+mb)*we(i, 2*mp1-2, k)
                                        cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-2, k) &
                                            +wb(i, np1+mb)*ve(i, 2*mp1-1, k)
                                        ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-1, k) &
                                            -wb(i, np1+mb)*ve(i, 2*mp1-2, k)
                                    end do
                                end do
                            end do

                            if (mlat /= 0) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        br(mp1, np1, k) = br(mp1, np1, k)+wb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                                        bi(mp1, np1, k) = bi(mp1, np1, k)-wb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                        cr(mp1, np1, k) = cr(mp1, np1, k)+wb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                                        ci(mp1, np1, k) = ci(mp1, np1, k)-wb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                    end do
                                end do

                            end if
                        end if

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*wo(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*wo(i, 2*mp1-2, k)
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*vo(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*vo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
                case (1)
                    !
                    !==> case ityp=1 ,  no symmetries but cr and ci equal zero
                    !
                    !    case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*ve(i, 1, k)
                            end do
                        end do
                    end do

                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*vo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 <= ndo1) then

                            do k=1, nt
                                do i=1, imm1
                                    do np1=mp1, ndo1, 2
                                        br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-2, k) &
                                            +wb(i, np1+mb)*we(i, 2*mp1-1, k)
                                        bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-1, k) &
                                            -wb(i, np1+mb)*we(i, 2*mp1-2, k)
                                    end do
                                end do
                            end do

                            if (mlat /= 0) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        br(mp1, np1, k) = br(mp1, np1, k)+wb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                                        bi(mp1, np1, k) = bi(mp1, np1, k)-wb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                    end do
                                end do
                            end if
                        end if

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*wo(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*wo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
                case (2)
                    !
                    !==> case ityp=2 ,  no symmetries but br and bi equal zero
                    !
                    !    case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*we(i, 1, k)
                            end do
                        end do
                    end do

                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*wo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 <= ndo1) then

                            do k=1, nt
                                do i=1, imm1
                                    do np1=mp1, ndo1, 2
                                        cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-2, k) &
                                            +wb(i, np1+mb)*ve(i, 2*mp1-1, k)
                                        ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-1, k) &
                                            -wb(i, np1+mb)*ve(i, 2*mp1-2, k)
                                    end do
                                end do
                            end do

                            if (mlat /= 0) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        cr(mp1, np1, k) = cr(mp1, np1, k)+wb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                                        ci(mp1, np1, k) = ci(mp1, np1, k)-wb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                    end do
                                end do
                            end if
                        end if

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*vo(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*vo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
                case (3)
                    !
                    !==> case ityp=3 ,  v even , w odd
                    !
                    !    case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*ve(i, 1, k)
                            end do
                        end do
                    end do

                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*wo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 <= ndo1) then

                            do k=1, nt
                                do i=1, imm1
                                    do np1=mp1, ndo1, 2
                                        cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-2, k) &
                                            +wb(i, np1+mb)*ve(i, 2*mp1-1, k)
                                        ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-1, k) &
                                            -wb(i, np1+mb)*ve(i, 2*mp1-2, k)
                                    end do
                                end do
                            end do

                            if (mlat /= 0) then

                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        cr(mp1, np1, k) = cr(mp1, np1, k)+wb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                                        ci(mp1, np1, k) = ci(mp1, np1, k)-wb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                    end do
                                end do
                            end if
                        end if

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*wo(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*wo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
                case (4)
                    !
                    !==> case ityp=4 ,  v even, w odd, and cr and ci equal 0.
                    !
                    !     case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*ve(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*wo(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*wo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
                case (5)
                    !
                    !==> case ityp=5   v even, w odd, and br and bi equal zero
                    !
                    !     case m=0
                    !
                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*wo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 > ndo1) exit

                        do  k=1, nt
                            do  i=1, imm1
                                do  np1=mp1, ndo1, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*ve(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*ve(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp1, ndo1, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)+wb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-wb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                            end do
                        end do
                    end do
                case (6)
                    !
                    !==> case ityp=6 ,  v odd , w even
                    !
                    !     case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*we(i, 1, k)
                            end do
                        end do
                    end do

                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*vo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 <= ndo1) then

                            do k=1, nt
                                do i=1, imm1
                                    do np1=mp1, ndo1, 2
                                        br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-2, k) &
                                            +wb(i, np1+mb)*we(i, 2*mp1-1, k)
                                        bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-1, k) &
                                            -wb(i, np1+mb)*we(i, 2*mp1-2, k)
                                    end do
                                end do
                            end do


                            if (mlat /= 0) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        br(mp1, np1, k) = br(mp1, np1, k)+wb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                                        bi(mp1, np1, k) = bi(mp1, np1, k)-wb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                    end do
                                end do
                            end if

                        end if

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*vo(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*vo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
                case (7)
                    !
                    !==> case ityp=7   v odd, w even, and cr and ci equal zero
                    !
                    !    case m=0
                    !
                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*vo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 > ndo1) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, ndo1, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*we(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*we(i, 2*mp1-2, k)
                                end do
                            end do
                        end do


                        if (mlat == 0) exit

                        do  k=1, nt
                            do np1=mp1, ndo1, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+wb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)-wb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                            end do
                        end do
                    end do
                case (8)
                    !
                    !==> case ityp=8   v odd, w even, and both br and bi equal zero
                    !
                    !     case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*we(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*vo(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*vo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do  k=1, nt
                            do  np1=mp2, ndo2, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
            end select

        end subroutine vhags1

    end subroutine vhags



    subroutine vhagsi(nlat, nlon, wvhags, lvhags, dwork, ldwork, ierror)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: wvhags(lvhags)
        integer(ip), intent(in)  :: lvhags
        real(wp),    intent(out) :: dwork(ldwork)
        integer(ip), intent(in)  :: ldwork
        integer(ip), intent(out) :: ierror
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip)    :: imid, lmn
        integer(ip)    :: workspace_indices(7)
        type(HFFTpack) :: hfft
        !----------------------------------------------------------------------

        imid = (nlat+1)/2
        lmn = (nlat*(nlat+1))/2

        !
        !==> Check validity of input arguments
        !
        if (nlat < 3) then
            ierror = 1
            return
        else if (nlon < 1) then
            ierror = 2
            return
        else if (lvhags < 2*(imid*lmn)+nlon+15) then
            ierror = 3
            return
        else if (ldwork < (nlat*(3*nlat+9)+2)/2) then
            ierror = 4
            return
        else
            ierror = 0
        end if

        !
        !==> Compute workspace indices
        !
        workspace_indices = get_workspace_indices(nlat, imid, lmn)

        associate( &
            jw1 => workspace_indices(1), &
            jw2 => workspace_indices(2), &
            jw3 => workspace_indices(3), &
            iw1 => workspace_indices(4), &
            iw2 => workspace_indices(5), &
            iw3 => workspace_indices(6), &
            iw4 => workspace_indices(7) &
            )

            call vhgai1(nlat, imid, wvhags(jw1), wvhags(jw2), &
                dwork(iw1), dwork(iw2), dwork(iw3), dwork(iw4))

            call hfft%initialize(nlon, wvhags(jw3))

        end associate

    contains

        pure function get_workspace_indices(nlat, imid, lmn) result (return_value)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: imid
            integer(ip), intent(in)  :: lmn
            integer(ip)               :: return_value(7)
            !----------------------------------------------------------------------

            associate( i => return_value )
                !
                !==> set pointers
                !
                i(1) = 1
                i(2) = i(1)+imid*lmn
                i(3) = i(2)+imid*lmn
                i(4) = 1
                i(5) = i(4)+nlat
                i(6) = i(5)+nlat
                i(7) = i(6)+3*imid*nlat

            end associate

        end function get_workspace_indices


        subroutine vhgai1(nlat, imid, vb, wb, dthet, dwts, dpbar, work)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: imid
            real(wp),    intent(out) :: vb(imid, *)
            real(wp),    intent(out) :: wb(imid, *)
            real(wp),    intent(out) :: dthet(nlat)
            real(wp),    intent(out) :: dwts(nlat)
            real(wp),    intent(out) :: dpbar(imid, nlat, 3)
            real(wp),    intent(out) :: work(*)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer(ip)         :: i, local_error_flag, id, ix, iy
            integer(ip)         :: m, mn, n, nm, np, nz
            real(wp)            :: abel, bbel, cbel, dcf
            integer(ip)         :: dummy_integer
            real(wp)            :: dummy_real
            type(SpherepackAux) :: sphere_aux
            !----------------------------------------------------------------------

            !
            !==> Compute gaussian grid
            !
            call compute_gaussian_latitudes_and_weights(nlat, dthet, dwts, dummy_real, dummy_integer, local_error_flag)

            !
            !==> Compute associated legendre functions
            !
            !    Set m=n=0 legendre polynomials for all theta(i)
            !
            dpbar(:,1,1) = cos(PI/4)
            vb(:,1) = 0.0_wp
            wb(:,1) = 0.0_wp

            !
            !==> main loop for remaining vb, and wb
            !
            do n=1, nlat-1
                nm = mod(n-2, 3)+1
                nz = mod(n-1, 3)+1
                np = mod(n, 3)+1
                !
                !==> compute dpbar for m=0
                !
                call sphere_aux%dnlfk(0, n, work)
                mn = indx(0, n, nlat)
                do i=1, imid
                    call sphere_aux%dnlft(0, n, dthet(i), work, dpbar(i, 1, np))
                end do
                !
                !==> compute dpbar for m=1
                !
                call sphere_aux%dnlfk(1, n, work)

                mn = indx(1, n, nlat)

                do i=1, imid
                    call sphere_aux%dnlft(1, n, dthet(i), work, dpbar(i, 2, np))
                end do

                !
                !==> compute and store dpbar for m=2, n
                !
                if (n >= 2) then
                    do m=2, n
                        abel = sqrt(real((2*n+1)*(m+n-2)*(m+n-3))/ &
                            real((2*n-3)*(m+n-1)*(m+n)))
                        bbel = sqrt(real((2*n+1)*(n-m-1)*(n-m))/ &
                            real((2*n-3)*(m+n-1)*(m+n)))
                        cbel = sqrt(real((n-m+1)*(n-m+2))/ &
                            real((m+n-1)*(m+n)))
                        id = indx(m, n, nlat)

                        if (m >= n-1) then
                            dpbar(1:imid, m+1, np) = &
                                abel*dpbar(1:imid, m-1, nm)-cbel*dpbar(1:imid, m-1, np)
                        else
                            dpbar(1:imid, m+1, np) = &
                                abel*dpbar(1:imid, m-1, nm)+bbel*dpbar(1:imid, m+1, nm) &
                                -cbel*dpbar(1:imid, m-1, np)
                        end if
                    end do
                end if
                !
                !     compute the derivative of the functions
                !
                ix = indx(0, n, nlat)
                iy = indx(n, n, nlat)
                vb(1:imid, ix) = -dpbar(1:imid, 2, np)*dwts(1:imid)
                vb(1:imid, iy) = dpbar(1:imid, n, np)/sqrt(real(2*(n+1), kind=wp))*dwts(1:imid)

                if (n==1) then
                    !
                    !==> compute the vector harmonic w(theta) = m*pbar/cos(theta)
                    !
                    !     set wb=0 for m=0
                    !
                    ix = indx(0, n, nlat)
                    wb(1:imid, ix) = 0.0_wp
                else
                    dcf = sqrt(real(4*n*(n+1), kind=wp))
                    do m=1, n-1
                        ix = indx(m, n, nlat)
                        abel = sqrt(real((n+m)*(n-m+1), kind=wp))/dcf
                        bbel = sqrt(real((n-m)*(n+m+1), kind=wp))/dcf
                        vb(1:imid, ix) = &
                            (abel*dpbar(1:imid, m, np)-bbel*dpbar(1:imid, m+2, np))&
                            * dwts(1:imid)
                    end do
                end if
                !
                !==> compute wb for m=1, n
                !
                dcf = sqrt(real(2*n+1, kind=wp)/real(4*n*(n+1)*(n+n-1), kind=wp))
                do m=1, n
                    ix = indx(m, n, nlat)
                    abel = dcf*sqrt(real((n+m)*(n+m-1), kind=wp))
                    bbel = dcf*sqrt(real((n-m)*(n-m-1), kind=wp))
                    if (m >= n-1) then
                        wb(1:imid, ix) = abel * dpbar(1:imid, m, nz) * dwts(1:imid)
                    else
                        wb(1:imid, ix) = &
                            (abel*dpbar(1:imid, m, nz) + bbel*dpbar(1:imid, m+2, nz))&
                            * dwts(1:imid)
                    end if
                end do
            end do

        end subroutine vhgai1



        pure function indx(m, n, nlat) result (return_value)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: nlat
            integer              :: return_value
            !----------------------------------------------------------------------

            return_value = m*nlat-(m*(m+1))/2+n+1

        end function indx

    end subroutine vhagsi

end module module_vhags
