!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                      SPHEREPACK                               *
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
! ... file vhaec.f
!
!     this file contains code and documentation for subroutines
!     vhaec and vhaeci
!
! ... files which must be loaded with vhaec.f
!
!     type_SpherepackAux.f, type_RealPeriodicTransform.f
!
!                                                                              
!     subroutine vhaec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, 
!    +                 mdab, ndab, wvhaec, lvhaec, work, lwork, ierror)
!
!     subroutine vhaec performs the vector spherical harmonic analysis
!     on the vector field (v, w) and stores the result in the arrays
!     br, bi, cr, and ci. v(i, j) and w(i, j) are the colatitudinal 
!     (measured from the north pole) and east longitudinal components
!     respectively, located at colatitude theta(i) = (i-1)*pi/(nlat-1)
!     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
!     representation of (v, w) is given at output parameters v, w in 
!     subroutine vhsec.  
!
!     input parameters
!
!     nlat   the number of colatitudes on the full sphere including the
!            poles. for example, nlat = 37 for a five degree grid.
!            nlat determines the grid increment in colatitude as
!            pi/(nlat-1).  if nlat is odd the equator is located at
!            grid point i=(nlat+1)/2. if nlat is even the equator is
!            located half way between points i=nlat/2 and i=nlat/2+1.
!            nlat must be at least 3. note: on the half sphere, the
!            number of grid points in the colatitudinal direction is
!            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
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
!     nt     the number of analyses.  in the program that calls vhaec, 
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
!            components at colatitude theta(i) = (i-1)*pi/(nlat-1)
!            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges
!            are defined above at the input parameter ityp.
!
!     idvw   the first dimension of the arrays v, w as it appears in
!            the program that calls vhaec. if ityp .le. 2 then idvw
!            must be at least nlat.  if ityp .gt. 2 and nlat is
!            even then idvw must be at least nlat/2. if ityp .gt. 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls vhaec. jdvw must be at least nlon.
!
!     mdab   the first dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhaec. mdab must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhaec. ndab must be at
!            least nlat.
!
!     wvhaec an array which must be initialized by subroutine vhaeci.
!            once initialized, wvhaec can be used repeatedly by vhaec
!            as long as nlon and nlat remain unchanged.  wvhaec must
!            not be altered between calls of vhaec.
!
!     lvhaec the dimension of the array wvhaec as it appears in the
!            program that calls vhaec. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhaec must be at least
!
!            4*nlat*l2+3*max(l1-2, 0)*(nlat+nlat-l1-1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhaec. define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if ityp .le. 2 then lwork must be at least
!
!                    nlat*(2*nt*nlon+max(6*l2, nlon))
!
!            if ityp .gt. 2 then lwork must be at least
!
!                    l2*(2*nt*nlon+max(6*nlat, nlon))
!
!     **************************************************************
!
!     output parameters
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!     cr, ci  that contain the vector spherical harmonic coefficients
!            in the spectral representation of v(i, j) and w(i, j) given 
!            in the discription of subroutine vhsec. br(mp1, np1), 
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
!            = 9  error in the specification of lvhaec
!            = 10 error in the specification of lwork
!
!
! *******************************************************************
!
!     subroutine vhaeci(nlat, nlon, wvhaec, lvhaec, dwork, ldwork, ierror)
!
!     subroutine vhaeci initializes the array wvhaec which can then be
!     used repeatedly by subroutine vhaec until nlat or nlon is changed.
!
!     input parameters
!
!     nlat   the number of colatitudes on the full sphere including the
!            poles. for example, nlat = 37 for a five degree grid.
!            nlat determines the grid increment in colatitude as
!            pi/(nlat-1).  if nlat is odd the equator is located at
!            grid point i=(nlat+1)/2. if nlat is even the equator is
!            located half way between points i=nlat/2 and i=nlat/2+1.
!            nlat must be at least 3. note: on the half sphere, the
!            number of grid points in the colatitudinal direction is
!            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     lvhaec the dimension of the array wvhaec as it appears in the
!            program that calls vhaec. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhaec must be at least
!
!            4*nlat*l2+3*max(l1-2, 0)*(nlat+nlat-l1-1)+nlon+15
!
!
!     dwork  a real work array that does not have to be saved.
!
!     ldwork the dimension of the array dwork as it appears in the
!            program that calls vhaec. ldwork must be at least
!            2*(nlat+2)
!
!
!     **************************************************************
!
!     output parameters
!
!     wvhaec an array which is initialized for use by subroutine vhaec.
!            once initialized, wvhaec can be used repeatedly by vhaec
!            as long as nlat or nlon remain unchanged.  wvhaec must not
!            be altered between calls of vhaec.
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lvhaec
!            = 4  error in the specification of ldwork
!
!
module vector_analysis_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI

    use type_SpherepackAux, only: &
        SpherepackAux

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: vhaec, vhaes, vhagc, vhags
    public :: vhaeci, vhaesi, vhagci, vhagsi

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: FOUR = 4.0_wp

    ! Declare interfaces for submodule implementation
    interface
        module subroutine vhaec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhaec, lvhaec, work, lwork, ierror)

            ! Dummy arguments
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
            real(wp),    intent(in)  :: wvhaec(lvhaec)
            integer(ip), intent(in)  :: lvhaec
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vhaec

        module subroutine vhaes(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhaes, lvhaes, work, lwork, ierror)

            ! Dummy arguments
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
            real(wp),    intent(in)  :: wvhaes(lvhaes)
            integer(ip), intent(in)  :: lvhaes
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vhaes

        module subroutine vhagc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhagc, lvhagc, work, lwork, ierror)

            ! Dummy arguments
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
            real(wp),    intent(in)  :: wvhagc(lvhagc)
            integer(ip), intent(in)  :: lvhagc
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vhagc

        module subroutine vhags(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhags, lvhags, work, lwork, ierror)

            ! Dummy arguments
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
        end subroutine vhags

        module subroutine vhaeci(nlat, nlon, wvhaec, lvhaec, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhaec(lvhaec)
            integer(ip), intent(in)  :: lvhaec
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine vhaeci

        module subroutine vhaesi(nlat, nlon, wvhaes, lvhaes, work, lwork, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhaes(lvhaes)
            integer(ip), intent(in)  :: lvhaes
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine vhaesi

        module subroutine vhagci(nlat, nlon, wvhagc, lvhagc, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhagc(lvhagc)
            integer(ip), intent(in)  :: lvhagc
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine vhagci

        module subroutine vhagsi(nlat, nlon, wvhags, lvhags, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhags(lvhags)
            integer(ip), intent(in)  :: lvhags
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine vhagsi
    end interface

    type, public :: VhaesAux
    contains
        ! Type-bound procedures
        procedure, nopass :: vhaes
        procedure, nopass :: vhaesi
        procedure, nopass :: get_lvhaes
        procedure, nopass :: get_lwork => get_vhaes_lwork
        procedure, nopass :: get_ldwork => get_vhaes_ldwork
        procedure, nopass :: get_legendre_workspace_size => &
            get_vhaes_legendre_workspace_size
    end type VhaesAux

    type, public :: VhagsAux
    contains
        ! Type-bound procedures
        procedure, nopass :: vhags
        procedure, nopass :: vhagsi
        procedure, nopass :: get_lvhags
        procedure, nopass :: get_ldwork
        procedure, nopass :: get_legendre_workspace_size
    end type VhagsAux

contains

    pure function get_lvhags(nlat, nlon) result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value

        ! Local variables
        integer(ip)         :: l1, l2
        type(SpherepackAux) :: sphere_aux

        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = ((nlat+1)**2)*nlat/2+nlon+15
        !return_value = max(3*nlat*(nlat+1)+2, l1*l2*(2*nlat-l1+1)+nlon+15+2*nlat)

    end function get_lvhags

    pure function get_ldwork(nlat) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip)               :: return_value


        return_value = (3*nlat*(nlat+3)+2)/2

    end function get_ldwork

    pure function get_legendre_workspace_size(nlat, nlon, nt, ityp) result (return_value)

        ! Dummy arguments

        integer(ip),           intent(in) :: nlat
        integer(ip),           intent(in) :: nlon
        integer(ip), optional, intent(in) :: nt
        integer(ip), optional, intent(in) :: ityp
        integer(ip)                        :: return_value

        ! Local variables

        integer(ip) :: nt_op, ityp_op, l2


        !
        !  Address optional arguments
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
        !  Compute workspace size
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

    pure function get_lvhaes(nlat, nlon) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value

        ! Local variables

        integer(ip)         :: l1, l2
        type(SpherepackAux) :: sphere_aux


        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = l1*l2*(2*nlat-l1+1)+nlon+15

    end function get_lvhaes

    pure function get_vhaes_lwork(nlat, nlon) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value

        ! Local variables

        integer(ip)         :: l1, l2
        type(SpherepackAux) :: sphere_aux


        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = 3*(max(l1-2, 0)*(2*nlat-l1-1))/2+5*l2*nlat

    end function get_vhaes_lwork

    pure function get_vhaes_ldwork(nlat) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip)               :: return_value


        return_value = 2*(nlat+1)

    end function get_vhaes_ldwork

    pure function get_vhaes_legendre_workspace_size(nlat, nlon, nt, ityp) result (return_value)

        ! Dummy arguments

        integer(ip),           intent(in) :: nlat
        integer(ip),           intent(in) :: nlon
        integer(ip), optional, intent(in) :: nt
        integer(ip), optional, intent(in) :: ityp
        integer(ip)                        :: return_value

        ! Local variables

        integer(ip) :: nt_op, ityp_op, l2


        !
        !  Address optional arguments
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
        !  Compute workspace size
        !
        if (ityp_op <= 2) then
            ! Set workspace size
            return_value = (2*nt_op+1)*nlat*nlon
        else
            ! Compute parity
            select case (mod(nlat, 2))
                case (0)
                    l2 = nlat/2
                case default
                    l2 = (nlat + 1)/2
            end select
            ! Set workspace size
            return_value = (2*nt_op+1)*l2*nlon
        end if

    end function get_vhaes_legendre_workspace_size

end module vector_analysis_routines
