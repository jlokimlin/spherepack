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
! ... file divec.f
!
!     this file includes documentation and code for
!     subroutine divec          i
!
! ... files which must be loaded with divec.f
!
!     type_SpherepackAux.f, type_HFFTpack.f, vhaec.f, shsec.f
!
!
!     subroutine divec(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, 
!    +                 wshsec, lshsec, work, lwork, ierror)
!
!     given the vector spherical harmonic coefficients br and bi, precomputed
!     by subroutine vhaec for a vector field (v, w), subroutine divec
!     computes the divergence of the vector field in the scalar array dv.
!     dv(i, j) is the divergence at the colatitude
!
!            theta(i) = (i-1)*pi/(nlat-1)
!
!     and east longitude
!
!            lambda(j) = (j-1)*2*pi/nlon
!
!     on the sphere.  i.e.
!
!            dv(i, j) = 1/sint*[ d(sint*v(i, j))/dtheta + d(w(i, j))/dlambda ]
!
!     where sint = sin(theta(i)).  w is the east longitudinal and v
!     is the colatitudinal component of the vector field from which
!     br, bi were precomputed.  required associated legendre polynomials
!     are recomputed rather than stored as they are in subroutine dives.
!
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
!            nlon = 72 for a five degree grid. nlon must be greater than
!            3.  the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!
!     isym   a parameter which determines whether the divergence is
!            computed on the full or half sphere as follows:
!
!      = 0
!
!            the symmetries/antsymmetries described in isym=1, 2 below
!            do not exist in (v, w) about the equator.  in this case the
!            divergence is neither symmetric nor antisymmetric about
!            the equator.  the divergence is computed on the entire
!            sphere.  i.e., in the array dv(i, j) for i=1, ..., nlat and
!            j=1, ..., nlon.
!
!      = 1
!
!            w is antisymmetric and v is symmetric about the equator.
!            in this case the divergence is antisymmetyric about
!            the equator and is computed for the northern hemisphere
!            only.  i.e., if nlat is odd the divergence is computed
!            in the array dv(i, j) for i=1, ..., (nlat+1)/2 and for
!            j=1, ..., nlon.  if nlat is even the divergence is computed
!            in the array dv(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!      = 2
!            w is symmetric and v is antisymmetric about the equator
!            in this case the divergence is symmetyric about the
!            equator and is computed for the northern hemisphere
!            only.  i.e., if nlat is odd the divergence is computed
!            in the array dv(i, j) for i=1, ..., (nlat+1)/2 and for
!            j=1, ..., nlon.  if nlat is even the divergence is computed
!            in the array dv(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!
!     nt     nt is the number of scalar and vector fields.  some
!            computational efficiency is obtained for multiple fields.
!            can be three dimensional corresponding to an indexed multiple
!            vector field.  in this case multiple scalar synthesis will
!            be performed to compute the divergence for each field.  the
!            third index is the synthesis index which assumes the values
!            k=1, ..., nt.  for a single synthesis set nt = 1.  the
!            description of the remaining parameters is simplified by
!            assuming that nt=1 or that all the arrays are two dimensional.
!
!     idv    the first dimension of the array dv as it appears in
!            the program that calls divec. if isym = 0 then idv
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then idv must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then idv must be at least (nlat+1)/2.
!
!     jdv    the second dimension of the array dv as it appears in
!            the program that calls divec. jdv must be at least nlon.
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!            that contain vector spherical harmonic coefficients
!            of the vector field (v, w) as computed by subroutine vhaec.
!     ***    br and bi must be computed by vhaec prior to calling
!            divec.
!
!     mdb    the first dimension of the arrays br and bi as it
!            appears in the program that calls divec. mdb must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndb    the second dimension of the arrays br and bi as it
!            appears in the program that calls divec. ndb must be at
!            least nlat.
!
!
!     wshsec an array which must be initialized by subroutine shseci.
!            once initialized, 
!            wshsec can be used repeatedly by divec as long as nlon
!            and nlat remain unchanged.  wshsec must not be altered
!            between calls of divec.
!
!
!     lshsec the dimension of the array wshsec as it appears in the
!            program that calls divec. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshsec must be at least
!
!            2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls divec. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2                    if nlat is even or
!               l2 = (nlat+1)/2                if nlat is odd
!
!            if isym is zero then lwork must be at least
!
!               nlat*(nt*nlon+max(3*l2, nlon)+2*nt*l1+1)
!
!            if isym is not zero then lwork must be at least
!
!               l2*(nt*nlon+max(3*nlat, nlon)) + nlat*(2*nt*l1+1)
!
!
!     **************************************************************
!
!     output parameters
!
!
!    dv     a two or three dimensional array (see input parameter nt)
!           that contains the divergence of the vector field (v, w)
!           whose coefficients br, bi where computed by subroutine
!           vhaec.  dv(i, j) is the divergence at the colatitude point
!           theta(i) = (i-1)*pi/(nlat-1) and longitude point
!           lambda(j) = (j-1)*2*pi/nlon. the index ranges are defined
!           above at the input parameter isym.
!
!
!    ierror = 0  no errors
!           = 1  error in the specification of nlat
!           = 2  error in the specification of nlon
!           = 3  error in the specification of isym
!           = 4  error in the specification of nt
!           = 5  error in the specification of idv
!           = 6  error in the specification of jdv
!           = 7  error in the specification of mdb
!           = 8  error in the specification of ndb
!           = 9  error in the specification of lshsec
!           = 10 error in the specification of lwork
!
module module_divec

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use scalar_synthesis_routines, only: &
        shsec

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: divec

contains

    subroutine divec(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
        wshsec, lshsec, work, lwork, ierror)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        integer(ip), intent(in)     :: isym
        integer(ip), intent(in)     :: nt
        real(wp),    intent(inout)  :: dv(idv, jdv,*)
        integer(ip), intent(in)     :: idv
        integer(ip), intent(in)     :: jdv
        real(wp),    intent(inout)  :: br(mdb, ndb,*)
        real(wp),    intent(inout)  :: bi(mdb, ndb,*)
        integer(ip), intent(in)     :: mdb
        integer(ip), intent(in)     :: ndb
        real(wp),    intent(inout)  :: wshsec(lshsec)
        integer(ip), intent(in)     :: lshsec
        real(wp),    intent(inout)  :: work(lwork)
        integer(ip), intent(in)     :: lwork
        integer(ip), intent(out)    :: ierror
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: l1, l2, ia, ib, mn, is, ls, mab, nln
        integer(ip) :: iwk, lwk, lzz1, labc, imid
        integer(ip) :: mmax, lwmin, lwkmin
        !----------------------------------------------------------------------

        imid = (nlat+1)/2
        mmax = min(nlat, (nlon+2)/2)
        imid = (nlat+1)/2
        lzz1 = 2*nlat*imid
        labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
        lwmin = lzz1+labc+nlon+15

        select case (isym)
            case (1:)
                ls = imid
            case default
                ls = nlat
        end select

        nln = nt*ls*nlon

        !
        !  Set first dimension for a, b (as requiried by shsec)
        !
        mab = min(nlat, nlon/2+1)
        mn = mab*nlat*nt
        l1 = min(nlat, (nlon+2)/2)
        l2 = (nlat+1)/2

        select case (isym)
            case (0)
                lwkmin =  nlat*(nt*nlon+max(3*l2, nlon)+2*nt*l1+1)
            case default
                lwkmin = l2*(nt*nlon+max(3*nlat, nlon)) + nlat*(2*nt*l1+1)
        end select

        !
        !  check validity of input parameters
        !

        ! Initialize error flag
        ierror = 0

        ! Check case 1
        if (nlat < 3) then
            ierror = 1
            return
        else if (nlon < 4) then
            ierror = 2
            return
        else if (isym < 0 .or. isym > 2) then
            ierror = 3
            return
        else if (nt < 0) then
            ierror = 4
            return
        else if ( &
            (isym == 0 .and. idv < nlat) &
            .or. &
            (isym > 0 .and. idv < imid) &
            ) then
            ierror = 5
            return
        else if (jdv < nlon) then
            ierror = 6
            return
        else if (mdb < min(nlat, (nlon+1)/2)) then
            ierror = 7
            return
        else if  (ndb < nlat) then
            ierror = 8
            return
        else if  (lshsec < lwmin) then
            ierror = 9
            return
        else if (lwork < lwkmin) then
            ierror = 10
            return
        else
            ierror = 0
        end if

        !
        !     set work space pointers
        !
        ia = 1
        ib = ia+mn
        is = ib+mn
        iwk = is+nlat
        lwk = lwork-2*mn-nlat

        call divec1(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
            work(ia), work(ib), mab, work(is), wshsec, lshsec, work(iwk), lwk, &
            ierror)

    contains

        subroutine divec1(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
            a, b, mab, sqnn, wshsec, lshsec, wk, lwk, ierror)
            real :: a
            real :: b
            real :: bi
            real :: br
            real :: dv
            real :: fn
            integer(ip) :: idv
            integer(ip) :: ierror
            integer(ip) :: isym
            integer(ip) :: jdv
            integer(ip) :: k
            integer(ip) :: lshsec
            integer(ip) :: lwk
            integer(ip) :: m
            integer(ip) :: mab
            integer(ip) :: mdb
            integer(ip) :: mmax
            integer(ip) :: n
            integer(ip) :: ndb
            integer(ip) :: nlat
            integer(ip) :: nlon
            integer(ip) :: nt
            real :: sqnn
            real :: wk
            real :: wshsec
            dimension dv(idv, jdv, nt), br(mdb, ndb, nt), bi(mdb, ndb, nt)
            dimension a(mab, nlat, nt), b(mab, nlat, nt), sqnn(nlat)
            dimension wshsec(lshsec), wk(lwk)
            !
            !     set coefficient multiplyers
            !
            do n=2, nlat
                fn = real(n - 1, kind=wp)
                sqnn(n) = sqrt(fn * (fn + 1.0_wp))
            end do
            !
            !     compute divergence scalar coefficients for each vector field
            !
            do k=1, nt
                a(1: mab, 1: nlat, k) = 0.0_wp
                b(1: mab, 1: nlat, k) = 0.0_wp
                !
                !     compute m=0 coefficients
                !
                do n=2, nlat
                    a(1, n, k) = -sqnn(n)*br(1, n, k)
                    b(1, n, k) = -sqnn(n)*bi(1, n, k)
                end do
                !
                !     compute m>0 coefficients using vector spherepack value for mmax
                !
                mmax = min(nlat, (nlon+1)/2)
                do m=2, mmax
                    do n=m, nlat
                        a(m, n, k) = -sqnn(n)*br(m, n, k)
                        b(m, n, k) = -sqnn(n)*bi(m, n, k)
                    end do
                end do
            end do
            !
            !     synthesize a, b into dv
            !
            call shsec(nlat, nlon, isym, nt, dv, idv, jdv, a, b, &
                mab, nlat, wshsec, lshsec, wk, lwk, ierror)

        end subroutine divec1


    end subroutine divec

end module module_divec
