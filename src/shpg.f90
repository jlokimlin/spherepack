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
!                           August 2003
!
! ... in file shpg.f
!
!     this file contains code and documentation for subroutines
!     shpgi and shpg.
!
! ... files which must be loaded with shpg.f
!
!     type_RealPeriodicTransform.f
!
!     shpgi initializes the arrays wshp and iwshp for subsequent 
!     use in subroutine shpg, which performs the harmonic projection 
!     which is equivalent to a harmonic analysis followed by 
!     harmonic synthesis but faster and with less memory.
!     (see description of subroutine shpg below).
!
!     subroutine shpgi(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, 
!    1 liwshp, work, lwork, ierror)
!
!     shpgi initializes arrays wshp and iwshp for repeated use
!     by subroutine shpg ....
!
!     input parameters
!
!     nlat   the number of colatitudes on the full sphere including the
!            poles. for example, nlat = 37 for a five degree grid.
!            nlat determines the grid increment in colatitude as
!            pi/(nlat-1).  if nlat is odd the equator is located at
!            grid point i=(nlat+1)/2. if nlat is even the equator is
!            located half way between points i=nlat/2 and i=nlat/2+1.
!            nlat must be at least 3.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than or equal to 4. the efficiency of the computation is
!            improved when nlon is a product of small prime numbers.
!            nlon must be at least 4.
!
!     isym   currently not used, no equatorial symmetries assumed, 
!            only whole sphere computations.    
!
!     mtrunc the highest longitudinal wave number retained in the
!            projection. It must be less than or equal to
!            the minimum of nlat-1 and nlon/2. The first wave
!            number is zero. For example, if wave numbers 0 and
!            1 are desired then mtrunc = 1.
!
!     lwshp  the dimension of the array wshp as it appears in the
!            program that calls shpgi. It must be at least
!            2*(nlat+1)**2+nlon+log2(nlon)
!
!     liwshp the dimension of the array iwshp as it appears in the
!            program that calls shpgi. It must be at least
!            4*(nlat+1).
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shpgi. It must be at least
!            1.25*(nlat+1)**2+7*nlat+8.
!
!     **************************************************************
!
!     output parameters
!
!     wshp   a single precision array that must be saved for
!            repeated use by subroutine shpg.        
!
!     iwshp  an integer array that must be saved for repeated
!            use by subroutine shpg.        
!
!     work   a real work array that does 
!            not have to be saved.
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of isym
!            = 4  error in the specification of mtrunc
!            = 5  error in the specification of lwshp
!            = 6  error in the specification of liwshp
!            = 7  error in the specification of lwork
!
module module_shpg

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, &
        MACHINE_EPSILON

    use type_SpherepackAux, only: &
        SpherepackAux

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: shpg
    public :: shpgi

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: THREE = 3.0_wp

contains

    subroutine shpgi(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, &
        liwshp, work, lwork, ierror)


        type(SpherepackAux) :: sphere_aux
        integer :: ierror
        integer :: isym
        integer :: iw1
        integer :: iw2
        integer :: iw3
        integer :: iw4
        integer :: iwshp
        integer :: jw1
        integer :: jw2
        integer :: jw3
        integer :: jw4
        integer :: ktot
        integer :: kw1
        integer :: kw10
        integer :: kw11
        integer :: kw2
        integer :: kw3
        integer :: kw4
        integer :: kw5
        integer :: kw6
        integer :: kw7
        integer :: kw8
        integer :: kw9
        integer :: liwshp
        integer :: log2n
        integer :: lw1
        integer :: lwork
        integer :: lwshp
        integer :: mlwk
        integer :: mmax
        integer :: mtrunc
        integer :: nlat
        integer :: nloc1
        integer :: nloc2
        integer :: nlon
        integer :: nte
        real(wp) :: wshp
        real work(*)
        dimension wshp(*), iwshp(*)
        !
        ierror = 1
        if (nlat<1) return
        ierror = 2
        if (nlon<1) return
        !      ierror = 3
        !      if (isym.lt.0_wp .or. isym.gt.2) return
        ierror = 4
        mmax = min(nlat-1, nlon/2)
        if (mtrunc<0 .or. mtrunc>mmax) return
        ierror = 5
        lw1 = 2*(nlat+1)**2
        log2n = log(real(nlon))/log(TWO)
        if (lwshp<lw1+nlon+log2n) return
        ierror = 6
        if (liwshp<4*(nlat+1)) return
        ierror = 7
        mlwk = 1.25*(nlat+1)**2+7*nlat+8
        if (lwork <mlwk) return
        ierror = 0
        !
        call sphere_aux%hfft%initialize(nlon, wshp(lw1+1))

        nte = (nlat+1)/2
        nloc1 = 2*nte*nte
        nloc2 = nlat+1
        iw1 = 1
        iw2 = iw1+nloc1
        iw3 = iw2+nloc1
        iw4 = iw3+nloc1
        jw1 = 1
        jw2 = jw1+nloc2
        jw3 = jw2+nloc2
        jw4 = jw3+nloc2
        kw1 = 1
        kw2 = kw1+nte
        kw3 = kw2+nte
        kw4 = kw3+2*nte
        kw5 = kw4+2*nte
        kw6 = kw5+nte
        kw7 = kw6+nte
        kw8 = kw7+4*nte
        kw9 = kw8+2*nte
        kw10 = kw9+nloc1
        kw11 = kw10+nloc1
        ktot = kw11+nte*nte

        call shpgi_lower_routine(nlat, nlon, isym, mtrunc, nte, ierror, wshp(iw1), wshp(iw2), &
            wshp(iw3), wshp(iw4), iwshp(jw1), iwshp(jw2), iwshp(jw3), &
            iwshp(jw4), work(kw1), work(kw2), work(kw3), work(kw4), work(kw5), &
            work(kw6), work(kw7), work(kw8), work(kw9), work(kw10), work(kw11))

    end subroutine shpgi

    subroutine shpgi_lower_routine(nlat, nlon, isym, mtrunc, idp, ierror, &
        pe, po, ze, zo, ipse, jzse, ipso, jzso, &
        cp, wx, thet, gwts, xx, z, a, b, ped, pod, u)

        real(wp) :: dfn
        real(wp) :: dmax
        integer :: i
        integer :: idp
        integer :: ierr
        integer :: ierror
        integer :: iip
        integer :: ipse
        integer :: ipso
        integer :: isym
        integer :: it
        integer :: j
        integer :: js
        integer :: jzse
        integer :: jzso
        integer :: k
        integer :: lock
        
        integer :: m
        integer :: modn
        integer :: mp1
        integer :: ms2
        integer :: mtrunc
        integer :: mxtr
        integer :: n
        integer :: nec
        integer :: nem
        integer :: nlat
        integer :: nlon
        integer :: nmx
        integer :: noc
        integer :: nom
        integer :: ns2
        integer :: nshe
        integer :: nsho
        integer :: nte
        integer :: nto
        real(wp) :: pe
        real(wp) :: po
        real(wp) :: sum1
        real(wp) :: toe
        real(wp) :: tusl
        real(wp) :: ze
        real(wp) :: zo
        real(wp) :: zort

        real(wp) :: summation, a1, b1, c1
        real(wp) :: cp(idp), wx(idp), &
            thet(nlat), gwts(nlat), xx(idp), z(idp), a(4*idp), &
            b(2*idp), ped(idp, idp, 2), pod(idp, idp, 2), u(idp, idp)

        dimension pe(idp, idp, 2), po(idp, idp, 2), ze(idp, idp, 2), &
            zo(idp, idp, 2), &
            ipse(idp, 2), jzse(idp, 2), ipso(idp, 2), jzso(idp, 2), &
            nshe(2), nsho(2)
        dimension zort(64, 64, 2)

        type(SpherepackAux) :: sphere_aux

        ns2 = nlat/2
        modn = nlat-2*ns2
        nte = (nlat+1)/2
        nto = nlat-nte
        tusl = ZERO
        toe = ZERO

        ! Compute gauss grid distribution
        call compute_gaussian_latitudes_and_weights(nlat, thet, gwts, ierr)

        gwts(1:nto) = TWO * gwts(1:nto)

        ! Compute n**2 basis (even functions)
        do n=1, 2*nlat-2
            dfn = n
            a(n) = sqrt(dfn * (dfn + ONE))
        end do

        do n=1, nlat-1
            dfn = n
            b(n) = sqrt((TWO * dfn + THREE)/(TWO * dfn - ONE))
        end do

        mxtr = min(nlat-1, nlon/2, mtrunc)
        iip = 2
        do 200 mp1=1, mxtr+1
            m = mp1-1
            iip = 3-iip
            ms2 = mp1/2
            nem = (nlat-m+1)/2
            nec = nte-nem

            ! Compute associated legendre functions
            if (m<=1) then
                do 205 j=1, nem
                    n = 2*j+m-2
                    call sphere_aux%compute_fourier_coefficients(m, n, cp)
                    do i=1, nte
                        call sphere_aux%compute_legendre_polys_from_fourier_coeff(m, n, thet(i), cp, ped(i, j+nec, iip))
                    end do
205             continue
            !
            else
                !
                do 207 j=1, nem
                    n = 2*j+m-2
                    if (m>1.and.n>mxtr) then
                        do i=1, nte
                            u(i, j+nec) = ped(i, j+nec, iip)
                        end do
                        goto 207
                    end if
                    a1 = b(n-1)*a(n+m-3)/a(n+m-1)
                    b1 = a(n-m+1)/a(n+m-1)
                    if (n-m<=1) then
                        do i=1, nte
                            u(i, j+nec) = a1*ped(i, j+nec-1, iip) &
                                - b1*ped(i, j+nec, iip)
                        end do
                    else
                        c1 = b(n-1)*a(n-m-1)/a(n+m-1)
                        do i=1, nte
                            u(i, j+nec) = a1*ped(i, j+nec-1, iip) &
                                - b1*ped(i, j+nec, iip) + c1*u(i, j+nec-1)
                        end do
                    end if
207             continue
                do j=1, nem
                    do i=1, nte
                        ped(i, j+nec, iip) = u(i, j+nec)
                    end do
                end do
            end if
            if (nec<=0) goto 200
            !
            !     generate orthogonal vector with
            !     random numbers
            call random_seed()
            call random_number(xx(1:nte))
            !
            it = 0
            201 do i=1, nte
                z(i) = ZERO
                wx(i) = gwts(i)*xx(i)
            end do
            do 220 j=1, nte
                if (j==nec) goto 220
                call accumulate_inner_products(nte, wx, ped(1, j, iip), z)
220         continue
            !
            do i=1, nte
                xx(i) = xx(i)-z(i)
            end do
            call compute_normal_gaussian_grid(nte, xx, idp, gwts)
            it = it+1
            if (it<=2) goto 201
            do i=1, nte
                ped(i, nec, iip) = xx(i)
            end do
200     continue
        !
        !     reorder if mtrunc is less than nlat-1
        !         case of even functions
        !
        nmx = nlat-mxtr
        if (modn==1) then
            nshe(1) = nmx/2
            nshe(2) = (nmx-1)/2
        else
            nshe(1) = (nmx-1)/2
            nshe(2) = nmx/2
        end if
        !
        do 210 mp1=1, 2
            do j=1, nte
                js = j+nshe(mp1)
                if (js>nte) js = js-nte
                do i=1, nte
                    u(i, js) = ped(i, j, mp1)
                end do
            end do
            do j=1, nte
                do i=1, nte
                    ped(i, j, mp1) = u(i, j)
                end do
            end do
210     continue
        !
        call truncate(0, nte, idp, ped(1, 1, 1), nte, ipse(1, 1))
        call truncate(0, nte, idp, ped(1, 1, 2), nte, ipse(1, 2))
        !
        ! Compute the analysis matrices
        !
        do 250 iip=1, 2
            do i=1, nte
                lock = 0
                do j=1, nte
                    summation = ped(j, i, iip)*gwts(j)
                    ze(j, i, iip) =  summation
                    pe(i, j, iip) = ped(i, j, iip)
                    if (abs(summation)>MACHINE_EPSILON .and. lock==0) then
                        lock = 1
                        jzse(i, iip) = j
                    end if
                end do
            end do
250     continue
        !
        !     check orthogonality of pe(i, j, mp1)  mp1=1, 2
        !
        do iip=1, 2
            dmax = ZERO
            do i=1, nte
                do j=1, nte
                    sum1 = ZERO
                    do k=1, nte
                        sum1 = sum1+ze(k, i, iip)*pe(k, j, iip)
                    end do
                    zo(i, j, iip) = sum1
                    if (i/=j) then
                        dmax = max(dmax, abs(sum1))
                    else
                        dmax = max(dmax, abs(sum1-ONE))
                    end if
                end do
            end do
        end do
        !
        ! Compute n**2 basis (odd functions)
        !
        iip = 2
        do 300 mp1=1, mxtr+1
            iip = 3-iip
            m = mp1-1
            ms2 = mp1/2
            nem = (nlat-m+1)/2
            nom = nlat-m-nem
            noc = nto-nom
            !
            ! Compute associated legendre functions
            !
            if (m<=1) then
                do 305 j=1, nom
                    n = 2*j+m-1
                    call sphere_aux%compute_fourier_coefficients(m, n, cp)
                    do i=1, nte
                        call sphere_aux%compute_legendre_polys_from_fourier_coeff(m, n, thet(i), cp, pod(i, j+noc, iip))
                    end do
                    if (modn>0) pod(nte, j+noc, iip) = ZERO
305             continue
            !
            else
                !
                do 307 j=1, nom
                    n = 2*j+m-1
                    if (m>1.and.n>mxtr) then
                        do i=1, nte
                            u(i, j+noc) = pod(i, j+noc, iip)
                        end do
                        goto 304
                    end if
                    a1 = b(n-1)*a(n+m-3)/a(n+m-1)
                    b1 = a(n-m+1)/a(n+m-1)
                    if (n-m<=1) then
                        do i=1, nte
                            u(i, j+noc) = a1*pod(i, j+noc-1, iip) &
                                - b1*pod(i, j+noc, iip)
                        end do
                    else
                        c1 = b(n-1)*a(n-m-1)/a(n+m-1)
                        do i=1, nte
                            u(i, j+noc) = a1*pod(i, j+noc-1, iip) &
                                - b1*pod(i, j+noc, iip) + c1*u(i, j+noc-1)
                        end do
                    end if
304                 if (modn==1) u(nte, j+noc) = ZERO
307             continue
                do j=1, nom
                    do i=1, nte
                        pod(i, j+noc, iip) = u(i, j+noc)
                    end do
                end do
            end if
            !
            if (noc<=0) goto 300
            call random_number(xx(1:nte))
            if (modn==1) xx(nte) = ZERO
            it = 0
            306 do i=1, nte
                z(i) = ZERO
                wx(i) = gwts(i)*xx(i)
            end do
            do 330 j=1, nto
                if (j==noc) goto 330
                call accumulate_inner_products(nte, wx, pod(1, j, iip), z(1))
330         continue
            !
            do i=1, nte
                xx(i) = xx(i)-z(i)
            end do
            call compute_normal_gaussian_grid(nte, xx, idp, gwts)
            it = it+1
            if (it<=2) goto 306
            do i=1, nte
                pod(i, noc, iip) = xx(i)
            end do
            if (modn==1) pod(nte, noc, iip) = ZERO
300     continue
        !
        nmx = nlat-mxtr
        if (modn==1) then
            nsho(1) = (nmx-1)/2
            nsho(2) = nmx/2
        else
            nsho(1) = nmx/2
            nsho(2) = (nmx-1)/2
        end if
        !
        do 310 mp1=1, 2
            do j=1, nto
                js = j+nsho(mp1)
                if (js>nto) js = js-nto
                do i=1, nte
                    u(i, js) = pod(i, j, mp1)
                end do
            end do
            do j=1, nto
                do i=1, nte
                    pod(i, j, mp1) = u(i, j)
                end do
            end do
310     continue
        !
        call truncate(0, nte, idp, pod(1, 1, 1), nto, ipso(1, 1))
        call truncate(0, nte, idp, pod(1, 1, 2), nto, ipso(1, 2))
        !
        ! Compute the analysis matrices (odd functions)
        !
        do iip=1, 2
            do i=1, nto
                lock = 0
                do j=1, nto
                    summation = pod(j, i, iip)*gwts(j)
                    zo(j, i, iip) = summation
                    po(i, j, iip) = pod(i, j, iip)
                    if (abs(summation)>MACHINE_EPSILON .and. lock==0) then
                        lock = 1
                        jzso(i, iip) = j
                    end if
                end do
            end do
        end do
        !
        !     check orthogonality of po(i, j, mp1)  mp1=1, 2
        !
        do iip=1, 2
            dmax = ZERO
            do i=1, nto
                do j=1, nto
                    sum1 = ZERO
                    do k=1, nto
                        sum1 = sum1+zo(k, i, iip)*po(k, j, iip)
                    end do
                    zort(i, j, iip) = sum1
                    if (i/=j) then
                        dmax = max(dmax, abs(sum1))
                    else
                        dmax = max(dmax, abs(sum1-ONE))
                    end if
                end do
            end do
        end do

    end subroutine shpgi_lower_routine

    !
    !
    ! ... file shpg.f
    !
    ! ... files which must be loaded with shpg.f
    !
    !     type_RealPeriodicTransform.f
    !
    !     shpg computes the harmonic projection, which is
    !     equivalent to a harmonic analysis (forward) followed
    !     by a harmonic synthesis (backward transform).
    !     shpg uses the n**2 projection or complement when appropriate
    !     as well as  odd/even factorization and zero truncation on an
    !     on a Gaussian distributed grid as defined in the JCP paper
    !     "Generalized discrete spherical harmonic transforms"
    !     by Paul N. Swarztrauber and William F. Spotz
    !     J. Comp. Phys., 159(2000) pp. 213-230.
    !
    !     subroutine shpg(nlat, nlon, isym, mtrunc, x, y, idxy,
    !    1        wshp, lwshp, iwshp, liwshp, work, lwork, ierror)
    !
    !     shpg projects the array x onto the set of functions represented
    !     by a discrete set of spherical harmonics.
    !
    !     input parameters
    !
    !     nlat   the number of colatitudes on the full sphere including the
    !            poles. for example, nlat = 37 for a five degree grid.
    !            nlat determines the grid increment in colatitude as
    !            pi/(nlat-1).  if nlat is odd the equator is located at
    !            grid point i=(nlat+1)/2. if nlat is even the equator is
    !            located half way between points i=nlat/2 and i=nlat/2+1.
    !            nlat must be at least 3.
    !
    !     nlon   the number of distinct londitude points.  nlon determines
    !            the grid increment in longitude as 2*pi/nlon. for example
    !            nlon = 72 for a five degree grid. nlon must be greater
    !            than or equal to 4. the efficiency of the computation is
    !            improved when nlon is a product of small prime numbers.
    !            nlon must be at least 4.
    !
    !     isym   currently not used.
    !
    !     mtrunc the highest longitudinal wave number retained in the
    !            projection. It must be less than or equal to
    !            the minimum of nlat-1 and nlon/2. The first wave
    !            number is zero. For example, if wave numbers 0 and
    !            1 are desired then mtrunc = 1.

    !            zero.
    !
    !     x      a two dimensional array that contains the the nlat
    !            by nlon array x(i, j) defined at the colatitude point
    !            theta(i) = (i-1)*pi/(nlat-1) and longitude point phi(j) =
    !            (j-1)*2*pi/nlon.
    !
    !     idxy   the first dimension of the arrays x and y as they
    !            appear in the program that calls shpg. It must be
    !            at least nlat.
    !
    !     wshp   a single precision array that must be saved for
    !            repeated use by subroutine shpg.
    !
    !     lwshp  the dimension of the array wshp as it appears in the
    !            program that calls shpgi. It must be at least
    !            2*(nlat+1)**2+nlon+log2(nlon)
    !
    !     iwshp  an integer array that must be saved for repeated
    !            use by subroutine shpg.
    !
    !
    !     liwshp the dimension of the array iwshp as it appears in the
    !            program that calls shpgi. It must be at least
    !            4*(nlat+1).
    !
    !     work   a single precision work array that does
    !            not have to be saved.
    !
    !     lwork  the dimension of the array work as it appears in the
    !            program that calls shpg. It must be at least
    !            max(nlat*nlon, 4*(nlat+1)).
    !
    !     **************************************************************
    !
    !     output parameters
    !
    !     y      an nlat by nlon single precision array that contains
    !            the projection of x onto the set of functions that
    !            can be represented by the discrete set of spherical
    !            harmonics. The arrays x(i, j) and y(i, j) are located
    !            at colatitude point theta(i) = (i-1)*pi/(nlat-1) and
    !            longitude point phi(j) = (j-1)*2*pi/nlon.
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of isym
    !            = 4  error in the specification of mtrunc
    !            = 5  error in the specification of lwshp
    !            = 6  error in the specification of liwshp
    !            = 7  error in the specification of lwork
    !
    subroutine shpg(nlat, nlon, isym, mtrunc, x, y, idxy, &
        wshp, lwshp, iwshp, liwshp, work, lwork, ierror)


        type(SpherepackAux) :: sphere_aux
        integer :: i
        integer :: idxy
        integer :: ierror
        integer :: isym
        integer :: iw1
        integer :: iw2
        integer :: iw3
        integer :: iw4
        integer :: iwshp
        integer :: j
        integer :: jw1
        integer :: jw2
        integer :: jw3
        integer :: jw4
        integer :: liwshp
        integer :: log2n
        integer :: lw1
        integer :: lwork
        integer :: lwshp
        integer :: mmax
        integer :: mtrunc
        integer :: mwrk
        integer :: nlat
        integer :: nloc1
        integer :: nloc2
        integer :: nlon
        integer :: nte
        
        real(wp) :: work
        real(wp) :: wshp
        real(wp) :: x
        real(wp) :: y
        !
        dimension wshp(*), iwshp(*), work(*), x(idxy, nlon), y(idxy, nlon)
        !
        ierror = 1
        if (nlat<1) return
        ierror = 2
        if (nlon<1) return
        !      ierror = 3
        !      if (isym.lt.0_wp .or. isym.gt.2) return
        ierror = 4
        mmax = min(nlat-1, nlon/2)
        if (mtrunc<0 .or. mtrunc>mmax) return
        ierror = 5
        log2n = log(real(nlon))/log(2.0_wp)
        lw1 = 2*(nlat+1)**2
        if (lwshp<lw1+nlon+log2n) return
        ierror = 6
        if (liwshp<4*(nlat+1)) return
        ierror = 7
        mwrk = max(nlat*nlon, 4*(nlat+1))
        if (lwork <mwrk) return
        ierror = 0
        !
        do j=1, nlon
            do i=1, nlat
                y(i, j) = x(i, j)
            end do
        end do
        call sphere_aux%hfft%forward(nlat, nlon, y, idxy, wshp(lw1+1), work)
        !
        nte = (nlat+1)/2
        nloc1 = 2*nte*nte
        nloc2 = nlat+1
        iw1 = 1
        iw2 = iw1+nloc1
        iw3 = iw2+nloc1
        iw4 = iw3+nloc1
        jw1 = 1
        jw2 = jw1+nloc2
        jw3 = jw2+nloc2
        jw4 = jw3+nloc2
        !
        call shpg_lower_routine(nlat, nlon, isym, mtrunc, y, y, idxy, ierror, &
            nte, wshp(iw1), wshp(iw2), wshp(iw3), wshp(iw4), iwshp(jw1), &
            iwshp(jw2), iwshp(jw3), iwshp(jw4), work(jw1), &
            work(jw2), work(jw3), work(jw4))

        call sphere_aux%hfft%backward(nlat, nlon, y, idxy, wshp(lw1+1), work)

        y(1: nlat,:) = y(1:nlat,:)/nlon

    end subroutine shpg

    subroutine shpg_lower_routine(nlat, nlon, isym, mtrunc, sx, sy, idxy, ierror, &
        idp, pe, po, ze, zo, ipse, jzse, ipso, jzso, xe, xo, ye, yo)

        integer :: i
        integer :: idp
        integer :: idxy
        integer :: ierror
        integer :: iip
        integer :: ipse
        integer :: ipso
        integer :: isym
        integer :: j
        integer :: js
        integer :: jzse
        integer :: jzso
        integer :: lag
        integer :: m
        integer :: modn
        integer :: mp1
        integer :: mpm
        integer :: ms2
        integer :: mtrunc
        integer :: mxtr
        integer :: nec
        integer :: nem
        integer :: nlat
        integer :: nlon
        integer :: nmx
        integer :: noc
        integer :: nom
        integer :: ns2
        integer :: nshe
        integer :: nsho
        integer :: nte
        integer :: nto
        real(wp) :: pe
        real(wp) :: po
        real(wp) :: sx
        real(wp) :: sy
        real(wp) :: xe
        real(wp) :: xo
        real(wp) :: ye
        real(wp) :: yo
        real(wp) :: ze
        real(wp) :: zo
        !
        dimension sx(idxy, nlon), sy(idxy, nlon), nshe(2), nsho(2), &
            pe(idp, idp, 2), po(idp, idp, 2), ze(idp, idp, 2), zo(idp, idp, 2), &
            ipse(idp, 2), jzse(idp, 2), ipso(idp, 2), jzso(idp, 2), &
            xe(idp, 2), xo(idp, 2), ye(idp, 2), yo(idp, 2)
        !
        ns2 = nlat/2
        modn = nlat-ns2-ns2
        nte = (nlat+1)/2
        nto = nlat-nte
        !
        mxtr = min(nlat-1, nlon/2, mtrunc)
        nmx = nlat-mxtr
        if (modn==1) then
            nshe(1) = nmx/2
            nshe(2) = (nmx-1)/2
            nsho(1) = (nmx-1)/2
            nsho(2) = nmx/2
        else
            nshe(1) = (nmx-1)/2
            nshe(2) = nmx/2
            nsho(1) = nmx/2
            nsho(2) = (nmx-1)/2
        end if
        !
        iip = 2
        do 100 mp1=1, mxtr+1
            iip = 3-iip
            if (mxtr==nlat-1.and.mp1==1) then
                do i=1, nlat
                    sy(i, mp1) = sx(i, mp1)
                end do
                goto 100
            end if
            m = mp1-1
            mpm = max(1, m+m)
            ms2 = mp1/2
            nem = (nlat-m+1)/2-nshe(iip)
            nom = (nlat-m)/2-nsho(iip)
            nec = nte-nem
            noc = nto-nom

            do i=1, nte
                xe(i, 1) = HALF * (sx(i, mpm)+sx(nlat+1-i, mpm))
                xo(i, 1) = HALF * (sx(i, mpm)-sx(nlat+1-i, mpm))
            end do

            if (mpm<nlon) then
                do i=1, nte
                    xe(i, 2) = HALF * (sx(i, mpm+1)+sx(nlat+1-i, mpm+1))
                    xo(i, 2) = HALF * (sx(i, mpm+1)-sx(nlat+1-i, mpm+1))
                end do
            end if

            lag = 0
            if (m==0.or.mpm==nlon) lag = 1
            if (3*nec<2*nem.or.nem==0) then
                call matrix_multiply(lag, nte, nec, idp, pe(1, 1, iip), nte, idp, &
                    ze(1, 1, iip), xe, ye, ipse(1, iip), jzse(1, iip))
                do i=1, nte
                    ye(i, 1) = xe(i, 1)-ye(i, 1)
                end do
                if (mpm<nlon.and.m/=0) then
                    do i=1, nte
                        ye(i, 2) = xe(i, 2)-ye(i, 2)
                    end do
                end if
            else
                call matrix_multiply(lag, nte, nem, idp, pe(1, nec+1, iip), nte, idp, &
                    ze(1, nec+1, iip), xe, ye, ipse(nec+1, iip), jzse(nec+1, iip))
            end if
            if (3*noc<2*nom.or.nom==0) then
                call matrix_multiply(lag, nto, noc, idp, po(1, 1, iip), nto, idp, &
                    zo(1, 1, iip), xo, yo, ipso(1, iip), jzso(1, iip))
                do i=1, nto
                    yo(i, 1) = xo(i, 1)-yo(i, 1)
                end do
                if (mpm<nlon.and.m/=0) then
                    do i=1, nto
                        yo(i, 2) = xo(i, 2)-yo(i, 2)
                    end do
                end if
            else
                call matrix_multiply(lag, nto, nom, idp, po(1, noc+1, iip), nto, idp, &
                    zo(1, noc+1, iip), xo, yo, ipso(noc+1, iip), jzso(noc+1, iip))
            end if
            do i=1, nto
                sy(i, mpm) = ye(i, 1)+yo(i, 1)
                sy(nlat+1-i, mpm) = ye(i, 1)-yo(i, 1)
            end do
            if (nte>nto) sy(nte, mpm) = ye(nte, 1)
            if (mpm<nlon.and.m/=0) then
                do i=1, nto
                    sy(i, mpm+1) = ye(i, 2)+yo(i, 2)
                    sy(nlat+1-i, mpm+1) = ye(i, 2)-yo(i, 2)
                end do
                if (nte>nto) sy(nte, mpm+1) = ye(nte, 2)
            end if
100     continue
        !
        js = mxtr+mxtr+2
        do j=js, nlon
            do i=1, nlat
                sy(i, j) = ZERO
            end do
        end do

    end subroutine shpg_lower_routine

    subroutine matrix_multiply(lag, lr, lc, ld, a, mc, md, b, x, y, is, js)

        real(wp) :: a
        real(wp) :: b
        integer :: i
        integer :: is
        integer :: j
        integer :: js
        integer :: k
        integer :: kmx
        integer :: lag
        integer :: lc
        integer :: ld
        integer :: lr
        integer :: mc
        integer :: md
        real(wp) :: sum1
        real(wp) :: sum2
        real(wp) :: x
        real(wp) :: y
        dimension a(ld, *), b(md, *), x(ld, 2), y(ld, 2), &
            is(*), js(*)

        kmx = min(lr+1, ld)

        select case (lag)
            case(1)
                y(1: kmx, 1) = ZERO
                if (lc > 0) then
                    do i=1, lc
                        sum1 = ZERO
                        do j=js(i), mc
                            sum1 = sum1 + b(j, i)*x(j, 1)
                        end do
                        do k=is(i), lr
                            y(k, 1) = y(k, 1)+sum1*a(k, i)
                        end do
                    end do
                end if
            case default
                y(1: kmx, :) = ZERO
                if (lc > 0) then
                    do i=1, lc
                        sum1 = ZERO
                        sum2 = ZERO
                        do j=js(i), mc
                            sum1 = sum1 + b(j, i)*x(j, 1)
                            sum2 = sum2 + b(j, i)*x(j, 2)
                        end do
                        do k=is(i), lr
                            y(k, 1) = y(k, 1)+sum1*a(k, i)
                            y(k, 2) = y(k, 2)+sum2*a(k, i)
                        end do
                    end do
                end if
        end select

    end subroutine matrix_multiply

    subroutine truncate(irc, n, idp, a, nrc, ijs)

        integer(ip) :: i
        integer(ip) :: idp
        integer(ip) :: ijs(n)
        integer(ip) :: irc
        integer(ip) :: j
        integer(ip) :: n
        integer(ip) :: nrc
        real(wp) :: a(idp, *)

        !     irc = 0 for columns , or irc = 1 for rows
        select case (irc)
            case(0)
                outer_loop: do j=1, nrc
                    do i=1, n
                        ijs(j) = i
                        if (abs(a(i, j)) > MACHINE_EPSILON) cycle outer_loop
                    end do
                end do outer_loop
            case default
                default_outer_loop: do i=1, nrc
                    do j=1, n
                        ijs(i) = j
                        if (abs(a(i, j)) > MACHINE_EPSILON) cycle default_outer_loop
                    end do
                end do default_outer_loop
        end select

    end subroutine truncate

    ! Purpose:
    !
    ! Accumulate inner products of x with respect to y.
    !
    subroutine accumulate_inner_products(n, x, y, z)

        ! Dummy arguments
        integer(ip), intent(in)    :: n
        real(wp),    intent(in)    :: x(n)
        real(wp),    intent(in)    :: y(n)
        real(wp),    intent(inout) :: z(n)

        !  Let the intrinsic function dot_product take care of optimization.
        z = z + dot_product(x,y) * y

    end subroutine accumulate_inner_products

    subroutine compute_normal_gaussian_grid(n, x, id, q)

        ! Dummy arguments
        integer(ip), intent(in)    :: n
        real(wp),    intent(inout) :: x(n)
        integer(ip), intent(in)    :: id
        real(wp),    intent(in)    :: q(n)

        ! Local variables
        integer(ip) :: i
        real(wp)    :: sqs

        ! Normalize x
        sqs = ZERO
        do i=1, n
            sqs = sqs+q(i)*(x(i)**2)
        end do

        x = x/sqrt(sqs)

    end subroutine compute_normal_gaussian_grid

end module module_shpg
