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
!                            August 2003
!
! ... file shpe.f
!
!     this file contains code and documentation for subroutines
!     shpei and shpe.
!
! ... files which must be loaded with shpe.f
!
!     type_SpherepackAux.f, type_RealPeriodicTransform.f
!
!     subroutine shpei initializes arrays wshp and iwshp for
!     subsequent repeated use by subroutine shpe, which
!     performs the harmonic projection equivalent to a
!     harmonic analysis followed by harmonic synthesis
!     but faster and with less memory. (see description of
!     subroutine shpe below)
!
!     subroutine shpei(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, 
!    1 liwshp, work, lwork, ierror)
!
!     shpei initializes arrays wshp and iwshp for repeated use
!     by subroutine shpe ....
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
!            nlon must beat least 4.
!
!     isym   currently not used.    
!
!     mtrunc the highest longitudinal wave number retained in the
!            projection. It must be less than or equal to
!            the minimum of nlat-1 and nlon/2. The first wave
!            number is zero. For example, if wave numbers 0 and
!            1 are desired then mtrunc = 1.
!
!     lwshp  the dimension of the array wshp as it appears in the
!            program that calls shpei. It must be at least
!            2*(nlat+1)**2+nlon+log2(nlon)
!
!     liwshp the dimension of the array iwshp as it appears in the
!            program that calls shpei. It must be at least
!            4*(nlat+1).
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shpei. It must be at least
!            1.25*(nlat+1)**2+7*nlat+8.
!
!     **************************************************************
!
!     output parameters
!
!     wshp   a single precision array that must be saved for
!            repeated use by subroutine shpe.        
!
!     iwshp  an integer array that must be saved for repeated
!            use by subroutine shpe.        
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
module module_shpe

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, &
        MACHINE_EPSILON

    use type_SpherepackAux, only: &
        SpherepackAux

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: shpe
    public :: shpei

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp

contains

    subroutine shpei(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, &
        liwshp, work, lwork, ierror)

        type(SpherepackAux) :: sphere_aux
        integer(ip) :: ierror
        integer(ip) :: isym
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip) :: iw3
        integer(ip) :: iw4
        integer(ip) :: iwshp
        integer(ip) :: jw1
        integer(ip) :: jw2
        integer(ip) :: jw3
        integer(ip) :: jw4
        integer(ip) :: kw1
        integer(ip) :: kw10
        integer(ip) :: kw11
        integer(ip) :: kw12
        integer(ip) :: kw13
        integer(ip) :: kw2
        integer(ip) :: kw3
        integer(ip) :: kw4
        integer(ip) :: kw5
        integer(ip) :: kw6
        integer(ip) :: kw7
        integer(ip) :: kw8
        integer(ip) :: kw9
        integer(ip) :: liwshp
        integer(ip) :: log2n
        integer(ip) :: lw1
        integer(ip) :: lwork
        integer(ip) :: lwshp
        integer(ip) :: mlwk
        integer(ip) :: mmax
        integer(ip) :: mtrunc
        integer(ip) :: nlat
        integer(ip) :: nloc1
        integer(ip) :: nloc2
        integer(ip) :: nlon
        integer(ip) :: nte
        real(wp) :: wshp
        real work(*)
        dimension wshp(*), iwshp(*)
        !
        ierror = 1
        if (nlat < 3) return
        ierror = 2
        if (nlon < 4) return
        !      ierror = 3
        !      if (isym.lt.0  .or.  isym.gt.2) return
        ierror = 4
        mmax = min(nlat-1, nlon/2)
        if (mtrunc<0  .or.  mtrunc>mmax) return
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
        kw1 = 1
        kw2 = kw1+nte
        kw3 = kw2+nte
        kw4 = kw3+nte
        kw5 = kw4+nte+1
        kw6 = kw5+nte
        kw7 = kw6+nte
        kw8 = kw7+nte
        kw9 = kw8+nte
        kw10 = kw9+nloc2+nloc2
        kw11 = kw10+nloc2

        kw12 = kw11+nloc1
        kw13 = kw12+nloc1
        !
        call shpei_lower_routine(nlat, nlon, isym, mtrunc, nte, ierror, wshp(iw1), wshp(iw2), &
            wshp(iw3), wshp(iw4), iwshp(jw1), iwshp(jw2), iwshp(jw3), &
            iwshp(jw4), work(kw1), work(kw2), work(kw3), work(kw4), work(kw5), &
            work(kw6), work(kw7), work(kw8), work(kw9), work(kw10), work(kw11), &
            work(kw12), work(kw11), work(kw12), work(kw13))

    end subroutine shpei

    subroutine shpei_lower_routine(nlat, nlon, isym, mtrunc, idp, ierror, &
        pe, po, ze, zo, ipse, jzse, ipso, jzso, &
        cp, work, wx, s, e, thet, xx, z, a, b, we, ped, wo, pod, u)

        real(wp) :: dfn
        integer(ip) :: i
        integer(ip) :: idp
        integer(ip) :: ierror
        integer(ip) :: info
        integer(ip) :: iip
        integer(ip) :: ipse
        integer(ip) :: ipso
        integer(ip) :: isym
        integer(ip) :: it
        integer(ip) :: j
        integer(ip) :: js
        integer(ip) :: jzse
        integer(ip) :: jzso
        integer(ip) :: k
        integer(ip) :: lock
        integer(ip) :: m
        integer(ip) :: modn
        integer(ip) :: mp1
        integer(ip) :: mrank
        integer(ip) :: ms2
        integer(ip) :: mtrunc
        integer(ip) :: mxtr
        integer(ip) :: n
        integer(ip) :: nem
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nom
        integer(ip) :: nrank
        integer(ip) :: ns2
        integer(ip) :: nshe
        integer(ip) :: nsho
        integer(ip) :: nte
        integer(ip) :: nto
        real(wp) :: pe
        real(wp) :: po
        real(wp) :: toe
        real(wp) :: tusl
        real(wp) :: ze
        real(wp) :: zo
        real (wp) :: summation, dthet, v(1,1), a1, b1, c1
        real cp(idp), work(idp), wx(idp), s(idp+1), &
            e(idp), thet(idp), xx(idp), z(idp), u(idp, idp), &
            we(idp, idp, 2), ped(idp, idp, 2), a(4*idp), b(2*idp), &
            wo(idp, idp, 2), pod(idp, idp, 2)

        dimension pe(idp, idp, 2), po(idp, idp, 2), ze(idp, idp, 2), &
            zo(idp, idp, 2), &
            ipse(idp, 2), jzse(idp, 2), ipso(idp, 2), jzso(idp, 2), &
            nshe(2), nsho(2)

        type(SpherepackAux) :: sphere_aux

        ns2 = nlat/2
        modn = nlat-ns2-ns2
        nte = (nlat+1)/2
        nto = nlat-nte
        tusl = ZERO
        toe = ZERO
        !
        !     compute grid distribution
        !
        dthet = pi/(nlat-1)
        do i=1, nte
            thet(i) = (i-1)*dthet
        end do
        !
        !     compute weight matrices for even functions
        !
        do 40 mp1=1, 2
            m = mp1-1
            mrank = nlat-m-m
            nem = (mrank+1)/2
            do j=1, nem
                n = 2*j+m-2
                call sphere_aux%compute_fourier_coefficients(m, n, cp)
                do i=1, nte
                    call sphere_aux%compute_legendre_polys_from_fourier_coeff(m, n, thet(i), cp, ped(i, j, mp1))
                end do
                if (m>0) ped(1, j, mp1) = ZERO
            end do
            call dsvdc(ped(m+1, 1, mp1), idp, nem, nem, s, e, u, &
                idp, v(1,1), idp, work, 10, info)

            do j=1, nem
                s(j) = ONE/(s(j)*s(j))
            end do
            !
            !     compute weight matrix as u  s sup -2 u transpose
            !
            do j=1, nte
                do i=1, nte
                    we(i, j, mp1) = ZERO
                end do
            end do
            do i=1, nem
                do j=1, nem
                    summation = ZERO
                    do k=1, nem
                        summation = summation+s(k)*u(i, k)*u(j, k)
                    end do
                    we(i+m, j+m, mp1) = summation
                end do
            end do
40      continue
        we(1, 1, 2) = ONE
        !
        !     compute n**2 basis (even functions)
        !
        do n=1, nlat+nlat-2
            dfn = n
            a(n) = sqrt(dfn*(dfn+ONE))
        end do
        do n=1, nlat-1
            dfn = n
            b(n) = sqrt((dfn+dfn+3.0)/(dfn+dfn-ONE))
        end do
        !
        mxtr = min(nlat-1, nlon/2, mtrunc)
        iip = 2
        do 200 mp1=1, mxtr+1
            m = mp1-1
            iip = 3-iip
            ms2 = mp1/2
            nrank = ms2+ms2
            mrank = nlat-nrank
            nem = (mrank+1)/2
            !
            !     compute associated legendre functions
            !
            if (m<=1) then
                do 205 j=1, nem
                    n = 2*j+m-2
                    call sphere_aux%compute_fourier_coefficients(m, n, cp)
                    do i=1, nte
                        call sphere_aux%compute_legendre_polys_from_fourier_coeff(m, n, thet(i), cp, ped(i, j+ms2, iip))
                    end do
202                 if (m>0) ped(1, j+ms2, iip) = ZERO
205             continue
            !
            else
                !
                do 207 j=1, nem
                    n = 2*j+m-2
                    if (m>1 .and. n>mxtr) then
                        do i=1, nte
                            u(i, j+ms2) = ped(i, j+ms2, iip)
                        end do
                        goto 207
                    end if
                    a1 = b(n-1)*a(n+m-3)/a(n+m-1)
                    b1 = a(n-m+1)/a(n+m-1)
                    if (n-m<=1) then
                        do i=1, nte
                            u(i, j+ms2) = a1*ped(i, j+ms2-1, iip) &
                                - b1*ped(i, j+ms2, iip)
                        end do
                    else
                        c1 = b(n-1)*a(n-m-1)/a(n+m-1)
                        do i=1, nte
                            u(i, j+ms2) = a1*ped(i, j+ms2-1, iip) &
                                - b1*ped(i, j+ms2, iip) + c1*u(i, j+ms2-1)
                        end do
                    end if
207             continue
                do j=1, nem
                    do i=1, nte
                        ped(i, j+ms2, iip) = u(i, j+ms2)
                    end do
                end do
            end if
            !
            if (ms2<=0. .or. ms2>=nte) goto 200
            !
            ! initialize array with random numbers using
            ! Fortran90 intrinsics RANDOM_{SEED, NUMBER}
            !
            ! old code commented out
            !     do i=1, nte
            !     xx(i) = rand()
            !     end do
            !
            ! replacement code
            !
            call random_seed()
            call random_number(xx(1:nte))
            it = 0
            201 do i=1, nte
                z(i) = ZERO
                wx(i) = ZERO
                do j=1, nte
                    wx(i) = wx(i)+we(i, j, iip)*xx(j)
                end do
            end do
            do 220 j=1, nte
                if (j==ms2) goto 220
                call accumulate_inner_products(nte, wx, ped(1, j, iip), z)
220         continue
            !
            do i=1, nte
                xx(i) = xx(i)-z(i)
            end do
            call normal(nte, xx, idp, we(1, 1, iip))
            it = it+1
            if (it<=2) goto 201
            do i=1, nte
                ped(i, ms2, iip) = xx(i)
            end do
200     continue
        !
        !     reorder if mtrunc is less than nlat-1
        !         case of even functions
        !
        if (modn==0) then
            nshe(1) = (nlat-mtrunc-1)/2
            nshe(2) = (nlat-mtrunc-2)/2
        else
            nshe(1) = (nlat-mtrunc)/2
            nshe(2) = (nlat-mtrunc-1)/2
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
        call trunc(0, nte, idp, ped(1, 1, 1), nte, ipse(1, 1))
        call trunc(0, nte, idp, ped(1, 1, 2), nte, ipse(1, 2))
        !
        !     compute the analysis matrices
        !
        do 250 iip=1, 2
            do i=1, nte
                lock = 0
                do j=1, nte
                    summation = ZERO
                    do k=1, nte
                        summation = summation+ped(k, i, iip)*we(k, j, iip)
                    end do
                    pe(i, j, iip) = ped(i, j, iip)
                    ze(j, i, iip) =  summation
                    if (abs(summation)>MACHINE_EPSILON .and. lock==0) then
                        lock = 1
                        jzse(i, iip) = j
                    end if
                end do
            end do
250     continue
        !
        !     compute weight matrices for odd functions
        !
        do mp1=1, 2
            m = mp1-1
            mrank = nlat-m-m
            nem = (mrank+1)/2
            nom = mrank-nem
            do j=1, nom
                n = 2*j+m-1
                call sphere_aux%compute_fourier_coefficients(m, n, cp)
                do i=1, nte
                    call sphere_aux%compute_legendre_polys_from_fourier_coeff(m, n, thet(i), cp, pod(i, j, mp1))
                end do
                if (modn==1) pod(nte, j, mp1) = ZERO
            end do
            call dsvdc(pod(m+1, 1, mp1), idp, nom, nom, s, e, u, &
                idp, v(1,1), idp, work, 10, info)
            !
            do j=1, nom
                s(j) = ONE/(s(j)*s(j))
            end do
            !
            !     compute weight matrix as u  s sup -2 u transpose
            !
            do j=1, nte
                do i=1, nte
                    wo(i, j, mp1) = ZERO
                end do
            end do
            do i=1, nom
                do j=1, nom
                    summation = ZERO
                    do k=1, nom
                        summation = summation+s(k)*u(i, k)*u(j, k)
                    end do
                    wo(i+m, j+m, mp1) = summation
                end do
            end do
        end do

        wo(1, 1, 2) = ONE
        if (modn==1) then
            wo(nte, nte, 1) = ONE
            wo(nte, nte, 2) = ONE
        end if
        !
        !     compute n**2 basis (odd functions)
        !
        iip = 2
        do 300 mp1=1, mxtr+1
            iip = 3-iip
            m = mp1-1
            ms2 = mp1/2
            nrank = ms2+ms2
            mrank = nlat-nrank
            nem = (mrank+1)/2
            nom = mrank-nem
            !
            !     compute associated legendre functions
            !
            if (m<=1) then
                do 305 j=1, nom
                    n = 2*j+m-1
                    call sphere_aux%compute_fourier_coefficients(m, n, cp)
                    do i=1, nte
                        call sphere_aux%compute_legendre_polys_from_fourier_coeff(m, n, thet(i), cp, pod(i, j+ms2, iip))
                    end do
302                 if (modn==1) pod(nte, j+ms2, iip) = ZERO
                    if (m>0) pod(1, j+ms2, iip) = ZERO
305             continue
            !
            else
                !
                do 307 j=1, nom
                    n = 2*j+m-1
                    if (m>1 .and. n>mxtr) then
                        do i=1, nte
                            u(i, j+ms2) = pod(i, j+ms2, iip)
                        end do
                        goto 304
                    end if
                    a1 = b(n-1)*a(n+m-3)/a(n+m-1)
                    b1 = a(n-m+1)/a(n+m-1)
                    if (n-m<=1) then
                        do i=1, nte
                            u(i, j+ms2) = a1*pod(i, j+ms2-1, iip) &
                                - b1*pod(i, j+ms2, iip)
                        end do
                    else
                        c1 = b(n-1)*a(n-m-1)/a(n+m-1)
                        do i=1, nte
                            u(i, j+ms2) = a1*pod(i, j+ms2-1, iip) &
                                - b1*pod(i, j+ms2, iip) + c1*u(i, j+ms2-1)
                        end do
                    end if
304                 if (modn==1) u(nte, j+ms2) = ZERO
307             continue
                do j=1, nom
                    do i=1, nte
                        pod(i, j+ms2, iip) = u(i, j+ms2)
                    end do
                end do
            end if
            !
            if (ms2<=0. .or. ms2>=nto) goto 300
            !
            ! initialize array with random numbers using
            ! Fortran 90 intrinsics random_seed and random_number
            !
            ! old code commented out
            !
            !     do i=1, nte
            !     xx(i) = rand()
            !     end do
            ! replacement code
            !
            call random_number(xx(1:nte))
            if (modn==1) xx(nte) = ZERO
            it = 0
            306 do i=1, nte
                z(i) = ZERO
                wx(i) = ZERO
                do j=1, nto
                    wx(i) = wx(i)+wo(i, j, iip)*xx(j)
                end do
            end do

            do j=1, nto
                if (j==ms2) cycle
                call accumulate_inner_products(nte, wx, pod(1, j, iip), z(1))
            end do

            do i=1, nte
                xx(i) = xx(i)-z(i)
            end do

            call normal(nte, xx, idp, wo(1, 1, iip))

            it = it+1
            if (it<=2) goto 306

            do i=1, nte
                pod(i, ms2, iip) = xx(i)
            end do

            if (modn==1) pod(nte, ms2, iip) = ZERO
300     continue
        !
        !     reorder if mtrunc is less than nlat-1
        !        case of odd functions
        !
        if (modn==0) then
            nsho(1) = (nlat-mtrunc)/2
            nsho(2) = (nlat-mtrunc-1)/2
        else
            nsho(1) = (nlat-mtrunc-1)/2
            nsho(2) = (nlat-mtrunc-2)/2
        end if

        do mp1=1, 2
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
        end do

        call trunc(0, nte, idp, pod(1, 1, 1), nto, ipso(1, 1))
        call trunc(0, nte, idp, pod(1, 1, 2), nto, ipso(1, 2))
        !
        !     compute the analysis matrices (odd functions)
        !
        do iip=1, 2
            do i=1, nto
                lock = 0
                do j=1, nto
                    summation = ZERO
                    do k=1, nte
                        summation = summation+pod(k, i, iip)*wo(k, j, iip)
                    end do
                    po(i, j, iip) = pod(i, j, iip)
                    zo(j, i, iip) = summation
                    if (abs(summation)>MACHINE_EPSILON .and. lock==0) then
                        lock = 1
                        jzso(i, iip) = j
                    end if
                end do
            end do
        end do

    end subroutine shpei_lower_routine
    !
    !
    ! ... file shpe.f
    !
    ! ... files which must be loaded with shpe.f
    !
    !     type_SpherepackAux.f, type_RealPeriodicTransform.f
    !
    !     the n**2 projection with complement, odd/even
    !     factorization and zero truncation on an
    !     equally spaced grid as defined in the JCP paper
    !     "Generalized discrete spherical harmonic transforms"
    !     by Paul N. Swarztrauber and William F. Spotz
    !     It is equivalent to a harmonic analysis followed
    !     by a synthesis except faster and requires less memory.
    !
    !     subroutine shpe(nlat, nlon, isym, mtrunc, x, y, idxy,
    !    1        wshp, lwshp, iwshp, liwshp, work, lwork, ierror)
    !
    !     shpe projects the array x onto the set of functions represented
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
    !            nlon must beat least 4.
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
    !            appear in the program that calls shpe. It must be
    !            at least nlat.
    !
    !     wshp   a single precision array that must be saved for
    !            repeated use by subroutine shpe.
    !
    !     lwshp  the dimension of the array wshp as it appears in the
    !            program that calls shpei. It must be at least
    !            2*(nlat+1)**2+nlon+log2(nlon)
    !
    !     iwshp  an integer array that must be saved for repeated
    !            use by subroutine shpe.
    !
    !
    !     liwshp the dimension of the array iwshp as it appears in the
    !            program that calls shpei. It must be at least
    !            4*(nlat+1).
    !
    !     work   a single precision work array that does
    !            not have to be saved.
    !
    !     lwork  the dimension of the array work as it appears in the
    !            program that calls shpe. It must be at least
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
    subroutine shpe(nlat, nlon, isym, mtrunc, x, y, idxy, &
        wshp, lwshp, iwshp, liwshp, work, lwork, ierror)

        type(SpherepackAux) :: sphere_aux
        
        integer(ip) :: idxy
        integer(ip) :: ierror
        integer(ip) :: isym
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip) :: iw3
        integer(ip) :: iw4
        integer(ip) :: iwshp(liwshp)
        integer(ip) :: jw1
        integer(ip) :: jw2
        integer(ip) :: jw3
        integer(ip) :: jw4
        integer(ip) :: liwshp
        integer(ip) :: log2n
        integer(ip) :: lw1
        integer(ip) :: lwork
        integer(ip) :: lwshp
        integer(ip) :: mmax
        integer(ip) :: mtrunc
        integer(ip) :: mwrk
        integer(ip) :: nlat
        integer(ip) :: nloc1
        integer(ip) :: nloc2
        integer(ip) :: nlon
        integer(ip) :: nte
        real(wp), intent(out) :: work(lwork)
        real(wp) :: wshp(lwshp)
        real(wp) :: x(idxy, nlon), y(idxy, nlon)

        ! Check input arguments
        ierror = 1
        if (nlat < 3) return
        ierror = 2
        if (nlon < 4) return
        !      ierror = 3
        !      if (isym.lt.0  .or.  isym.gt.2) return
        ierror = 4
        mmax = min(nlat-1, nlon/2)
        if (mtrunc<0  .or.  mtrunc>mmax) return
        ierror = 5
        log2n = log(real(nlon))/log(TWO)
        lw1 = 2*(nlat+1)**2
        if (lwshp<lw1+nlon+log2n) return
        ierror = 6
        if (liwshp<4*(nlat+1)) return
        ierror = 7
        mwrk = max(nlat*nlon, 4*(nlat+1))
        if (lwork <mwrk) return
        ierror = 0

        y(1:nlat,:) = x(1:nlat,:)

        call sphere_aux%hfft%forward(nlat, nlon, y, idxy, wshp(lw1+1), work)


        nte = (nlat+1)/2
        nloc1 = 2*nte*nte
        nloc2 = nlat+1

        ! Set workspace index pointers
        iw1 = 1
        iw2 = iw1+nloc1
        iw3 = iw2+nloc1
        iw4 = iw3+nloc1
        jw1 = 1
        jw2 = jw1+nloc2
        jw3 = jw2+nloc2
        jw4 = jw3+nloc2

        call shpe_lower_routine(nlat, nlon, isym, mtrunc, y, y, idxy, ierror, &
            nte, wshp(iw1), wshp(iw2), wshp(iw3), wshp(iw4), iwshp(jw1), &
            iwshp(jw2), iwshp(jw3), iwshp(jw4), work(jw1), &
            work(jw2), work(jw3), work(jw4))

        call sphere_aux%hfft%backward(nlat, nlon, y, idxy, wshp(lw1+1), work)

        y(1:nlat,:) = y(1:nlat,:)/nlon

    end subroutine shpe

    subroutine shpe_lower_routine(nlat, nlon, isym, mtrunc, sx, sy, idxy, ierror, &
        idp, pe, po, ze, zo, ipse, jzse, ipso, jzso, xe, xo, ye, yo)

        integer(ip) :: i
        integer(ip) :: idp
        integer(ip) :: idxy
        integer(ip) :: ierror
        integer(ip) :: iip
        integer(ip) :: ipse
        integer(ip) :: ipso
        integer(ip) :: isym
        integer(ip) :: j
        integer(ip) :: js
        integer(ip) :: jzse
        integer(ip) :: jzso
        integer(ip) :: m
        integer(ip) :: modn
        integer(ip) :: mp1
        integer(ip) :: mpm
        integer(ip) :: mrank
        integer(ip) :: ms2
        integer(ip) :: mtrunc
        integer(ip) :: mxtr
        integer(ip) :: nec
        integer(ip) :: nem
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: noc
        integer(ip) :: nom
        integer(ip) :: nrank
        integer(ip) :: ns2
        integer(ip) :: nshe
        integer(ip) :: nsho
        integer(ip) :: nte
        integer(ip) :: nto
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
        if (modn==0) then
            nshe(1) = (nlat-mtrunc-1)/2
            nshe(2) = (nlat-mtrunc-2)/2
            nsho(1) = (nlat-mtrunc)/2
            nsho(2) = (nlat-mtrunc-1)/2
        else
            nshe(1) = (nlat-mtrunc)/2
            nshe(2) = (nlat-mtrunc-1)/2
            nsho(1) = (nlat-mtrunc-1)/2
            nsho(2) = (nlat-mtrunc-2)/2
        end if
        mxtr = min(nlat-1, nlon/2, mtrunc)
        iip = 2

        main_loop: do mp1=1, mxtr+1
            iip = 3-iip
            if (mxtr==nlat-1 .and. mp1<=2) then
                do i=1, nlat
                    sy(i, mp1) = sx(i, mp1)
                end do
                if (mp1==2) then
                    sy(1, 2) = ZERO
                    sy(nlat, 2) = ZERO
                end if
                if (nlon>=3) then
                    sy(1, 3) = ZERO
                    sy(nlat, 3) = ZERO
                    do i=2, nlat-1
                        sy(i, 3) = sx(i, 3)
                    end do
                end if
                cycle main_loop
            end if
            m = mp1-1
            mpm = max(1, m+m)
            ms2 = mp1/2
            mrank = min(nlat-m, nlat-ms2-ms2)
            !      mrank = mxtr+1-ms2-ms2
            nrank = nlat-mrank
            nem = (mrank+1)/2-nshe(iip)
            nom = mrank-(mrank+1)/2-nsho(iip)
            nec = nte-nem
            noc = nto-nom
            !
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
            if (3*nec<2*nem .or. nem==0) then
                call tmxmx(nte, nec, idp, pe(1, 1, iip), nte, idp, &
                    ze(1, 1, iip), xe, ye, ipse(1, iip), jzse(1, iip))
                do i=1, nte
                    ye(i, 1) = xe(i, 1)-ye(i, 1)
                end do
                if (mpm<nlon .and. m/=0) then
                    do i=1, nte
                        ye(i, 2) = xe(i, 2)-ye(i, 2)
                    end do
                end if
            else
                call tmxmx(nte, nem, idp, pe(1, nec+1, iip), nte, idp, &
                    ze(1, nec+1, iip), xe, ye, ipse(nec+1, iip), jzse(nec+1, iip))
            end if
            if (3*noc<2*nom .or. nom==0) then
                call tmxmx(nto, noc, idp, po(1, 1, iip), nto, idp, &
                    zo(1, 1, iip), xo, yo, ipso(1, iip), jzso(1, iip))
                do i=1, nte
                    yo(i, 1) = xo(i, 1)-yo(i, 1)
                end do
                if (mpm<nlon .and. m/=0) then
                    do i=1, nte
                        yo(i, 2) = xo(i, 2)-yo(i, 2)
                    end do
                end if
            else
                call tmxmx(nto, nom, idp, po(1, noc+1, iip), nto, idp, &
                    zo(1, noc+1, iip), xo, yo, ipso(noc+1, iip), jzso(noc+1, iip))
            end if
            do i=1, nte
                sy(i, mpm) = ye(i, 1)+yo(i, 1)
                sy(nlat+1-i, mpm) = ye(i, 1)-yo(i, 1)
            end do
            if (mpm<nlon .and. m/=0) then
                do i=1, nte
                    sy(i, mpm+1) = ye(i, 2)+yo(i, 2)
                    sy(nlat+1-i, mpm+1) = ye(i, 2)-yo(i, 2)
                end do
            end if

        end do main_loop

        js = mxtr+mxtr+2

        do j=js, nlon
            do i=1, nlat
                sy(i, j) = ZERO
            end do
        end do

    end subroutine shpe_lower_routine

    subroutine mxm(lr, lc, ld, a, mc, md, b, nd, c)

        integer(ip) :: i
        integer(ip) :: j
        integer(ip) :: k
        integer(ip) :: lc
        integer(ip) :: ld
        integer(ip) :: lr
        integer(ip) :: mc
        integer(ip) :: md
        integer(ip) :: nd
        real a(ld, *), b(md, *), c(nd, *)
        do i=1, lr
            do j=1, mc
                c(i, j) = ZERO
                do k=1, lc
                    c(i, j) = c(i, j)+a(i, k)*b(k, j)
                end do
            end do
        end do

    end subroutine mxm

    subroutine smxm(lr, lc, ld, a, mc, md, b, nd, c)

        real(wp) :: a
        real(wp) :: b
        real(wp) :: c
        integer(ip) :: i
        integer(ip) :: j
        integer(ip) :: k
        integer(ip) :: lc
        integer(ip) :: ld
        integer(ip) :: lr
        integer(ip) :: mc
        integer(ip) :: md
        integer(ip) :: nd
        dimension a(ld, *), b(md, *), c(nd, *)
        do i=1, lr
            do j=1, mc
                c(i, j) = ZERO
                do k=1, lc
                    c(i, j) = c(i, j)+a(i, k)*b(k, j)
                end do
            end do
        end do

    end subroutine smxm

    subroutine mxmx(lr, lc, ld, a, mc, md, b, x, y)

        real(wp) :: a
        real(wp) :: b
        integer(ip) :: i
        integer(ip) :: j
        integer(ip) :: k
        integer(ip) :: lc
        integer(ip) :: ld
        integer(ip) :: lr
        integer(ip) :: mc
        integer(ip) :: md
        real(wp) :: sum1
        real(wp) :: sum2
        real(wp) :: x
        real(wp) :: y
        dimension a(ld, *), b(md, *), x(ld, 2), y(ld, 2)
        do k=1, lr
            y(k, 1) = ZERO
            y(k, 2) = ZERO
        end do
        !
        if (lc <= 0) return
        do i=1, lc
            sum1 = ZERO
            sum2 = ZERO
            do j=1, mc
                sum1 = sum1 + b(i, j)*x(j, 1)
                sum2 = sum2 + b(i, j)*x(j, 2)
            end do
            do k=1, lr
                y(k, 1) = y(k, 1)+sum1*a(k, i)
                y(k, 2) = y(k, 2)+sum2*a(k, i)
            end do
        end do

    end subroutine mxmx

    subroutine dmxmx(lr, lc, ld, a, mc, md, b, x, y)

        integer(ip) :: i
        integer(ip) :: j
        integer(ip) :: k
        integer(ip) :: lc
        integer(ip) :: ld
        integer(ip) :: lr
        integer(ip) :: mc
        integer(ip) :: md
        real a(ld, *), b(md, *), x(ld, 2), y(ld, 2), &
            sum1, sum2
        do k=1, lr
            y(k, 1) = ZERO
            y(k, 2) = ZERO
        end do
        !
        if (lc <= 0) return
        do i=1, lc
            sum1 = ZERO
            sum2 = ZERO
            do j=1, mc
                sum1 = sum1 + b(i, j)*x(j, 1)
                sum2 = sum2 + b(i, j)*x(j, 2)
            end do
            do k=1, lr
                y(k, 1) = y(k, 1)+sum1*a(k, i)
                y(k, 2) = y(k, 2)+sum2*a(k, i)
            end do
        end do

    end subroutine dmxmx

    subroutine tmxmx(lr, lc, ld, a, mc, md, b, x, y, is, js)

        real(wp) :: a
        real(wp) :: b
        integer(ip) :: i
        integer(ip) :: is
        integer(ip) :: j
        integer(ip) :: js
        integer(ip) :: k
        integer(ip) :: kmx
        integer(ip) :: lc
        integer(ip) :: ld
        integer(ip) :: lr
        integer(ip) :: mc
        integer(ip) :: md
        real(wp) :: sum1
        real(wp) :: sum2
        real(wp) :: x
        real(wp) :: y
        dimension a(ld, *), b(md, *), x(ld, 2), y(ld, 2), &
            is(*), js(*)

        kmx = min(lr+1, ld)
        do k=1, kmx
            y(k, 1) = ZERO
            y(k, 2) = ZERO
        end do

        if (lc <= 0) return

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

    end subroutine tmxmx

    subroutine trunc(irc, n, idp, a, nrc, ijs)

        integer(ip) :: i
        integer(ip) :: idp
        integer(ip) :: ijs(n)
        integer(ip) :: irc
        integer(ip) :: j
        integer(ip) :: n
        integer(ip) :: nrc
        real(wp) :: a(idp, *)
        !
        !     irc = 0 for columns , or irc = 1 for rows
        !
        select case (irc)
            case(0)
                loop_20: do j=1, nrc
                    do i=1, n
                        ijs(j) = i
                        if (abs(a(i, j)) > MACHINE_EPSILON) cycle loop_20
                    end do
                end do loop_20
            case default
                default_outer_loop: do i=1, nrc
                    do j=1, n
                        ijs(i) = j
                        if (abs(a(i, j)) > MACHINE_EPSILON) cycle default_outer_loop
                    end do
                end do default_outer_loop
        end select

    end subroutine trunc

    ! Purpose:
    !
    ! Accumulate inner products of x with respect to y.
    !
    subroutine accumulate_inner_products(n, x, y, z)

        ! Dummy arguments
        integer(ip), intent(in)  :: n
        real(wp),    intent(in)  :: x(n)
        real(wp),    intent(in)  :: y(n)
        real(wp),    intent(out) :: z(n)

        ! Local variables
        integer(ip) :: i
        real(wp)    :: summation

        summation = dot_product(x, y)

        do i=1, n
            z(i) = z(i)+summation*y(i)
        end do

    end subroutine accumulate_inner_products

    subroutine normal(n, x, id, q)

        integer(ip) :: i
        integer(ip) :: id
        integer(ip) :: j
        integer(ip) :: n
        dimension x(n), q(id, n)
        real x, q, summation, sqs
        !
        !     normalize x
        !
        sqs = ZERO
        do i=1, n
            summation = ZERO
            do j=1, n
                summation = summation+q(i, j)*x(j)
            end do
            sqs = sqs+summation*x(i)
        end do

        sqs = sqrt(sqs)
        do i=1, n
            x(i) = x(i)/sqs
        end do

    end subroutine normal

    subroutine coe(moe, n, x, dmax)

        integer(ip) :: i
        integer(ip) :: moe
        integer(ip) :: n
        integer(ip) :: nh
        real(wp) :: x(n), dmax

        nh = (n+1)/2
        dmax = ZERO

        select case (moe)
            case(0)
                do i=1, nh
                    dmax = max(dmax, abs(x(i)-x(n-i+1)))
                    x(i) = HALF * (x(i)+x(n-i+1))
                    x(n-i+1) = x(i)
                end do
            case default
                do i=1, nh
                    dmax = max(dmax, abs(x(i)+x(n-i+1)))
                    x(i) = HALF * (x(i)-x(n-i+1))
                    x(n-i+1) = -x(i)
                end do
                if (mod(n, 2)/=0) x(nh) = ZERO
        end select

    end subroutine coe

    subroutine dsvdc(x, ldx, n, p, s, e, u, ldu, v, ldv, work, job, info)

        integer ldx, n, p, ldu, ldv, job, info
        real x(ldx, 1), s(1), e(1), u(ldu, 1), v(ldv, 1), work(1)
        !
        !
        !     dsvdc is a subroutine to reduce a real nxp matrix x
        !     by orthogonal transformations u and v to diagonal form.  the
        !     diagonal elements s(i) are the singular values of x.  the
        !     columns of u are the corresponding left singular vectors,
        !     and the columns of v the right singular vectors.
        !
        !     on entry
        !
        !         x         real(ldx, p), where ldx.ge.n.
        !                   x contains the matrix whose singular value
        !                   decomposition is to be computed.  x is
        !                   destroyed by dsvdc.
        !
        !         ldx       integer.
        !                   ldx is the leading dimension of the array x.
        !
        !         n         integer.
        !                   n is the number of rows of the matrix x.
        !
        !         p         integer.
        !                   p is the number of columns of the matrix x.
        !
        !         ldu       integer.
        !                   ldu is the leading dimension of the array u.
        !                   (see below).
        !
        !         ldv       integer.
        !                   ldv is the leading dimension of the array v.
        !                   (see below).
        !
        !         work      real(n).
        !                   work is a scratch array.
        !
        !         job       integer.
        !                   job controls the computation of the singular
        !                   vectors.  it has the decimal expansion ab
        !                   with the following meaning
        !
        !                        a.eq.0    do not compute the left singular
        !                                  vectors.
        !                        a.eq.1    return the n left singular vectors
        !                                  in u.
        !                        a.ge.2    return the first min(n, p) singular
        !                                  vectors in u.
        !                        b.eq.0    do not compute the right singular
        !                                  vectors.
        !                        b.eq.1    return the right singular vectors
        !                                  in v.
        !
        !     on return
        !
        !         s         real(mm), where mm=min(n+1, p).
        !                   the first min(n, p) entries of s contain the
        !                   singular values of x arranged in descending
        !                   order of magnitude.
        !
        !         e         real(p),
        !                   e ordinarily contains zeros.  however see the
        !                   discussion of info for exceptions.
        !
        !         u         real(ldu, k), where ldu.ge.n.  if
        !                                   joba.eq.1 then k.eq.n, if joba.ge.2
        !                                   then k.eq.min(n, p).
        !                   u contains the matrix of left singular vectors.
        !                   u is not referenced if joba.eq.ZERO  if n.le.p
        !                   or if joba.eq.2, then u may be identified with x
        !                   in the subroutine call.
        !
        !         v         real(ldv, p), where ldv.ge.p.
        !                   v contains the matrix of right singular vectors.
        !                   v is not referenced if job.eq.ZERO  if p.le.n,
        !                   then v may be identified with x in the
        !                   subroutine call.
        !
        !         info      integer.
        !                   the singular values (and their corresponding
        !                   singular vectors) s(info+1), s(info+2), ..., s(m)
        !                   are correct (here m=min(n, p)).  thus if
        !                   info.eq.0, all the singular values and their
        !                   vectors are correct.  in any event, the matrix
        !                   b = trans(u)*x*v is the bidiagonal matrix
        !                   with the elements of s on its diagonal and the
        !                   elements of e on its super-diagonal (trans(u)
        !                   is the transpose of u).  thus the singular
        !                   values of x and b are the same.
        !
        !     linpack. this version dated 08/14/78 .
        !              correction made to shift 2/84.
        !     g.w. stewart, university of maryland, argonne national lab.
        !
        !     dsvdc uses the following functions and subprograms.
        !
        !     external drot
        !     blas daxpy, ddot, dscal, dswap, dnrm2, drotg
        !     fortran dabs, dmax1, max0, min0, mod, dsqrt
        !
        !     internal variables
        !
        integer i, iter, j, jobu, k, kase, kk
        integer l, ll, lls, lm1, lp1, ls, lu
        integer  m, maxit, mm, mm1, mp1, nct, nctp1, ncu, nrt, nrtp1
        real t, b, c, cs, el, emm1, f, g, scale_rename
        real shift, sl, sm, sn, smm1, t1, test, ztest
        logical wantu, wantv
        !
        !
        !     set the maximum number of iterations.
        !
        maxit = 30
        !
        !     determine what is to be computed.
        !
        wantu = .false.
        wantv = .false.
        jobu = mod(job, 100)/10
        ncu = n
        if (jobu > 1) ncu = min(n, p)
        if (jobu /= 0) wantu = .true.
        if (mod(job, 10) /= 0) wantv = .true.
        !
        !     reduce x to bidiagonal form, storing the diagonal elements
        !     in s and the super-diagonal elements in e.
        !
        info = 0
        nct = min(n-1, p)
        nrt = max(0, min(p-2, n))
        lu = max(nct, nrt)
        if (lu < 1) goto 170
        do 160 l = 1, lu
            lp1 = l + 1
            if (l > nct) goto 20
            !
            !           compute the transformation for the l-th column and
            !           place the l-th diagonal in s(l).
            !
            s(l) = get_norm2(n-l+1, x(l, l), 1)
            if (s(l) == ZERO) goto 10
            if (x(l, l) /= ZERO) s(l) = sign(s(l), x(l, l))
            call scale_vector_by_constant(n-l+1, ONE/s(l), x(l, l), 1)
            x(l, l) = ONE + x(l, l)
10      continue
        s(l) = -s(l)
20  continue
    if (p < lp1) goto 50
    do 40 j = lp1, p
        if (l > nct) goto 30
        if (s(l) == ZERO) goto 30
        !
        !              apply the transformation.
        !
        t = -get_dot_product(n-l+1, x(l, l), 1, x(l, j), 1)/x(l, l)
        call daxpy(n-l+1, t, x(l, l), 1, x(l, j), 1)
30  continue
    !
    !           place the l-th row of x into  e for the
    !           subsequent calculation of the row transformation.
    !
    e(j) = x(l, j)
40 continue
50 continue
   if (.not.wantu  .or.  l > nct) goto 70
   !
   !           place the transformation in u for subsequent back
   !           multiplication.
   !
   do 60 i = l, n
       u(i, l) = x(i, l)
60 continue
70 continue
   if (l > nrt) goto 150
   !
   !           compute the l-th row transformation and place the
   !           l-th super-diagonal in e(l).
   !
   e(l) = get_norm2(p-l, e(lp1), 1)
   if (e(l) == ZERO) goto 80
   if (e(lp1) /= ZERO) e(l) = sign(e(l), e(lp1))
   call scale_vector_by_constant(p-l, ONE/e(l), e(lp1), 1)
   e(lp1) = ONE + e(lp1)
80 continue
   e(l) = -e(l)
   if (lp1 > n  .or.  e(l) == ZERO) goto 120
   !
   !              apply the transformation.
   !
   do 90 i = lp1, n
       work(i) = ZERO
90 continue
   do 100 j = lp1, p
       call daxpy(n-l, e(j), x(lp1, j), 1, work(lp1), 1)
100 continue
    do 110 j = lp1, p
        call daxpy(n-l, -e(j)/e(lp1), work(lp1), 1, x(lp1, j), 1)
110 continue
120 continue
    if (.not.wantv) goto 140
    !
    !              place the transformation in v for subsequent
    !              back multiplication.
    !
    do 130 i = lp1, p
        v(i, l) = e(i)
130 continue
140 continue
150 continue
160 continue
170 continue
    !
    !     set up the final bidiagonal matrix or order m.
    !
    m = min(p, n+1)
    nctp1 = nct + 1
    nrtp1 = nrt + 1
    if (nct < p) s(nctp1) = x(nctp1, nctp1)
    if (n < m) s(m) = ZERO
    if (nrtp1 < m) e(nrtp1) = x(nrtp1, m)
    e(m) = ZERO
    !
    !     if required, generate u.
    !
    if (.not.wantu) goto 300
    if (ncu < nctp1) goto 200
    do 190 j = nctp1, ncu
        do 180 i = 1, n
            u(i, j) = ZERO
180     continue
        u(j, j) = ONE
190 continue
200 continue
    if (nct < 1) goto 290
    do 280 ll = 1, nct
        l = nct - ll + 1
        if (s(l) == ZERO) goto 250
        lp1 = l + 1
        if (ncu < lp1) goto 220
        do 210 j = lp1, ncu
            t = -get_dot_product(n-l+1, u(l, l), 1, u(l, j), 1)/u(l, l)
            call daxpy(n-l+1, t, u(l, l), 1, u(l, j), 1)
210     continue
220 continue
    call scale_vector_by_constant(n-l+1, -ONE, u(l, l), 1)
    u(l, l) = ONE + u(l, l)
    lm1 = l - 1
    if (lm1 < 1) goto 240
    do 230 i = 1, lm1
        u(i, l) = ZERO
230 continue
240 continue
    goto 270
250 continue
    do 260 i = 1, n
        u(i, l) = ZERO
260 continue
    u(l, l) = ONE
270 continue
280 continue
290 continue
300 continue
    !
    !     if it is required, generate v.
    !
    if (.not.wantv) goto 350
    do 340 ll = 1, p
        l = p - ll + 1
        lp1 = l + 1
        if (l > nrt) goto 320
        if (e(l) == ZERO) goto 320
        do 310 j = lp1, p
            t = -get_dot_product(p-l, v(lp1, l), 1, v(lp1, j), 1)/v(lp1, l)
            call daxpy(p-l, t, v(lp1, l), 1, v(lp1, j), 1)
310     continue
320 continue
    do 330 i = 1, p
        v(i, l) = ZERO
330 continue
    v(l, l) = ONE
340 continue
350 continue
    !
    !     main iteration loop for the singular values.
    !
    mm = m
    iter = 0
360 continue
    !
    !        quit if all the singular values have been found.
    !
    !     ...exit
    if (m == 0) goto 620
    !
    !        if too many iterations have been performed, set
    !        flag and return.
    !
    if (iter < maxit) goto 370
    info = m
    !     ......exit
    goto 620
370 continue
    !
    !        this section of the program inspects for
    !        negligible elements in the s and e arrays.  on
    !        completion the variables kase and l are set as follows.
    !
    !           kase = 1     if s(m) and e(l-1) are negligible and l.lt.m
    !           kase = 2     if s(l) is negligible and l.lt.m
    !           kase = 3     if e(l-1) is negligible, l.lt.m, and
    !                        s(l), ..., s(m) are not negligible (qr step).
    !           kase = 4     if e(m-1) is negligible (convergence).
    !
    do 390 ll = 1, m
        l = m - ll
        !        ...exit
        if (l == 0) goto 400
        test = abs(s(l)) + abs(s(l+1))
        ztest = test + abs(e(l))
        if (ztest /= test) goto 380
        e(l) = ZERO
        !        ......exit
        goto 400
380 continue
390 continue
400 continue
    if (l /= m - 1) goto 410
    kase = 4
    goto 480
410 continue
    lp1 = l + 1
    mp1 = m + 1
    do 430 lls = lp1, mp1
        ls = m - lls + lp1
        !           ...exit
        if (ls == l) goto 440
        test = ZERO
        if (ls /= m) test = test + abs(e(ls))
        if (ls /= l + 1) test = test + abs(e(ls-1))
        ztest = test + abs(s(ls))
        if (ztest /= test) goto 420
        s(ls) = ZERO
        !           ......exit
        goto 440
420 continue
430 continue
440 continue
    if (ls /= l) goto 450
    kase = 3
    goto 470
450 continue
    if (ls /= m) goto 460
    kase = 1
    goto 470
460 continue
    kase = 2
    l = ls
470 continue
480 continue
    l = l + 1
    !
    !        perform the task indicated by kase.
    !
    goto (490, 520, 540, 570), kase
!
!        deflate negligible s(m).
!
490 continue
    mm1 = m - 1
    f = e(m-1)
    e(m-1) = ZERO
    do 510 kk = l, mm1
        k = mm1 - kk + l
        t1 = s(k)
        call construct_given_plane_rotation(t1, f, cs, sn)
        s(k) = t1
        if (k == l) goto 500
        f = -sn*e(k-1)
        e(k-1) = cs*e(k-1)
500 continue
    if (wantv) call apply_plane_rotation(p, v(1, k), 1, v(1, m), 1, cs, sn)
510 continue
    goto 610
!
!        split at negligible s(l).
!
520 continue
    f = e(l-1)
    e(l-1) = ZERO
    do 530 k = l, m
        t1 = s(k)
        call construct_given_plane_rotation(t1, f, cs, sn)
        s(k) = t1
        f = -sn*e(k)
        e(k) = cs*e(k)
        if (wantu) call apply_plane_rotation(n, u(1, k), 1, u(1, l-1), 1, cs, sn)
530 continue
    goto 610
!
!        perform one qr step.
!
540 continue
    !
    !           calculate the shift.
    !
    scale_rename = dmax1(abs(s(m)), abs(s(m-1)), abs(e(m-1)), &
        abs(s(l)), abs(e(l)))
    sm = s(m)/scale_rename
    smm1 = s(m-1)/scale_rename
    emm1 = e(m-1)/scale_rename
    sl = s(l)/scale_rename
    el = e(l)/scale_rename
    b = ((smm1 + sm)*(smm1 - sm) + emm1**2)/TWO
    c = (sm*emm1)**2
    shift = ZERO
    if (b == ZERO .and. c == ZERO) goto 550
    shift = sqrt(b**2+c)
    if (b < ZERO) shift = -shift
    shift = c/(b + shift)
550 continue
    f = (sl + sm)*(sl - sm) + shift
    g = sl*el
    !
    !           chase zeros.
    !
    mm1 = m - 1
    do 560 k = l, mm1
        call construct_given_plane_rotation(f, g, cs, sn)
        if (k /= l) e(k-1) = f
        f = cs*s(k) + sn*e(k)
        e(k) = cs*e(k) - sn*s(k)
        g = sn*s(k+1)
        s(k+1) = cs*s(k+1)
        if (wantv) call apply_plane_rotation(p, v(1, k), 1, v(1, k+1), 1, cs, sn)
        call construct_given_plane_rotation(f, g, cs, sn)
        s(k) = f
        f = cs*e(k) + sn*s(k+1)
        s(k+1) = -sn*e(k) + cs*s(k+1)
        g = sn*e(k+1)
        e(k+1) = cs*e(k+1)
        if (wantu .and. k < n) &
            call apply_plane_rotation(n, u(1, k), 1, u(1, k+1), 1, cs, sn)
560 continue
    e(m-1) = f
    iter = iter + 1
    goto 610
!
!        convergence.
!
570 continue
    !
    !           make the singular value  positive.
    !
    if (s(l) >= ZERO) goto 580
    s(l) = -s(l)
    if (wantv) call scale_vector_by_constant(p, -ONE, v(1, l), 1)
580 continue
    !
    !           order the singular value.
    !
590 if (l == mm) goto 600
    !           ...exit
    if (s(l) >= s(l+1)) goto 600
    t = s(l)
    s(l) = s(l+1)
    s(l+1) = t
    if (wantv .and. l < p) &
        call swap_vectors(p, v(1, l), 1, v(1, l+1), 1)
    if (wantu .and. l < n) &
        call swap_vectors(n, u(1, l), 1, u(1, l+1), 1)
    l = l + 1
    goto 590
600 continue
    iter = 0
    m = m - 1
610 continue
    goto 360
620 continue

end subroutine dsvdc

! Purpose:
!
! Computes constant times a vector plus a vector.
!
! Jack dongarra, linpack, 3/11/78.
! Modified 12/3/93, array(1) declarations changed to array(*)
! Modified 01/27/17, whole array operations to aid compiler optimization
!
subroutine daxpy(n, da, dx, incx, dy, incy)

    ! Dummy arguments
    integer(ip), intent(in)    :: n
    real(wp),    intent(in)    :: da
    real(wp),    intent(in)    :: dx(*)
    integer(ip), intent(in)    :: incx
    real(wp),    intent(inout) :: dy(*)
    integer(ip), intent(in)    :: incy

    associate( &
        x => dx(1:n:incx), &
        y => dy(1:n:incy) &
        )
        y = y + da * x
    end associate

end subroutine daxpy

! Purpose:
!
! Forms the dot product of two vectors.
! Jack dongarra, linpack, 3/11/78.
! Modified 12/3/93, array(1) declarations changed to array(*)
! Modified 01/27/17, the function now wraps around the intrinsic dot_product
!
pure function get_dot_product(n, dx, incx, dy, incy) &
    result (return_value)

    ! Dummy arguments
    integer(ip), intent(in) :: n
    real(wp),    intent(in) :: dx(*)
    integer(ip), intent(in) :: incx
    real(wp),    intent(in) :: dy(*)
    integer(ip), intent(in) :: incy
    real(wp)                :: return_value

    return_value = dot_product(dx(:n:incx), dy(:n:incy))

end function get_dot_product

! Purpose:
!
!  Returns the euclidean norm of a vector via the function
!  name, so that
!
!     return_value := sqrt( transpose(x) * x )
!
! This version written on 25-October-1982.
! Modified on 14-October-1993 to inline the call to DLASSQ.
! Sven Hammarling, Nag Ltd.
!
! Modified 01/27/17, function now wraps around the intrinsic norm2
!
pure function get_norm2(n, x, incx) &
    result (return_value)

    ! Dummy arguments
    integer(ip), intent(in) :: n
    real(wp),    intent(in) :: x(*)
    integer(ip), intent(in) :: incx
    real(wp)                :: return_value

    return_value = norm2(x(1:n:incx))

end function get_norm2

! Purpose:
!
! Applies a plane rotation.
! Jack dongarra, linpack, 3/11/78.
! Modified 12/3/93, array(1) declarations changed to array(*)
!
subroutine  apply_plane_rotation(n, dx, incx, dy, incy, c, s)

    ! Dummy arguments
    integer(ip), intent(in)    :: n
    real(wp),    intent(inout) :: dx(*)
    integer(ip), intent(in)    :: incx
    real(wp),    intent(inout) :: dy(*)
    integer(ip), intent(in)    :: incy
    real(wp),    intent(in)    :: c
    real(wp),    intent(in)    :: s

    associate( &
        x => dx(1:n:incx), &
        y => dy(1:n:incy) &
        )

        block
            real(wp) :: temp(size(x))

            temp = c * x + s * y
            y = c * y - s * x
            x = temp
        end block
    end associate

end subroutine apply_plane_rotation

! Purpose:
!
! Construct givens plane rotation.
! jack dongarra, linpack, 3/11/78.
!
subroutine construct_given_plane_rotation(da, db, c, s)

    real(wp), intent(inout) :: da, db
    real(wp), intent(out)   :: c, s

    ! Local variables
    real(wp) :: roe, scale_factor, r, z

    if ( abs(da) > abs(db) ) then
        roe = da
    else
        roe = db
    end if

    scale_factor = abs(da) + abs(db)

    if (scale_factor == ZERO) then
        c = ONE
        s = ZERO
        r = ZERO
        z = ZERO
    else
        r = scale_factor * hypot((da/scale_factor),(db/scale_factor))
        r = sign(ONE, roe) * r
        c = da/r
        s = db/r
        z = ONE
        if ( abs(da) > abs(db) ) z = s
        if ( abs(da) <= abs(db) .and. c /= ZERO ) z = ONE/c
    end if

    da = r
    db = z

end subroutine construct_given_plane_rotation

! Purpose:
!
! Scales a vector by a constant.
!
! Jack dongarra, linpack, 3/11/78.
! Modified 3/93 to return if incx <= 0.0
! Modified 12/3/93, array(1) declarations changed to array(*)
! Modified 01/27/17, uses array operations to aid compiler optimization
!
subroutine scale_vector_by_constant(n, da, dx, incx)

    ! Dummy arguments
    integer(ip), intent(in)    :: n
    real(wp),    intent(in)    :: da
    real(wp),    intent(inout) :: dx(*)
    integer(ip), intent(in)    :: incx

    dx(1:n:incx) = da * dx(1:n:incx)

end subroutine  scale_vector_by_constant

! Purpose:
!
! Interchanges two vectors.
! uses unrolled loops for increments equal one.
! Jack dongarra, linpack, 3/11/78.
! Modified 12/3/93, array(1) declarations changed to array(*)
! Modified 01/27/17, whole array operations to improve compiler optimization
!
subroutine swap_vectors(n, dx, incx, dy, incy)

    ! Dummy arguments
    integer(ip), intent(in)    :: n
    real(wp),    intent(inout) :: dx(*)
    integer(ip), intent(in)    :: incx
    real(wp),    intent(inout) :: dy(*)
    integer(ip), intent(in)    :: incy

    associate( &
        x => dx(1:n:incx), &
        y => dy(1:n:incy) &
        )

        block
            real(wp) :: temp(size(x))

            temp = x
            x = y
            y = temp
        end block
    end associate

end subroutine swap_vectors

end module module_shpe
