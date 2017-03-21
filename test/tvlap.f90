!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                         Spherepack                            *
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
!     1/976
!
!     a program for testing all vector laplacian and its inverse
!
!     (1) first set the vector function (rotation about equator) using
!         v = cos(phi) and w = -cos(theta)*sin(phi)
!
!     (2) set vector laplacian ananlytically
!         vlap = -2.*cos(phi)=-2.*v, wlap = -2.*w
!         (i.e., L(v, w) = -2.*(v, w) so (v, w) is an eigenfunction for the
!         vector Laplacian with eigenvalue -2.
!
!     (3) compute the coefficients br, bi, cr, ci of (v, w) using vector analysis
!
!     (3) compute the vector laplacian of (v, w) using vlapec, vlapes, vlapgc, vlapgs
!
!     (4) compare (3) with (2)
!
!     (5) invert (4) and compare with (v, w)
!
program tvlap

    use spherepack

    ! Explicit typing only
    implicit none

    real(wp) :: cosp
    real(wp) :: cost
    real(wp) :: dlat
    real(wp) :: dphi
    real(wp) :: err2v
    real(wp) :: err2w
    integer(ip) :: i
    integer(ip) :: solver_case
    
    integer(ip) :: error_flag
    integer(ip) :: j
    integer(ip) :: k
    
    real(wp) :: phi

    real(wp) :: sinp
    real(wp) :: sint
    real(wp) :: theta

    ! Set dimensions with parameter statements
    integer(ip), parameter :: nlat = 29, nlon = 16, nt = 1
    integer(ip), parameter :: mdbc = min(nlat, nlon/2)
    integer(ip), parameter :: isym = 0, ityp = 0
    real(wp), allocatable  :: wavetable(:)
    real(wp), dimension(mdbc, nlat, nt) :: br, bi, cr, ci
    real(wp), dimension(nlat)           :: gaussian_latitudes, gaussian_weights
    real(wp), dimension(nlat, nlon, nt) :: v, w, ve, we, vlap, wlap, vlape, wlape
    real(wp), parameter                 :: ZERO = 0.0_wp, TWO = 2.0_wp

    call iout(nlat, "nlat")
    call iout(nlon, "nlon")
    call iout(nt, "  nt")

    ! Set equally spaced colatitude and longitude increments
    dphi = TWO_PI/nlon
    dlat = PI/(nlat-1)

    call compute_gaussian_latitudes_and_weights(nlat, &
        gaussian_latitudes, gaussian_weights, error_flag)

    call name("compute_gaussian_latitudes_and_weights")
    call iout(error_flag, " error_flag")
    call vecout(gaussian_latitudes, "gaussian_latitudes", nlat)
    !
    !     test all divergence and inverse divergence subroutines
    !
    do solver_case=1, 4
        call name("****")
        call name("****")
        call iout(solver_case, "icas")

        ! Set vector field v, w
        do k=1, nt
            do j=1, nlon
                phi = real(j - 1, kind=wp) * dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    select case (solver_case)
                        case (0:2)
                            theta = real(i - 1, kind=wp) * dlat
                        case default
                            theta = gaussian_latitudes(i)
                    end select
                    cost = cos(theta)
                    sint = sin(theta)
                    if (k==1) then
                        v(i, j, k) = cosp
                        w(i, j, k) = -cost*sinp
                        vlape(i, j, k) = -TWO * v(i, j, k)
                        wlape(i, j, k) = -TWO * w(i, j, k)
                    end if
                end do
            end do
        end do

        ! Store exact solution
        ve = v
        we = w

        select case (solver_case)
            case (1)

                call name("**ec")
                !
                !     analyze vector field
                !
                call initialize_vhaec(nlat, nlon, wavetable, error_flag)
                call name("vhai")
                call iout(error_flag, "error_flag = ")

                call vhaec(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, mdbc, &
                    nlat, wavetable, error_flag)
                call name("vha ")
                call iout(error_flag, "error_flag = ")

                !     compute vector laplacian
                call initialize_vhsec(nlat, nlon, wavetable, error_flag)
                call name("vhsi")
                call iout(error_flag, "error_flag = ")

                call vlapec(nlat, nlon, ityp, nt, vlap, wlap, nlat, nlon, br, bi, &
                    cr, ci, mdbc, nlat, wavetable, error_flag)
                call name("vlap")
                call iout(error_flag, "error_flag = ")


            case (2)

                call name("**es")
                !
                !     analyze vector field
                !
                call initialize_vhaes(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhaes")
                call iout(error_flag, "error_flag = ")

                call vhaes(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, cr, ci, mdbc, &
                    nlat, wavetable, error_flag)
                call name("vhaes")
                call iout(error_flag, "error_flag = ")

                call initialize_vhses(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhses")
                call iout(error_flag, "error_flag = ")

                call vlapes(nlat, nlon, isym, nt, vlap, wlap, nlat, nlon, br, bi, &
                    cr, ci, mdbc, nlat, wavetable, error_flag)
                call name("vlapes")
                call iout(error_flag, "error_flag = ")

            case (3)

                call name("**gc")
                !
                !     analyze vector field
                !
                call initialize_vhagc(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhagc")
                call iout(error_flag, "error_flag = ")

                call vhagc(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, mdbc, &
                    nlat, wavetable, error_flag)
                call name("vhagc")
                call iout(error_flag, "error_flag = ")

                !     compute vector laplacian
                !

                call initialize_vhsgc(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhsgc")
                call iout(error_flag, "error_flag = ")

                call vlapgc(nlat, nlon, ityp, nt, vlap, wlap, nlat, nlon, br, bi, &
                    cr, ci, mdbc, nlat, wavetable, error_flag)
                call name("vlapgc")
                call iout(error_flag, "error_flag = ")

            case (4)

                call name("**gs")
                !
                !     analyze vector field
                !
                call initialize_vhags(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhags")
                call iout(error_flag, "error_flag = ")

                call vhags(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, mdbc, &
                    nlat, wavetable, error_flag)
                call name("vhags")
                call iout(error_flag, "error_flag = ")

                !     compute vector laplacian
                call initialize_vhsgs(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhsgs")
                call iout(error_flag, "error_flag = ")

                call vlapgs(nlat, nlon, ityp, nt, vlap, wlap, nlat, nlon, br, bi, &
                    cr, ci, mdbc, nlat, wavetable, error_flag)
                call name("vlapgs")
                call iout(error_flag, "error_flag = ")

        end select

        ! Set and print least squares discretization error in vlap, wlap
        err2v = norm2(vlape - vlap)
        err2w = norm2(wlape - wlap)
        call vout(err2v, "errv")
        call vout(err2w, "errw")

        ! Now recompute (v, w) inverting (vlap, wlap) ivlap codes
        v = ZERO
        w = ZERO

        select case (solver_case)
            case (1)
                call name("analyze vector field (vlap, wlap)")
                !
                !     analyze vector field (vlap, wlap)
                !
                call initialize_vhaec(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhaec")
                call iout(error_flag, "error_flag = ")

                call vhaec(nlat, nlon, ityp, nt, vlap, wlap, nlat, nlon, &
                    br, bi, cr, ci, mdbc, nlat, wavetable, error_flag)
                call name("vhaec ")
                call iout(error_flag, "error_flag = ")

                call initialize_vhsec(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhsec")
                call iout(error_flag, "error_flag = ")

                call ivlapec(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, &
                    cr, ci, mdbc, nlat, wavetable, error_flag)
                call name("ivlapec")
                call iout(error_flag, "error_flag = ")

            case (2)

                call name("analyze vector field (vlap, wlap)")
                !
                !     analyze vector field (vlap, wlap)
                !
                call initialize_vhaes(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhaes")
                call iout(error_flag, "error_flag = ")

                call vhaes(nlat, nlon, isym, nt, vlap, wlap, nlat, nlon, &
                    br, bi, cr, ci, mdbc, nlat, wavetable, error_flag)
                call name("vhaes")
                call iout(error_flag, "error_flag = ")

                call initialize_vhses(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhses")
                call iout(error_flag, "error_flag = ")

                call ivlapes(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, &
                    cr, ci, mdbc, nlat, wavetable, error_flag)
                call name("ivlapes")
                call iout(error_flag, "error_flag = ")

            case (3)

                call name("analyze vector field (vlap, wlap)")

                !
                !     analyze vector field (vlap, wlap)
                !
                call initialize_vhagc(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhagc")
                call iout(error_flag, "error_flag = ")

                call vhagc(nlat, nlon, ityp, nt, vlap, wlap, nlat, nlon, &
                    br, bi, cr, ci, mdbc, nlat, wavetable, error_flag)
                call name("vhagc")
                call iout(error_flag, "error_flag = ")

                call initialize_vhsgc(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhsgc")
                call iout(error_flag, "error_flag = ")

                call ivlapgc(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, &
                    cr, ci, mdbc, nlat, wavetable, error_flag)
                call name("ivlapgc")
                call iout(error_flag, "error_flag = ")

            case default

                call name("analyze vector field (vlap, wlap)")

                !
                !     analyze vector field (vlap, wlap)
                !
                call initialize_vhags(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhags")
                call iout(error_flag, "error_flag = ")

                call vhags(nlat, nlon, ityp, nt, vlap, wlap, nlat, nlon, &
                    br, bi, cr, ci, mdbc, nlat, wavetable, error_flag)
                call name("vhags")
                call iout(error_flag, "error_flag = ")

                call initialize_vhsgs(nlat, nlon, wavetable, error_flag)
                call name("initialize_vhsgs")
                call iout(error_flag, "error_flag = ")

                call ivlapgs(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, &
                    cr, ci, mdbc, nlat, wavetable, error_flag)
                call name("ivlapgs")
                call iout(error_flag, "error_flag = ")

        end select

        ! Compare this v, w with original
        err2v = norm2(ve - v)
        err2w = norm2(we - w)
        call vout(err2v, "errv")
        call vout(err2w, "errw")
    end do

    ! Release memory
    deallocate (wavetable)

end program tvlap
