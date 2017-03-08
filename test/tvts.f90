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
!
!     4/97
!
!     a program for testing all theta derivative subroutines
!     vtses, vtsec, vtsgs, vtsgc
!
!
!     (1) first set a valid vector field (v, w) in terms of x, y, z
!         cartesian coordinates
!
!     (2) analytically compute (vt, wt) from (1)
!
!     (3) compute (vt, wt) using vtses, vtsec, vtsgs, vtsgc and compare with (2)
!
program tvts

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use spherepack

    ! Explicit typing only
    implicit none

    real(wp) :: cosp
    real(wp) :: cost
    real(wp) :: dlat
    real(wp) :: dphi
    real(wp) :: dxdt
    real(wp) :: dydt
    real(wp) :: dzdt
    real(wp) :: emz
    real(wp) :: err2v
    real(wp) :: err2w
    real(wp) :: ex
    real(wp) :: ey
    real(wp) :: ez
    integer(ip) :: i
    integer(ip) :: icase
    
    integer(ip) :: error_flag
    integer(ip) :: j
    integer(ip) :: k
    real(wp) :: phi
    real(wp) :: sinp
    real(wp) :: sint
    real(wp) :: theta
    real(wp) :: x
    real(wp) :: y
    real(wp) :: z
    real(wp), allocatable :: wsave(:)
    !
    !     set dimensions with parameter statements
    !
    integer(ip), parameter              :: NLAT= 25, NLON= 19, NT = 3, ITYP = 0
    real(wp), dimension(NLAT, NLAT, NT) :: br, bi, cr, ci
    real(wp), dimension(NLAT)           :: gaussian_latitudes, gaussian_weights
    real(wp), dimension(NLAT, NLON, NT) :: v, w, vt, wt, exact_vt, exact_wt


    write (stdout, '(/a/)') '     tvts *** TEST RUN *** '

    !  Print dimension variables
    call iout(NLAT, "nlat")
    call iout(NLON, "nlon")
    call iout(NT, "  nt")

    ! Set equally spaced colatitude and longitude grid increments
    dphi = TWO_PI/NLON
    dlat = PI/(NLAT-1)

    ! Compute nlat-many gaussian latitudinal points
    call compute_gaussian_latitudes_and_weights(NLAT, gaussian_latitudes, gaussian_weights, error_flag)
    call name("compute_gaussian_latitudes_and_weights")
    call iout(error_flag, "error_flag = ")
    call vecout(gaussian_latitudes, "gaussian_latitudes", NLAT)
    !
    !     test all theta derivative subroutines
    !
    do icase=1, 4

        call name("*****************************************")

        ! Set vector field v, w and compute theta derivatives in (exact_vt, exact_wt)
        do k=1, NT
            do j=1, NLON
                phi = real(j - 1, kind=wp) * dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, NLAT
                    select case (icase)
                        case(0:2)
                            theta = real(i - 1, kind=wp) * dlat
                        case default
                            theta = gaussian_latitudes(i)
                    end select
                    cost = cos(theta)
                    sint = sin(theta)
                       !
                       !    set x, y, z and their theta derivatives at colatitude theta and longitude p
                       !
                    x = sint*cosp
                    dxdt = cost*cosp
                    y = sint*sinp
                    dydt = cost*sinp
                    z = cost
                    dzdt = -sint
                    !
                    !     set (v, w) field corresponding to stream function
                    !     S = exp(y)+exp(-z) and velocity potential function
                    !     P = exp(x)+exp(z)
                    !
                    ex = exp(x)
                    ey = exp(y)
                    ez = exp(z)
                    emz = exp(-z)
                    w(i, j, k) =-ex*sinp+emz*sint+ey*cost*sinp
                    v(i, j, k) =-ey*cosp-ez*sint+ex*cost*cosp
                    !
                    !     set theta derivatives differentiating w, v above
                    exact_wt(i, j, k) = -ex*dxdt*sinp+emz*(-dzdt*sint+cost) &
                        +ey*sinp*(dydt*cost-sint)
                    exact_vt(i, j, k) = -ey*dydt*cosp-ez*(dzdt*sint+cost) &
                        +ex*cosp*(dxdt*cost-sint)
                end do
            end do
        end do

        select case (icase)
            case (1)

                call name("testing vtsec and initialize_vtsec")

                call initialize_vhaec(NLAT, NLON, wsave, error_flag)
                call name("initialize_vhaec")
                call iout(error_flag, "error_flag = ")

                call vhaec(NLAT, NLON, ITYP, NT, v, w, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wsave, error_flag)
                call name("vhaec")
                call iout(error_flag, "error_flag = ")

                ! Now compute theta derivatives of v, w
                call initialize_vtsec(NLAT, NLON, wsave, error_flag)
                call name("initialize_vtsec")
                call iout(error_flag, "error_flag = ")
                call vtsec(NLAT, NLON, ITYP, NT, vt, wt, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wsave, error_flag)
                call name("vtsec")
                call iout(error_flag, "error_flag = ")
            case (2)

                call name("testing vtses and initialize_vtses")

                call initialize_vhaes(NLAT, NLON, wsave, error_flag)
                call name("initialize_vhaes")
                call iout(error_flag, "error_flag = ")

                call vhaes(NLAT, NLON, ITYP, NT, v, w, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wsave, error_flag)
                call name("vhaes")
                call iout(error_flag, "error_flag = ")

                ! Now compute theta derivatives of v, w
                call initialize_vtses(NLAT, NLON, wsave, error_flag)
                call name("initialize_vtses")
                call iout(error_flag, "error_flag = ")

                call vtses(NLAT, NLON, ITYP, NT, vt, wt, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wsave, error_flag)
                call name("vtses")
                call iout(error_flag, "error_flag = ")
            case (3)


                call name("testing vtsgc and initialize_vtsgc")

                call initialize_vhagc(NLAT, NLON, wsave, error_flag)
                call name("initialize_vhagc")
                call iout(error_flag, "error_flag = ")

                call vhagc(NLAT, NLON, ITYP, NT, v, w, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wsave, error_flag)
                call name("vhagc")
                call iout(error_flag, "error_flag = ")

                ! Now synthesize v, w from br, bi, cr, ci and compare with original
                call initialize_vtsgc(NLAT, NLON, wsave, error_flag)
                call name("initialize_vtsgc")
                call iout(error_flag, "error_flag = ")

                call vtsgc(NLAT, NLON, ITYP, NT, vt, wt, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wsave, error_flag)
                call name("vtsgc")
                call iout(error_flag, "error_flag = ")
            case (4)

                call name("testing vtsgs and initialize_vtsgs")

                call initialize_vhags(NLAT, NLON, wsave, error_flag)
                call name("initialize_vhags")
                call iout(error_flag, "error_flag = ")

                call vhags(NLAT, NLON, ITYP, NT, v, w, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wsave, error_flag)
                call name("vhags")
                call iout(error_flag, "error_flag = ")

                call initialize_vtsgs(NLAT, NLON, wsave, error_flag)
                call name("initialize_vtsgs")
                call iout(error_flag, "error_flag = ")

                call vtsgs(NLAT, NLON, ITYP, NT, vt, wt, NLAT, NLON, br, bi, cr, ci, NLAT, &
                    NLAT, wsave, error_flag)
                call name("vtsgs")
                call iout(error_flag, "error_flag = ")
        end select

        ! Compute discretization error in vt, wt
        err2v = norm2(vt- exact_vt)
        err2w = norm2(wt - exact_wt)

        ! Set and print least squares error in v, w
        call vout(err2v, "errv")
        call vout(err2w, "errw")
    end do

    ! Release memory
    deallocate (wsave)

end program tvts
