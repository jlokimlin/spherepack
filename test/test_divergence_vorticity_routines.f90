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
!     3/97
!
!     a program for testing all divergence, vorticity and idvt(ec, es, gc, gs) routines
!
!     (1) first set a valid vector field by setting a stream function sf and velocity
!         potential function sv as polys in x, y, z restricted to the sphere.  Then
!         derive (v, w) and dv, vt from sf and sv analytically by differentiation.
!         (see tvha.f)
!
!     (2) compute the coefficients br, bi, cr, ci of (v, w) using vector analysis
!
!     (3) compute the divergence and vorticity of (v, w) using div, vrt (es, ec, gc, gs)
!
!     (4) compare with divergence and vorticity from (1)
!
!     (5) invert dv, vt with idvt(ec, es, gc, gs) and compare with vector field from (1)
!
program tidvt

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use spherepack, only: &
        ip, & ! Integer precision
        wp, & ! Working precision
        Sphere, &
        Regularsphere, &
        GaussianSphere

    ! Explicit typing only
    implicit none

    ! Dictionary
    class(Sphere), allocatable :: solver

    ! Test gaussian grid
    allocate (GaussianSphere :: solver)
    call test_divergence_vorticity_routines(solver)
    deallocate (solver)

    ! Test regular grid
    allocate (RegularSphere :: solver)
    call test_divergence_vorticity_routines(solver)
    deallocate (solver)

contains

    subroutine test_divergence_vorticity_routines(sphere_type)

        ! Dummy arguments
        class(Sphere), intent(inout)  :: sphere_type

        ! Local variables
        integer(ip), parameter        :: NLONS = 16
        integer(ip), parameter        :: NLATS = 25
        integer(ip), parameter        :: NSYNTHS = 3
        integer(ip)                   :: i, j, k ! Counters
        real(wp)                      :: exact_polar_component(NLATS, NLONS, NSYNTHS)
        real(wp)                      :: exact_azimuthal_component(NLATS, NLONS, NSYNTHS)
        real(wp)                      :: exact_vorticity(NLATS, NLONS, NSYNTHS)
        real(wp)                      :: exact_divergence(NLATS, NLONS, NSYNTHS)
        real(wp)                      :: approximate_polar_component(NLATS, NLONS, NSYNTHS)
        real(wp)                      :: approximate_azimuthal_component(NLATS, NLONS, NSYNTHS)
        real(wp)                      :: approximate_vorticity(NLATS, NLONS, NSYNTHS)
        real(wp)                      :: approximate_divergence(NLATS, NLONS, NSYNTHS)
        character(len=:), allocatable :: previous_vorticity_error, previous_divergence_error
        character(len=:), allocatable :: previous_polar_inversion_error
        character(len=:), allocatable :: previous_azimuthal_inversion_error


        !  Set up workspace arrays
        select type(sphere_type)
            !
            !  For gaussian sphere
            !
            class is (GaussianSphere)

            !
            !  Initialize gaussian sphere object
            !
            call sphere_type%create(nlat=NLATS, nlon=NLONS)

            ! Allocate known error from previous platform
            allocate (previous_vorticity_error, source='     vorticity error     = 5.351275e-14')
            allocate (previous_divergence_error, source='     divergence error    = 7.149836e-14')
            allocate (previous_polar_inversion_error, source='     polar inversion error     = 1.776357e-15')
            allocate (previous_azimuthal_inversion_error, source='     azimuthal inversion error = 1.776357e-15')
            !
            !  For regular sphere
            !
            class is (RegularSphere)

            !
            !  Initialize regular sphere object
            !
            call sphere_type%create(nlat=NLATS, nlon=NLONS)

            ! Allocate known error from previous platform
            allocate (previous_vorticity_error, source='     vorticity error     = 3.019807e-14')
            allocate (previous_divergence_error, source='     divergence error    = 4.141132e-14')
            allocate (previous_polar_inversion_error, source='     polar inversion error     = 1.554312e-15')
            allocate (previous_azimuthal_inversion_error, source='     azimuthal inversion error = 1.110223e-15')
        end select

        !
        !  set scalar stream and velocity potential fields as polys in x, y, z
        !    and then set v, w from st, sv scalar fields
        !
        associate (&
            ve => exact_polar_component, &
            we => exact_azimuthal_component, &
            vte => exact_vorticity, &
            dve => exact_divergence, &
            radial => sphere_type%unit_vectors%radial, &
            theta => sphere_type%unit_vectors%polar, &
            phi => sphere_type%unit_vectors%azimuthal &
           )
            do k=1, NSYNTHS
                do j=1, NLONS
                    do i=1, NLATS
                        associate (&
                            x => radial(i, j)%x, & !sint*cosp
                            y => radial(i, j)%y, & !sint*sinp
                            z => radial(i, j)%z, & !cost
                            dxdt => theta(i, j)%x, &! cost*cosp
                            dxdp => -radial(i, j)%y, & !-sint*sinp
                            dydt => theta(i, j)%y, &! cost*sinp
                            dydp => radial(i, j)%x, & ! sint*cosp
                            dzdt => theta(i, j)%z, & ! -sint
                            dzdp => phi(i, j)%z, & ! 0.0
                            cosp => phi(i, j)%y, &
                            sinp => -phi(i, j)%x &
                           )
                            select case (k)
                                case (1)
                                    !
                                    !   sf = x
                                    !   sv = y
                                    !
                                    !   v = -1/sint*dstdp + dsvdt
                                    !
                                    !   w =  1/sint*dsvdp + dstdt
                                    !
                                    !   dv = 1/sint*[d(sint*v)/dt + dwdp]  = dvdt + ct/st*v + 1/st*dwdp
                                    !
                                    !   vt = 1/sint*[-dv/dp + d(sint*w)/dt) = dwdt + ct/st*w - 1/st*dvdp
                                    !
                                    ve(i, j, k) = sinp + dydt ! sinp + cost*sinp
                                    we(i, j, k) = cosp + dxdt !cosp + cost*cosp
                                    dve(i, j, k) = -2.0_wp * y !-2.0*sint*sinp
                                    vte(i, j, k) = -2.0_wp * x ! -2.0*sint*cosp
                                case (2)
                                    !
                                    !   sf = y
                                    !   sv = z
                                    !
                                    ve(i, j, k) = -cosp + dzdt !-cosp-sint
                                    we(i, j, k) = dydt ! cost*sinp
                                    vte(i, j, k) = -2.0_wp * y ! -2.*sint*sinp
                                    dve(i, j, k) = -2.0_wp * z !-2.*cost
                                case (3)
                                    !
                                    !   st = x
                                    !   sv = z
                                    !
                                    ve(i, j, k) = sinp + dzdt ! sinp - sint
                                    we(i, j, k) = dxdt ! cost*cosp
                                    vte(i, j, k) = -2.0_wp * x !-2.*sint*cosp
                                    dve(i, j, k) = -2.0_wp * z !-2.*cost
                            end select
                        end associate
                    end do
                end do
            end do
        end associate

        !
        !  Compute vt and dv from (ve, we)
        !
        do k=1, NSYNTHS
            associate (&
                ve => exact_polar_component(:, :, k), &
                we => exact_azimuthal_component(:, :, k), &
                vt => approximate_vorticity(:, :, k), &
                dv => approximate_divergence(:, :, k) &
               )
                !
                !  Compute vorticity and divergence
                !
                call sphere_type%get_vorticity_and_divergence_from_velocities(ve, we, vt, dv)
            end associate
        end do

        !
        !  Compute "error" in dv, vt
        !
        associate (&
            vt => approximate_vorticity, &
            dv => approximate_divergence, &
            vte => exact_vorticity,  &
            dve => exact_divergence &
           )
            associate (&
                err2vt => maxval(abs(vt-vte)), &
                err2dv => maxval(abs(dv-dve)) &
               )
                !
                !  Print earlier output from platform with 64-bit floating point
                !    arithmetic followed by the output from this computer
                !
                write( stdout, '(a)') ''
                write( stdout, '(a)') '     tidvt *** TEST RUN *** '
                write( stdout, '(a)') ''
                write( stdout, '(a)') '     grid type = '//sphere_type%grid%grid_type
                write( stdout, '(a)') '     Testing vorticity'
                write( stdout, '(2(a, i3))') '     nlat = ', NLATS, ' nlon = ', NLONS
                write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
                write( stdout, '(a)') previous_vorticity_error
                write( stdout, '(a)') previous_divergence_error
                write( stdout, '(a)') '     The output from your computer is: '
                write( stdout, '(a, 1pe15.6)') '     vorticity error     = ', err2vt
                write( stdout, '(a, 1pe15.6)') '     divergence error    = ', err2dv
                write( stdout, '(a)') ''
            end associate
        end associate

        !
        !  Now compute (v, w) inverting vte, dve
        !
        do k=1, NSYNTHS
            associate (&
                v => approximate_polar_component(:, :, k), &
                w => approximate_azimuthal_component(:, :, k), &
                vte => exact_vorticity(:, :, k), &
                dve => exact_divergence(:, :, k) &
               )
                call sphere_type%get_velocities_from_vorticity_and_divergence(vte, dve, v, w)
            end associate
        end do
        !
        !  compare this v, w with original
        !
        associate (&
            ve => exact_polar_component, &
            we => exact_azimuthal_component, &
            v => approximate_polar_component, &
            w => approximate_azimuthal_component &
           )
            associate (&
                err2v => maxval(abs(v-ve)), &
                err2w => maxval(abs(w-we)) &
               )
                !
                !  Print earlier output from platform with 64-bit floating point
                !    arithmetic followed by the output from this computer
                !
                write( stdout, '(a)') ''
                write( stdout, '(a)') '     test divergence and vorticity *** TEST RUN *** '
                write( stdout, '(a)') ''
                write( stdout, '(a)') '     grid type = '//sphere_type%grid%grid_type
                write( stdout, '(a)') '     Testing vorticity inversion'
                write( stdout, '(2(a, i3))') '     nlat = ', NLATS, ' nlon = ', NLONS
                write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
                write( stdout, '(a)') previous_polar_inversion_error
                write( stdout, '(a)') previous_azimuthal_inversion_error
                write( stdout, '(a)') '     The output from your computer is: '
                write( stdout, '(a, 1pe15.6)') '     polar inversion error     = ', err2v
                write( stdout, '(a, 1pe15.6)') '     azimuthal inversion error = ', err2w
                write( stdout, '(a)') ''
            end associate
        end associate

        !  Release memory
        deallocate (previous_vorticity_error)
        deallocate (previous_divergence_error)
        deallocate (previous_polar_inversion_error)
        deallocate (previous_azimuthal_inversion_error)
        call sphere_type%destroy()

    end subroutine test_divergence_vorticity_routines

end program tidvt
