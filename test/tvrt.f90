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
!     *       A Package of Fortran77 Subroutines and Programs         *
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
!     6/98
!
!     a program for testing all vorticity and inverse vorticity routines
!
!
!     (1) first set a stream function and velocity potential scalar fields as
!         polys in x,y,z restricted to the sphere
!
!     (2) derive a vector field (v,w) from (1)
!
!     (3) compute the vorticity vt of (2) and compare with the vorticity
!         computed analytically
!
!     (4) compute vector field (ve,we) using br,bi,cr,ci from (v,w) with
!         br=bi=0.0
!
!     (5) invert the vorticity in (3) and compare with (4)
!

program tvrt

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT

    use modern_spherepack_library, only: &
        Sphere, &
        Regularsphere, &
        GaussianSphere

    ! Explicit typing only
    implicit none

    !----------------------------------------------------------------------
    ! Dictionary
    !----------------------------------------------------------------------
    type (GaussianSphere) :: gaussian_sphere
    type (RegularSphere)  :: regular_sphere
    !----------------------------------------------------------------------

    call test_vorticity_routines(gaussian_sphere)
    call test_vorticity_routines(regular_sphere)


contains

    subroutine test_vorticity_routines(sphere_type )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: sphere_type
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip), parameter        :: NLONS = 14
        integer (ip), parameter        :: NLATS = 24
        integer (ip), parameter        :: NSYNTHS = 3
        integer (ip)                   :: i, j, k !! Counters
        real (wp)                      :: exact_polar_component(NLATS, NLONS, NSYNTHS)
        real (wp)                      :: exact_azimuthal_component(NLATS, NLONS, NSYNTHS)
        real (wp)                      :: exact_vorticity(NLATS, NLONS, NSYNTHS)
        real (wp)                      :: approximate_polar_component(NLATS, NLONS, NSYNTHS)
        real (wp)                      :: approximate_azimuthal_component(NLATS, NLONS, NSYNTHS)
        real (wp)                      :: approximate_vorticity(NLATS, NLONS, NSYNTHS)
        real (wp)                      :: vorticity_error, polar_inversion_error, azimuthal_inversion_error
        character (len=:), allocatable :: previous_vorticity_error
        character (len=:), allocatable :: previous_polar_inversion_error
        character (len=:), allocatable :: previous_azimuthal_inversion_error
        !----------------------------------------------------------------------

        !
        !==> Set up workspace arrays
        !
        select type (sphere_type)
            !
            !==> For gaussian sphere
            !
            class is (GaussianSphere)

            !
            !==> Initialize gaussian sphere object
            !
            call sphere_type%create(nlat=NLATS, nlon=NLONS)

            ! Allocate known error from previous platform
            allocate( previous_vorticity_error, source='     vorticity error     = 8.246182e-14' )
            allocate( previous_polar_inversion_error, source='     polar inversion error     = 1.998401e-15' )
            allocate( previous_azimuthal_inversion_error, source='     azimuthal inversion error = 4.107825e-15' )
            !
            !==> For regular sphere
            !
            class is (RegularSphere)

            !
            !==> Initialize regular sphere object
            !
            call sphere_type%create(nlat=NLATS, nlon=NLONS)

            ! Allocate known error from previous platform
            allocate( previous_vorticity_error, source='     vorticity error     = 2.420286e-14' )
            allocate( previous_polar_inversion_error, source='     polar inversion error     = 1.110223e-15' )
            allocate( previous_azimuthal_inversion_error, source='     azimuthal inversion error = 1.887379e-15' )
        end select

        !
        !==> set scalar stream and velocity potential fields as polys in x,y,z
        !    and then set v,w from st,sv scalar fields
        !
        associate( &
            ve => exact_polar_component, &
            we => exact_azimuthal_component, &
            vte => exact_vorticity, &
            radial => sphere_type%unit_vectors%radial, &
            theta => sphere_type%unit_vectors%polar, &
            phi => sphere_type%unit_vectors%azimuthal &
            )
            do k=1, NSYNTHS
                do j=1, NLONS
                    do i=1, NLATS
                        associate( &
                            x => radial(i,j)%x, & !sint*cosp
                            y => radial(i,j)%y, & !sint*sinp
                            z => radial(i,j)%z, & !cost
                            dxdt => theta(i,j)%x, &! cost*cosp
                            dxdp => -radial(i,j)%y, & !-sint*sinp
                            dydt => theta(i,j)%y, &! cost*sinp
                            dydp => radial(i,j)%x, & ! sint*cosp
                            dzdt => theta(i,j)%z, & ! -sint
                            dzdp => phi(i,j)%z,& ! 0.0
                            cosp => phi(i,j)%y, &
                            sinp => -phi(i,j)%x &
                            )
                            select case (k)
                                case (1)
                                    ve(i,j,k) = sinp + dydt !sinp + cost*sinp
                                    we(i,j,k) = cosp + dxdt !cosp + cost*cosp
                                    vte(i,j,k) = -2.0_wp * x !-2.0*sint*cosp
                                case (2)
                                    ve(i,j,k) = -cosp + dzdt !-cosp-sint
                                    we(i,j,k) = dydt !cost*sinp
                                    vte(i,j,k) = -2.0_wp * y !-2.*sint*sinp
                                case (3)
                                    ve(i,j,k) = sinp + dzdt !sinp - sint
                                    we(i,j,k) = dxdt !cost*cosp
                                    vte(i,j,k) = -2.0_wp * x !-2.*sint*cosp
                            end select
                        end associate
                    end do
                end do
            end do
        end associate

        !
        !==> Compute vt from (ve,we)
        !
        do k=1, NSYNTHS
            associate( &
                ve => exact_polar_component(:,:,k), &
                we => exact_azimuthal_component(:,:,k), &
                vt => approximate_vorticity(:,:,k) &
                )
                call sphere_type%get_vorticity(ve, we, vt)
            end associate
        end do

        !
        !==> Compute vorticity error
        !
        associate( &
            err2 => vorticity_error, &
            vt => approximate_vorticity, &
            vte => exact_vorticity &
            )
            err2 = maxval(abs(vt-vte)) !norm2(vt - vte)/size(vt)
        end associate

        !
        !==> Print earlier output from platform with 64-bit floating point
        !    arithmetic followed by the output from this computer
        !
        write( stdout, '(A)') ''
        write( stdout, '(A)') '     tvrt *** TEST RUN *** '
        write( stdout, '(A)') ''
        write( stdout, '(A)') '     grid type = '//sphere_type%grid%grid_type
        write( stdout, '(A)') '     Testing vorticity'
        write( stdout, '(2(A,I2))') '     nlat = ', NLATS,' nlon = ', NLONS
        write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') previous_vorticity_error
        write( stdout, '(A)') '     The output from your computer is: '
        write( stdout, '(A,1pe15.6)') '     vorticity error     = ', vorticity_error
        write( stdout, '(A)' ) ''

        !
        !==> Now recompute (v,w) inverting vte
        !
        do k=1, NSYNTHS
            associate( &
                v => approximate_polar_component(:,:,k), &
                w => approximate_azimuthal_component(:,:,k), &
                vte => exact_vorticity(:,:,k) &
                )
                call sphere_type%invert_vorticity(vte, v, w)
            end associate
        end do

        !
        !==> compute vector field (ve,we) using br,bi,cr,ci from (ve,we) with
        !    br = bi = 0.0
        !
        do k = 1, NSYNTHS
            associate( &
                ve => exact_polar_component(:,:,k), &
                we => exact_azimuthal_component(:,:,k) &
                )
                !
                !==> Get polar coefficients (br, bi) and azimuthal coefficients (cr, ci) from (ve, we)
                !
                call sphere_type%vector_analysis_from_spherical_components(ve, we)
                !
                !==> Set polar coefficients to zero
                !
                associate( &
                    br => sphere_type%workspace%real_polar_harmonic_coefficients, &
                    bi => sphere_type%workspace%imaginary_polar_harmonic_coefficients &
                    )
                    br = 0.0_wp; bi = 0.0_wp
                end associate
                !
                !==> Synthesize azimuthal coefficients (cr, ci) into (ve, we)
                !
                call sphere_type%perform_vector_synthesis(ve, we)
            end associate
        end do
        !
        !==> compare this v,w with original
        !
        associate( &
            err2v => polar_inversion_error, &
            err2w => azimuthal_inversion_error, &
            ve => exact_polar_component, &
            we => exact_azimuthal_component, &
            v => approximate_polar_component, &
            w => approximate_azimuthal_component &
            )
            err2v = maxval( abs(v-ve) ) !norm2(v - ve)/size(v)
            err2w = maxval( abs(w-we) ) !norm2(w - we)/size(w)
        end associate

        !
        !==> Print earlier output from platform with 64-bit floating point
        !    arithmetic followed by the output from this computer
        !
        write( stdout, '(A)') ''
        write( stdout, '(A)') '     tvrt *** TEST RUN *** '
        write( stdout, '(A)') ''
        write( stdout, '(A)') '     grid type = '//sphere_type%grid%grid_type
        write( stdout, '(A)') '     Testing vorticity inversion'
        write( stdout, '(2(A,I2))') '     nlat = ', NLATS,' nlon = ', NLONS
        write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') previous_polar_inversion_error
        write( stdout, '(A)') previous_azimuthal_inversion_error
        write( stdout, '(A)') '     The output from your computer is: '
        write( stdout, '(A,1pe15.6)') '     polar inversion error     = ', polar_inversion_error
        write( stdout, '(A,1pe15.6)') '     azimuthal inversion error = ', azimuthal_inversion_error
        write( stdout, '(A)' ) ''

        !
        !==> Release memory
        !
        deallocate( previous_vorticity_error )
        deallocate( previous_polar_inversion_error )
        deallocate( previous_azimuthal_inversion_error )
        call sphere_type%destroy()

    end subroutine test_vorticity_routines

end program tvrt
