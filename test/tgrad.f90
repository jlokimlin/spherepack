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
!     12/96
!
!     a program for testing all gradient and inverse gradient routines
!
!     (1) first a scalar field is set in st by restricting a poly in x,y,z
!         to the sphere surface
!
!     (2) a scalar analysis is used to compute the coefs a,b of st
!
!     (3) a,b are input to the various gradient routines to compute a vector field
!         (v,w)
!
!     (4) the vector field (v,w) is compared with the gradient of st obtained
!         analytically
!
!     (5) the inverse gradient of (v,w) is computed and compared with (1)
!
program tgrad

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT

    use spherepack_library, only: &
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

    call test_gradient_routines(gaussian_sphere)
    call test_gradient_routines(regular_sphere)


contains

    subroutine test_gradient_routines(sphere_type )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: sphere_type
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip), parameter        :: NLATS = 33
        integer (ip), parameter        :: NLONS = 18
        integer (ip), parameter        :: NSYNTHS = 4
        integer (ip)                   :: i, j, k !! Counters
        real (wp)                      :: original_polar_gradient_component(NLATS,NLONS,NSYNTHS)
        real (wp)                      :: original_azimuthal_gradient_component(NLATS,NLONS,NSYNTHS)
        real (wp)                      :: original_scalar_function(NLATS,NLONS,NSYNTHS)
        real (wp)                      :: approximate_polar_gradient_component(NLATS,NLONS,NSYNTHS)
        real (wp)                      :: approximate_azimuthal_gradient_component(NLATS,NLONS,NSYNTHS)
        real (wp)                      :: approximate_scalar_function(NLATS,NLONS,NSYNTHS)
        character (len=:), allocatable :: previous_gradient_inversion_error
        character (len=:), allocatable :: previous_polar_gradient_error
        character (len=:), allocatable :: previous_azimuthal_gradient_error
        !----------------------------------------------------------------------

        !
        !==> Set up workspace arrays
        !
        select type (sphere_type)
            !
            !==> For gaussian sphere
            !
            class is (GaussianSphere)

            !  Initialize gaussian sphere object
            call sphere_type%create(nlat=NLATS, nlon=NLONS)

            ! Allocate known error from previous platform
            allocate( previous_polar_gradient_error, source='     polar inversion error     = 2.953193e-14' )
            allocate( previous_azimuthal_gradient_error, source='     azimuthal inversion error = 1.720846e-14' )
            allocate( previous_gradient_inversion_error, source='     gradient error     = 6.106227e-16' )
            !
            !==> For regular sphere
            !
            class is (RegularSphere)

            ! Initialize regular sphere
            call sphere_type%create(nlat=NLATS, nlon=NLONS)

            ! Allocate known error from previous platform

            allocate( previous_polar_gradient_error, source='     polar inversion error     = 2.126077e-14' )
            allocate( previous_azimuthal_gradient_error, source='     azimuthal inversion error = 3.663736e-15' )
            allocate( previous_gradient_inversion_error, source='     gradient error     = 6.661338e-16' )
        end select

        !
        !==> set scalar stream and velocity potential fields as polys in x,y,z
        !    and then set v,w from st,sv scalar fields
        !
        associate( &
            ve => original_polar_gradient_component, &
            we => original_azimuthal_gradient_component, &
            sfe => original_scalar_function, &
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
                                    associate( &
                                        dsfdt => x*dydt+y*dxdt, &
                                        dsfdp => x*dydp+y*dxdp &
                                        )
                                        sfe(i,j,k) = x*y
                                        ve(i,j,k) = dsfdt
                                        we(i,j,k) = (cosp*dydp+sinp*dxdp)
                                    end associate
                                case (2)
                                    associate( &
                                        dsfdp => x*dzdp+z*dxdp, &
                                        dsfdt => x*dzdt+z*dxdt &
                                        )
                                        sfe(i,j,k) = x*z
                                        ve(i,j,k) = dsfdt
                                        we(i,j,k) = cosp*dzdp-z*sinp
                                    end associate
                                case (3)
                                    associate( &
                                        dsfdt => y*dzdt + z*dydt, &
                                        dsfdp => y*dzdp+ z*dydp &
                                        )
                                        sfe(i,j,k) = y*z
                                        ve(i,j,k) = dsfdt
                                        we(i,j,k) = sinp*dzdp + z*cosp
                                    end associate
                                case (4)
                                    associate( &
                                        dsfdt => x*y*dzdt + x*z*dydt + y*z*dxdt, &
                                        dsfdp => x*y*dzdp + x*z*dydp + y*z*dxdp &
                                        )
                                        sfe(i,j,k) = x*y*z
                                        ve(i,j,k) = dsfdt
                                        we(i,j,k) = cosp*y*dzdp+cosp*z*dydp+sinp*z*dxdp
                                    end associate
                            end select
                        end associate
                    end do
                end do
            end do
        end associate

        !
        !==> Compute gradient
        !
        do k=1, NSYNTHS
            associate( &
                v => approximate_polar_gradient_component(:,:,k), &
                w => approximate_azimuthal_gradient_component(:,:,k), &
                sfe => original_scalar_function(:,:,k) &
                )
                call sphere_type%get_gradient(sfe, v, w)
            end associate
        end do

        !
        !==> compare this v,w with original
        !
        associate( &
            ve => original_polar_gradient_component, &
            we => original_azimuthal_gradient_component, &
            v => approximate_polar_gradient_component, &
            w => approximate_azimuthal_gradient_component &
            )
            associate( &
                err2v => maxval(abs(v-ve)), &
                err2w => maxval(abs(w-we)) &
                )
                !
                !==> Print earlier output from platform with 64-bit floating point
                !    arithmetic followed by the output from this computer
                !
                write( stdout, '(A)') ''
                write( stdout, '(A)') '     tgrad *** TEST RUN *** '
                write( stdout, '(A)') ''
                write( stdout, '(A)') '     grid type = '//sphere_type%grid%grid_type
                write( stdout, '(A)') '     Testing gradient'
                write( stdout, '(2(A,I2))') '     nlat = ', NLATS,' nlon = ', NLONS
                write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
                write( stdout, '(A)') previous_polar_gradient_error
                write( stdout, '(A)') previous_azimuthal_gradient_error
                write( stdout, '(A)') '     The output from your computer is: '
                write( stdout, '(A,1pe15.6)') '     polar gradient error     = ', err2v
                write( stdout, '(A,1pe15.6)') '     azimuthal gradient error = ', err2w
                write( stdout, '(A)' ) ''
            end associate
        end associate

        !
        !==> Now recompute sf inverting (ve, we)
        !
        do k=1, NSYNTHS
            associate( &
                ve => original_polar_gradient_component(:,:,k), &
                we => original_azimuthal_gradient_component(:,:,k), &
                sf => approximate_scalar_function(:,:,k) &
                )
                call sphere_type%invert_gradient(ve, we, sf)
            end associate
        end do
        !
        !==> Compute gradient inversion error
        !
        associate( &
            sf => approximate_scalar_function, &
            sfe => original_scalar_function &
            )
            associate( err2 => maxval(abs(sf-sfe)) )
                !
                !==> Print earlier output from platform with 64-bit floating point
                !    arithmetic followed by the output from this computer
                !
                write( stdout, '(A)') ''
                write( stdout, '(A)') '     tgrad *** TEST RUN *** '
                write( stdout, '(A)') ''
                write( stdout, '(A)') '     grid type = '//sphere_type%grid%grid_type
                write( stdout, '(A)') '     Testing gradient inversion'
                write( stdout, '(2(A,I2))') '     nlat = ', NLATS,' nlon = ', NLONS
                write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
                write( stdout, '(A)') previous_gradient_inversion_error
                write( stdout, '(A)') '     The output from your computer is: '
                write( stdout, '(A,1pe15.6)') '     gradient inversion error     = ', err2
                write( stdout, '(A)' ) ''
            end associate
        end associate
        !
        !==> Release memory
        !
        deallocate( previous_gradient_inversion_error )
        deallocate( previous_polar_gradient_error )
        deallocate( previous_azimuthal_gradient_error )
        call sphere_type%destroy()

    end subroutine test_gradient_routines

end program tgrad
