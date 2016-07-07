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
!
!     11/96
!
!     a program for testing all vector analysis and synthesis subroutines
!
!     (1) first a scalar stream function and a scalar velocity potential function
!         are set in st,sv by restricting polys in x,y,z to the sphere surface
!
!     (2) the vector vield (v,w) is set by analytically differenting the scalar fields in (1)
!         using the standard formula relating a vector field and the stream and velocity
!         potential scalar fields in colatitude X longitude spherical coordinates
!
!          v = -1/sin(theta)*d(st)/dphi + d(sv)/dtheta
!
!          w =  1/sin(theta)*d(sv)/dphi + d(st)/dtheta
!
!     (3) a vector analysis is performed on (v,w)
!
!     (4) a vector synthesis is performed using coeffs from (3)
!
!     (5) the synthesized vector field from (4) is compared with the vector field from (2)
!
!     note:  vhaec,vhaes,vhagc,vhags,vhsec,vhses,vhsgc,vhsgs are all tested!!
!
program tvha

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
    class (Sphere), allocatable :: sphere_dat
    !----------------------------------------------------------------------

    !
    !==> Test gaussian case
    !
    allocate( GaussianSphere :: sphere_dat )

    call test_vector_analysis_and_synthesis_routines(sphere_dat)

    deallocate( sphere_dat )

    !
    !==> Test regular case
    !
    allocate( RegularSphere :: sphere_dat )

    call test_vector_analysis_and_synthesis_routines(sphere_dat)

    deallocate( sphere_dat )



contains




    subroutine test_vector_analysis_and_synthesis_routines(sphere_type)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: sphere_type
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip), parameter        :: NLONS = 19
        integer (ip), parameter        :: NLATS = 25
        integer (ip), parameter        :: NSYNTHS = 2
        integer (ip)                   :: i, j, k !! Counters
        real (wp)                      :: polar_component(NLATS,NLONS,NSYNTHS)
        real (wp)                      :: azimuthal_component(NLATS,NLONS,NSYNTHS)
        real (wp)                      :: synthesized_polar(NLATS,NLONS,NSYNTHS)
        real (wp)                      :: synthesized_azimuthal(NLATS,NLONS,NSYNTHS)
        character (len=:), allocatable :: previous_polar_error, previous_azimuthal_error
        !----------------------------------------------------------------------

        !
        !==> Set up workspace arrays
        !
        select type (sphere_type)
            type is (GaussianSphere)

            !  Initialize gaussian sphere object
            sphere_type = GaussianSphere(NLATS, NLONS)

            ! Allocate known error from previous platform
            allocate( previous_polar_error, source='     polar error     = 5.107026e-15' )
            allocate( previous_azimuthal_error, source='     azimuthal error = 9.325873e-15' )

            type is (RegularSphere)

            ! Initialize regular sphere
            sphere_type = RegularSphere(NLATS, NLONS)

            ! Allocate known error from previous platform
            allocate( previous_polar_error, source='     polar error     = 5.329071e-15' )
            allocate( previous_azimuthal_error, source='     azimuthal error = 7.771561e-15' )
        end select

        !
        !==> Set scalar stream and velocity potential fields as polys in x,y,z
        !    and then set v,w from st,sv scalar fields
        !
        associate( &
            ve => polar_component, &
            we => azimuthal_component, &
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
                                        dstdt => x*dydt+y*dxdt, &
                                        dsvdt => y*dzdt+z*dydt &
                                        )
                                        ve(i,j,k) = -(cosp*dydp+sinp*dxdp) + dsvdt
                                        we(i,j,k) = sinp*dzdp + dxdt + dstdt
                                    end associate
                                case (2)
                                    associate( &
                                        dstdt => x*dzdt+z*dxdt, &
                                        dsvdt => x*dydt+y*dxdt &
                                        )
                                        !
                                        !          v = -1/sin(theta)*d(st)/dphi + d(sv)/dtheta
                                        !
                                        !          w =  1/sin(theta)*d(sv)/dphi + d(st)/dtheta
                                        !
                                        ve(i,j,k) = z*sinp + dsvdt
                                        we(i,j,k) = cosp*dydp+ sinp*dxdp + dstdt
                                    end associate
                            end select
                        end associate
                    end do
                end do
            end do
        end associate


        !
        !==> Perform analysis then synthesis
        !
        do k = 1, NSYNTHS
            associate( &
                ve => polar_component(:,:,k), &
                we => azimuthal_component(:,:,k), &
                v => synthesized_polar(:,:,k), &
                w => synthesized_azimuthal(:,:,k) &
                )

                ! Analyse function into (real) coefficients
                call sphere_type%vector_analysis_from_spherical_components(ve, we)

                ! Synthesize function from (real) coefficients
                call sphere_type%perform_vector_synthesis(v, w)
            end associate
        end do


        !
        !==> Compute discretization error
        !
        associate( &
            ve => polar_component, &
            we => azimuthal_component, &
            v => synthesized_polar, &
            w => synthesized_azimuthal &
            )
            associate( &
                err2v => maxval(abs(v - ve)), &
                err2w => maxval(abs(w - we)) &
                )
                !
                !==> Print earlier output from platform with 64-bit floating point
                !    arithmetic followed by the output from this computer
                !
                write( stdout, '(a)') ''
                write( stdout, '(a)') '     tvha *** TEST RUN *** '
                write( stdout, '(a)') ''
                write( stdout, '(a)') '     grid type = '//sphere_type%grid%grid_type
                write( stdout, '(a)') '     Testing vector analysis and synthesis'
                write( stdout, '(2(A,I2))') '     nlat = ', NLATS,' nlon = ', NLONS
                write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
                write( stdout, '(a)') previous_polar_error
                write( stdout, '(a)') previous_azimuthal_error
                write( stdout, '(a)') '     The output from your computer is: '
                write( stdout, '(a,1pe15.6)') '     polar error     = ', err2v
                write( stdout, '(a,1pe15.6)') '     azimuthal error = ', err2w
                write( stdout, '(a)' ) ''
            end associate
        end associate
        !
        !==> Release memory
        !
        call sphere_type%destroy()
        deallocate( previous_polar_error )
        deallocate( previous_azimuthal_error )

    end subroutine test_vector_analysis_and_synthesis_routines

end program tvha
