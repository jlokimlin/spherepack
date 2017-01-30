module type_ShallowWaterSolver

    use spherepack_library, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, &
        Regularsphere

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: ShallowWaterSolver

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    
    type, extends(RegularSphere), public :: ShallowWaterSolver
    contains
        ! Type-bound procedures
        procedure         :: spec_to_grid
        procedure         :: grid_to_spec
        procedure         :: get_vrtdivspec
        procedure         :: get_uv
        procedure, nopass :: atanxy
        procedure, nopass :: get_initial_velocity
        procedure, nopass :: sine_transform
        procedure, nopass :: cosine_transform
    end type ShallowWaterSolver

contains

    ! Purpose:
    !
    ! Computes the initial unrotated longitudinal velocity
    ! See: section 3.3
    pure function get_initial_velocity(amp, thetad) result (return_value)

        ! Dummy arguments
        real(wp), intent(in) :: amp
        real(wp), intent(in) :: thetad
        real(wp)             :: return_value

        ! Local variables
        real(wp)            :: thetab, thetae, xe, x

        thetab = -PI/6
        thetae = PI/2
        xe = 3.0e-1_wp

        x = xe*(thetad-thetab)/(thetae-thetab)

        return_value = ZERO

        if (x <= ZERO .or. x >= xe) return

        return_value = amp*exp(-ONE /x-ONE/(xe-x)+4.0_wp/xe)

    end function get_initial_velocity

    pure function atanxy(x,y) &
        result (return_value)

        ! Dummy arguments
        real(wp), intent(in) :: x
        real(wp), intent(in) :: y
        real(wp)             :: return_value

        return_value = ZERO

        if (x == ZERO .and. y == ZERO) return

        return_value = atan2(y,x)

    end function atanxy

    subroutine sine_transform(x)

        ! Dummy arguments
        real(wp), intent(inout) :: x(:)

        ! Local variables
        integer(ip) :: i, j
        real(wp)    :: arg

        associate( n => size(x) )

            arg = PI/(n+1)

            block
                real(wp) :: w(n)

                do j=1, n
                    w(j) = ZERO
                    do i=1, n
                        w(j) = w(j)+x(i)*sin(real(i*j, kind=wp)*arg)
                    end do
                end do

                x = TWO * w
            end block
        end associate

    end subroutine sine_transform

    pure function cosine_transform(theta, cf) &
        result (return_value)

        ! Dummy arguments
        real(wp), intent(in) :: theta
        real(wp), intent(in) :: cf(:)
        real(wp)             :: return_value

        ! Local variables
        integer(ip) ::  i

        return_value = ZERO

        do i=1, size(cf)
            return_value = return_value+cf(i)*cos(real(i, kind=wp)*theta)
        end do

    end function cosine_transform

    ! Purpose:
    !
    ! Converts complex spectral coefficients (dataspec) to
    ! gridded data array (datagrid).
    subroutine grid_to_spec(self, datagrid, dataspec)

        ! Dummy arguments
        class(ShallowWaterSolver), intent(inout)  :: self
        real(wp),                  intent(in)     :: datagrid(:,:)
        complex(wp),               intent(out)    :: dataspec(:)

        ! Local variables
        integer(ip) :: n, m

        associate( &
            nlat => size(datagrid, dim=2), &
            nlon => size(datagrid, dim=1)  &
            )
            block
                real(wp) :: temp(nlat,nlon)

                !  Transpose data
                temp = transpose(datagrid)

                !  spherical harmonic analysis
                call self%perform_scalar_analysis(temp)

                associate( &
                    ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT, &
                    a => self%workspace%real_harmonic_coefficients, &
                    b => self%workspace%imaginary_harmonic_coefficients &
                    )

                    !  Fill complex array dataspec with result.
                    dataspec = HALF * cmplx( &
                        [((a(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                        [((b(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                        kind=wp)

                    !  Reset constants
                    a = ZERO
                    b = ZERO
                end associate
            end block
        end associate

    end subroutine grid_to_spec

    ! Purpose:
    !
    ! Converts complex spectral coefficients (dataspec) to
    ! gridded data array (datagrid).
    !
    subroutine spec_to_grid(self, dataspec, datagrid)

        ! Dummy arguments
        class(ShallowWaterSolver), intent(inout)  :: self
        real(wp),                  intent(inout)  :: datagrid(:,:)
        complex(wp),               intent(in)     :: dataspec(:)

        ! Local variables
        integer(ip) :: n, m, nm, i

        associate( &
            nlat => size(datagrid, dim=2), &
            nlon => size(datagrid, dim=1)  &
            )
            block
                real(wp) :: temp(nlat,nlon)

                !  fill two real arrays (a, b) with contents of dataspec.
                associate( &
                    indxn => self%INDEX_DEGREE_N, &
                    indxm => self%INDEX_ORDER_M, &
                    a => self%workspace%real_harmonic_coefficients, &
                    b => self%workspace%imaginary_harmonic_coefficients &
                    )

                    a = ZERO
                    b = ZERO
                    do nm=1, size(dataspec)
                        n = indxn(nm) ! Set degree n
                        m = indxm(nm) ! Set order m
                        a(m + 1, n + 1) = TWO * real(dataspec(nm))
                        b(m + 1, n + 1) = TWO * aimag(dataspec(nm))
                    end do

                    !  Perform spherical harmonic synthesis
                    call self%perform_scalar_synthesis(temp)

                    !  Reset coefficients
                    a = ZERO
                    b = ZERO
                end associate

                !  Transpose data
                datagrid = transpose(temp)
            end block
        end associate

    end subroutine spec_to_grid

    ! Purpose:
    !
    ! Calculate spectral coefficients of vorticity and divergence
    ! (vrtspec, divspec) given input gridded winds (ugrid, vgrid).
    !
    subroutine get_vrtdivspec(self, ugrid, vgrid, vrtspec, divspec)

        ! Dummy arguments
        class(ShallowWaterSolver), intent(inout)  :: self
        real(wp),                  intent(in)     :: ugrid(:,:)
        real(wp),                  intent(in)     :: vgrid(:,:)
        complex(wp),               intent(out)    :: vrtspec(:)
        complex(wp),               intent(out)    :: divspec(:)

        ! Local variables
        real(wp)    :: fn
        integer(ip) :: n, m ! Counters

        associate( &
            nlat => size(ugrid, dim=2), &
            nlon => size(ugrid, dim=1)  &
            )
            block
                real(wp) :: v(nlat,nlon)
                real(wp) :: w(nlat,nlon)
                real(wp) :: sqnn(nlat)

                !  Transpose data.
                !  minus sign to account for difference between
                !  mathematical and geophysical spherical coords
                v = -transpose(vgrid)
                w = transpose(ugrid)

                !  Calculate vector spherical harmonic analysis.
                call self%vector_analysis_from_spherical_components(v, w)

                !  Multiply vector harmonic coefficients of winds by
                !  appropriate factors to convert into vorticity and
                !  divergence coefficients.
                associate( &
                    ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT, &
                    nlat => self%NUMBER_OF_LATITUDES, &
                    rsphere => self%RADIUS_OF_SPHERE, &
                    a => self%workspace%real_harmonic_coefficients, &
                    b => self%workspace%imaginary_harmonic_coefficients, &
                    br => self%workspace%real_polar_harmonic_coefficients, &
                    bi => self%workspace%imaginary_polar_harmonic_coefficients, &
                    cr => self%workspace%real_azimuthal_harmonic_coefficients, &
                    ci => self%workspace%imaginary_azimuthal_harmonic_coefficients &
                    )

                    do n=1, nlat
                        fn = real(n - 1, kind=wp)
                        sqnn(n) = sqrt(fn * (fn + ONE))
                    end do

                    a = ZERO
                    b = ZERO
                    do n=1, nlat
                        a(:,n) = -(sqnn(n)/rsphere)*br(:,n)
                        b(:,n) = -(sqnn(n)/rsphere)*bi(:,n)
                    end do

                    divspec = HALF * cmplx( &
                        [((a(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                        [((b(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                        kind=wp)

                    a = ZERO
                    b = ZERO
                    do n=1, nlat
                        a(:,n) = (sqnn(n)/rsphere)*cr(:,n)
                        b(:,n) = (sqnn(n)/rsphere)*ci(:,n)
                    end do

                    vrtspec = HALF * cmplx( &
                        [((a(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                        [((b(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                        kind=wp)

                    !  Reset coefficients
                    a = ZERO
                    b = ZERO
                    br = ZERO
                    bi = ZERO
                    cr = ZERO
                    ci = ZERO
                end associate
            end block
        end associate

    end subroutine get_vrtdivspec

    ! Purpose:
    !
    ! Given spectral coefficients of vorticity and divergence
    ! (vrtspec, divspec) compute gridded winds (ugrid, vgrid)
    !
    subroutine get_uv(self, vrtspec, divspec, ugrid, vgrid)

        ! Dummy arguments
        class(ShallowWaterSolver), intent(inout) :: self
        complex(wp),               intent(in)    :: vrtspec(:)
        complex(wp),               intent(in)    :: divspec(:)
        real(wp),                  intent(out)   :: ugrid(:,:)
        real(wp),                  intent(out)   :: vgrid(:,:)

        ! Local variables
        real(wp)    :: fn
        integer(ip) :: n, m, nm, i ! Counters

        associate( &
            nlat => size(ugrid, dim=2),&
            nlon => size(ugrid, dim=1) &
            )
            block
                real(wp) :: v(nlat,nlon)
                real(wp) :: w(nlat,nlon)
                real(wp) :: isqnn(nlat)

                ! Multiply spectral coefficients of vorticity and divergence
                ! by appropriate factors to convert them into vector harmonic
                ! coefficients of winds.
                associate( &
                    indxn => self%INDEX_DEGREE_N, &
                    indxm => self%INDEX_ORDER_M, &
                    nlat => self%NUMBER_OF_LATITUDES, &
                    rsphere => self%RADIUS_OF_SPHERE, &
                    a => self%workspace%real_harmonic_coefficients, &
                    b => self%workspace%imaginary_harmonic_coefficients, &
                    br => self%workspace%real_polar_harmonic_coefficients, &
                    bi => self%workspace%imaginary_polar_harmonic_coefficients, &
                    cr => self%workspace%real_azimuthal_harmonic_coefficients, &
                    ci => self%workspace%imaginary_azimuthal_harmonic_coefficients &
                    )

                    isqnn(1) = ZERO
                    do n=2, nlat
                        fn = real(n - 1, kind=wp)
                        isqnn(n) = rsphere/sqrt(fn*(fn+ONE))
                    end do

                    ! Preset real harmonic coefficients to 0.0
                    a = ZERO
                    b = ZERO

                    ! Preset polar coefficients to 0.0
                    br = ZERO
                    bi = ZERO

                    do nm=1, size(divspec)
                        n = indxn(nm) ! Set degree n
                        m = indxm(nm) ! Set order m
                        a(m + 1, n + 1) = -TWO * real(divspec(nm))
                        b(m + 1, n + 1) = -TWO * aimag(divspec(nm))
                    end do

                    do n=1, nlat
                        br(:,n) = isqnn(n)*a(:,n)
                        bi(:,n) = isqnn(n)*b(:,n)
                    end do

                    ! Preset azimuthal coefficients to 0.0
                    cr = ZERO
                    ci = ZERO

                    do nm=1, size(divspec)
                        n = indxn(nm) ! Set degree n
                        m = indxm(nm) ! Set order m
                        a(m + 1, n + 1) = TWO * real(vrtspec(nm))
                        b(m + 1, n + 1) = TWO * aimag(vrtspec(nm))
                    end do

                    do n=1, nlat
                        cr(:,n) = isqnn(n)*a(:,n)
                        ci(:,n) = isqnn(n)*b(:,n)
                    end do

                    ! Compute vector harmonic synthesis to get winds on grid.
                    call self%perform_vector_synthesis(v, w)

                    !  Reset coefficients
                    a = ZERO
                    b = ZERO
                    br = ZERO
                    bi = ZERO
                    cr = ZERO
                    ci = ZERO
                end associate

                ! Transpose data
                ! minus sign to account for differences
                ! between mathematical and geophysical spherical coords.
                vgrid = -transpose(v)
                ugrid = transpose(w)
            end block
        end associate

    end subroutine get_uv

end module type_ShallowWaterSolver
