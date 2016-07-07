module type_ShallowWaterSolver

    use spherepack_library, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, &
        Regularsphere

    ! Explicit typing only
    implicit none

    ! Declare derived data type
    type, extends (RegularSphere) :: ShallowWaterSolver
        !--------------------------------------------------------------
        ! Class variables
        !--------------------------------------------------------------
    contains
        !--------------------------------------------------------------
        ! Class methods
        !--------------------------------------------------------------
        procedure         :: spec_to_grid
        procedure         :: grid_to_spec
        procedure         :: get_vrtdivspec
        procedure         :: get_uv
        procedure, nopass :: atanxy
        procedure, nopass :: get_initial_velocity
        procedure, nopass :: sine_transform
        procedure, nopass :: cosine_transform
        !--------------------------------------------------------------
    end type ShallowWaterSolver



contains


    pure function get_initial_velocity(amp, thetad) result (return_value)
        !
        !     computes the initial unrotated longitudinal velocity
        !     see section 3.3
        !
        !------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------
        real (wp), intent (in) :: amp
        real (wp), intent (in) :: thetad
        real (wp)              :: return_value
        !------------------------------------------------------
        ! Dictionary: local variables
        !------------------------------------------------------
        real (wp)            :: thetab, thetae, xe, x
        !------------------------------------------------------

        thetab = -PI/6
        thetae = PI/2
        xe = 3.0e-1_wp

        x = xe*(thetad-thetab)/(thetae-thetab)

        return_value = 0.0_wp

        if (x <= 0.0_wp .or. x >= xe) return

        return_value = amp*exp(-1.0_wp /x-1.0_wp/(xe-x)+4.0_wp/xe)

    end function get_initial_velocity



    pure function atanxy(x,y) result (return_value)
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        real (wp), intent (in) :: x
        real (wp), intent (in) :: y
        real                   :: return_value
        !--------------------------------------------------------------

        return_value = 0.0_wp

        if (x == 0.0_wp .and. y == 0.0_wp) return

        return_value = atan2(y,x)

    end function atanxy



    subroutine sine_transform(n, x)
        !
        ! Purpose:
        !
        ! Computes the sine transform
        !
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)     :: n
        real (wp),    intent (in out) :: x(n)
        !--------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------
        integer (ip)           :: i, j
        real (wp)              :: arg
        real (wp), allocatable :: w(:)
        !--------------------------------------------------------------

        arg = PI/(n+1)

        ! Allocate memory
        allocate( w(n) )

        do j=1, n
            w(j) = 0.0_wp
            do i=1, n
                w(j) = w(j)+x(i)*sin(real(i*j, kind=wp)*arg)
            end do
        end do

        x = 2.0_wp * w

        ! Release memory
        deallocate( w )

    end subroutine sine_transform



    pure function cosine_transform(theta, n, cf) result (return_value)
        !
        ! Purpose:
        !
        ! Returns the cosine transform
        !
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        real (wp),    intent (in) :: theta
        integer (ip), intent (in) :: n
        real (wp),    intent (in) :: cf(n)
        real (wp)                 :: return_value
        !--------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------
        integer (ip) ::  i
        !--------------------------------------------------------------

        return_value = 0.0_wp

        do i=1, n
            return_value = return_value+cf(i)*cos(real(i, kind=wp)*theta)
        end do

    end function cosine_transform



    subroutine grid_to_spec(this, datagrid, dataspec)
        !
        ! Purpose:
        !
        ! Converts complex spectral coefficients (dataspec) to
        ! gridded data array (datagrid).
        !
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        class (ShallowWaterSolver), intent (in out) :: this
        real (wp),                  intent (in)     :: datagrid(:,:)
        complex (wp),               intent (out)    :: dataspec(:)
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        integer (ip)           :: n, m
        real (wp), allocatable :: temp(:,:)
        !--------------------------------------------------------------

        associate( &
            nlat => size(datagrid, dim=2),&
            nlon => size(datagrid, dim=1) &
            )
            !
            !==> Allocate memory
            !
            allocate( temp(nlat,nlon) )

        end associate

        !
        !==> Transpose data
        !
        temp = transpose(datagrid)

        !
        !==> spherical harmonic analysis
        !
        call this%perform_scalar_analysis(temp)

        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients &
            )

            !
            !==> Fill complex array dataspec with result.
            !
            dataspec = 0.5_wp * cmplx( &
                [((a(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                [((b(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                kind=wp)
            !
            !==> Reset constants
            !
            a = 0.0_wp
            b = 0.0_wp
        end associate
        !
        !==> Release memory
        !
        deallocate( temp )

    end subroutine grid_to_spec



    subroutine spec_to_grid(this, dataspec, datagrid)
        !
        ! Purpose:
        !
        ! Converts complex spectral coefficients (dataspec) to
        ! gridded data array (datagrid).
        !
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        class (ShallowWaterSolver), intent (in out) :: this
        real (wp),                  intent (in out) :: datagrid(:,:)
        complex (wp),               intent (in)     :: dataspec(:)
        !--------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------
        integer (ip)           :: n, m, nm, i
        real (wp), allocatable :: temp(:,:)
        !--------------------------------------------------------------

        associate( &
            nlat => size(datagrid, dim=2),&
            nlon => size(datagrid, dim=1) &
            )
            !
            !==> Allocate memory
            !
            allocate( temp(nlat,nlon) )

        end associate


        !
        !==> fill two real arrays (a, b) with contents of dataspec.
        !
        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients &
            )

            a = 0.0_wp
            b = 0.0_wp
            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    nm = sum([(i, i=ntrunc+1, ntrunc-m+3, -1)])+n-m+1
                    a(m, n) = 2.0_wp * real(dataspec(nm))
                    b(m, n) = 2.0_wp * aimag(dataspec(nm))
                end do
            end do

            !
            !==> Perform spherical harmonic synthesis
            !
            call this%perform_scalar_synthesis(temp)

            !
            !==> Reset coefficients
            !
            a = 0.0_wp
            b = 0.0_wp
        end associate
        !
        !==> Transpose data
        !
        datagrid = transpose(temp)

        !
        !==> Release memory
        !
        deallocate( temp )

    end subroutine spec_to_grid



    subroutine get_vrtdivspec(this, ugrid, vgrid, vrtspec, divspec)
        !
        ! Purpose:
        !
        ! Calculate spectral coefficients of vorticity and divergence
        ! (vrtspec, divspec) given input gridded winds (ugrid, vgrid).
        !
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        class (ShallowWaterSolver), intent (in out) :: this
        real (wp),                  intent (in)     :: ugrid(:,:)
        real (wp),                  intent (in)     :: vgrid(:,:)
        complex (wp),               intent (out)    :: vrtspec(:)
        complex (wp),               intent (out)    :: divspec(:)
        !--------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------
        real (wp)              :: fn
        real (wp), allocatable :: v(:,:), w(:,:), sqnn(:)
        integer (ip)           :: n, m !! Counters
        !--------------------------------------------------------------

        associate( &
            nlat => size(ugrid, dim=2),&
            nlon => size(ugrid, dim=1) &
            )
            !
            !==> Allocate memory
            !
            allocate( v(nlat,nlon) )
            allocate( w(nlat,nlon) )
            allocate( sqnn(nlat) )

        end associate

        !
        !==> Transpose data.
        !    minus sign to account for difference between
        !    mathematical and geophysical spherical coords
        !
        v = -transpose(vgrid)
        w = transpose(ugrid)

        !
        !==> Calculate vector spherical harmonic analysis.
        !
        call this%vector_analysis_from_spherical_components(v, w)

        !
        !==> Multiply vector harmonic coefficients of winds by
        !    appropriate factors to convert into vorticity and
        !    divergence coefficients.
        !
        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            nlat => this%NUMBER_OF_LATITUDES, &
            rsphere => this%RADIUS_OF_SPHERE, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            br => this%workspace%real_polar_harmonic_coefficients, &
            bi => this%workspace%imaginary_polar_harmonic_coefficients, &
            cr => this%workspace%real_azimuthal_harmonic_coefficients, &
            ci => this%workspace%imaginary_azimuthal_harmonic_coefficients &
            )

            do n=1, nlat
                fn = real(n - 1, kind=wp)
                sqnn(n) = sqrt(fn * (fn + 1.0_wp))
            end do

            a = 0.0_wp
            b = 0.0_wp
            do n=1, nlat
                a(:,n) = -(sqnn(n)/rsphere)*br(:,n)
                b(:,n) = -(sqnn(n)/rsphere)*bi(:,n)
            end do

            divspec = 0.5_wp * cmplx( &
                [((a(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                [((b(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                kind=wp)

            a = 0.0_wp
            b = 0.0_wp
            do n=1, nlat
                a(:,n) = (sqnn(n)/rsphere)*cr(:,n)
                b(:,n) = (sqnn(n)/rsphere)*ci(:,n)
            end do

            vrtspec = 0.5_wp * cmplx( &
                [((a(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                [((b(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                kind=wp)

            !
            !==> Reset coefficients
            !
            a = 0.0_wp
            b = 0.0_wp
            br = 0.0_wp
            bi = 0.0_wp
            cr = 0.0_wp
            ci = 0.0_wp
        end associate

        !
        !==> Release memory
        !
        deallocate( v )
        deallocate( w )
        deallocate( sqnn )

    end subroutine get_vrtdivspec



    subroutine get_uv(this, vrtspec, divspec, ugrid, vgrid)
        !
        ! Purpose:
        !
        ! Given spectral coefficients of vorticity and divergence
        ! (vrtspec, divspec) compute gridded winds (ugrid, vgrid)
        !
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        class (ShallowWaterSolver), intent (in out) :: this
        complex (wp),               intent (in)     :: vrtspec(:)
        complex (wp),               intent (in)     :: divspec(:)
        real (wp),                  intent (out)    :: ugrid(:,:)
        real (wp),                  intent (out)    :: vgrid(:,:)
        !--------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------
        real (wp)              :: fn
        real (wp), allocatable :: v(:,:), w(:,:), isqnn(:)
        integer (ip)           :: n, m, nm, i !! Counters
        !--------------------------------------------------------------

        associate( &
            nlat => size(ugrid, dim=2),&
            nlon => size(ugrid, dim=1) &
            )
            !
            !==> Allocate memory
            !
            allocate( v(nlat,nlon) )
            allocate( w(nlat,nlon) )
            allocate( isqnn(nlat) )

        end associate

        !
        !==> multiply spectral coefficients of vorticity and divergence
        !    by appropriate factors to convert them into vector harmonic
        !    coefficients of winds.
        !
        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            nlat => this%NUMBER_OF_LATITUDES, &
            rsphere => this%RADIUS_OF_SPHERE, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            br => this%workspace%real_polar_harmonic_coefficients, &
            bi => this%workspace%imaginary_polar_harmonic_coefficients, &
            cr => this%workspace%real_azimuthal_harmonic_coefficients, &
            ci => this%workspace%imaginary_azimuthal_harmonic_coefficients &
            )

            isqnn(1) = 0.0_wp
            do n=2, nlat
                fn = real(n - 1, wp)
                isqnn(n) = rsphere/sqrt(fn*(fn+1.0_wp))
            end do

            a = 0.0_wp
            b = 0.0_wp
            br = 0.0_wp
            bi = 0.0_wp

            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    nm = sum([(i, i=ntrunc+1, ntrunc-m+3, -1)])+n-m+1
                    a(m, n) = -2.0_wp * real(divspec(nm))
                    b(m, n) = -2.0_wp * aimag(divspec(nm))
                end do
            end do

            do n=1, nlat
                br(:,n) = isqnn(n)*a(:,n)
                bi(:,n) = isqnn(n)*b(:,n)
            end do

            cr = 0.0_wp
            ci = 0.0_wp
            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    nm = sum([(i, i=ntrunc+1, ntrunc-m+3, -1)])+n-m+1
                    a(m, n) = 2.0_wp * real(vrtspec(nm))
                    b(m, n) = 2.0_wp * aimag(vrtspec(nm))
                end do
            end do

            do n=1, nlat
                cr(:,n) = isqnn(n)*a(:,n)
                ci(:,n) = isqnn(n)*b(:,n)
            end do



            !
            !==> compute vector harmonic synthesis to get winds on grid.
            !
            call this%perform_vector_synthesis(v, w)

            !
            !==> Reset coefficients
            !
            a = 0.0_wp
            b = 0.0_wp
            br = 0.0_wp
            bi = 0.0_wp
            cr = 0.0_wp
            ci = 0.0_wp

        end associate
        !
        !==> Transpose data
        !    minus sign to account for differences
        !    between mathematical and geophysical spherical coords.
        !
        vgrid = -transpose(v)
        ugrid = transpose(w)

        !
        !==> Release memory
        !
        deallocate( v )
        deallocate( w )
        deallocate( isqnn )

    end subroutine get_uv



end module type_ShallowWaterSolver
