program test
    !
    ! Purpose:

    use, intrinsic :: iso_fortran_env, only: &
        REAL64, &
        INT32

    use type_vector_mod

    use type_sphere_mod

    ! Explicit typing only
    implicit none

    !--------------------------------------------------------------------------------
    ! Dictionary
    !--------------------------------------------------------------------------------
    integer, parameter    :: WP     = REAL64  !! 64 bit real
    integer, parameter    :: IP     = INT32    !! 32 bit integer

    !--------------------------------------------------------------------------------

    ! Test all the procedures
    call Test_all()

contains
    !
    !o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
    !
    subroutine Test_all()
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP), parameter    :: nlon = 128
        integer (IP), parameter    :: nlat = nlon/2 + 1
        type (sphere_t)            :: this
        !--------------------------------------------------------------------------------

        ! Create sphere object
        call this%Create( nlat, nlon )

        print *, 'Spherepack wrapper validation tests'
        print *, ' '
        print *, 'nlat = ', nlat, 'nlon = ', nlon

        ! Test all the subroutines
        call Test_scalar_analysis_and_synthesis( this )
        call Test_vector_analysis_and_synthesis( this )
        call Test_compute_surface_integral( this )
        call Test_invert_helmholtz( this )
        call Test_get_gradient( this )
        call Test_get_vorticity( this )
        call Test_get_rotation_operator( this )

        ! Destroy sphere object
        call this%Destroy()

    end subroutine Test_all
    !
    !o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
    !
    subroutine Test_scalar_analysis_and_synthesis( this )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)    :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                            :: nlat, nlon, k, l
        real (WP)                               :: real_error, complex_error
        type (vector_t)                         :: u
        real (WP), dimension (this%nlat, this%nlon) :: original_function
        real (WP), dimension (this%nlat, this%nlon) :: real_synthesized_function
        real (WP), dimension (this%nlat, this%nlon) :: complex_synthesized_function
        !--------------------------------------------------------------------------------

        ! Constants
        nlat = this%nlat
        nlon = this%nlon

        ! initialize arrays
        original_function = 0.0_WP

        do l = 1, nlon
            do k = 1, nlat

                u = this%radial_unit_vector(:,k,l)

                original_function(k,l) = exp(u%x + u%y + u%z)

            end do
        end do

        ! real case
        call this%Perform_scalar_analysis( original_function )

        call this%Perform_scalar_synthesis( real_synthesized_function )

        real_error = maxval(abs(original_function &
            - real_synthesized_function))

        ! complex case
        call this%Perform_complex_analysis( original_function )

        call this% Perform_complex_synthesis( complex_synthesized_function )

        complex_error = maxval(abs(original_function &
            - complex_synthesized_function))

        ! print errors to console
        print *, " "
        print *, "Test_scalar_analysis_and_synthesis"
        ! print real error to console
        print *, "real error"
        write (*,101) real_error
        ! print complex error to console
        print *, "complex error"
        write (*,101) complex_error
101     format (1x, "l-infinity error = ", es23.16 )
        print *, " "
        print *, "******************************************** "

    end subroutine Test_scalar_analysis_and_synthesis
    !
    !o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
    !
    subroutine Test_vector_analysis_and_synthesis( this )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)     :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)       :: nlat, nlon, k, l
        type (vector_t) :: u, phi, theta, vector_field
        real (WP)      :: polar_error, azimuthal_error
        real (WP), dimension (3,this%nlat, this%nlon) :: &
            vector_function
        real (WP), dimension (this%nlat, this%nlon) :: &
            original_polar_component, original_azimuthal_component, &
            approximate_polar_component, approximate_azimuthal_component
        !--------------------------------------------------------------------------------

        nlat = this%nlat
        nlon = this%nlon

        ! initialize arrays
        vector_function = 0.0_WP
        original_polar_component = 0.0_WP
        original_azimuthal_component = 0.0_WP

        ! compute the vector field that gives rise original components
        do l = 1, nlon
            do k = 1, nlat

                u = this%radial_unit_vector(:,k,l)
                theta = this%polar_unit_vector(:,k,l)
                phi = this%azimuthal_unit_vector(:,k,l)

                vector_field = [1.0E+3_WP, 1.0E+2_WP, 1.0E+1_WP]

                vector_function(:, k, l) = vector_field

                original_polar_component(k,l) = &
                    theta.dot.vector_field

                original_azimuthal_component(k,l) = &
                    phi.dot.vector_field
            end do
        end do

        ! analyze the vector function
        call this%Perform_vector_analysis( vector_function )

        ! synthesize the function from the coefficients
        call this%Perform_vector_synthesis( &
            approximate_polar_component, &
            approximate_azimuthal_component)

        ! set errors
        polar_error = &
            maxval(abs(original_polar_component &
            - approximate_polar_component))

        azimuthal_error = &
            maxval(abs(original_azimuthal_component &
            - approximate_azimuthal_component))

        ! print the errors to console
        print *, "Test_real_vector_analysis_and_synthesis"
        ! polar error
        print *, " "
        print *, "polar error"
        write (*,101) polar_error
        ! azimuthal error
        print *, " "
        print *, "azimuthal error"
        write (*,101) azimuthal_error
101     format (1x, "l-infinity error = ", es23.16 )
        print *, " "
        print *, "******************************************** "

    end subroutine Test_vector_analysis_and_synthesis
    !
    !o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
    !
    subroutine Test_compute_surface_integral( this )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)    :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        real (WP),dimension (this%nlat, this%nlon) :: scalar_function, constant_function
        integer (IP)                           :: nlat, nlon, k,l
        real (WP)                              :: LAMBDA, exact_value, approximate_value
        real (WP)                              :: constant_exact_value, constant_approximate_value
        type (vector_t)                        :: u, jhat
        real (WP), parameter                   :: FOUR_PI = 16.0_WP * atan( 1.0_WP ) !! To calculate integral
        !--------------------------------------------------------------------------------

        nlat = this%nlat
        nlon = this%nlon

        ! set the constant function
        constant_function = 1.0_WP

        ! initialize the array
        scalar_function = 0.0_WP

        LAMBDA = exp(1.0_WP)

        ! set the scalar function
        do k = 1, nlat
            do l = 1, nlon

                u = this%radial_unit_vector(:,k,l)
                jhat = [ 0.0_WP, 1.0_WP, 0.0_WP ]
                scalar_function(k,l) = exp( u.dot.jhat)

            end do
        end do

        ! compute the surface integral
        approximate_value = this%Compute_surface_integral( scalar_function )
        constant_approximate_value = this%Compute_surface_integral( constant_function )

        ! set the exact valuei
        constant_exact_value =  FOUR_PI

        exact_value = 14.768_WP

        ! print the error to the console
        write (*,101)  abs(approximate_value - exact_value)
        print *, " "
        write (*,101)  abs(constant_approximate_value - constant_exact_value)
101     format (1x, "surface_integral error = ", es23.16 )
        print *, " "
        print *, "******************************************** "


    end subroutine Test_Compute_surface_integral
    !
    !o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
    !
    subroutine Test_invert_helmholtz( this )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)    :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                            :: nlat, nlon, k, l
        type (vector_t)                         :: u
        real (WP)                               :: error
        real (WP), dimension (this%nlat, this%nlon) :: exact_solution
        real (WP), dimension (this%nlat, this%nlon) :: source_term
        real (WP), dimension (this%nlat, this%nlon) :: approximate_solution
        !--------------------------------------------------------------------------------

        nlat = this%nlat
        nlon = this%nlon

        ! initialize the scalar function and the exact solution
        do l = 1, nlon
            do k = 1, nlat

                u = this%radial_unit_vector(:,k,l)

                exact_solution(k,l) = (1.0_WP + u%x * u%y) * exp(u%z)

                source_term(k,l) = &
                    -(u%x * u%y*(u%z *u%z + 6.0_WP * (u%z + 1.0_WP)) &
                    + u%z * ( u%z + 2.0_WP)) * exp(u%z)

            end do

        end do

        ! solve helmholtz's equation
        call this%Invert_helmholtz( 1.0_WP, source_term, approximate_solution)

        ! set error
        error = &
            maxval(abs(exact_solution - approximate_solution))

        print *, " "
        print *, "Test_invert_helmholtz"
        write (*,101) error
101     format (1x, "l-infinity error = ", es23.16 )
        print *, " "
        print *, "******************************************** "

    end subroutine Test_invert_helmholtz
    !
    !o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
    !
    subroutine Test_get_gradient( this )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)    :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)      :: nlat, nlon, k, l
        real (WP)         :: csc, sint, cost, sinp, cosp, &
            polar_error, azimuthal_error
        type (vector_t)   :: u
        real (WP),dimension (this%nlat, this%nlon)  :: &
            scalar_function, &
            exact_polar_component, exact_azimuthal_component, &
            approximate_polar_component, approximate_azimuthal_component
        !--------------------------------------------------------------------------------

        nlat = this%nlat
        nlon = this%nlon

        ! initialize arrays
        scalar_function = 0.0_WP
        exact_polar_component = 0.0_WP
        exact_azimuthal_component = 0.0_WP

        do l = 1, nlon
            do k = 1, nlat

                sint = this%sint(k)
                cost = this%cost(k)

                sinp = this%sinp(l)
                cosp = this%cosp(l)

                u = this%radial_unit_vector(:,k,l)

                scalar_function(k,l) = exp( u%x + u%y + u%z )

                csc = 1.0_WP /sint

                exact_polar_component(k,l) = &
                    scalar_function(k,l) * (cost * cosp &
                    - sint + cost * sinp)

                exact_azimuthal_component(k,l) = &
                    scalar_function(k,l) * csc * (u%x - u%y)
            end do

        end do

        ! Calculate the gradient
        call this%Get_gradient( &
            scalar_function, &
            approximate_polar_component, &
            approximate_azimuthal_component )

        ! set errors
        polar_error = maxval(abs(exact_polar_component &
            - approximate_polar_component))

        azimuthal_error = maxval(abs(exact_azimuthal_component&
            - approximate_azimuthal_component))

        print *, " "
        print *, "Test_get_gradient"
        ! print the polar error
        print *, " "
        print *, "polar component"
        write (*,101) polar_error
        ! print the azimuthal error
        print *, " "
        print *, "azimuthal component"
        write (*,101) azimuthal_error
101     format (1x, "l-infinity error = ", es23.16 )
        print *, " "
        print *, "******************************************** "

    end subroutine Test_get_gradient
    !
    !o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
    !
    subroutine Test_get_vorticity( this )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)  :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)       :: nlat, nlon, k, l
        type (vector_t) :: u, omega, rotation_operator
        real (WP)      :: scalar_function, sint, cost, sinp, cosp, &
            cot, phi_derivative, theta_derivative, error
        real (WP), dimension (3,this%nlat, this%nlon) :: &
            vector_function
        real (WP),dimension (this%nlat, this%nlon)    :: &
            exact_solution, approximate_solution
        !--------------------------------------------------------------------------------

        nlat = this%nlat
        nlon = this%nlon

        ! initialize arrays
        vector_function = 0.0_WP
        exact_solution = 0.0_WP

        omega =  [ 1.0E+1_WP, 1.0E+2_WP, 1.0E+3_WP ]

        do l = 1, nlon
            do k = 1, nlat

                u = this%radial_unit_vector(:,k,l)

                sint = this%sint(k)
                cost = this%cost(k)

                sinp = this%sinp(l)
                cosp = this%cosp(l)

                scalar_function = exp(u%x + u%y + u%z)

                vector_function(:,k,l) = omega * scalar_function

                theta_derivative = &
                    (cost * cosp  + cost * sinp - sint ) &
                    * scalar_function

                phi_derivative = (u%x - u%y) * scalar_function

                cot = &
                    1.0_WP/ tan(this%grid%latitudes(k))

                rotation_operator = &
                    [ -sinp * theta_derivative - cosp * cot * phi_derivative, &
                    cosp * theta_derivative - sinp * cot * phi_derivative, &
                    phi_derivative ]

                exact_solution(k,l) = rotation_operator.dot.omega
            end do

        end do

        ! check if the work space arrays are set up for the (real) vector transform
        call this%Get_vorticity( vector_function, approximate_solution )

        ! set error
        error = maxval(abs(exact_solution - approximate_solution))

        ! print error to console
        print *, " "
        print *, "Test_get_vorticity"
        write (*,101)  error
101     format (1x, "l-infinity error = ", es23.16 )
        print *, " "
        print *, "******************************************** "

    end subroutine Test_get_vorticity
    !
    !o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
    !
    subroutine Test_get_rotation_operator( this )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)    :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                                   :: nlat, nlon, k,l
        real (WP)                                      :: sint, cost, sinp, cosp
        real (WP)                                      :: cot, theta, phi_derivative, theta_derivative
        real (WP)                                      :: error
        type (vector_t)                                :: u
        real (WP), dimension (3, this%nlat, this%nlon) :: exact_solution, approximate_solution
        real (WP),dimension (this%nlat, this%nlon)     :: scalar_function
        !--------------------------------------------------------------------------------

        ! Set constants
        nlat = this%nlat
        nlon = this%nlon

        ! initialize arrays
        scalar_function = 0.0_WP
        exact_solution = 0.0_WP

        do l = 1, nlon
            do k = 1, nlat

                u = this%radial_unit_vector(:,k,l)

                scalar_function(k,l) = exp(u%x + u%y + u%z)

                sint = this%sint(k)
                cost = this%cost(k)

                sinp = this%sinp(l)
                cosp = this%cosp(l)

                theta_derivative = (cost * cosp &
                    + cost*sinp - sint) * exp(u%x + u%y + u%z)

                phi_derivative = (u%x - u%y) * exp(u%x + u%y + u%z)

                theta = this%grid%latitudes(k)

                cot = 1.0_WP/ tan(theta)

                exact_solution(1, k, l) = &
                    -sinp * theta_derivative &
                    - cosp * cot * phi_derivative

                exact_solution(2,k,l) = &
                    cosp * theta_derivative &
                    - sinp * cot * phi_derivative

                exact_solution(3,k,l) = &
                    phi_derivative
            end do

        end do

        ! compute the rotation operator applied to the scalar function
        call this%Get_rotation_operator( scalar_function, approximate_solution)

        ! set error
        error = &
            maxval( &
            abs( exact_solution &
            -approximate_solution ) )

        ! print error to console
        print *, " "
        print *, "Test_get_rotation_operator"
        write (*,101)  error
101     format (1x, "l-infinity error = ", es23.16 )
        print *, " "
        print *, "******************************************** "

    end subroutine Test_get_rotation_operator
    !
    !o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
    !
end program test
