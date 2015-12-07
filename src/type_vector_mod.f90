!*****************************************************************************************
!
!  Purpose:
!
!  Defines a derived data type for 3-dimensional cartesian vector calculations
!
!*****************************************************************************************
!
module type_vector_mod

    use, intrinsic :: iso_fortran_env, only: &
        REAL64, &
        INT32

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: vector_t
    public :: vector_ptr
    public :: assignment(=)
    public :: operator(*)

    !---------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !---------------------------------------------------------------------------------
    integer, parameter    :: WP   = REAL64  !! 64 bit real
    integer, parameter    :: IP   = INT32   !! 32 bit integer
    !---------------------------------------------------------------------------------

    ! Declare derived data type
    type :: vector_t

        ! All components are public unless stated otherwise
        !---------------------------------------------------------------------------------
        ! Real constants
        !---------------------------------------------------------------------------------
        real (WP)          :: x = 0.0_WP
        real (WP)          :: y = 0.0_WP
        real (WP)          :: z = 0.0_WP
        !---------------------------------------------------------------------------------

    contains

        ! All methods are private unless stated otherwise
        private

        !---------------------------------------------------------------------------------
        ! Private methods
        !---------------------------------------------------------------------------------
        procedure         :: Add_vectors
        procedure         :: Subtract_vectors
        procedure         :: Get_vector_divide_real
        procedure         :: Get_vector_divide_integer
        procedure         :: Get_dot_product
        final             :: Finalize
        !---------------------------------------------------------------------------------
        ! Public generic methods
        !---------------------------------------------------------------------------------
        generic, public   :: operator (.dot.) => Get_dot_product
        generic, public   :: operator (+)     => Add_vectors
        generic, public   :: operator (-)     => Subtract_vectors
        generic, public   :: operator (/)     => Get_vector_divide_real, Get_vector_divide_integer
        !---------------------------------------------------------------------------------

    end type vector_t

    ! declare interface operators
    interface assignment (=)

        module procedure Convert_array_to_vector
        module procedure Convert_vector_to_array

    end interface

    interface operator (*)

        module procedure Get_vector_times_real
        module procedure Get_real_times_vector
        module procedure Get_vector_times_integer
        module procedure Get_integer_times_vector
        module procedure Get_cross_product

    end interface

    ! Pointer of "vector_t" for creating array of pointers of "vector_t".
    type :: vector_ptr

        type(vector_t), pointer :: p => null()

    end type vector_ptr

contains
    !
    !*****************************************************************************************
    !
    subroutine Convert_array_to_vector( this, array )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (vector_t), intent (out)         :: this
        real (WP), dimension (3), intent (in)  :: array
        !--------------------------------------------------------------------------------

        this%x = array(1)
        this%y = array(2)
        this%z = array(3)

    end subroutine Convert_array_to_vector
    !
    !*****************************************************************************************
    !
    subroutine Convert_vector_to_array( array, this )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (WP), dimension (3), intent (out) :: array
        class (vector_t), intent (in)          :: this
        !--------------------------------------------------------------------------------

        array(1) = this%x
        array(2) = this%y
        array(3) = this%z

    end subroutine Convert_vector_to_array
    !
    !*****************************************************************************************
    !
    function Add_vectors( vec_1, vec_2 ) result ( return_value )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)               :: return_value
        class (vector_t), intent (in) :: vec_1
        class (vector_t), intent (in) :: vec_2
        !--------------------------------------------------------------------------------

        return_value%x = vec_1%x + vec_2%x
        return_value%y = vec_1%y + vec_2%y
        return_value%z = vec_1%z + vec_2%z

    end function Add_vectors
    !
    !*****************************************************************************************
    !
    function Subtract_vectors( vec_1, vec_2 ) result ( return_value )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)               :: return_value
        class (vector_t), intent (in) :: vec_1
        class (vector_t), intent (in) :: vec_2
        !--------------------------------------------------------------------------------

        return_value%x = vec_1%x - vec_2%x
        return_value%y = vec_1%y - vec_2%y
        return_value%z = vec_1%z - vec_2%z

    end function Subtract_vectors
    !
    !*****************************************************************************************
    !
    function Get_vector_times_real( vec_1, real_2 ) result ( return_value )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)               :: return_value
        class (vector_t), intent (in) :: vec_1
        real (WP), intent (in)         :: real_2
        !--------------------------------------------------------------------------------

        return_value%x = vec_1%x * real_2
        return_value%y = vec_1%y * real_2
        return_value%z = vec_1%z * real_2

    end function Get_vector_times_real
    !
    !*****************************************************************************************
    !
    function Get_real_times_vector(real_1, vec_2 ) result ( return_value )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)               :: return_value
        real (WP), intent (in)        :: real_1
        class (vector_t), intent (in) :: vec_2
        !--------------------------------------------------------------------------------

        return_value%x = real_1 * vec_2%x
        return_value%y = real_1 * vec_2%y
        return_value%z = real_1 * vec_2%z

    end function Get_real_times_vector
    !
    !*****************************************************************************************
    !
    function Get_vector_times_integer( vec_1, int_2 ) result ( return_value )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)                :: return_value
        class (vector_t), intent (in)  :: vec_1
        integer (IP), intent (in)      :: int_2
        !--------------------------------------------------------------------------------

        return_value%x = vec_1%x * real( int_2, WP)
        return_value%y = vec_1%y * real( int_2, WP)
        return_value%z = vec_1%z * real( int_2, WP)

    end function Get_vector_times_integer
    !
    !*****************************************************************************************
    !
    function Get_integer_times_vector( int_1, vec_2 ) result ( return_value )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)               :: return_value
        integer (IP), intent (in)     :: int_1
        class (vector_t), intent (in) :: vec_2
        !--------------------------------------------------------------------------------

        return_value%x = real( int_1, WP) * vec_2%x
        return_value%y = real( int_1, WP) * vec_2%y
        return_value%z = real( int_1, WP) * vec_2%z

    end function Get_integer_times_vector
    !
    !*****************************************************************************************
    !
    function Get_vector_divide_real( vec_1, real_2 ) result ( return_value )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)                :: return_value
        class (vector_t), intent (in)  :: vec_1
        real (WP), intent (in)         :: real_2
        !--------------------------------------------------------------------------------

        return_value%x = vec_1%x / real_2
        return_value%y = vec_1%y / real_2
        return_value%z = vec_1%z / real_2

    end function Get_vector_divide_real
    !
    !*****************************************************************************************
    !
    function Get_vector_divide_integer( vec_1, int_2 ) result ( return_value )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)                :: return_value
        class (vector_t), intent (in)  :: vec_1
        integer (IP), intent (in)      :: int_2
        !--------------------------------------------------------------------------------

        return_value%x = vec_1%x / int_2
        return_value%y = vec_1%y / int_2
        return_value%z = vec_1%z / int_2

    end function Get_vector_divide_integer
    !
    !*****************************************************************************************
    !
    function Get_dot_product( vec_1, vec_2 ) result ( return_value )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (WP)                     :: return_value
        class (vector_t), intent (in) :: vec_1
        class (vector_t), intent (in) :: vec_2
        !--------------------------------------------------------------------------------

        return_value = &
            vec_1%x*vec_2%x &
            + vec_1%y*vec_2%y &
            + vec_1%z*vec_2%z 

    end function Get_dot_product
    !
    !*****************************************************************************************
    !
    function Get_cross_product( vec_1, vec_2 ) result ( return_value )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)               :: return_value
        class (vector_t), intent (in) :: vec_1
        class (vector_t), intent (in) :: vec_2
        !--------------------------------------------------------------------------------

        return_value%x = vec_1%y*vec_2%z - vec_1%z*vec_2%y
        return_value%y = vec_1%z*vec_2%x - vec_1%x*vec_2%z
        return_value%z = vec_1%x*vec_2%y - vec_1%y*vec_2%x

    end function Get_cross_product
    !
    !*****************************************************************************************
    !
    subroutine Finalize( this )
        !
        ! Purpose:
        !< Finalize object
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t), intent (in out) :: this
        !--------------------------------------------------------------------------------

        ! Reset constants
        this%x = 0.0_WP
        this%y = 0.0_WP
        this%z = 0.0_WP

    end subroutine Finalize
    !
    !*****************************************************************************************
    !
end module type_vector_mod
