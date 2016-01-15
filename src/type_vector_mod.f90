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
        WP => REAL64, &
        IP => INT32

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: vector_t
    public :: vector_ptr
    public :: assignment(=)
    public :: operator(*)

    ! Declare derived data type
    type, public :: vector_t

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
        procedure          :: Add_vectors
        procedure          :: Subtract_vectors
        procedure          :: Get_vector_divide_real
        procedure          :: Get_vector_divide_integer
        procedure          :: Get_dot_product
        procedure          :: Convert_array_to_vector
        procedure, nopass  :: Convert_vector_to_array
        procedure          :: Copy_vector_to_vector
        procedure          :: Get_vector_times_real
        procedure, nopass  :: Get_real_times_vector
        procedure          :: Get_vector_times_integer
        procedure, nopass  :: Get_integer_times_vector
        procedure          :: Get_cross_product
        !---------------------------------------------------------------------------------
        ! Public methods
        !---------------------------------------------------------------------------------
        procedure, public   :: Get_norm
        generic,   public   :: operator (.dot.) => Get_dot_product
        generic,   public   :: operator (+)     => Add_vectors
        generic,   public   :: operator (-)     => Subtract_vectors
        generic,   public   :: operator (/)     => Get_vector_divide_real, Get_vector_divide_integer
        !---------------------------------------------------------------------------------
        ! Finalizer
        !---------------------------------------------------------------------------------
        final              :: Finalize_vector
        !---------------------------------------------------------------------------------

    end type vector_t

    ! declare interface operators
    interface assignment (=)

        module procedure Convert_array_to_vector
        module procedure Convert_vector_to_array
        module procedure Copy_vector_to_vector

    end interface

    interface operator (*)

        module procedure Get_vector_times_real
        module procedure Get_real_times_vector
        module procedure Get_vector_times_integer
        module procedure Get_integer_times_vector
        module procedure Get_cross_product

    end interface

    ! To create array of pointers of TYPE (vector_t).
    type, public :: vector_ptr

        ! All components are public unless stated otherwise
        type (vector_t), pointer :: p => null()

    contains

        ! All methods are private unless stated otherwise
        private

        !---------------------------------------------------------------------------------
        ! Private methods
        !---------------------------------------------------------------------------------
        final             :: Finalize_vector_ptr
        !---------------------------------------------------------------------------------

    end type vector_ptr

contains
    !
    !*****************************************************************************************
    !
    subroutine Convert_array_to_vector( this, array )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (vector_t),   intent (out) :: this
        real (WP),          intent (in)  :: array(:)
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
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (WP),        intent (out) :: array(:)
        class (vector_t), intent (in)  :: this
        !--------------------------------------------------------------------------------

        array(1) = this%x
        array(2) = this%y
        array(3) = this%z

    end subroutine Convert_vector_to_array
    !
    !*****************************************************************************************
    !
    subroutine Copy_vector_to_vector( this, vector_to_be_copied )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (vector_t), intent (out) :: this
        class (vector_t), intent (in)  :: vector_to_be_copied
        !--------------------------------------------------------------------------------

        this%x = vector_to_be_copied%x
        this%y = vector_to_be_copied%y
        this%z = vector_to_be_copied%z

    end subroutine Copy_vector_to_vector
    !
    !*****************************************************************************************
    !
    function Add_vectors( vec_1, vec_2 ) result ( return_value )
        !
        !< Purpose:
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
        !< Purpose:
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
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)               :: return_value
        class (vector_t), intent (in) :: vec_1
        real (WP),        intent (in) :: real_2
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
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)               :: return_value
        real (WP),        intent (in) :: real_1
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
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)                :: return_value
        class (vector_t), intent (in)  :: vec_1
        integer (IP),     intent (in)  :: int_2
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
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)               :: return_value
        integer (IP),     intent (in) :: int_1
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
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)                :: return_value
        class (vector_t), intent (in)  :: vec_1
        real (WP),        intent (in)  :: real_2
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
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_t)                :: return_value
        class (vector_t), intent (in)  :: vec_1
        integer (IP),     intent (in)  :: int_2
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
        !< Purpose:
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
        !< Purpose:
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
    function Get_norm( this ) result ( return_value )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (WP)                         :: return_value
        class (vector_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        real (WP)                         :: array(3)
        !--------------------------------------------------------------------------------

        call Convert_vector_to_array( array, this )

        return_value = norm2( array )

    end function Get_norm
    !
    !*****************************************************************************************
    !
    elemental subroutine Finalize_vector( this )
        !
        !< Purpose:
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

    end subroutine Finalize_vector
    !
    !*****************************************************************************************
    !
    elemental subroutine Finalize_vector_ptr( this )
        !
        !< Purpose:
        !< Finalize object
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (vector_ptr), intent (in out) :: this
        !--------------------------------------------------------------------------------

        ! Check if pointer is associated
        if ( associated( this%p ) ) then

            ! Nullify pointer
            nullify( this%p )

        end if


    end subroutine Finalize_vector_ptr
    !
    !*****************************************************************************************
    !
end module type_vector_mod
