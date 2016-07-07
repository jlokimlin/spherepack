!
!  Purpose:
!
!  Defines a class for 3-dimensional cartesian vector calculations
!
module type_ThreeDimensionalVector

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: ThreeDimensionalVector
    public :: ThreeDimensionalVectorPointer
    public :: assignment(=)
    public :: operator(*)

    ! Declare derived data type
    type, public :: ThreeDimensionalVector
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        real (wp), public :: x = 0.0_wp
        real (wp), public :: y = 0.0_wp
        real (wp), public :: z = 0.0_wp
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure,          private  :: add_vectors
        procedure,          private  :: subtract_vectors
        procedure,          private  :: divide_vector_by_float
        procedure,          private  :: divide_vector_by_integer
        procedure,          private  :: get_dot_product
        procedure,          private  :: assign_vector_from_int_array
        procedure,          private  :: assign_vector_from_float_array
        procedure, nopass,  private  :: assign_array_from_vector
        procedure,          private  :: copy_vector
        procedure,          private  :: multiply_vector_times_float
        procedure, nopass,  private  :: multiply_real_times_vector
        procedure,          private  :: multiply_vector_times_integer
        procedure, nopass,  private  :: multiply_integer_times_vector
        procedure,          private  :: get_cross_product
        procedure,          public   :: get_norm
        generic,            public   :: operator (.dot.) => get_dot_product
        generic,            public   :: operator (.cross.) => get_cross_product
        generic,            public   :: operator (+) => add_vectors
        generic,            public   :: operator (-) => subtract_vectors
        generic,            public   :: operator (/) => &
            divide_vector_by_float,&
            divide_vector_by_integer
        final                        :: finalize_three_dimensional_vector
        !----------------------------------------------------------------------
    end type ThreeDimensionalVector


    ! Interface for assignment operator
    interface assignment (=)
        module procedure assign_vector_from_int_array
        module procedure assign_vector_from_float_array
        module procedure assign_array_from_vector
        module procedure copy_vector
    end interface

    ! Interface for multiplication operator
    interface  operator (*)
        module procedure multiply_vector_times_float
        module procedure multiply_real_times_vector
        module procedure multiply_vector_times_integer
        module procedure multiply_integer_times_vector
        module procedure get_cross_product
    end interface


    ! To create array of pointers
    type, public :: ThreeDimensionalVectorPointer
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        type (ThreeDimensionalVector), pointer :: p => null()
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        final :: finalize_three_dimensional_vector_pointer
        !----------------------------------------------------------------------
    end type ThreeDimensionalVectorPointer



contains


    subroutine assign_vector_from_int_array(this, array)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (ThreeDimensionalVector),  intent (out) :: this
        integer (ip),                    intent (in)  :: array(:)
        !----------------------------------------------------------------------

        select type (this)
            type is (ThreeDimensionalVector)
            this%x = real(array(1), kind=wp)
            this%y = real(array(2), kind=wp)
            this%z = real(array(3), kind=wp)
        end select

    end subroutine assign_vector_from_int_array



    subroutine assign_vector_from_float_array(this, array)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (ThreeDimensionalVector),  intent (out) :: this
        real (wp),                       intent (in)  :: array(:)
        !----------------------------------------------------------------------

        select type (this)
            type is (ThreeDimensionalVector)
            this%x = array(1)
            this%y = array(2)
            this%z = array(3)
        end select

    end subroutine assign_vector_from_float_array


    subroutine assign_array_from_vector( array, this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real (wp),                      intent (out)  :: array(:)
        class (ThreeDimensionalVector), intent (in)   :: this
        !----------------------------------------------------------------------

        select type (this)
            type is (ThreeDimensionalVector)
            array(1) = this%x
            array(2) = this%y
            array(3) = this%z
        end select

    end subroutine assign_array_from_vector


    subroutine copy_vector(this, vector_to_be_copied )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (ThreeDimensionalVector), intent (out) :: this
        class (ThreeDimensionalVector), intent (in)  :: vector_to_be_copied
        !----------------------------------------------------------------------

        this%x = vector_to_be_copied%x
        this%y = vector_to_be_copied%y
        this%z = vector_to_be_copied%z

    end subroutine copy_vector


    function add_vectors(this, vec) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (ThreeDimensionalVector)               :: return_value
        class (ThreeDimensionalVector), intent (in) :: this
        class (ThreeDimensionalVector), intent (in) :: vec
        !----------------------------------------------------------------------

        return_value%x = this%x + vec%x
        return_value%y = this%y + vec%y
        return_value%z = this%z + vec%z

    end function add_vectors


    function subtract_vectors(this, vec) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (ThreeDimensionalVector), intent (in) :: this
        class (ThreeDimensionalVector), intent (in) :: vec
        type (ThreeDimensionalVector)               :: return_value
        !----------------------------------------------------------------------

        return_value%x = this%x - vec%x
        return_value%y = this%y - vec%y
        return_value%z = this%z - vec%z

    end function subtract_vectors


    pure function multiply_vector_times_float(this, float) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (ThreeDimensionalVector), intent (in) :: this
        real (wp),                      intent (in) :: float
        type (ThreeDimensionalVector)               :: return_value
        !----------------------------------------------------------------------

        return_value%x = this%x * float
        return_value%y = this%y * float
        return_value%z = this%z * float

    end function multiply_vector_times_float


    function multiply_real_times_vector(float, vec) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real (wp),                      intent (in) :: float
        class (ThreeDimensionalVector), intent (in) :: vec
        type (ThreeDimensionalVector)               :: return_value
        !----------------------------------------------------------------------

        return_value%x = float * vec%x
        return_value%y = float * vec%y
        return_value%z = float * vec%z

    end function multiply_real_times_vector


    function multiply_vector_times_integer(this, int) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (ThreeDimensionalVector), intent (in)  :: this
        integer (ip),                   intent (in)  :: int
        type (ThreeDimensionalVector)                :: return_value
        !----------------------------------------------------------------------

        return_value%x = this%x * real(int, kind=wp)
        return_value%y = this%y * real(int, kind=wp)
        return_value%z = this%z * real(int, kind=wp)

    end function multiply_vector_times_integer


    pure function multiply_integer_times_vector(int, vec) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (ThreeDimensionalVector)               :: return_value
        integer (ip),                   intent (in) :: int
        class (ThreeDimensionalVector), intent (in) :: vec
        !----------------------------------------------------------------------

        return_value%x = real(int, kind=wp) * vec%x
        return_value%y = real(int, kind=wp) * vec%y
        return_value%z = real(int, kind=wp) * vec%z

    end function multiply_integer_times_vector


    pure function divide_vector_by_float(this, float) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (ThreeDimensionalVector), intent (in)  :: this
        real (wp),                      intent (in)  :: float
        type (ThreeDimensionalVector)                :: return_value
        !----------------------------------------------------------------------

        return_value%x = this%x / float
        return_value%y = this%y / float
        return_value%z = this%z / float

    end function divide_vector_by_float



    pure function divide_vector_by_integer(this, int) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (ThreeDimensionalVector), intent (in)  :: this
        integer (ip),                   intent (in)  :: int
        type (ThreeDimensionalVector)                :: return_value
        !----------------------------------------------------------------------

        return_value%x = this%x / int
        return_value%y = this%y / int
        return_value%z = this%z / int

    end function divide_vector_by_integer



    function get_dot_product(this, vec) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real (wp)                                   :: return_value
        class (ThreeDimensionalVector), intent (in) :: this
        class (ThreeDimensionalVector), intent (in) :: vec
        !----------------------------------------------------------------------

        associate( &
            a => [ this%x, this%y, this%z ], &
            b => [ vec%x, vec%y, vec%z ] &
            )
            return_value = dot_product(a, b)
        end associate

    end function get_dot_product


    function get_cross_product(this, vec) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (ThreeDimensionalVector), intent (in) :: this
        class (ThreeDimensionalVector), intent (in) :: vec
        type (ThreeDimensionalVector)               :: return_value
        !----------------------------------------------------------------------

        return_value%x = this%y*vec%z - this%z*vec%y
        return_value%y = this%z*vec%x - this%x*vec%z
        return_value%z = this%x*vec%y - this%y*vec%x

    end function get_cross_product


    function get_norm(this) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (ThreeDimensionalVector), intent (in out) :: this
        real (wp)                                       :: return_value
        !----------------------------------------------------------------------

        associate( array => [ this%x, this%y, this%z ] )
            return_value = norm2(array)
        end associate

    end function get_norm


    elemental subroutine finalize_three_dimensional_vector(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (ThreeDimensionalVector), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Reset constants
        this%x = 0.0_wp
        this%y = 0.0_wp
        this%z = 0.0_wp

    end subroutine finalize_three_dimensional_vector


    elemental subroutine finalize_three_dimensional_vector_pointer(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (ThreeDimensionalVectorPointer), intent (in out) :: this
        !----------------------------------------------------------------------
        ! Check if pointer is associated
        if ( associated(this%p ) ) then
            nullify(this%p )
        end if

    end subroutine finalize_three_dimensional_vector_pointer


end module type_ThreeDimensionalVector
