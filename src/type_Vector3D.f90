!
!  Purpose:
!
!  Defines a class for 3-dimensional cartesian vector calculations
!
module type_Vector3D

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: Vector3D
    public :: Vector3DPointer
    public :: assignment(=)
    public :: operator(*)


    ! Declare derived data type
    type, public :: Vector3D
        !----------------------------------------------------------------------
        ! Type components
        !----------------------------------------------------------------------
        real(wp), public :: x = 0.0_wp
        real(wp), public :: y = 0.0_wp
        real(wp), public :: z = 0.0_wp
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Type-bound procedures
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
        final                        :: finalize_vector3d
        !----------------------------------------------------------------------
        ! Generic type-bound procedures
        !----------------------------------------------------------------------
        generic, public :: operator (.dot.) => get_dot_product
        generic, public :: operator (.cross.) => get_cross_product
        generic, public :: operator (+) => add_vectors
        generic, public :: operator (-) => subtract_vectors
        generic, public :: operator (/) => &
            divide_vector_by_float,&
            divide_vector_by_integer
        !----------------------------------------------------------------------
    end type Vector3D

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

    ! Declare derived data type to create array of pointers
    type, public :: Vector3DPointer
        !----------------------------------------------------------------------
        ! Type components
        !----------------------------------------------------------------------
        type(Vector3D), pointer :: ptr => null()
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Type-bound procedures
        !----------------------------------------------------------------------
        final :: finalize_vector3d_pointer
        !----------------------------------------------------------------------
    end type Vector3DPointer

contains

    subroutine assign_vector_from_int_array(self, array)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Vector3D),  intent(out) :: self
        integer(ip),      intent(in)  :: array(:)
        !----------------------------------------------------------------------

        select type(self)
            type is (Vector3D)
            self%x = real(array(1), kind=wp)
            self%y = real(array(2), kind=wp)
            self%z = real(array(3), kind=wp)
        end select

    end subroutine assign_vector_from_int_array

    subroutine assign_vector_from_float_array(self, array)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Vector3D),  intent(out) :: self
        real(wp),         intent(in)  :: array(:)
        !----------------------------------------------------------------------

        select type(self)
            type is (Vector3D)
            self%x = array(1)
            self%y = array(2)
            self%z = array(3)
        end select

    end subroutine assign_vector_from_float_array

    subroutine assign_array_from_vector(array, self)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        real(wp),        intent(out)  :: array(:)
        class(Vector3D), intent(in)   :: self
        !----------------------------------------------------------------------

        select type(self)
            type is (Vector3D)
            array(1) = self%x
            array(2) = self%y
            array(3) = self%z
        end select

    end subroutine assign_array_from_vector

    subroutine copy_vector(self, vector_to_be_copied)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Vector3D), intent(out) :: self
        class(Vector3D), intent(in)  :: vector_to_be_copied
        !----------------------------------------------------------------------

        self%x = vector_to_be_copied%x
        self%y = vector_to_be_copied%y
        self%z = vector_to_be_copied%z

    end subroutine copy_vector

    function add_vectors(self, vec) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        type(Vector3D)               :: return_value
        class(Vector3D), intent(in) :: self
        class(Vector3D), intent(in) :: vec
        !----------------------------------------------------------------------

        return_value%x = self%x + vec%x
        return_value%y = self%y + vec%y
        return_value%z = self%z + vec%z

    end function add_vectors

    function subtract_vectors(self, vec) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Vector3D), intent(in) :: self
        class(Vector3D), intent(in) :: vec
        type(Vector3D)               :: return_value
        !----------------------------------------------------------------------

        return_value%x = self%x - vec%x
        return_value%y = self%y - vec%y
        return_value%z = self%z - vec%z

    end function subtract_vectors

    pure function multiply_vector_times_float(self, float) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Vector3D), intent(in) :: self
        real(wp),        intent(in) :: float
        type(Vector3D)               :: return_value
        !----------------------------------------------------------------------

        return_value%x = self%x * float
        return_value%y = self%y * float
        return_value%z = self%z * float

    end function multiply_vector_times_float

    function multiply_real_times_vector(float, vec) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        real(wp),        intent(in) :: float
        class(Vector3D), intent(in) :: vec
        type(Vector3D)               :: return_value
        !----------------------------------------------------------------------

        return_value%x = float * vec%x
        return_value%y = float * vec%y
        return_value%z = float * vec%z

    end function multiply_real_times_vector

    function multiply_vector_times_integer(self, int) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Vector3D), intent(in)  :: self
        integer(ip),     intent(in)  :: int
        type(Vector3D)                :: return_value
        !----------------------------------------------------------------------

        return_value%x = self%x * real(int, kind=wp)
        return_value%y = self%y * real(int, kind=wp)
        return_value%z = self%z * real(int, kind=wp)

    end function multiply_vector_times_integer

    pure function multiply_integer_times_vector(int, vec) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip),     intent(in) :: int
        class(Vector3D), intent(in) :: vec
        type(Vector3D)               :: return_value
        !----------------------------------------------------------------------

        return_value%x = real(int, kind=wp) * vec%x
        return_value%y = real(int, kind=wp) * vec%y
        return_value%z = real(int, kind=wp) * vec%z

    end function multiply_integer_times_vector

    pure function divide_vector_by_float(self, float) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Vector3D), intent(in)  :: self
        real(wp),                      intent(in)  :: float
        type(Vector3D)                :: return_value
        !----------------------------------------------------------------------

        return_value%x = self%x / float
        return_value%y = self%y / float
        return_value%z = self%z / float

    end function divide_vector_by_float

    pure function divide_vector_by_integer(self, int) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Vector3D), intent(in)  :: self
        integer(ip),                   intent(in)  :: int
        type(Vector3D)                :: return_value
        !----------------------------------------------------------------------

        return_value%x = self%x / int
        return_value%y = self%y / int
        return_value%z = self%z / int

    end function divide_vector_by_integer

    function get_dot_product(self, vec) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        real(wp)                    :: return_value
        class(Vector3D), intent(in) :: self
        class(Vector3D), intent(in) :: vec
        !----------------------------------------------------------------------

        associate( &
            a => [ self%x, self%y, self%z ], &
            b => [ vec%x, vec%y, vec%z ] &
            )

            return_value = dot_product(a, b)

        end associate

    end function get_dot_product

    function get_cross_product(self, vec) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Vector3D), intent(in) :: self
        class(Vector3D), intent(in) :: vec
        type(Vector3D)              :: return_value
        !----------------------------------------------------------------------

        return_value%x = self%y*vec%z - self%z*vec%y
        return_value%y = self%z*vec%x - self%x*vec%z
        return_value%z = self%x*vec%y - self%y*vec%x

    end function get_cross_product

    function get_norm(self) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Vector3D), intent(inout)  :: self
        real(wp)                        :: return_value
        !----------------------------------------------------------------------

        associate( array => [ self%x, self%y, self%z ] )
            return_value = norm2(array)
        end associate

    end function get_norm

    elemental subroutine finalize_vector3d(self)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        type(Vector3D), intent(inout)  :: self
        !----------------------------------------------------------------------

        ! Reset constants
        self%x = 0.0_wp
        self%y = 0.0_wp
        self%z = 0.0_wp

    end subroutine finalize_vector3d

    elemental subroutine finalize_vector3d_pointer(self)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        type(Vector3DPointer), intent(inout)  :: self
        !----------------------------------------------------------------------

        ! Check if pointer is associated
        if (associated(self%ptr)) nullify(self%ptr )

    end subroutine finalize_vector3d_pointer

end module type_Vector3D
