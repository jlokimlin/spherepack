!*****************************************************************************************
!
! Purpose:
!
! Defines the derived data type "workspace_t" required to
! invoke the spherepack library
!
module type_workspace_mod

    use, intrinsic :: iso_fortran_env, only: &
        REAL64, &
        INT32

    ! Explicit typing only
    implicit none

    ! Everything is private except the derived data type itself
    private

    ! Local variables confined to the module
    integer, parameter    :: WP     = REAL64                    !! 64 bit real
    integer, parameter    :: IP     = INT32                     !! 32 bit integer
    character(200)        :: error_message                      !! Probably long enough
    integer (IP)          :: allocate_status                    !! To check allocation status
    integer (IP)          :: deallocate_status                  !! To check deallocation status

    ! Declare derived data type
    type, public :: workspace_t

        ! Workspace arrays for Legendre transform
        real (WP), dimension (:), allocatable :: work
        real (WP), dimension (:), allocatable :: dwork

        ! Workspace arrays for scalar transform - GAU
        real (WP), dimension (:), allocatable :: wshags
        real (WP), dimension (:), allocatable :: wshsgs

        ! Workspace arrays for scalar transform - REG
        real (WP), dimension (:), allocatable :: wshaes
        real (WP), dimension (:), allocatable :: wshses

        ! Workspace arrays for vector transform - GAU
        real (WP), dimension (:), allocatable :: wvhags
        real (WP), dimension (:), allocatable :: wvhsgs

        ! Workspace arrays for vector transform - GAU
        real (WP), dimension (:), allocatable :: wvhaes
        real (WP), dimension (:), allocatable :: wvhses

        ! Scalar transform coefficients
        real (WP), dimension (:,:), pointer :: a => null()
        real (WP), dimension (:,:), pointer :: b => null()

        ! Vector transform coefficients
        real (WP), dimension (:,:), pointer :: br => null()
        real (WP), dimension (:,:), pointer :: bi => null()
        real (WP), dimension (:,:), pointer :: cr => null()
        real (WP), dimension (:,:), pointer :: ci => null()

        logical, private :: initialized = .false.

    contains

        ! SPHEREPACK 3.2 methods
        !        procedure :: Set_up_scalar_analysis
        !        procedure :: Set_up_scalar_synthesis
        !        procedure :: Set_up_vector_analysis
        !        procedure :: Set_up_vector_synthesis

        procedure, non_overridable :: Create
        procedure, non_overridable :: Destroy

    end type workspace_t

contains
    !
    !*****************************************************************************************
    !
    subroutine Create( me, nlat, nlon, grid_type )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: me
        integer (IP), intent (in)            :: nlat, nlon
        character (*), intent (in), optional :: grid_type
        !--------------------------------------------------------------------------------

        ! Check status
        if ( me%initialized ) then
            print *, 'ERROR: You must destroy "workspace" before re-instantiating'
            return
        end if

        ! Check if optional grid type argument is present
        if ( present( grid_type ) ) then

            ! Check if the grid type is equally spaced
            if ( grid_type .eq. 'REG' ) then

                ! Set up transforms for regular grids
                call Initialize_scalar_transform_reg( me, nlat, nlon )
                call Initialize_vector_transform_reg( me, nlat, nlon )

            ! Check if the grid type is gaussian
            else if ( grid_type .eq. 'GAU' ) then

                ! Set up transforms for gaussian grids
                call Initialize_scalar_transform_gau( me, nlat, nlon )
                call Initialize_vector_transform_gau( me, nlat, nlon )

            else

                ! Handle invalid grid type
                print *, 'ERROR: optional argument grid_type = ', grid_type
                print *, 'must be either REG or GAU (default GAU)'
                stop

            end if
        else

            ! Set up default transforms for gaussian grids
            call Initialize_scalar_transform_gau( me, nlat, nlon )
            call Initialize_vector_transform_gau( me, nlat, nlon )

        end if


        ! Set status
        me%initialized = .true.

    end subroutine Create
    !
    !*****************************************************************************************
    !
    subroutine Destroy( me )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: me
        !--------------------------------------------------------------------------------

        ! Check status
        if ( .not. me%initialized ) return

        ! Deallocate workspace arrays for Legendre transform
        if ( allocated( me%work ) ) then
            deallocate ( &
                me%work, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for "work":', &
                    trim( error_message )
            end if
        end if

        if ( allocated( me%dwork ) ) then
            deallocate ( &
                me%dwork, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for "dwork":', &
                    trim( error_message )
            end if
        end if

        ! Deallocate workspace arrays for scalar transform - GAU
        if ( allocated( me%wshags ) ) then
            deallocate ( &
                me%wshags, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for "wshags":', &
                    trim( error_message )
            end if
        end if

        if ( allocated( me%wshsgs ) ) then
            deallocate ( &
                me%wshsgs, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for "wshsgs":', &
                    trim( error_message )
            end if
        end if

        ! Deallocate workspace arrays for scalar transform - REG
        if ( allocated( me%wshaes ) ) then
            deallocate ( &
                me%wshaes, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for "wshaes":', &
                    trim( error_message )
            end if
        end if

        if ( allocated( me%wshses ) ) then
            deallocate ( &
                me%wshses, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for "wshses":', &
                    trim( error_message )
            end if
        end if

        ! Deallocate workspace arrays for vector transform - GAU
        if ( allocated( me%wvhags ) ) then
            deallocate ( &
                me%wvhags, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for "wvhags":',&
                    trim( error_message )
            end if
        end if

        if ( allocated( me%wvhsgs ) ) then
            deallocate ( &
                me%wvhsgs, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for "wvhsgs":', &
                    trim( error_message )
            end if
        end if

        ! Deallocate workspace arrays for vector transform - REG
        if ( allocated( me%wvhaes ) ) then
            deallocate ( &
                me%wvhaes, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for "wvhaes":',&
                    trim( error_message )
            end if
        end if

        if ( allocated( me%wvhses ) ) then
            deallocate ( &
                me%wvhses, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for "wvhses":', &
                    trim( error_message )
            end if
        end if

        ! Deallocate scalar transform pointers
        if ( associated(me%a) ) then

            ! Deallocate pointer
            deallocate( &
                me%a, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for a:', &
                    trim( error_message )
            end if

            ! Nullify pointer
            nullify( me%a )

        end if

        if ( associated(me%b) ) then

            ! Deallocate pointer
            deallocate( &
                me%b, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for b:', &
                    trim( error_message )
            end if

            ! Nullify pointer
            nullify( me%b )
        end if

        ! Deallocate vector transform pointers
        if ( associated(me%br) ) then

            ! Deallocate pointer
            deallocate( &
                me%br, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed:', trim( error_message )
            end if

            ! Nullify pointer
            nullify( me%br )

        end if

        if ( associated(me%bi) ) then

            ! Deallocate pointer
            deallocate( &
                me%bi, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed:', trim( error_message )
            end if

            ! Nullify pointer
            nullify( me%bi )

        end if

        if ( associated(me%cr) ) then

            ! Deallocate pointer
            deallocate( &
                me%cr, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for cr:', &
                    trim( error_message )
            end if

            ! Nullify pointer
            nullify( me%cr )

        end if

        if ( associated(me%ci) ) then

            ! Deallocate pointer
            deallocate( &
                me%ci, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed for ci:', &
                    trim( error_message )
            end if

            ! Nullify pointer
            nullify( me%ci )

        end if

        ! Reset status
        me%initialized = .false.

    end subroutine Destroy
    !
    !*****************************************************************************************
    !
    subroutine Initialize_scalar_transform_gau( me, nlat, nlon )
        !
        ! Purpose:
        ! Set the various workspace arrays and pointer to perform the
        ! (real) scalar harmonic transform.
        !
        ! Remark:
        ! This subroutine must be called first before calling
        ! Initialize_vector_transform_gau
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: me
        integer (IP), intent (in)            :: nlat, nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: ierror
        integer (IP) :: LWORK, LDWORK, LSHAGS, LSHSGS
        !--------------------------------------------------------------------------------

        ! Check status
        if ( me%initialized ) then
            print *, 'ERROR: You must destroy "workspace" before re-instantiating'
            return
        end if

        ! (STEP 1) set up scalar analysis

        ! Compute LWORK
        LWORK = Get_lwork( nlat, nlon )

        ! Compute LDWORK
        LDWORK = Get_ldwork( nlat )

        ! Compute LSHAGS
        LSHAGS = Get_lshags( nlat, nlon )

        ! Allocate arrays
        allocate ( &
            me%work( 1:LWORK ), &
            me%dwork( 1:LDWORK ), &
            me%wshags( LSHAGS ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed in '&
                &//'Initialize_scalar_transform_gau:', &
                trim( error_message )
            return
        end if

        ! Compute workspace arrays for analysis

        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#shags.html
        call Shagsi( &
            nlat, nlon, &
            me%wshags, LSHAGS, &
            me%work, LWORK, &
            me%dwork, LDWORK, ierror)

        ! Check the error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Shagsi'
            return
        end if

        ! (STEP 2) Set up scalar synthesis

        ! Compute LSHSGS
        LSHSGS = Get_lshsgs( nlat, nlon )

        ! Allocate array
        allocate ( &
            me%wshsgs( 1:LSHSGS ), &
            stat = allocate_status, &
            errmsg = error_message )

                ! Check allocation status
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed in '&
                &//'Initialize_scalar_transform_gau:', &
                trim( error_message )
            return
        end if

        ! Compute workspace arrays for synthesis
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#shsgs.html
        call Shsgsi( &
            nlat, nlon, &
            me%wshsgs, LSHSGS, &
            me%work, LWORK, &
            me%dwork, LDWORK, ierror)

        ! check error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Shsgsi'
            return
        end if

        ! (STEP 3) Allocate pointers for the (real) scalar transform
        allocate ( &
            me%a( 1:nlat, 1:nlat ), &
            me%b( 1:nlat, 1:nlat ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed in '&
                &//'Initialize_scalar_transform_gau:', &
                trim( error_message )
            return
        end if

    end subroutine Initialize_scalar_transform_gau
    !
    !*****************************************************************************************
    !
    subroutine Initialize_vector_transform_gau( me, nlat, nlon )
        !
        ! Purpose:
        ! Sets the various workspace arrays and pointers
        ! required for the (real) vector transform
        !
        ! Remark:
        ! The subroutine "Set_work_space_arrays_for_scalar_transform"
        ! must be called first to initialize the "dwork" component
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: me
        integer (IP), intent (in)            :: nlat, nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    :: ierror
        integer (IP)    :: LVHAGS, LDWORK, LVHSGS
        !--------------------------------------------------------------------------------

        ! Check status
        if ( me%initialized ) then
            print *, 'ERROR: You must destroy "workspace" before re-instantiating'
            return
        end if

        ! (Step 1) set up vector analysis

        ! Compute LVHAGS
        LVHAGS = Get_lvhags( nlat, nlon )

        ! Remark:
        ! me%dwork was already initialized in the previous call to
        ! in Set_work_space_arrays_for_scalar_transform
        LDWORK = Get_ldwork( nlat )

        ! Allocate array
        allocate ( &
            me%wvhags( 1:LVHAGS ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed in '&
                &//'Initialize_vector_transform_gau: ', &
                trim( error_message )
            return
        end if

        ! Compute workspace arrays for vector analysis

        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#vhags.html
        call Vhagsi( &
            nlat, nlon, &
            me%wvhags, LVHAGS, &
            me%dwork, LDWORK, ierror)

        ! Check error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vhagsi'
            return
        end if

        ! (Step 2) set up vector synthesis

        ! Compute LVHSGS
        LVHSGS = Get_lvhsgs( nlat, nlon )

        ! Allocate array
        allocate ( &
            me%wvhsgs( 1:LVHSGS ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed in '&
                &//'Initialize_vector_transform_gau:', &
                trim( error_message )
            return
        end if

        ! Compute workspace arrays for vector synthesis
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#vhsgs.html
        call Vhsgsi( &
            nlat, nlon, &
            me%wvhsgs, LVHSGS, &
            me%dwork, LDWORK, ierror )

        ! check the error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vhsgsi'
            return
        end if

        ! (Step 3) Allocate pointers for the vector transform coefficients
        allocate ( &
            me%br( 1:nlat, 1:nlat ), &
            me%bi( 1:nlat, 1:nlat ), &
            me%cr( 1:nlat, 1:nlat ), &
            me%ci( 1:nlat, 1:nlat ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed in '&
                &//'Initialize_vector_transform_gau:', &
                trim( error_message )
            return
        end if

    end subroutine Initialize_vector_transform_gau
    !
    !*****************************************************************************************
    !
    subroutine Initialize_scalar_transform_reg( me, nlat, nlon )
        !
        ! Purpose:
        ! Set the various workspace arrays and pointer to perform the
        ! (real) scalar harmonic transform.
        !
        ! Remark:
        ! This subroutine must be called first before calling
        ! Initialize_vector_transform_reg
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: me
        integer (IP), intent (in)            :: nlat, nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: ierror
        integer (IP) :: LWORK, LDWORK, LSHAES, LSHSES
        !--------------------------------------------------------------------------------

        ! Check status
        if ( me%initialized ) then
            print *, 'ERROR: You must destroy "workspace" before re-instantiating'
            return
        end if

        ! (STEP 1) set up scalar analysis

        ! Compute LWORK
        LWORK = Get_lwork( nlat, nlon )

        ! Compute LDWORK
        LDWORK = Get_ldwork( nlat )

        ! Compute LSHAGS
        LSHAES = Get_lshaes( nlat, nlon )

        ! Allocate arrays
        allocate ( &
            me%work( 1:LWORK ), &
            me%dwork( 1:LDWORK ), &
            me%wshaes( LSHAES ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed in '&
                &//'Initialize_scalar_transform_reg: ',&
                trim( error_message )
            return
        end if

        ! Compute workspace arrays for analysis

        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#shags.html
        call Shaesi( &
            nlat, nlon, &
            me%wshaes, LSHAES, &
            me%work, LWORK, &
            me%dwork, LDWORK, ierror)

        ! Check the error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Shaesi'
            return
        end if

        ! (STEP 2) Set up scalar synthesis

        ! Compute LSHSGS
        LSHSES = Get_lshses( nlat, nlon )

        ! Allocate array
        allocate ( &
            me%wshses( 1:LSHSES ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed in '&
                &//'Initialize_scalar_transform_reg: ',&
                trim( error_message )
            return
        end if

        ! Compute workspace arrays for synthesis
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#shsgs.html
        call Shsesi( &
            nlat, nlon, &
            me%wshses, LSHSES, &
            me%work, LWORK, &
            me%dwork, LDWORK, ierror)

        ! check error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Shsesi'
            return
        end if

        ! (STEP 3) Allocate pointers for the (real) scalar transform
        allocate ( &
            me%a( 1:nlat, 1:nlat ), &
            me%b( 1:nlat, 1:nlat ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed in '&
                &//'Initialize_scalar_transform_reg: ',&
                trim( error_message )
            return
        end if

    end subroutine Initialize_scalar_transform_reg
    !
    !*****************************************************************************************
    !
    subroutine Initialize_vector_transform_reg( me, nlat, nlon )
        !
        ! Purpose:
        ! Sets the various workspace arrays and pointers
        ! required for the (real) vector transform
        !
        ! Remark:
        ! Initialize_scalar_transform_reg
        ! must be called first to initialize the dwork component
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: me
        integer (IP), intent (in)            :: nlat, nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    :: ierror
        integer (IP)    :: LVHAES, LDWORK, LVHSES
        !--------------------------------------------------------------------------------

        ! Check status
        if ( me%initialized ) then
            print *, 'ERROR: You must destroy "workspace" before re-instantiating'
            return
        end if

        ! (Step 1) set up vector analysis

        ! Compute LVHAGS
        LVHAES = Get_lvhaes( nlat, nlon )

        ! Remark:
        ! me%dwork was already initialized in the previous call to
        ! in Set_work_space_arrays_for_scalar_transform
        LDWORK = Get_ldwork( nlat )

        ! Allocate array
        allocate ( &
            me%wvhaes( 1:LVHAES ), &
            stat = allocate_status, &
            errmsg = error_message )
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed:', trim( error_message )
            return
        end if

        ! Compute workspace arrays for vector analysis

        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#vhags.html
        call Vhaesi( &
            nlat, nlon, &
            me%wvhaes, LVHAES, &
            me%dwork, LDWORK, ierror)

        ! Check error status
        if ( ierror /= 0 ) then
            print *, 'Error ', ierror, ' in Vhaesi'
            return
        end if

        ! (Step 2) set up vector synthesis

        ! Compute LVHSGS
        LVHSES = Get_lvhses( nlat, nlon )

        ! Allocate array
        allocate ( &
            me%wvhses( 1:LVHSES ), &
            stat = allocate_status, &
            errmsg = error_message )
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed:', trim( error_message )
            return
        end if

        ! Compute workspace arrays for vector synthesis
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#vhsgs.html
        call Vhsesi( &
            nlat, nlon, &
            me%wvhses, LVHSES, &
            me%dwork, LDWORK, ierror )

        ! check the error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vhsesi'
            return
        end if

        ! (Step 3) Allocate pointers for the vector transform coefficients
        allocate ( &
            me%br( 1:nlat, 1:nlat ), &
            me%bi( 1:nlat, 1:nlat ), &
            me%cr( 1:nlat, 1:nlat ), &
            me%ci( 1:nlat, 1:nlat ), &
            stat = allocate_status, &
            errmsg = error_message )
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed in '&
                &//'Initialize_vector_transform_reg:', &
                trim( error_message )
            return
        end if

    end subroutine Initialize_vector_transform_reg
    !
    !*****************************************************************************************
    !
    pure function Get_lwork( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)              :: return_value
        integer (IP), intent (in) :: nlat
        integer (IP), intent (in) :: nlon
        !--------------------------------------------------------------------------------

        return_value = (4 * nlon + 2) * nlat

    end function Get_lwork
    !
    !*****************************************************************************************
    !
    pure function Get_ldwork( nlat ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)              :: return_value
        integer (IP), intent (in) :: nlat
        !--------------------------------------------------------------------------------

        return_value = (3 * nlat * (nlat + 3) + 2)/2

    end function Get_ldwork
    !
    !*****************************************************************************************
    !
    pure function Get_lshags( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)               :: return_value
        integer (IP), intent (in)  :: nlat
        integer (IP), intent (in)  :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)  ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            nlat * (3 * (l1 + l2) - 2) &
            +(l1 - 1) * (l2 * (2 * nlat - l1) - 3 * l1)/2 &
            + nlon + 15

    end function Get_lshags
    !
    !*****************************************************************************************
    !
    pure function Get_lshaes( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)               :: return_value
        integer (IP), intent (in)  :: nlat
        integer (IP), intent (in)  :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)  ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            (l1 * l2 * (nlat + nlat - l1 + 1))/2 + nlon + 15

    end function Get_lshaes
    !
    !*****************************************************************************************
    !
    pure function Get_lshsgs( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)              :: return_value
        integer (IP), intent (in) :: nlat
        integer (IP), intent (in) :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            nlat * (3 * (l1 + l2) - 2) &
            +(l1 - 1) * (l2 * (2 * nlat - l1) - 3 * l1)/2 &
            + nlon + 15

    end function Get_lshsgs
    !
    !*****************************************************************************************
    !
    pure function Get_lshses( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)              :: return_value
        integer (IP), intent (in) :: nlat
        integer (IP), intent (in) :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            (l1 * l2 * (nlat + nlat - l1 + 1))/2 + nlon + 15

    end function Get_lshses
    !
    !*****************************************************************************************
    !
    pure function Get_lvhags( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)              :: return_value
        integer (IP), intent (in) :: nlat
        integer (IP), intent (in) :: nlon
        !--------------------------------------------------------------------------------

        return_value = &
            (nlat + 1) * (nlat + 1) * nlat / 2 &
            + nlon + 15

    end function Get_lvhags
    !
    !*****************************************************************************************
    !
    pure function Get_lvhaes( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)               :: return_value
        integer (IP), intent (in)  :: nlat
        integer (IP), intent (in)  :: nlon
        !--------------------------------------------------------------------------------
        integer (IP) ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            l1 * l2 * (nlat + nlat - l1 + 1) + nlon + 15

    end function Get_lvhaes
    !
    !*****************************************************************************************
    !
    pure function Get_lvhsgs( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)              :: return_value
        integer (IP), intent (in) :: nlat
        integer (IP), intent (in) :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)  ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, nlon/2 )
        else
            l1 = min( nlat,(nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value =  &
            l1 * l2 * (nlat + nlat - l1 + 1) &
            + nlon + 15 + 2 * nlat

    end function Get_lvhsgs
    !
    !*****************************************************************************************
    !
    pure function Get_lvhses( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)               :: return_value
        integer (IP), intent (in)  :: nlat
        integer (IP), intent (in)  :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, nlon/2 )
        else
            l1 = min( nlat,(nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value =  &
            l1 * l2 * (nlat + nlat - l1 + 1) + nlon + 15

    end function Get_lvhses
    !
    !*****************************************************************************************
    !
end module type_workspace_mod
