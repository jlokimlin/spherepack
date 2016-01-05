module spherepack_wrapper_mod
    
    use type_vector_mod, only: &
        vector_t, &
        vector_ptr, &
        assignment(=), &
        operator(*)
    
    use type_grid_mod, only:&
        grid_t
        
    use type_sphere_mod, only: &
        sphere_t
    
    use type_workspace_mod, only: &
        workspace_t

    ! Explicit typing only
    implicit none
    
    ! Everything is private unless stated otherwise
    private
    public :: vector_t
    public :: vector_ptr
    public :: assignment(=)
    public :: operator(*)
    public :: grid_t
    public :: sphere_t
    public :: workspace_t

end module spherepack_wrapper_mod
