module spherepack_wrapper_library
    
    use type_ThreeDimensionalVector, only: &
        ThreeDimensionalVector, &
        ThreeDimensionalVectorPointer, &
        assignment(=), &
        operator(*)
    
    use type_GaussianGrid, only:&
        GaussianGrid
        
    use type_SpherepackWrapper, only: &
        SpherepackWrapper
    
    use type_SpherepackWorkspace, only: &
        SpherepackWorkspace

    ! Explicit typing only
    implicit none
    
    ! Everything is private unless stated otherwise
    private
    public :: ThreeDimensionalVector
    public :: ThreeDimensionalVectorPointer
    public :: assignment(=)
    public :: operator(*)
    public :: GaussianGrid
    public :: SpherepackWrapper
    public :: SpherepackWorkspace

end module spherepack_wrapper_library
