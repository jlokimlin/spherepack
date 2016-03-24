module spherepack_wrapper_library
    
    use type_ThreeDimensionalVector, only: &
        ThreeDimensionalVector, &
        ThreeDimensionalVectorPointer, &
        assignment(=), &
        operator(*)
    
    use type_GaussianGrid, only:&
        GaussianGrid
        
    use type_GaussianSphere, only: &
        GaussianSphere
    
    ! Explicit typing only
    implicit none
    
    ! Everything is private unless stated otherwise
    private
    public :: ThreeDimensionalVector
    public :: ThreeDimensionalVectorPointer
    public :: assignment(=)
    public :: operator(*)
    public :: GaussianGrid
    public :: GaussianSphere

end module spherepack_wrapper_library
