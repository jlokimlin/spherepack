module spherepack_wrapper_library
    
    use type_ThreeDimensionalVector, only: &
        ThreeDimensionalVector, &
        ThreeDimensionalVectorPointer, &
        assignment(=), &
        operator(*)
    
    use type_Sphere, only: &
        Sphere

    use type_GaussianGrid, only:&
        GaussianGrid
        
    use type_GaussianSphere, only: &
        GaussianSphere
    
    use type_RegularGrid, only: &
        RegularGrid

    use type_RegularSphere, only: &
        RegularSphere

    ! Explicit typing only
    implicit none
    
    ! Everything is private unless stated otherwise
    private
    public :: ThreeDimensionalVector
    public :: ThreeDimensionalVectorPointer
    public :: assignment(=)
    public :: operator(*)
    public :: Sphere
    public :: GaussianGrid
    public :: GaussianSphere
    public :: RegularGrid
    public :: RegularSphere

end module spherepack_wrapper_library
