program test_symmetric_sparse_norms
    use stdlib_linalg_constants, only: dp
    use stdlib_sparse
    implicit none

    integer, parameter :: n = 3
    real(dp) :: A(n,n)
    type(COO_dp_type) :: COO_lower, COO_upper, COO_full
    real(dp) :: norm_fro_lower, norm_fro_upper, norm_fro_full
    real(dp) :: norm_one_lower, norm_one_upper, norm_one_full  
    real(dp) :: norm_inf_lower, norm_inf_upper, norm_inf_full
    real(dp) :: expected_fro, expected_one, expected_inf
    real(dp), parameter :: tol = 1.0e-12_dp
    
    ! Create a symmetric test matrix:
    ! A = [2, 1, 0]
    !     [1, 3, 4]
    !     [0, 4, 5]
    A = 0.0_dp
    A(1,1) = 2.0_dp; A(1,2) = 1.0_dp
    A(2,1) = 1.0_dp; A(2,2) = 3.0_dp; A(2,3) = 4.0_dp  
    A(3,2) = 4.0_dp; A(3,3) = 5.0_dp
    
    print *, "Original symmetric matrix:"
    print *, A(1,:)
    print *, A(2,:)
    print *, A(3,:)
    
    ! Expected norms for the full symmetric matrix:
    ! Frobenius: sqrt(2^2 + 1^2 + 1^2 + 3^2 + 4^2 + 4^2 + 5^2) = sqrt(72)
    ! 1-norm: max column sum = max(3, 8, 9) = 9
    ! Infinity norm: max row sum = max(3, 8, 9) = 9
    expected_fro = sqrt(72.0_dp)
    expected_one = 9.0_dp
    expected_inf = 9.0_dp
    
    ! Convert full matrix to COO format
    call dense2coo(A, COO_full)
    COO_full%storage = sparse_full
    
    ! Create lower triangular storage (store only lower triangle)
    call dense2coo(A, COO_lower)
    COO_lower%storage = sparse_lower
    ! Remove upper triangular elements
    call remove_upper_triangle(COO_lower)
    
    ! Create upper triangular storage (store only upper triangle)  
    call dense2coo(A, COO_upper)
    COO_upper%storage = sparse_upper
    ! Remove lower triangular elements
    call remove_lower_triangle(COO_upper)
    
    print *, ""
    print *, "Lower triangular storage has", COO_lower%nnz, "non-zeros"
    print *, "Upper triangular storage has", COO_upper%nnz, "non-zeros"
    print *, "Full storage has", COO_full%nnz, "non-zeros"
    
    ! Test norms for full storage
    norm_fro_full = sparse_norm(COO_full, 'frobenius')
    norm_one_full = sparse_norm(COO_full, 1)
    norm_inf_full = sparse_norm(COO_full, 'inf')
    
    ! Test norms for lower triangular storage
    norm_fro_lower = sparse_norm(COO_lower, 'frobenius')
    norm_one_lower = sparse_norm(COO_lower, 1)
    norm_inf_lower = sparse_norm(COO_lower, 'inf')
    
    ! Test norms for upper triangular storage
    norm_fro_upper = sparse_norm(COO_upper, 'frobenius')
    norm_one_upper = sparse_norm(COO_upper, 1)
    norm_inf_upper = sparse_norm(COO_upper, 'inf')
    
    print *, ""
    print *, "Full storage norms:"
    print *, "Frobenius norm:", norm_fro_full
    print *, "1-norm:        ", norm_one_full
    print *, "Infinity norm: ", norm_inf_full
    
    print *, ""
    print *, "Lower triangular storage norms:"
    print *, "Frobenius norm:", norm_fro_lower
    print *, "1-norm:        ", norm_one_lower
    print *, "Infinity norm: ", norm_inf_lower
    
    print *, ""
    print *, "Upper triangular storage norms:"
    print *, "Frobenius norm:", norm_fro_upper
    print *, "1-norm:        ", norm_one_upper
    print *, "Infinity norm: ", norm_inf_upper
    
    print *, ""
    print *, "Expected norms:"
    print *, "Frobenius norm:", expected_fro
    print *, "1-norm:        ", expected_one
    print *, "Infinity norm: ", expected_inf
    
    ! Check if all storage types give the same result
    if (abs(norm_fro_full - expected_fro) > tol .or. &
        abs(norm_fro_lower - expected_fro) > tol .or. &
        abs(norm_fro_upper - expected_fro) > tol) then
        print *, "ERROR: Frobenius norms don't match!"
        stop 1
    end if
    
    if (abs(norm_one_full - expected_one) > tol .or. &
        abs(norm_one_lower - expected_one) > tol .or. &
        abs(norm_one_upper - expected_one) > tol) then
        print *, "ERROR: 1-norms don't match!"
        stop 1
    end if
    
    if (abs(norm_inf_full - expected_inf) > tol .or. &
        abs(norm_inf_lower - expected_inf) > tol .or. &
        abs(norm_inf_upper - expected_inf) > tol) then
        print *, "ERROR: Infinity norms don't match!"
        stop 1
    end if
    
    print *, ""
    print *, "SUCCESS: All symmetric storage norms match!"

contains

    subroutine remove_upper_triangle(matrix)
        type(COO_dp_type), intent(inout) :: matrix
        integer :: i, j, new_nnz
        integer, allocatable :: new_index(:,:)
        real(dp), allocatable :: new_data(:)
        
        ! Count lower triangular elements (i >= j)
        new_nnz = 0
        do i = 1, matrix%nnz
            if (matrix%index(1,i) >= matrix%index(2,i)) then
                new_nnz = new_nnz + 1
            end if
        end do
        
        ! Allocate new arrays
        allocate(new_index(2, new_nnz))
        allocate(new_data(new_nnz))
        
        ! Copy lower triangular elements
        j = 0
        do i = 1, matrix%nnz
            if (matrix%index(1,i) >= matrix%index(2,i)) then
                j = j + 1
                new_index(:,j) = matrix%index(:,i)
                new_data(j) = matrix%data(i)
            end if
        end do
        
        ! Update matrix
        deallocate(matrix%index, matrix%data)
        matrix%index = new_index
        matrix%data = new_data
        matrix%nnz = new_nnz
    end subroutine

    subroutine remove_lower_triangle(matrix)
        type(COO_dp_type), intent(inout) :: matrix
        integer :: i, j, new_nnz
        integer, allocatable :: new_index(:,:)
        real(dp), allocatable :: new_data(:)
        
        ! Count upper triangular elements (i <= j)
        new_nnz = 0
        do i = 1, matrix%nnz
            if (matrix%index(1,i) <= matrix%index(2,i)) then
                new_nnz = new_nnz + 1
            end if
        end do
        
        ! Allocate new arrays
        allocate(new_index(2, new_nnz))
        allocate(new_data(new_nnz))
        
        ! Copy upper triangular elements
        j = 0
        do i = 1, matrix%nnz
            if (matrix%index(1,i) <= matrix%index(2,i)) then
                j = j + 1
                new_index(:,j) = matrix%index(:,i)
                new_data(j) = matrix%data(i)
            end if
        end do
        
        ! Update matrix
        deallocate(matrix%index, matrix%data)
        matrix%index = new_index
        matrix%data = new_data
        matrix%nnz = new_nnz
    end subroutine

end program test_symmetric_sparse_norms