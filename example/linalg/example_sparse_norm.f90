program example_sparse_norm
    use stdlib_linalg_constants, only: dp
    use stdlib_sparse
    implicit none

    integer, parameter :: m = 4, n = 4
    real(dp) :: A(m,n)
    type(COO_dp_type) :: COO
    type(CSR_dp_type) :: CSR  
    real(dp) :: norm_frobenius, norm_one, norm_inf
    
    ! Create a simple test matrix
    ! A = [2, 0, 1, 0]
    !     [0, 3, 0, 0]  
    !     [0, 0, 0, 4]
    !     [1, 0, 0, 5]
    A = 0.0_dp
    A(1,1) = 2.0_dp; A(1,3) = 1.0_dp
    A(2,2) = 3.0_dp
    A(3,4) = 4.0_dp
    A(4,1) = 1.0_dp; A(4,4) = 5.0_dp
    
    print *, "Original dense matrix:"
    print *, A(1,:)
    print *, A(2,:)  
    print *, A(3,:)
    print *, A(4,:)
    
    ! Convert to COO format
    call dense2coo(A, COO)
    
    ! Convert COO to CSR
    call coo2csr(COO, CSR)
    
    ! Compute different norms using COO format
    norm_frobenius = sparse_norm(COO, 'frobenius')
    norm_one = sparse_norm(COO, 1)
    norm_inf = sparse_norm(COO, 'inf')
    
    print *, ""
    print *, "COO format norms:"
    print *, "Frobenius norm:", norm_frobenius  
    print *, "1-norm:        ", norm_one
    print *, "Infinity norm: ", norm_inf
    
    ! Compute norms using CSR format 
    norm_frobenius = sparse_norm(CSR, 'frobenius')
    norm_one = sparse_norm(CSR, 1)
    norm_inf = sparse_norm(CSR, 'inf')
    
    print *, ""
    print *, "CSR format norms:"
    print *, "Frobenius norm:", norm_frobenius
    print *, "1-norm:        ", norm_one
    print *, "Infinity norm: ", norm_inf
    
    ! Compare with expected values
    ! Frobenius norm should be sqrt(2^2 + 3^2 + 1^2 + 4^2 + 1^2 + 5^2) = sqrt(56) ≈ 7.48
    ! 1-norm should be max column sum = max(3, 3, 1, 9) = 9
    ! Infinity norm should be max row sum = max(3, 3, 4, 6) = 6
    
    print *, ""
    print *, "Expected values:"
    print *, "Frobenius norm: ~7.48 (sqrt(56))"
    print *, "1-norm:         9.0"
    print *, "Infinity norm:  6.0"
    
end program example_sparse_norm