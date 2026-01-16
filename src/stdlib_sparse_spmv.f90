!! The `stdlib_sparse_spmv` submodule provides matrix-vector product kernels.
!!
! This code was modified from https://github.com/jalvesz/FSPARSE by its author: Alves Jose
module stdlib_sparse_spmv
    use stdlib_sparse_constants
    use stdlib_sparse_kinds
    use stdlib_linalg_blas, only: gemv, gemm
    implicit none
    private

    !! Version experimental
    !!
    !! Apply the sparse matrix-vector product $$y = \alpha * op(M) * x + \beta * y $$
    !! [Specifications](../page/specs/stdlib_sparse.html#spmv)
    interface spmv
        module procedure :: spmv_coo_1d_sp
        module procedure :: spmv_bsr_1d_sp
        module procedure :: spmv_csr_1d_sp
        module procedure :: spmv_csc_1d_sp
        module procedure :: spmv_ell_1d_sp
        module procedure :: spmv_coo_2d_sp
        module procedure :: spmv_bsr_2d_sp
        module procedure :: spmv_csr_2d_sp
        module procedure :: spmv_csc_2d_sp
        module procedure :: spmv_ell_2d_sp
        module procedure :: spmv_sellc_sp
        module procedure :: spmv_coo_1d_dp
        module procedure :: spmv_bsr_1d_dp
        module procedure :: spmv_csr_1d_dp
        module procedure :: spmv_csc_1d_dp
        module procedure :: spmv_ell_1d_dp
        module procedure :: spmv_coo_2d_dp
        module procedure :: spmv_bsr_2d_dp
        module procedure :: spmv_csr_2d_dp
        module procedure :: spmv_csc_2d_dp
        module procedure :: spmv_ell_2d_dp
        module procedure :: spmv_sellc_dp
        module procedure :: spmv_coo_1d_csp
        module procedure :: spmv_bsr_1d_csp
        module procedure :: spmv_csr_1d_csp
        module procedure :: spmv_csc_1d_csp
        module procedure :: spmv_ell_1d_csp
        module procedure :: spmv_coo_2d_csp
        module procedure :: spmv_bsr_2d_csp
        module procedure :: spmv_csr_2d_csp
        module procedure :: spmv_csc_2d_csp
        module procedure :: spmv_ell_2d_csp
        module procedure :: spmv_sellc_csp
        module procedure :: spmv_coo_1d_cdp
        module procedure :: spmv_bsr_1d_cdp
        module procedure :: spmv_csr_1d_cdp
        module procedure :: spmv_csc_1d_cdp
        module procedure :: spmv_ell_1d_cdp
        module procedure :: spmv_coo_2d_cdp
        module procedure :: spmv_bsr_2d_cdp
        module procedure :: spmv_csr_2d_cdp
        module procedure :: spmv_csc_2d_cdp
        module procedure :: spmv_ell_2d_cdp
        module procedure :: spmv_sellc_cdp
    end interface
    public :: spmv

contains

    !! spmv_coo
    subroutine spmv_coo_1d_sp(matrix,vec_x,vec_y,alpha,beta,op)
        type(COO_sp_type), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: k, ik, jk

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_sp
        endif

        associate( data => matrix%data, index => matrix%index, storage => matrix%storage, nnz => matrix%nnz )
        select case(op_)
        case(sparse_op_none)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + alpha_*data(k) * vec_x(ik)
                end do

            end if
        case(sparse_op_transpose)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + alpha_*data(k) * vec_x(ik)
                end do

            end if
        end select
        end associate
    end subroutine

    subroutine spmv_coo_2d_sp(matrix,vec_x,vec_y,alpha,beta,op)
        type(COO_sp_type), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:,:)
        real(sp), intent(inout) :: vec_y(:,:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: k, ik, jk

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_sp
        endif

        associate( data => matrix%data, index => matrix%index, storage => matrix%storage, nnz => matrix%nnz )
        select case(op_)
        case(sparse_op_none)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + alpha_*data(k) * vec_x(:,ik)
                end do

            end if
        case(sparse_op_transpose)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + alpha_*data(k) * vec_x(:,ik)
                end do

            end if
        end select
        end associate
    end subroutine

    subroutine spmv_coo_1d_dp(matrix,vec_x,vec_y,alpha,beta,op)
        type(COO_dp_type), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: k, ik, jk

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_dp
        endif

        associate( data => matrix%data, index => matrix%index, storage => matrix%storage, nnz => matrix%nnz )
        select case(op_)
        case(sparse_op_none)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + alpha_*data(k) * vec_x(ik)
                end do

            end if
        case(sparse_op_transpose)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + alpha_*data(k) * vec_x(ik)
                end do

            end if
        end select
        end associate
    end subroutine

    subroutine spmv_coo_2d_dp(matrix,vec_x,vec_y,alpha,beta,op)
        type(COO_dp_type), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:,:)
        real(dp), intent(inout) :: vec_y(:,:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: k, ik, jk

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_dp
        endif

        associate( data => matrix%data, index => matrix%index, storage => matrix%storage, nnz => matrix%nnz )
        select case(op_)
        case(sparse_op_none)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + alpha_*data(k) * vec_x(:,ik)
                end do

            end if
        case(sparse_op_transpose)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + alpha_*data(k) * vec_x(:,ik)
                end do

            end if
        end select
        end associate
    end subroutine

    subroutine spmv_coo_1d_csp(matrix,vec_x,vec_y,alpha,beta,op)
        type(COO_csp_type), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: k, ik, jk

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_csp
        endif

        associate( data => matrix%data, index => matrix%index, storage => matrix%storage, nnz => matrix%nnz )
        select case(op_)
        case(sparse_op_none)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + alpha_*data(k) * vec_x(ik)
                end do

            end if
        case(sparse_op_transpose)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + alpha_*data(k) * vec_x(ik)
                end do

            end if
        case(sparse_op_hermitian)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*conjg(data(k)) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*conjg(data(k)) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + alpha_*conjg(data(k)) * vec_x(ik)
                end do

            end if
        end select
        end associate
    end subroutine

    subroutine spmv_coo_2d_csp(matrix,vec_x,vec_y,alpha,beta,op)
        type(COO_csp_type), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:,:)
        complex(sp), intent(inout) :: vec_y(:,:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: k, ik, jk

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_csp
        endif

        associate( data => matrix%data, index => matrix%index, storage => matrix%storage, nnz => matrix%nnz )
        select case(op_)
        case(sparse_op_none)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + alpha_*data(k) * vec_x(:,ik)
                end do

            end if
        case(sparse_op_transpose)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + alpha_*data(k) * vec_x(:,ik)
                end do

            end if
        case(sparse_op_hermitian)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*conjg(data(k)) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*conjg(data(k)) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + alpha_*conjg(data(k)) * vec_x(:,ik)
                end do

            end if
        end select
        end associate
    end subroutine

    subroutine spmv_coo_1d_cdp(matrix,vec_x,vec_y,alpha,beta,op)
        type(COO_cdp_type), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: k, ik, jk

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_cdp
        endif

        associate( data => matrix%data, index => matrix%index, storage => matrix%storage, nnz => matrix%nnz )
        select case(op_)
        case(sparse_op_none)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + alpha_*data(k) * vec_x(ik)
                end do

            end if
        case(sparse_op_transpose)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + alpha_*data(k) * vec_x(ik)
                end do

            end if
        case(sparse_op_hermitian)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*conjg(data(k)) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(ik) = vec_y(ik) + alpha_*conjg(data(k)) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + alpha_*conjg(data(k)) * vec_x(ik)
                end do

            end if
        end select
        end associate
    end subroutine

    subroutine spmv_coo_2d_cdp(matrix,vec_x,vec_y,alpha,beta,op)
        type(COO_cdp_type), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:,:)
        complex(dp), intent(inout) :: vec_y(:,:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: k, ik, jk

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_cdp
        endif

        associate( data => matrix%data, index => matrix%index, storage => matrix%storage, nnz => matrix%nnz )
        select case(op_)
        case(sparse_op_none)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + alpha_*data(k) * vec_x(:,ik)
                end do

            end if
        case(sparse_op_transpose)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + alpha_*data(k) * vec_x(:,ik)
                end do

            end if
        case(sparse_op_hermitian)
            if(storage == sparse_full) then
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*conjg(data(k)) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    jk = index(1,k)
                    ik = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + alpha_*conjg(data(k)) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + alpha_*conjg(data(k)) * vec_x(:,ik)
                end do

            end if
        end select
        end associate
    end subroutine


    !! spmv_csr
    subroutine spmv_csr_1d_sp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSR_sp_type), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        real(sp) :: aux, aux2
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_sp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
    
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do i = 1, nrows
                    aux = zero_sp
                    do j = rowptr(i), rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(col(j))
                    end do
                    vec_y(i) = vec_y(i) + alpha_ * aux
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do i = 1, nrows
                    aux = alpha_ * vec_x(i)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux
                    end do
                end do
                
            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = zero_sp
                    aux2 = alpha_ * vec_x(i)
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux2
                    end do
                    aux = alpha_ * aux + data(j) * aux2
                    vec_y(i) = vec_y(i) + aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = vec_x(i) * data(rowptr(i))
                    aux2 = alpha_ * vec_x(i)
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux2
                    end do
                    vec_y(i) = vec_y(i) + alpha_ * aux
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine spmv_csr_2d_sp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSR_sp_type), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:,:)
        real(sp), intent(inout) :: vec_y(:,:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        real(sp) :: aux(size(vec_x,dim=1)), aux2(size(vec_x,dim=1))
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_sp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
    
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do i = 1, nrows
                    aux = zero_sp
                    do j = rowptr(i), rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(:,col(j))
                    end do
                    vec_y(:,i) = vec_y(:,i) + alpha_ * aux
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do i = 1, nrows
                    aux = alpha_ * vec_x(:,i)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * aux
                    end do
                end do
                
            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = zero_sp
                    aux2 = alpha_ * vec_x(:,i)
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * aux2
                    end do
                    aux = alpha_ * aux + data(j) * aux2
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = vec_x(:,i) * data(rowptr(i))
                    aux2 = alpha_ * vec_x(:,i)
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * aux2
                    end do
                    vec_y(:,i) = vec_y(:,i) + alpha_ * aux
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine spmv_csr_1d_dp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSR_dp_type), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        real(dp) :: aux, aux2
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_dp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
    
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do i = 1, nrows
                    aux = zero_dp
                    do j = rowptr(i), rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(col(j))
                    end do
                    vec_y(i) = vec_y(i) + alpha_ * aux
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do i = 1, nrows
                    aux = alpha_ * vec_x(i)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux
                    end do
                end do
                
            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = zero_dp
                    aux2 = alpha_ * vec_x(i)
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux2
                    end do
                    aux = alpha_ * aux + data(j) * aux2
                    vec_y(i) = vec_y(i) + aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = vec_x(i) * data(rowptr(i))
                    aux2 = alpha_ * vec_x(i)
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux2
                    end do
                    vec_y(i) = vec_y(i) + alpha_ * aux
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine spmv_csr_2d_dp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSR_dp_type), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:,:)
        real(dp), intent(inout) :: vec_y(:,:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        real(dp) :: aux(size(vec_x,dim=1)), aux2(size(vec_x,dim=1))
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_dp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
    
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do i = 1, nrows
                    aux = zero_dp
                    do j = rowptr(i), rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(:,col(j))
                    end do
                    vec_y(:,i) = vec_y(:,i) + alpha_ * aux
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do i = 1, nrows
                    aux = alpha_ * vec_x(:,i)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * aux
                    end do
                end do
                
            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = zero_dp
                    aux2 = alpha_ * vec_x(:,i)
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * aux2
                    end do
                    aux = alpha_ * aux + data(j) * aux2
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = vec_x(:,i) * data(rowptr(i))
                    aux2 = alpha_ * vec_x(:,i)
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * aux2
                    end do
                    vec_y(:,i) = vec_y(:,i) + alpha_ * aux
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine spmv_csr_1d_csp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSR_csp_type), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        complex(sp) :: aux, aux2
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_csp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
    
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do i = 1, nrows
                    aux = zero_sp
                    do j = rowptr(i), rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(col(j))
                    end do
                    vec_y(i) = vec_y(i) + alpha_ * aux
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do i = 1, nrows
                    aux = alpha_ * vec_x(i)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux
                    end do
                end do
                
            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = zero_csp
                    aux2 = alpha_ * vec_x(i)
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux2
                    end do
                    aux = alpha_ * aux + data(j) * aux2
                    vec_y(i) = vec_y(i) + aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = vec_x(i) * data(rowptr(i))
                    aux2 = alpha_ * vec_x(i)
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux2
                    end do
                    vec_y(i) = vec_y(i) + alpha_ * aux
                end do
                
            else if( storage == sparse_full .and. op_==sparse_op_hermitian) then
                do i = 1, nrows
                    aux = alpha_ * vec_x(i)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(col(j)) = vec_y(col(j)) + conjg(data(j)) * aux
                    end do
                end do

            else if( storage == sparse_lower .and. op_==sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = zero_csp
                    aux2 = alpha_ * vec_x(i)
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + conjg(data(j)) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + conjg(data(j)) * aux2
                    end do
                    aux = alpha_ * aux + conjg(data(j)) * aux2
                    vec_y(i) = vec_y(i) + aux
                end do

            else if( storage == sparse_upper .and. op_==sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = vec_x(i) * conjg(data(rowptr(i)))
                    aux2 = alpha_ * vec_x(i)
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + conjg(data(j)) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + conjg(data(j)) * aux2
                    end do
                    vec_y(i) = vec_y(i) + alpha_ * aux
                end do
            end if
        end associate
    end subroutine
    
    subroutine spmv_csr_2d_csp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSR_csp_type), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:,:)
        complex(sp), intent(inout) :: vec_y(:,:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        complex(sp) :: aux(size(vec_x,dim=1)), aux2(size(vec_x,dim=1))
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_csp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
    
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do i = 1, nrows
                    aux = zero_sp
                    do j = rowptr(i), rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(:,col(j))
                    end do
                    vec_y(:,i) = vec_y(:,i) + alpha_ * aux
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do i = 1, nrows
                    aux = alpha_ * vec_x(:,i)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * aux
                    end do
                end do
                
            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = zero_csp
                    aux2 = alpha_ * vec_x(:,i)
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * aux2
                    end do
                    aux = alpha_ * aux + data(j) * aux2
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = vec_x(:,i) * data(rowptr(i))
                    aux2 = alpha_ * vec_x(:,i)
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * aux2
                    end do
                    vec_y(:,i) = vec_y(:,i) + alpha_ * aux
                end do
                
            else if( storage == sparse_full .and. op_==sparse_op_hermitian) then
                do i = 1, nrows
                    aux = alpha_ * vec_x(:,i)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(:,col(j)) = vec_y(:,col(j)) + conjg(data(j)) * aux
                    end do
                end do

            else if( storage == sparse_lower .and. op_==sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = zero_csp
                    aux2 = alpha_ * vec_x(:,i)
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + conjg(data(j)) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + conjg(data(j)) * aux2
                    end do
                    aux = alpha_ * aux + conjg(data(j)) * aux2
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            else if( storage == sparse_upper .and. op_==sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = vec_x(:,i) * conjg(data(rowptr(i)))
                    aux2 = alpha_ * vec_x(:,i)
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + conjg(data(j)) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + conjg(data(j)) * aux2
                    end do
                    vec_y(:,i) = vec_y(:,i) + alpha_ * aux
                end do
            end if
        end associate
    end subroutine
    
    subroutine spmv_csr_1d_cdp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSR_cdp_type), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        complex(dp) :: aux, aux2
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_cdp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
    
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do i = 1, nrows
                    aux = zero_dp
                    do j = rowptr(i), rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(col(j))
                    end do
                    vec_y(i) = vec_y(i) + alpha_ * aux
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do i = 1, nrows
                    aux = alpha_ * vec_x(i)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux
                    end do
                end do
                
            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = zero_cdp
                    aux2 = alpha_ * vec_x(i)
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux2
                    end do
                    aux = alpha_ * aux + data(j) * aux2
                    vec_y(i) = vec_y(i) + aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = vec_x(i) * data(rowptr(i))
                    aux2 = alpha_ * vec_x(i)
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux2
                    end do
                    vec_y(i) = vec_y(i) + alpha_ * aux
                end do
                
            else if( storage == sparse_full .and. op_==sparse_op_hermitian) then
                do i = 1, nrows
                    aux = alpha_ * vec_x(i)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(col(j)) = vec_y(col(j)) + conjg(data(j)) * aux
                    end do
                end do

            else if( storage == sparse_lower .and. op_==sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = zero_cdp
                    aux2 = alpha_ * vec_x(i)
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + conjg(data(j)) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + conjg(data(j)) * aux2
                    end do
                    aux = alpha_ * aux + conjg(data(j)) * aux2
                    vec_y(i) = vec_y(i) + aux
                end do

            else if( storage == sparse_upper .and. op_==sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = vec_x(i) * conjg(data(rowptr(i)))
                    aux2 = alpha_ * vec_x(i)
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + conjg(data(j)) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + conjg(data(j)) * aux2
                    end do
                    vec_y(i) = vec_y(i) + alpha_ * aux
                end do
            end if
        end associate
    end subroutine
    
    subroutine spmv_csr_2d_cdp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSR_cdp_type), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:,:)
        complex(dp), intent(inout) :: vec_y(:,:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        complex(dp) :: aux(size(vec_x,dim=1)), aux2(size(vec_x,dim=1))
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_cdp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
    
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do i = 1, nrows
                    aux = zero_dp
                    do j = rowptr(i), rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(:,col(j))
                    end do
                    vec_y(:,i) = vec_y(:,i) + alpha_ * aux
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do i = 1, nrows
                    aux = alpha_ * vec_x(:,i)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * aux
                    end do
                end do
                
            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = zero_cdp
                    aux2 = alpha_ * vec_x(:,i)
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * aux2
                    end do
                    aux = alpha_ * aux + data(j) * aux2
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = vec_x(:,i) * data(rowptr(i))
                    aux2 = alpha_ * vec_x(:,i)
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * aux2
                    end do
                    vec_y(:,i) = vec_y(:,i) + alpha_ * aux
                end do
                
            else if( storage == sparse_full .and. op_==sparse_op_hermitian) then
                do i = 1, nrows
                    aux = alpha_ * vec_x(:,i)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(:,col(j)) = vec_y(:,col(j)) + conjg(data(j)) * aux
                    end do
                end do

            else if( storage == sparse_lower .and. op_==sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = zero_cdp
                    aux2 = alpha_ * vec_x(:,i)
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + conjg(data(j)) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + conjg(data(j)) * aux2
                    end do
                    aux = alpha_ * aux + conjg(data(j)) * aux2
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            else if( storage == sparse_upper .and. op_==sparse_op_hermitian )then
                do i = 1 , nrows
                    aux  = vec_x(:,i) * conjg(data(rowptr(i)))
                    aux2 = alpha_ * vec_x(:,i)
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + conjg(data(j)) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + conjg(data(j)) * aux2
                    end do
                    vec_y(:,i) = vec_y(:,i) + alpha_ * aux
                end do
            end if
        end associate
    end subroutine
    

    !! spmv_bsr
    subroutine spmv_bsr_1d_sp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_sp_type), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(sp) :: alpha_, beta_
        character(1) :: op_
        character(1) :: trans_blas
        integer(ilp) :: nbrow, bi, bj, p
        integer(ilp) :: row_i_start, col_j_start, i_max_i, j_max_j
        integer(ilp) :: row_j_start, col_i_start, i_max_j, j_max_i


        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp; if(present(alpha)) alpha_ = alpha
        beta_ = zero_sp; if(present(beta)) beta_ = beta
        if(present(beta)) then
            vec_y = beta_ * vec_y
        else
            vec_y = zero_sp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage, &
            br => matrix%block_shape(1), bc => matrix%block_shape(2) )

            nbrow = size(rowptr) - 1

            select case(op_)
            case(sparse_op_none)
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemv('N', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_j_start:), incx=1, beta=one_sp, &
                            y=vec_y(row_i_start:), incy=1)

                        if (storage /= sparse_full .and. bi /= bj) then
                            row_j_start = (bj - 1) * br + 1
                            if (row_j_start > nrows) cycle
                            i_max_j = min(br, nrows - row_j_start + 1)
                            col_i_start = (bi - 1) * bc + 1
                            if (col_i_start > ncols) cycle
                            j_max_i = min(bc, ncols - col_i_start + 1)

                            call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                x=vec_x(col_i_start:), incx=1, beta=one_sp, &
                                y=vec_y(row_j_start:), incy=1)
                        end if
                    end do
                end do

            case(sparse_op_transpose)
                ! storage == sparse_full here (otherwise op is mapped to 'N')
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(row_i_start:), incx=1, beta=one_sp, &
                            y=vec_y(col_j_start:), incy=1)
                    end do
                end do

            end select
        end associate
    end subroutine

    subroutine spmv_bsr_2d_sp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_sp_type), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:,:)
        real(sp), intent(inout) :: vec_y(:,:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(sp) :: alpha_, beta_
        character(1) :: op_
        character(1) :: trans_blas
        integer(ilp) :: nbrow, bi, bj, p
        integer(ilp) :: row_i_start, col_j_start, i_max_i, j_max_j
        integer(ilp) :: row_j_start, col_i_start, i_max_j, j_max_i

        integer(ilp) :: nrhs

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp; if(present(alpha)) alpha_ = alpha
        beta_ = zero_sp; if(present(beta)) beta_ = beta
        if(present(beta)) then
            vec_y = beta_ * vec_y
        else
            vec_y = zero_sp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage, &
            br => matrix%block_shape(1), bc => matrix%block_shape(2) )

            nbrow = size(rowptr) - 1
            nrhs = size(vec_x,dim=2)

            select case(op_)
            case(sparse_op_none)
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemm('N','N', m=i_max_i, n=nrhs, k=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_j_start:, :), ldb=j_max_j, beta=one_sp, &
                            c=vec_y(row_i_start:, :), ldc=i_max_i)

                        if (storage /= sparse_full .and. bi /= bj) then
                            row_j_start = (bj - 1) * br + 1
                            if (row_j_start > nrows) cycle
                            i_max_j = min(br, nrows - row_j_start + 1)
                            col_i_start = (bi - 1) * bc + 1
                            if (col_i_start > ncols) cycle
                            j_max_i = min(bc, ncols - col_i_start + 1)

                            call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_sp, &
                                c=vec_y(row_j_start:, :), ldc=i_max_j)
                        end if
                    end do
                end do

            case(sparse_op_transpose)
                ! storage == sparse_full here (otherwise op is mapped to 'N')
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(row_i_start:, :), ldb=i_max_i, beta=one_sp, &
                            c=vec_y(col_j_start:, :), ldc=j_max_j)
                    end do
                end do

            end select
        end associate
    end subroutine

    subroutine spmv_bsr_1d_dp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_dp_type), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(dp) :: alpha_, beta_
        character(1) :: op_
        character(1) :: trans_blas
        integer(ilp) :: nbrow, bi, bj, p
        integer(ilp) :: row_i_start, col_j_start, i_max_i, j_max_j
        integer(ilp) :: row_j_start, col_i_start, i_max_j, j_max_i


        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp; if(present(alpha)) alpha_ = alpha
        beta_ = zero_dp; if(present(beta)) beta_ = beta
        if(present(beta)) then
            vec_y = beta_ * vec_y
        else
            vec_y = zero_dp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage, &
            br => matrix%block_shape(1), bc => matrix%block_shape(2) )

            nbrow = size(rowptr) - 1

            select case(op_)
            case(sparse_op_none)
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemv('N', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_j_start:), incx=1, beta=one_dp, &
                            y=vec_y(row_i_start:), incy=1)

                        if (storage /= sparse_full .and. bi /= bj) then
                            row_j_start = (bj - 1) * br + 1
                            if (row_j_start > nrows) cycle
                            i_max_j = min(br, nrows - row_j_start + 1)
                            col_i_start = (bi - 1) * bc + 1
                            if (col_i_start > ncols) cycle
                            j_max_i = min(bc, ncols - col_i_start + 1)

                            call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                x=vec_x(col_i_start:), incx=1, beta=one_dp, &
                                y=vec_y(row_j_start:), incy=1)
                        end if
                    end do
                end do

            case(sparse_op_transpose)
                ! storage == sparse_full here (otherwise op is mapped to 'N')
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(row_i_start:), incx=1, beta=one_dp, &
                            y=vec_y(col_j_start:), incy=1)
                    end do
                end do

            end select
        end associate
    end subroutine

    subroutine spmv_bsr_2d_dp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_dp_type), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:,:)
        real(dp), intent(inout) :: vec_y(:,:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(dp) :: alpha_, beta_
        character(1) :: op_
        character(1) :: trans_blas
        integer(ilp) :: nbrow, bi, bj, p
        integer(ilp) :: row_i_start, col_j_start, i_max_i, j_max_j
        integer(ilp) :: row_j_start, col_i_start, i_max_j, j_max_i

        integer(ilp) :: nrhs

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp; if(present(alpha)) alpha_ = alpha
        beta_ = zero_dp; if(present(beta)) beta_ = beta
        if(present(beta)) then
            vec_y = beta_ * vec_y
        else
            vec_y = zero_dp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage, &
            br => matrix%block_shape(1), bc => matrix%block_shape(2) )

            nbrow = size(rowptr) - 1
            nrhs = size(vec_x,dim=2)

            select case(op_)
            case(sparse_op_none)
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemm('N','N', m=i_max_i, n=nrhs, k=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_j_start:, :), ldb=j_max_j, beta=one_dp, &
                            c=vec_y(row_i_start:, :), ldc=i_max_i)

                        if (storage /= sparse_full .and. bi /= bj) then
                            row_j_start = (bj - 1) * br + 1
                            if (row_j_start > nrows) cycle
                            i_max_j = min(br, nrows - row_j_start + 1)
                            col_i_start = (bi - 1) * bc + 1
                            if (col_i_start > ncols) cycle
                            j_max_i = min(bc, ncols - col_i_start + 1)

                            call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_dp, &
                                c=vec_y(row_j_start:, :), ldc=i_max_j)
                        end if
                    end do
                end do

            case(sparse_op_transpose)
                ! storage == sparse_full here (otherwise op is mapped to 'N')
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(row_i_start:, :), ldb=i_max_i, beta=one_dp, &
                            c=vec_y(col_j_start:, :), ldc=j_max_j)
                    end do
                end do

            end select
        end associate
    end subroutine

    subroutine spmv_bsr_1d_csp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_csp_type), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(sp) :: alpha_, beta_
        character(1) :: op_
        character(1) :: trans_blas
        integer(ilp) :: nbrow, bi, bj, p
        integer(ilp) :: row_i_start, col_j_start, i_max_i, j_max_j
        integer(ilp) :: row_j_start, col_i_start, i_max_j, j_max_i

        complex(sp), allocatable :: xwork_bc(:)
        complex(sp), allocatable :: ywork_br(:)

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_csp; if(present(alpha)) alpha_ = alpha
        beta_ = zero_csp; if(present(beta)) beta_ = beta
        if(present(beta)) then
            vec_y = beta_ * vec_y
        else
            vec_y = zero_csp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage, &
            br => matrix%block_shape(1), bc => matrix%block_shape(2) )

            nbrow = size(rowptr) - 1

            select case(op_)
            case(sparse_op_none)
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemv('N', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_j_start:), incx=1, beta=one_csp, &
                            y=vec_y(row_i_start:), incy=1)

                        if (storage /= sparse_full .and. bi /= bj) then
                            row_j_start = (bj - 1) * br + 1
                            if (row_j_start > nrows) cycle
                            i_max_j = min(br, nrows - row_j_start + 1)
                            col_i_start = (bi - 1) * bc + 1
                            if (col_i_start > ncols) cycle
                            j_max_i = min(bc, ncols - col_i_start + 1)

                            call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                x=vec_x(col_i_start:), incx=1, beta=one_csp, &
                                y=vec_y(row_j_start:), incy=1)
                        end if
                    end do
                end do

            case(sparse_op_transpose)
                ! storage == sparse_full here (otherwise op is mapped to 'N')
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(row_i_start:), incx=1, beta=one_csp, &
                            y=vec_y(col_j_start:), incy=1)
                    end do
                end do

            case(sparse_op_hermitian)
                if (storage == sparse_full) then
                    trans_blas = 'C'
                    do bi = 1, nbrow
                        row_i_start = (bi - 1) * br + 1
                        if (row_i_start > nrows) cycle
                        i_max_i = min(br, nrows - row_i_start + 1)
                        do p = rowptr(bi), rowptr(bi+1)-1
                            bj = col(p)
                            col_j_start = (bj - 1) * bc + 1
                            if (col_j_start > ncols) cycle
                            j_max_j = min(bc, ncols - col_j_start + 1)

                            call gemv(trans_blas, m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                x=vec_x(row_i_start:), incx=1, beta=one_csp, &
                                y=vec_y(col_j_start:), incy=1)
                        end do
                    end do
                else
                    ! symmetric-storage: full matrix is symmetric; compute y = alpha * A^H * x
                    if(.not.allocated(xwork_bc)) allocate(xwork_bc(bc))
                    if(.not.allocated(ywork_br)) allocate(ywork_br(br))
                    do bi = 1, nbrow
                        row_i_start = (bi - 1) * br + 1
                        if (row_i_start > nrows) cycle
                        i_max_i = min(br, nrows - row_i_start + 1)
                        do p = rowptr(bi), rowptr(bi+1)-1
                            bj = col(p)
                            col_j_start = (bj - 1) * bc + 1
                            if (col_j_start > ncols) cycle
                            j_max_j = min(bc, ncols - col_j_start + 1)

                            ! y_i += conjg(A_ij) * x_j  == conjg( A_ij * conjg(x_j) )
                            xwork_bc(1:j_max_j) = conjg(vec_x(col_j_start:col_j_start+j_max_j-1))
                            ywork_br(1:i_max_i) = zero_csp
                            call gemv('N', m=i_max_i, n=j_max_j, alpha=one_csp, a=data(1:,1:,p), lda=br, &
                                x=xwork_bc, incx=1, beta=zero_csp, y=ywork_br, incy=1)
                            vec_y(row_i_start:row_i_start+i_max_i-1) = vec_y(row_i_start:row_i_start+i_max_i-1) + &
                                alpha_ * conjg(ywork_br(1:i_max_i))

                            if (bi /= bj) then
                                row_j_start = (bj - 1) * br + 1
                                if (row_j_start > nrows) cycle
                                i_max_j = min(br, nrows - row_j_start + 1)
                                col_i_start = (bi - 1) * bc + 1
                                if (col_i_start > ncols) cycle
                                j_max_i = min(bc, ncols - col_i_start + 1)

                                ! y_j += A_ij^H * x_i
                                call gemv('C', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                    x=vec_x(col_i_start:), incx=1, beta=one_csp, &
                                    y=vec_y(row_j_start:), incy=1)
                            end if
                        end do
                    end do
                end if
            end select
        end associate
    end subroutine

    subroutine spmv_bsr_2d_csp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_csp_type), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:,:)
        complex(sp), intent(inout) :: vec_y(:,:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(sp) :: alpha_, beta_
        character(1) :: op_
        character(1) :: trans_blas
        integer(ilp) :: nbrow, bi, bj, p
        integer(ilp) :: row_i_start, col_j_start, i_max_i, j_max_j
        integer(ilp) :: row_j_start, col_i_start, i_max_j, j_max_i

        integer(ilp) :: nrhs
        complex(sp), allocatable :: xwork_bc(:,:)
        complex(sp), allocatable :: ywork_br(:,:)

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_csp; if(present(alpha)) alpha_ = alpha
        beta_ = zero_csp; if(present(beta)) beta_ = beta
        if(present(beta)) then
            vec_y = beta_ * vec_y
        else
            vec_y = zero_csp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage, &
            br => matrix%block_shape(1), bc => matrix%block_shape(2) )

            nbrow = size(rowptr) - 1
            nrhs = size(vec_x,dim=2)

            select case(op_)
            case(sparse_op_none)
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemm('N','N', m=i_max_i, n=nrhs, k=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_j_start:, :), ldb=j_max_j, beta=one_csp, &
                            c=vec_y(row_i_start:, :), ldc=i_max_i)

                        if (storage /= sparse_full .and. bi /= bj) then
                            row_j_start = (bj - 1) * br + 1
                            if (row_j_start > nrows) cycle
                            i_max_j = min(br, nrows - row_j_start + 1)
                            col_i_start = (bi - 1) * bc + 1
                            if (col_i_start > ncols) cycle
                            j_max_i = min(bc, ncols - col_i_start + 1)

                            call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_csp, &
                                c=vec_y(row_j_start:, :), ldc=i_max_j)
                        end if
                    end do
                end do

            case(sparse_op_transpose)
                ! storage == sparse_full here (otherwise op is mapped to 'N')
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(row_i_start:, :), ldb=i_max_i, beta=one_csp, &
                            c=vec_y(col_j_start:, :), ldc=j_max_j)
                    end do
                end do

            case(sparse_op_hermitian)
                if (storage == sparse_full) then
                    trans_blas = 'C'
                    do bi = 1, nbrow
                        row_i_start = (bi - 1) * br + 1
                        if (row_i_start > nrows) cycle
                        i_max_i = min(br, nrows - row_i_start + 1)
                        do p = rowptr(bi), rowptr(bi+1)-1
                            bj = col(p)
                            col_j_start = (bj - 1) * bc + 1
                            if (col_j_start > ncols) cycle
                            j_max_j = min(bc, ncols - col_j_start + 1)

                            call gemm('C','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                b=vec_x(row_i_start:, :), ldb=i_max_i, beta=one_csp, &
                                c=vec_y(col_j_start:, :), ldc=j_max_j)
                        end do
                    end do
                else
                    ! symmetric-storage: full matrix is symmetric; compute y = alpha * A^H * x
                    if(.not.allocated(xwork_bc)) allocate(xwork_bc(bc,nrhs))
                    if(.not.allocated(ywork_br)) allocate(ywork_br(br,nrhs))
                    do bi = 1, nbrow
                        row_i_start = (bi - 1) * br + 1
                        if (row_i_start > nrows) cycle
                        i_max_i = min(br, nrows - row_i_start + 1)
                        do p = rowptr(bi), rowptr(bi+1)-1
                            bj = col(p)
                            col_j_start = (bj - 1) * bc + 1
                            if (col_j_start > ncols) cycle
                            j_max_j = min(bc, ncols - col_j_start + 1)

                            ! y_i += conjg(A_ij) * x_j  == conjg( A_ij * conjg(x_j) )
                            xwork_bc(1:j_max_j, :) = conjg(vec_x(col_j_start:col_j_start+j_max_j-1, :))
                            ywork_br(1:i_max_i, :) = zero_csp
                            call gemm('N','N', m=i_max_i, n=nrhs, k=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                b=xwork_bc, ldb=bc, beta=zero_csp, c=ywork_br, ldc=br)
                            vec_y(row_i_start:row_i_start+i_max_i-1, :) = vec_y(row_i_start:row_i_start+i_max_i-1, :) + &
                                alpha_ * conjg(ywork_br(1:i_max_i, :))

                            if (bi /= bj) then
                                row_j_start = (bj - 1) * br + 1
                                if (row_j_start > nrows) cycle
                                i_max_j = min(br, nrows - row_j_start + 1)
                                col_i_start = (bi - 1) * bc + 1
                                if (col_i_start > ncols) cycle
                                j_max_i = min(bc, ncols - col_i_start + 1)

                                ! y_j += A_ij^H * x_i
                                call gemm('C','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                    b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_csp, &
                                    c=vec_y(row_j_start:, :), ldc=i_max_j)
                            end if
                        end do
                    end do
                end if
            end select
        end associate
    end subroutine

    subroutine spmv_bsr_1d_cdp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_cdp_type), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(dp) :: alpha_, beta_
        character(1) :: op_
        character(1) :: trans_blas
        integer(ilp) :: nbrow, bi, bj, p
        integer(ilp) :: row_i_start, col_j_start, i_max_i, j_max_j
        integer(ilp) :: row_j_start, col_i_start, i_max_j, j_max_i

        complex(dp), allocatable :: xwork_bc(:)
        complex(dp), allocatable :: ywork_br(:)

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_cdp; if(present(alpha)) alpha_ = alpha
        beta_ = zero_cdp; if(present(beta)) beta_ = beta
        if(present(beta)) then
            vec_y = beta_ * vec_y
        else
            vec_y = zero_cdp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage, &
            br => matrix%block_shape(1), bc => matrix%block_shape(2) )

            nbrow = size(rowptr) - 1

            select case(op_)
            case(sparse_op_none)
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemv('N', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_j_start:), incx=1, beta=one_cdp, &
                            y=vec_y(row_i_start:), incy=1)

                        if (storage /= sparse_full .and. bi /= bj) then
                            row_j_start = (bj - 1) * br + 1
                            if (row_j_start > nrows) cycle
                            i_max_j = min(br, nrows - row_j_start + 1)
                            col_i_start = (bi - 1) * bc + 1
                            if (col_i_start > ncols) cycle
                            j_max_i = min(bc, ncols - col_i_start + 1)

                            call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                x=vec_x(col_i_start:), incx=1, beta=one_cdp, &
                                y=vec_y(row_j_start:), incy=1)
                        end if
                    end do
                end do

            case(sparse_op_transpose)
                ! storage == sparse_full here (otherwise op is mapped to 'N')
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(row_i_start:), incx=1, beta=one_cdp, &
                            y=vec_y(col_j_start:), incy=1)
                    end do
                end do

            case(sparse_op_hermitian)
                if (storage == sparse_full) then
                    trans_blas = 'C'
                    do bi = 1, nbrow
                        row_i_start = (bi - 1) * br + 1
                        if (row_i_start > nrows) cycle
                        i_max_i = min(br, nrows - row_i_start + 1)
                        do p = rowptr(bi), rowptr(bi+1)-1
                            bj = col(p)
                            col_j_start = (bj - 1) * bc + 1
                            if (col_j_start > ncols) cycle
                            j_max_j = min(bc, ncols - col_j_start + 1)

                            call gemv(trans_blas, m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                x=vec_x(row_i_start:), incx=1, beta=one_cdp, &
                                y=vec_y(col_j_start:), incy=1)
                        end do
                    end do
                else
                    ! symmetric-storage: full matrix is symmetric; compute y = alpha * A^H * x
                    if(.not.allocated(xwork_bc)) allocate(xwork_bc(bc))
                    if(.not.allocated(ywork_br)) allocate(ywork_br(br))
                    do bi = 1, nbrow
                        row_i_start = (bi - 1) * br + 1
                        if (row_i_start > nrows) cycle
                        i_max_i = min(br, nrows - row_i_start + 1)
                        do p = rowptr(bi), rowptr(bi+1)-1
                            bj = col(p)
                            col_j_start = (bj - 1) * bc + 1
                            if (col_j_start > ncols) cycle
                            j_max_j = min(bc, ncols - col_j_start + 1)

                            ! y_i += conjg(A_ij) * x_j  == conjg( A_ij * conjg(x_j) )
                            xwork_bc(1:j_max_j) = conjg(vec_x(col_j_start:col_j_start+j_max_j-1))
                            ywork_br(1:i_max_i) = zero_cdp
                            call gemv('N', m=i_max_i, n=j_max_j, alpha=one_cdp, a=data(1:,1:,p), lda=br, &
                                x=xwork_bc, incx=1, beta=zero_cdp, y=ywork_br, incy=1)
                            vec_y(row_i_start:row_i_start+i_max_i-1) = vec_y(row_i_start:row_i_start+i_max_i-1) + &
                                alpha_ * conjg(ywork_br(1:i_max_i))

                            if (bi /= bj) then
                                row_j_start = (bj - 1) * br + 1
                                if (row_j_start > nrows) cycle
                                i_max_j = min(br, nrows - row_j_start + 1)
                                col_i_start = (bi - 1) * bc + 1
                                if (col_i_start > ncols) cycle
                                j_max_i = min(bc, ncols - col_i_start + 1)

                                ! y_j += A_ij^H * x_i
                                call gemv('C', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                    x=vec_x(col_i_start:), incx=1, beta=one_cdp, &
                                    y=vec_y(row_j_start:), incy=1)
                            end if
                        end do
                    end do
                end if
            end select
        end associate
    end subroutine

    subroutine spmv_bsr_2d_cdp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_cdp_type), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:,:)
        complex(dp), intent(inout) :: vec_y(:,:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(dp) :: alpha_, beta_
        character(1) :: op_
        character(1) :: trans_blas
        integer(ilp) :: nbrow, bi, bj, p
        integer(ilp) :: row_i_start, col_j_start, i_max_i, j_max_j
        integer(ilp) :: row_j_start, col_i_start, i_max_j, j_max_i

        integer(ilp) :: nrhs
        complex(dp), allocatable :: xwork_bc(:,:)
        complex(dp), allocatable :: ywork_br(:,:)

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_cdp; if(present(alpha)) alpha_ = alpha
        beta_ = zero_cdp; if(present(beta)) beta_ = beta
        if(present(beta)) then
            vec_y = beta_ * vec_y
        else
            vec_y = zero_cdp
        endif

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage, &
            br => matrix%block_shape(1), bc => matrix%block_shape(2) )

            nbrow = size(rowptr) - 1
            nrhs = size(vec_x,dim=2)

            select case(op_)
            case(sparse_op_none)
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemm('N','N', m=i_max_i, n=nrhs, k=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_j_start:, :), ldb=j_max_j, beta=one_cdp, &
                            c=vec_y(row_i_start:, :), ldc=i_max_i)

                        if (storage /= sparse_full .and. bi /= bj) then
                            row_j_start = (bj - 1) * br + 1
                            if (row_j_start > nrows) cycle
                            i_max_j = min(br, nrows - row_j_start + 1)
                            col_i_start = (bi - 1) * bc + 1
                            if (col_i_start > ncols) cycle
                            j_max_i = min(bc, ncols - col_i_start + 1)

                            call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_cdp, &
                                c=vec_y(row_j_start:, :), ldc=i_max_j)
                        end if
                    end do
                end do

            case(sparse_op_transpose)
                ! storage == sparse_full here (otherwise op is mapped to 'N')
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(row_i_start:, :), ldb=i_max_i, beta=one_cdp, &
                            c=vec_y(col_j_start:, :), ldc=j_max_j)
                    end do
                end do

            case(sparse_op_hermitian)
                if (storage == sparse_full) then
                    trans_blas = 'C'
                    do bi = 1, nbrow
                        row_i_start = (bi - 1) * br + 1
                        if (row_i_start > nrows) cycle
                        i_max_i = min(br, nrows - row_i_start + 1)
                        do p = rowptr(bi), rowptr(bi+1)-1
                            bj = col(p)
                            col_j_start = (bj - 1) * bc + 1
                            if (col_j_start > ncols) cycle
                            j_max_j = min(bc, ncols - col_j_start + 1)

                            call gemm('C','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                b=vec_x(row_i_start:, :), ldb=i_max_i, beta=one_cdp, &
                                c=vec_y(col_j_start:, :), ldc=j_max_j)
                        end do
                    end do
                else
                    ! symmetric-storage: full matrix is symmetric; compute y = alpha * A^H * x
                    if(.not.allocated(xwork_bc)) allocate(xwork_bc(bc,nrhs))
                    if(.not.allocated(ywork_br)) allocate(ywork_br(br,nrhs))
                    do bi = 1, nbrow
                        row_i_start = (bi - 1) * br + 1
                        if (row_i_start > nrows) cycle
                        i_max_i = min(br, nrows - row_i_start + 1)
                        do p = rowptr(bi), rowptr(bi+1)-1
                            bj = col(p)
                            col_j_start = (bj - 1) * bc + 1
                            if (col_j_start > ncols) cycle
                            j_max_j = min(bc, ncols - col_j_start + 1)

                            ! y_i += conjg(A_ij) * x_j  == conjg( A_ij * conjg(x_j) )
                            xwork_bc(1:j_max_j, :) = conjg(vec_x(col_j_start:col_j_start+j_max_j-1, :))
                            ywork_br(1:i_max_i, :) = zero_cdp
                            call gemm('N','N', m=i_max_i, n=nrhs, k=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                b=xwork_bc, ldb=bc, beta=zero_cdp, c=ywork_br, ldc=br)
                            vec_y(row_i_start:row_i_start+i_max_i-1, :) = vec_y(row_i_start:row_i_start+i_max_i-1, :) + &
                                alpha_ * conjg(ywork_br(1:i_max_i, :))

                            if (bi /= bj) then
                                row_j_start = (bj - 1) * br + 1
                                if (row_j_start > nrows) cycle
                                i_max_j = min(br, nrows - row_j_start + 1)
                                col_i_start = (bi - 1) * bc + 1
                                if (col_i_start > ncols) cycle
                                j_max_i = min(bc, ncols - col_i_start + 1)

                                ! y_j += A_ij^H * x_i
                                call gemm('C','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                                    b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_cdp, &
                                    c=vec_y(row_j_start:, :), ldc=i_max_j)
                            end if
                        end do
                    end do
                end if
            end select
        end associate
    end subroutine


    !! spmv_csc
    subroutine spmv_csc_1d_sp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSC_sp_type), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        real(sp) :: aux

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_sp
        endif

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent(j=1:ncols)
                    aux = alpha_ * vec_x(j)
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(row(i)) = vec_y(row(i)) + data(i) * aux
                    end do
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent(j=1:ncols)
                    aux = zero_sp
                    do i = colptr(j), colptr(j+1)-1
                        aux = aux + data(i) * vec_x(row(i))
                    end do
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = vec_x(j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(row(i))
                        vec_y(row(i)) = vec_y(row(i)) + alpha_ * data(i) * vec_x(j)
                    end do
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = zero_sp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(j)
                        vec_y(i) = vec_y(i) + alpha_ * data(i) * vec_x(row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(j)
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            end if
        end associate
    end subroutine
    
    subroutine spmv_csc_2d_sp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSC_sp_type), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:,:)
        real(sp), intent(inout) :: vec_y(:,:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        real(sp) :: aux(size(vec_x,dim=1))

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_sp
        endif

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent(j=1:ncols)
                    aux = alpha_ * vec_x(:,j)
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(:,row(i)) = vec_y(:,row(i)) + data(i) * aux
                    end do
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent(j=1:ncols)
                    aux = zero_sp
                    do i = colptr(j), colptr(j+1)-1
                        aux = aux + data(i) * vec_x(:,row(i))
                    end do
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = vec_x(:,j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(:,row(i))
                        vec_y(:,row(i)) = vec_y(:,row(i)) + alpha_ * data(i) * vec_x(:,j)
                    end do
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = zero_sp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(:,j)
                        vec_y(:,i) = vec_y(:,i) + alpha_ * data(i) * vec_x(:,row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(:,j)
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            end if
        end associate
    end subroutine
    
    subroutine spmv_csc_1d_dp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSC_dp_type), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        real(dp) :: aux

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_dp
        endif

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent(j=1:ncols)
                    aux = alpha_ * vec_x(j)
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(row(i)) = vec_y(row(i)) + data(i) * aux
                    end do
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent(j=1:ncols)
                    aux = zero_dp
                    do i = colptr(j), colptr(j+1)-1
                        aux = aux + data(i) * vec_x(row(i))
                    end do
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = vec_x(j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(row(i))
                        vec_y(row(i)) = vec_y(row(i)) + alpha_ * data(i) * vec_x(j)
                    end do
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = zero_dp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(j)
                        vec_y(i) = vec_y(i) + alpha_ * data(i) * vec_x(row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(j)
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            end if
        end associate
    end subroutine
    
    subroutine spmv_csc_2d_dp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSC_dp_type), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:,:)
        real(dp), intent(inout) :: vec_y(:,:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        real(dp) :: aux(size(vec_x,dim=1))

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_dp
        endif

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent(j=1:ncols)
                    aux = alpha_ * vec_x(:,j)
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(:,row(i)) = vec_y(:,row(i)) + data(i) * aux
                    end do
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent(j=1:ncols)
                    aux = zero_dp
                    do i = colptr(j), colptr(j+1)-1
                        aux = aux + data(i) * vec_x(:,row(i))
                    end do
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = vec_x(:,j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(:,row(i))
                        vec_y(:,row(i)) = vec_y(:,row(i)) + alpha_ * data(i) * vec_x(:,j)
                    end do
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = zero_dp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(:,j)
                        vec_y(:,i) = vec_y(:,i) + alpha_ * data(i) * vec_x(:,row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(:,j)
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            end if
        end associate
    end subroutine
    
    subroutine spmv_csc_1d_csp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSC_csp_type), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        complex(sp) :: aux

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_csp
        endif

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent(j=1:ncols)
                    aux = alpha_ * vec_x(j)
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(row(i)) = vec_y(row(i)) + data(i) * aux
                    end do
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent(j=1:ncols)
                    aux = zero_sp
                    do i = colptr(j), colptr(j+1)-1
                        aux = aux + data(i) * vec_x(row(i))
                    end do
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = vec_x(j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(row(i))
                        vec_y(row(i)) = vec_y(row(i)) + alpha_ * data(i) * vec_x(j)
                    end do
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = zero_csp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(j)
                        vec_y(i) = vec_y(i) + alpha_ * data(i) * vec_x(row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(j)
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_full .and. op_==sparse_op_hermitian ) then
                do concurrent(j=1:ncols)
                    aux = zero_sp
                    do i = colptr(j), colptr(j+1)-1
                        aux = aux + conjg(data(i)) * vec_x(row(i))
                    end do
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_lower .and. op_==sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = vec_x(j) * conjg(data(colptr(j)))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + conjg(data(i)) * vec_x(row(i))
                        vec_y(row(i)) = vec_y(row(i)) + alpha_ * conjg(data(i)) * vec_x(j)
                    end do
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_upper .and. op_==sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = zero_csp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + conjg(data(i)) * vec_x(j)
                        vec_y(i) = vec_y(i) + alpha_ * conjg(data(i)) * vec_x(row(i))
                    end do
                    aux = aux + conjg(data(colptr(j))) * vec_x(j)
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do
            end if
        end associate
    end subroutine
    
    subroutine spmv_csc_2d_csp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSC_csp_type), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:,:)
        complex(sp), intent(inout) :: vec_y(:,:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        complex(sp) :: aux(size(vec_x,dim=1))

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_csp
        endif

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent(j=1:ncols)
                    aux = alpha_ * vec_x(:,j)
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(:,row(i)) = vec_y(:,row(i)) + data(i) * aux
                    end do
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent(j=1:ncols)
                    aux = zero_sp
                    do i = colptr(j), colptr(j+1)-1
                        aux = aux + data(i) * vec_x(:,row(i))
                    end do
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = vec_x(:,j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(:,row(i))
                        vec_y(:,row(i)) = vec_y(:,row(i)) + alpha_ * data(i) * vec_x(:,j)
                    end do
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = zero_csp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(:,j)
                        vec_y(:,i) = vec_y(:,i) + alpha_ * data(i) * vec_x(:,row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(:,j)
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_full .and. op_==sparse_op_hermitian ) then
                do concurrent(j=1:ncols)
                    aux = zero_sp
                    do i = colptr(j), colptr(j+1)-1
                        aux = aux + conjg(data(i)) * vec_x(:,row(i))
                    end do
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_lower .and. op_==sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = vec_x(:,j) * conjg(data(colptr(j)))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + conjg(data(i)) * vec_x(:,row(i))
                        vec_y(:,row(i)) = vec_y(:,row(i)) + alpha_ * conjg(data(i)) * vec_x(:,j)
                    end do
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_upper .and. op_==sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = zero_csp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + conjg(data(i)) * vec_x(:,j)
                        vec_y(:,i) = vec_y(:,i) + alpha_ * conjg(data(i)) * vec_x(:,row(i))
                    end do
                    aux = aux + conjg(data(colptr(j))) * vec_x(:,j)
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do
            end if
        end associate
    end subroutine
    
    subroutine spmv_csc_1d_cdp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSC_cdp_type), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        complex(dp) :: aux

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_cdp
        endif

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent(j=1:ncols)
                    aux = alpha_ * vec_x(j)
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(row(i)) = vec_y(row(i)) + data(i) * aux
                    end do
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent(j=1:ncols)
                    aux = zero_dp
                    do i = colptr(j), colptr(j+1)-1
                        aux = aux + data(i) * vec_x(row(i))
                    end do
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = vec_x(j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(row(i))
                        vec_y(row(i)) = vec_y(row(i)) + alpha_ * data(i) * vec_x(j)
                    end do
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = zero_cdp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(j)
                        vec_y(i) = vec_y(i) + alpha_ * data(i) * vec_x(row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(j)
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_full .and. op_==sparse_op_hermitian ) then
                do concurrent(j=1:ncols)
                    aux = zero_dp
                    do i = colptr(j), colptr(j+1)-1
                        aux = aux + conjg(data(i)) * vec_x(row(i))
                    end do
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_lower .and. op_==sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = vec_x(j) * conjg(data(colptr(j)))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + conjg(data(i)) * vec_x(row(i))
                        vec_y(row(i)) = vec_y(row(i)) + alpha_ * conjg(data(i)) * vec_x(j)
                    end do
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do

            else if( storage == sparse_upper .and. op_==sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = zero_cdp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + conjg(data(i)) * vec_x(j)
                        vec_y(i) = vec_y(i) + alpha_ * conjg(data(i)) * vec_x(row(i))
                    end do
                    aux = aux + conjg(data(colptr(j))) * vec_x(j)
                    vec_y(j) = vec_y(j) + alpha_ * aux
                end do
            end if
        end associate
    end subroutine
    
    subroutine spmv_csc_2d_cdp(matrix,vec_x,vec_y,alpha,beta,op)
        type(CSC_cdp_type), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:,:)
        complex(dp), intent(inout) :: vec_y(:,:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j
        complex(dp) :: aux(size(vec_x,dim=1))

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_cdp
        endif

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent(j=1:ncols)
                    aux = alpha_ * vec_x(:,j)
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(:,row(i)) = vec_y(:,row(i)) + data(i) * aux
                    end do
                end do

            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent(j=1:ncols)
                    aux = zero_dp
                    do i = colptr(j), colptr(j+1)-1
                        aux = aux + data(i) * vec_x(:,row(i))
                    end do
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_lower .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = vec_x(:,j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(:,row(i))
                        vec_y(:,row(i)) = vec_y(:,row(i)) + alpha_ * data(i) * vec_x(:,j)
                    end do
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_upper .and. op_/=sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = zero_cdp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(:,j)
                        vec_y(:,i) = vec_y(:,i) + alpha_ * data(i) * vec_x(:,row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(:,j)
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_full .and. op_==sparse_op_hermitian ) then
                do concurrent(j=1:ncols)
                    aux = zero_dp
                    do i = colptr(j), colptr(j+1)-1
                        aux = aux + conjg(data(i)) * vec_x(:,row(i))
                    end do
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_lower .and. op_==sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = vec_x(:,j) * conjg(data(colptr(j)))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + conjg(data(i)) * vec_x(:,row(i))
                        vec_y(:,row(i)) = vec_y(:,row(i)) + alpha_ * conjg(data(i)) * vec_x(:,j)
                    end do
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do

            else if( storage == sparse_upper .and. op_==sparse_op_hermitian )then
                do j = 1 , ncols
                    aux  = zero_cdp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + conjg(data(i)) * vec_x(:,j)
                        vec_y(:,i) = vec_y(:,i) + alpha_ * conjg(data(i)) * vec_x(:,row(i))
                    end do
                    aux = aux + conjg(data(colptr(j))) * vec_x(:,j)
                    vec_y(:,j) = vec_y(:,j) + alpha_ * aux
                end do
            end if
        end associate
    end subroutine
    

    !! spmv_ell
    subroutine spmv_ell_1d_sp(matrix,vec_x,vec_y,alpha,beta,op)
        type(ELL_sp_type), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j, k
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_sp
        endif
        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(i) = vec_y(i) + alpha_*data(i,k) * vec_x(j)
                end do
            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(j) = vec_y(j) + alpha_*data(i,k) * vec_x(i)
                end do
            else if( storage /= sparse_full .and. op_/=sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) then
                        vec_y(i) = vec_y(i) + alpha_*data(i,k) * vec_x(j)
                        if(i==j) cycle 
                        vec_y(j) = vec_y(j) + alpha_*data(i,k) * vec_x(i)
                    end if
                end do
            end if
        end associate
    end subroutine
    
    subroutine spmv_ell_2d_sp(matrix,vec_x,vec_y,alpha,beta,op)
        type(ELL_sp_type), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:,:)
        real(sp), intent(inout) :: vec_y(:,:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j, k
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_sp
        endif
        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,i) = vec_y(:,i) + alpha_*data(i,k) * vec_x(:,j)
                end do
            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,j) = vec_y(:,j) + alpha_*data(i,k) * vec_x(:,i)
                end do
            else if( storage /= sparse_full .and. op_/=sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) then
                        vec_y(:,i) = vec_y(:,i) + alpha_*data(i,k) * vec_x(:,j)
                        if(i==j) cycle 
                        vec_y(:,j) = vec_y(:,j) + alpha_*data(i,k) * vec_x(:,i)
                    end if
                end do
            end if
        end associate
    end subroutine
    
    subroutine spmv_ell_1d_dp(matrix,vec_x,vec_y,alpha,beta,op)
        type(ELL_dp_type), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j, k
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_dp
        endif
        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(i) = vec_y(i) + alpha_*data(i,k) * vec_x(j)
                end do
            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(j) = vec_y(j) + alpha_*data(i,k) * vec_x(i)
                end do
            else if( storage /= sparse_full .and. op_/=sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) then
                        vec_y(i) = vec_y(i) + alpha_*data(i,k) * vec_x(j)
                        if(i==j) cycle 
                        vec_y(j) = vec_y(j) + alpha_*data(i,k) * vec_x(i)
                    end if
                end do
            end if
        end associate
    end subroutine
    
    subroutine spmv_ell_2d_dp(matrix,vec_x,vec_y,alpha,beta,op)
        type(ELL_dp_type), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:,:)
        real(dp), intent(inout) :: vec_y(:,:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j, k
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_dp
        endif
        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,i) = vec_y(:,i) + alpha_*data(i,k) * vec_x(:,j)
                end do
            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,j) = vec_y(:,j) + alpha_*data(i,k) * vec_x(:,i)
                end do
            else if( storage /= sparse_full .and. op_/=sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) then
                        vec_y(:,i) = vec_y(:,i) + alpha_*data(i,k) * vec_x(:,j)
                        if(i==j) cycle 
                        vec_y(:,j) = vec_y(:,j) + alpha_*data(i,k) * vec_x(:,i)
                    end if
                end do
            end if
        end associate
    end subroutine
    
    subroutine spmv_ell_1d_csp(matrix,vec_x,vec_y,alpha,beta,op)
        type(ELL_csp_type), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j, k
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_csp
        endif
        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(i) = vec_y(i) + alpha_*data(i,k) * vec_x(j)
                end do
            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(j) = vec_y(j) + alpha_*data(i,k) * vec_x(i)
                end do
            else if( storage /= sparse_full .and. op_/=sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) then
                        vec_y(i) = vec_y(i) + alpha_*data(i,k) * vec_x(j)
                        if(i==j) cycle 
                        vec_y(j) = vec_y(j) + alpha_*data(i,k) * vec_x(i)
                    end if
                end do
            else if( storage == sparse_full .and. op_==sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(j) = vec_y(j) + alpha_*conjg(data(i,k)) * vec_x(i)
                end do
            else if( storage /= sparse_full .and. op_==sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) then
                        vec_y(i) = vec_y(i) + alpha_*conjg(data(i,k)) * vec_x(j)
                        if(i==j) cycle 
                        vec_y(j) = vec_y(j) + alpha_*conjg(data(i,k)) * vec_x(i)
                    end if
                end do
            end if
        end associate
    end subroutine
    
    subroutine spmv_ell_2d_csp(matrix,vec_x,vec_y,alpha,beta,op)
        type(ELL_csp_type), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:,:)
        complex(sp), intent(inout) :: vec_y(:,:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j, k
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_csp
        endif
        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,i) = vec_y(:,i) + alpha_*data(i,k) * vec_x(:,j)
                end do
            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,j) = vec_y(:,j) + alpha_*data(i,k) * vec_x(:,i)
                end do
            else if( storage /= sparse_full .and. op_/=sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) then
                        vec_y(:,i) = vec_y(:,i) + alpha_*data(i,k) * vec_x(:,j)
                        if(i==j) cycle 
                        vec_y(:,j) = vec_y(:,j) + alpha_*data(i,k) * vec_x(:,i)
                    end if
                end do
            else if( storage == sparse_full .and. op_==sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,j) = vec_y(:,j) + alpha_*conjg(data(i,k)) * vec_x(:,i)
                end do
            else if( storage /= sparse_full .and. op_==sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) then
                        vec_y(:,i) = vec_y(:,i) + alpha_*conjg(data(i,k)) * vec_x(:,j)
                        if(i==j) cycle 
                        vec_y(:,j) = vec_y(:,j) + alpha_*conjg(data(i,k)) * vec_x(:,i)
                    end if
                end do
            end if
        end associate
    end subroutine
    
    subroutine spmv_ell_1d_cdp(matrix,vec_x,vec_y,alpha,beta,op)
        type(ELL_cdp_type), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j, k
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_cdp
        endif
        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(i) = vec_y(i) + alpha_*data(i,k) * vec_x(j)
                end do
            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(j) = vec_y(j) + alpha_*data(i,k) * vec_x(i)
                end do
            else if( storage /= sparse_full .and. op_/=sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) then
                        vec_y(i) = vec_y(i) + alpha_*data(i,k) * vec_x(j)
                        if(i==j) cycle 
                        vec_y(j) = vec_y(j) + alpha_*data(i,k) * vec_x(i)
                    end if
                end do
            else if( storage == sparse_full .and. op_==sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(j) = vec_y(j) + alpha_*conjg(data(i,k)) * vec_x(i)
                end do
            else if( storage /= sparse_full .and. op_==sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) then
                        vec_y(i) = vec_y(i) + alpha_*conjg(data(i,k)) * vec_x(j)
                        if(i==j) cycle 
                        vec_y(j) = vec_y(j) + alpha_*conjg(data(i,k)) * vec_x(i)
                    end if
                end do
            end if
        end associate
    end subroutine
    
    subroutine spmv_ell_2d_cdp(matrix,vec_x,vec_y,alpha,beta,op)
        type(ELL_cdp_type), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:,:)
        complex(dp), intent(inout) :: vec_y(:,:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, j, k
        
        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_cdp
        endif
        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage )
            if( storage == sparse_full .and. op_==sparse_op_none ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,i) = vec_y(:,i) + alpha_*data(i,k) * vec_x(:,j)
                end do
            else if( storage == sparse_full .and. op_==sparse_op_transpose ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,j) = vec_y(:,j) + alpha_*data(i,k) * vec_x(:,i)
                end do
            else if( storage /= sparse_full .and. op_/=sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) then
                        vec_y(:,i) = vec_y(:,i) + alpha_*data(i,k) * vec_x(:,j)
                        if(i==j) cycle 
                        vec_y(:,j) = vec_y(:,j) + alpha_*data(i,k) * vec_x(:,i)
                    end if
                end do
            else if( storage == sparse_full .and. op_==sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,j) = vec_y(:,j) + alpha_*conjg(data(i,k)) * vec_x(:,i)
                end do
            else if( storage /= sparse_full .and. op_==sparse_op_hermitian ) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) then
                        vec_y(:,i) = vec_y(:,i) + alpha_*conjg(data(i,k)) * vec_x(:,j)
                        if(i==j) cycle 
                        vec_y(:,j) = vec_y(:,j) + alpha_*conjg(data(i,k)) * vec_x(:,i)
                    end if
                end do
            end if
        end associate
    end subroutine
    

    !! spmv_sellc
    subroutine spmv_sellc_sp(matrix,vec_x,vec_y,alpha,beta,op)
        !! This algorithm was gracefully provided by Ivan Privec and adapted by Jose Alves
        type(SELLC_sp_type), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, nz, rowidx, num_chunks, rm

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_sp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_sp
        endif

        associate( data => matrix%data, ia => matrix%rowptr , ja => matrix%col, cs => matrix%chunk_size, &
        &   nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage  )

        if( .not.any( [4, 8, 16] == cs ) ) then
            print *, "error: sellc chunk size not supported."
            return
        end if

        num_chunks = nrows / cs
        rm = nrows - num_chunks * cs
        if( storage == sparse_full .and. op_==sparse_op_none ) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_rm(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
            end if
            
        else if( storage == sparse_full .and. op_==sparse_op_transpose ) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_trans_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_trans_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_trans_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_rm_trans(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
            end if
        
        else
            print *, "error: sellc format for spmv operation not yet supported."
            return
        end if
        end associate

    contains
        pure subroutine chunk_kernel_4(n,a,col,x,y)
            integer, value      :: n
            real(sp), intent(in)  :: a(4,n), x(:)
            integer(ilp), intent(in) :: col(4,n)
            real(sp), intent(inout) :: y(4)
            integer :: j
            do j = 1, n
                y(:) = y(:) + alpha_ * a(:,j) * x(col(:,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_trans_4(n,a,col,x,y)
            integer, value      :: n
            real(sp), intent(in)  :: a(4,n), x(4)
            integer(ilp), intent(in) :: col(4,n)
            real(sp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 4
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_8(n,a,col,x,y)
            integer, value      :: n
            real(sp), intent(in)  :: a(8,n), x(:)
            integer(ilp), intent(in) :: col(8,n)
            real(sp), intent(inout) :: y(8)
            integer :: j
            do j = 1, n
                y(:) = y(:) + alpha_ * a(:,j) * x(col(:,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_trans_8(n,a,col,x,y)
            integer, value      :: n
            real(sp), intent(in)  :: a(8,n), x(8)
            integer(ilp), intent(in) :: col(8,n)
            real(sp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 8
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_16(n,a,col,x,y)
            integer, value      :: n
            real(sp), intent(in)  :: a(16,n), x(:)
            integer(ilp), intent(in) :: col(16,n)
            real(sp), intent(inout) :: y(16)
            integer :: j
            do j = 1, n
                y(:) = y(:) + alpha_ * a(:,j) * x(col(:,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_trans_16(n,a,col,x,y)
            integer, value      :: n
            real(sp), intent(in)  :: a(16,n), x(16)
            integer(ilp), intent(in) :: col(16,n)
            real(sp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 16
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
    
        pure subroutine chunk_kernel_rm(n,cs,r,a,col,x,y)
            integer, value      :: n, cs, r
            real(sp), intent(in)  :: a(cs,n), x(:)
            integer(ilp), intent(in) :: col(cs,n)
            real(sp), intent(inout) :: y(r)
            integer :: j
            do j = 1, n
                y(1:r) = y(1:r) + alpha_ * a(1:r,j) * x(col(1:r,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_rm_trans(n,cs,r,a,col,x,y)
            integer, value      :: n, cs, r
            real(sp), intent(in)  :: a(cs,n), x(r)
            integer(ilp), intent(in) :: col(cs,n)
            real(sp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, r
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine

    end subroutine
    
    subroutine spmv_sellc_dp(matrix,vec_x,vec_y,alpha,beta,op)
        !! This algorithm was gracefully provided by Ivan Privec and adapted by Jose Alves
        type(SELLC_dp_type), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, nz, rowidx, num_chunks, rm

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_dp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_dp
        endif

        associate( data => matrix%data, ia => matrix%rowptr , ja => matrix%col, cs => matrix%chunk_size, &
        &   nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage  )

        if( .not.any( [4, 8, 16] == cs ) ) then
            print *, "error: sellc chunk size not supported."
            return
        end if

        num_chunks = nrows / cs
        rm = nrows - num_chunks * cs
        if( storage == sparse_full .and. op_==sparse_op_none ) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_rm(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
            end if
            
        else if( storage == sparse_full .and. op_==sparse_op_transpose ) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_trans_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_trans_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_trans_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_rm_trans(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
            end if
        
        else
            print *, "error: sellc format for spmv operation not yet supported."
            return
        end if
        end associate

    contains
        pure subroutine chunk_kernel_4(n,a,col,x,y)
            integer, value      :: n
            real(dp), intent(in)  :: a(4,n), x(:)
            integer(ilp), intent(in) :: col(4,n)
            real(dp), intent(inout) :: y(4)
            integer :: j
            do j = 1, n
                y(:) = y(:) + alpha_ * a(:,j) * x(col(:,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_trans_4(n,a,col,x,y)
            integer, value      :: n
            real(dp), intent(in)  :: a(4,n), x(4)
            integer(ilp), intent(in) :: col(4,n)
            real(dp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 4
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_8(n,a,col,x,y)
            integer, value      :: n
            real(dp), intent(in)  :: a(8,n), x(:)
            integer(ilp), intent(in) :: col(8,n)
            real(dp), intent(inout) :: y(8)
            integer :: j
            do j = 1, n
                y(:) = y(:) + alpha_ * a(:,j) * x(col(:,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_trans_8(n,a,col,x,y)
            integer, value      :: n
            real(dp), intent(in)  :: a(8,n), x(8)
            integer(ilp), intent(in) :: col(8,n)
            real(dp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 8
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_16(n,a,col,x,y)
            integer, value      :: n
            real(dp), intent(in)  :: a(16,n), x(:)
            integer(ilp), intent(in) :: col(16,n)
            real(dp), intent(inout) :: y(16)
            integer :: j
            do j = 1, n
                y(:) = y(:) + alpha_ * a(:,j) * x(col(:,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_trans_16(n,a,col,x,y)
            integer, value      :: n
            real(dp), intent(in)  :: a(16,n), x(16)
            integer(ilp), intent(in) :: col(16,n)
            real(dp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 16
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
    
        pure subroutine chunk_kernel_rm(n,cs,r,a,col,x,y)
            integer, value      :: n, cs, r
            real(dp), intent(in)  :: a(cs,n), x(:)
            integer(ilp), intent(in) :: col(cs,n)
            real(dp), intent(inout) :: y(r)
            integer :: j
            do j = 1, n
                y(1:r) = y(1:r) + alpha_ * a(1:r,j) * x(col(1:r,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_rm_trans(n,cs,r,a,col,x,y)
            integer, value      :: n, cs, r
            real(dp), intent(in)  :: a(cs,n), x(r)
            integer(ilp), intent(in) :: col(cs,n)
            real(dp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, r
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine

    end subroutine
    
    subroutine spmv_sellc_csp(matrix,vec_x,vec_y,alpha,beta,op)
        !! This algorithm was gracefully provided by Ivan Privec and adapted by Jose Alves
        type(SELLC_csp_type), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(sp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, nz, rowidx, num_chunks, rm

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_csp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_csp
        endif

        associate( data => matrix%data, ia => matrix%rowptr , ja => matrix%col, cs => matrix%chunk_size, &
        &   nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage  )

        if( .not.any( [4, 8, 16] == cs ) ) then
            print *, "error: sellc chunk size not supported."
            return
        end if

        num_chunks = nrows / cs
        rm = nrows - num_chunks * cs
        if( storage == sparse_full .and. op_==sparse_op_none ) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_rm(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
            end if
            
        else if( storage == sparse_full .and. op_==sparse_op_transpose ) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_trans_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_trans_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_trans_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_rm_trans(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
            end if
        
        else if( storage == sparse_full .and. op_==sparse_op_hermitian ) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_herm_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_herm_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_herm_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_rm_herm(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
            end if
        else
            print *, "error: sellc format for spmv operation not yet supported."
            return
        end if
        end associate

    contains
        pure subroutine chunk_kernel_4(n,a,col,x,y)
            integer, value      :: n
            complex(sp), intent(in)  :: a(4,n), x(:)
            integer(ilp), intent(in) :: col(4,n)
            complex(sp), intent(inout) :: y(4)
            integer :: j
            do j = 1, n
                y(:) = y(:) + alpha_ * a(:,j) * x(col(:,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_trans_4(n,a,col,x,y)
            integer, value      :: n
            complex(sp), intent(in)  :: a(4,n), x(4)
            integer(ilp), intent(in) :: col(4,n)
            complex(sp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 4
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_herm_4(n,a,col,x,y)
            integer, value      :: n
            complex(sp), intent(in)  :: a(4,n), x(4)
            integer(ilp), intent(in) :: col(4,n)
            complex(sp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 4
                    y(col(k,j)) = y(col(k,j)) + alpha_ * conjg(a(k,j)) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_8(n,a,col,x,y)
            integer, value      :: n
            complex(sp), intent(in)  :: a(8,n), x(:)
            integer(ilp), intent(in) :: col(8,n)
            complex(sp), intent(inout) :: y(8)
            integer :: j
            do j = 1, n
                y(:) = y(:) + alpha_ * a(:,j) * x(col(:,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_trans_8(n,a,col,x,y)
            integer, value      :: n
            complex(sp), intent(in)  :: a(8,n), x(8)
            integer(ilp), intent(in) :: col(8,n)
            complex(sp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 8
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_herm_8(n,a,col,x,y)
            integer, value      :: n
            complex(sp), intent(in)  :: a(8,n), x(8)
            integer(ilp), intent(in) :: col(8,n)
            complex(sp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 8
                    y(col(k,j)) = y(col(k,j)) + alpha_ * conjg(a(k,j)) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_16(n,a,col,x,y)
            integer, value      :: n
            complex(sp), intent(in)  :: a(16,n), x(:)
            integer(ilp), intent(in) :: col(16,n)
            complex(sp), intent(inout) :: y(16)
            integer :: j
            do j = 1, n
                y(:) = y(:) + alpha_ * a(:,j) * x(col(:,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_trans_16(n,a,col,x,y)
            integer, value      :: n
            complex(sp), intent(in)  :: a(16,n), x(16)
            integer(ilp), intent(in) :: col(16,n)
            complex(sp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 16
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_herm_16(n,a,col,x,y)
            integer, value      :: n
            complex(sp), intent(in)  :: a(16,n), x(16)
            integer(ilp), intent(in) :: col(16,n)
            complex(sp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 16
                    y(col(k,j)) = y(col(k,j)) + alpha_ * conjg(a(k,j)) * x(k)
                end do
            end do
        end subroutine
    
        pure subroutine chunk_kernel_rm(n,cs,r,a,col,x,y)
            integer, value      :: n, cs, r
            complex(sp), intent(in)  :: a(cs,n), x(:)
            integer(ilp), intent(in) :: col(cs,n)
            complex(sp), intent(inout) :: y(r)
            integer :: j
            do j = 1, n
                y(1:r) = y(1:r) + alpha_ * a(1:r,j) * x(col(1:r,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_rm_trans(n,cs,r,a,col,x,y)
            integer, value      :: n, cs, r
            complex(sp), intent(in)  :: a(cs,n), x(r)
            integer(ilp), intent(in) :: col(cs,n)
            complex(sp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, r
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_rm_herm(n,cs,r,a,col,x,y)
            integer, value      :: n, cs, r
            complex(sp), intent(in)  :: a(cs,n), x(r)
            integer(ilp), intent(in) :: col(cs,n)
            complex(sp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, r
                    y(col(k,j)) = y(col(k,j)) + alpha_ * conjg(a(k,j)) * x(k)
                end do
            end do
        end subroutine

    end subroutine
    
    subroutine spmv_sellc_cdp(matrix,vec_x,vec_y,alpha,beta,op)
        !! This algorithm was gracefully provided by Ivan Privec and adapted by Jose Alves
        type(SELLC_cdp_type), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(dp) :: alpha_
        character(1) :: op_
        integer(ilp) :: i, nz, rowidx, num_chunks, rm

        op_ = sparse_op_none; if(present(op)) op_ = op
        alpha_ = one_cdp
        if(present(alpha)) alpha_ = alpha
        if(present(beta)) then
            vec_y = beta * vec_y
        else 
            vec_y = zero_cdp
        endif

        associate( data => matrix%data, ia => matrix%rowptr , ja => matrix%col, cs => matrix%chunk_size, &
        &   nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, storage => matrix%storage  )

        if( .not.any( [4, 8, 16] == cs ) ) then
            print *, "error: sellc chunk size not supported."
            return
        end if

        num_chunks = nrows / cs
        rm = nrows - num_chunks * cs
        if( storage == sparse_full .and. op_==sparse_op_none ) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_rm(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
            end if
            
        else if( storage == sparse_full .and. op_==sparse_op_transpose ) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_trans_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_trans_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_trans_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_rm_trans(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
            end if
        
        else if( storage == sparse_full .and. op_==sparse_op_hermitian ) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_herm_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_herm_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_herm_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
                end do
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_rm_herm(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x(rowidx:),vec_y)
            end if
        else
            print *, "error: sellc format for spmv operation not yet supported."
            return
        end if
        end associate

    contains
        pure subroutine chunk_kernel_4(n,a,col,x,y)
            integer, value      :: n
            complex(dp), intent(in)  :: a(4,n), x(:)
            integer(ilp), intent(in) :: col(4,n)
            complex(dp), intent(inout) :: y(4)
            integer :: j
            do j = 1, n
                y(:) = y(:) + alpha_ * a(:,j) * x(col(:,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_trans_4(n,a,col,x,y)
            integer, value      :: n
            complex(dp), intent(in)  :: a(4,n), x(4)
            integer(ilp), intent(in) :: col(4,n)
            complex(dp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 4
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_herm_4(n,a,col,x,y)
            integer, value      :: n
            complex(dp), intent(in)  :: a(4,n), x(4)
            integer(ilp), intent(in) :: col(4,n)
            complex(dp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 4
                    y(col(k,j)) = y(col(k,j)) + alpha_ * conjg(a(k,j)) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_8(n,a,col,x,y)
            integer, value      :: n
            complex(dp), intent(in)  :: a(8,n), x(:)
            integer(ilp), intent(in) :: col(8,n)
            complex(dp), intent(inout) :: y(8)
            integer :: j
            do j = 1, n
                y(:) = y(:) + alpha_ * a(:,j) * x(col(:,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_trans_8(n,a,col,x,y)
            integer, value      :: n
            complex(dp), intent(in)  :: a(8,n), x(8)
            integer(ilp), intent(in) :: col(8,n)
            complex(dp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 8
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_herm_8(n,a,col,x,y)
            integer, value      :: n
            complex(dp), intent(in)  :: a(8,n), x(8)
            integer(ilp), intent(in) :: col(8,n)
            complex(dp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 8
                    y(col(k,j)) = y(col(k,j)) + alpha_ * conjg(a(k,j)) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_16(n,a,col,x,y)
            integer, value      :: n
            complex(dp), intent(in)  :: a(16,n), x(:)
            integer(ilp), intent(in) :: col(16,n)
            complex(dp), intent(inout) :: y(16)
            integer :: j
            do j = 1, n
                y(:) = y(:) + alpha_ * a(:,j) * x(col(:,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_trans_16(n,a,col,x,y)
            integer, value      :: n
            complex(dp), intent(in)  :: a(16,n), x(16)
            integer(ilp), intent(in) :: col(16,n)
            complex(dp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 16
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_herm_16(n,a,col,x,y)
            integer, value      :: n
            complex(dp), intent(in)  :: a(16,n), x(16)
            integer(ilp), intent(in) :: col(16,n)
            complex(dp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, 16
                    y(col(k,j)) = y(col(k,j)) + alpha_ * conjg(a(k,j)) * x(k)
                end do
            end do
        end subroutine
    
        pure subroutine chunk_kernel_rm(n,cs,r,a,col,x,y)
            integer, value      :: n, cs, r
            complex(dp), intent(in)  :: a(cs,n), x(:)
            integer(ilp), intent(in) :: col(cs,n)
            complex(dp), intent(inout) :: y(r)
            integer :: j
            do j = 1, n
                y(1:r) = y(1:r) + alpha_ * a(1:r,j) * x(col(1:r,j))
            end do
        end subroutine
        pure subroutine chunk_kernel_rm_trans(n,cs,r,a,col,x,y)
            integer, value      :: n, cs, r
            complex(dp), intent(in)  :: a(cs,n), x(r)
            integer(ilp), intent(in) :: col(cs,n)
            complex(dp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, r
                    y(col(k,j)) = y(col(k,j)) + alpha_ * a(k,j) * x(k)
                end do
            end do
        end subroutine
        pure subroutine chunk_kernel_rm_herm(n,cs,r,a,col,x,y)
            integer, value      :: n, cs, r
            complex(dp), intent(in)  :: a(cs,n), x(r)
            integer(ilp), intent(in) :: col(cs,n)
            complex(dp), intent(inout) :: y(:)
            integer :: j, k
            do j = 1, n
                do k = 1, r
                    y(col(k,j)) = y(col(k,j)) + alpha_ * conjg(a(k,j)) * x(k)
                end do
            end do
        end subroutine

    end subroutine
    

end module
