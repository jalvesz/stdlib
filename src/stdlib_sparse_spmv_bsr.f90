submodule (stdlib_sparse_spmv) stdlib_sparse_spmv_bsr
    use stdlib_linalg_blas, only: gemv, gemm
contains

    !! spmv_bsr
    module subroutine spmv_bsr_1d_sp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_sp_type), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(sp) :: alpha_, beta_
        character(1) :: op_
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

            nbrow = nrows / br

            if( storage == sparse_full .and. op_ == sparse_op_none ) then
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
                    end do
                end do

            else if( storage == sparse_full .and. op_ == sparse_op_transpose ) then
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

            else if( storage == sparse_lower .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_i_start:), incx=1, beta=one_sp, &
                            y=vec_y(row_j_start:), incy=1)
                    end do
                end do

            else if( storage == sparse_upper .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_i_start:), incx=1, beta=one_sp, &
                            y=vec_y(row_j_start:), incy=1)
                    end do
                end do

            end if
        end associate
    end subroutine

    module subroutine spmv_bsr_2d_sp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_sp_type), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:,:)
        real(sp), intent(inout) :: vec_y(:,:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(sp) :: alpha_, beta_
        character(1) :: op_
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

            nbrow = nrows / br
            nrhs = size(vec_x,dim=2)

            if( storage == sparse_full .and. op_ == sparse_op_none ) then
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
                    end do
                end do

            else if( storage == sparse_full .and. op_ == sparse_op_transpose ) then
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

            else if( storage == sparse_lower .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_sp, &
                            c=vec_y(row_j_start:, :), ldc=i_max_j)
                    end do
                end do

            else if( storage == sparse_upper .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_sp, &
                            c=vec_y(row_j_start:, :), ldc=i_max_j)
                    end do
                end do

            end if
        end associate
    end subroutine

    module subroutine spmv_bsr_1d_dp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_dp_type), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(dp) :: alpha_, beta_
        character(1) :: op_
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

            nbrow = nrows / br

            if( storage == sparse_full .and. op_ == sparse_op_none ) then
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
                    end do
                end do

            else if( storage == sparse_full .and. op_ == sparse_op_transpose ) then
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

            else if( storage == sparse_lower .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_i_start:), incx=1, beta=one_dp, &
                            y=vec_y(row_j_start:), incy=1)
                    end do
                end do

            else if( storage == sparse_upper .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_i_start:), incx=1, beta=one_dp, &
                            y=vec_y(row_j_start:), incy=1)
                    end do
                end do

            end if
        end associate
    end subroutine

    module subroutine spmv_bsr_2d_dp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_dp_type), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:,:)
        real(dp), intent(inout) :: vec_y(:,:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        real(dp) :: alpha_, beta_
        character(1) :: op_
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

            nbrow = nrows / br
            nrhs = size(vec_x,dim=2)

            if( storage == sparse_full .and. op_ == sparse_op_none ) then
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
                    end do
                end do

            else if( storage == sparse_full .and. op_ == sparse_op_transpose ) then
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

            else if( storage == sparse_lower .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_dp, &
                            c=vec_y(row_j_start:, :), ldc=i_max_j)
                    end do
                end do

            else if( storage == sparse_upper .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_dp, &
                            c=vec_y(row_j_start:, :), ldc=i_max_j)
                    end do
                end do

            end if
        end associate
    end subroutine

    module subroutine spmv_bsr_1d_csp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_csp_type), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(sp) :: alpha_, beta_
        character(1) :: op_
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

            nbrow = nrows / br

            if( storage == sparse_full .and. op_ == sparse_op_none ) then
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
                    end do
                end do

            else if( storage == sparse_full .and. op_ == sparse_op_transpose ) then
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

            else if( storage == sparse_lower .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_i_start:), incx=1, beta=one_csp, &
                            y=vec_y(row_j_start:), incy=1)
                    end do
                end do

            else if( storage == sparse_upper .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_i_start:), incx=1, beta=one_csp, &
                            y=vec_y(row_j_start:), incy=1)
                    end do
                end do

            else if( storage == sparse_full .and. op_ == sparse_op_hermitian ) then
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemv('C', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(row_i_start:), incx=1, beta=one_csp, &
                            y=vec_y(col_j_start:), incy=1)
                    end do
                end do

            else if( storage == sparse_lower .and. op_ == sparse_op_hermitian ) then
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

                        xwork_bc(1:j_max_j) = conjg(vec_x(col_j_start:col_j_start+j_max_j-1))
                        ywork_br(1:i_max_i) = zero_csp
                        call gemv('N', m=i_max_i, n=j_max_j, alpha=one_csp, a=data(1:,1:,p), lda=br, &
                            x=xwork_bc, incx=1, beta=zero_csp, y=ywork_br, incy=1)
                        vec_y(row_i_start:row_i_start+i_max_i-1) = vec_y(row_i_start:row_i_start+i_max_i-1) + &
                            alpha_ * conjg(ywork_br(1:i_max_i))

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemv('C', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_i_start:), incx=1, beta=one_csp, &
                            y=vec_y(row_j_start:), incy=1)
                    end do
                end do

            else if( storage == sparse_upper .and. op_ == sparse_op_hermitian ) then
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

                        xwork_bc(1:j_max_j) = conjg(vec_x(col_j_start:col_j_start+j_max_j-1))
                        ywork_br(1:i_max_i) = zero_csp
                        call gemv('N', m=i_max_i, n=j_max_j, alpha=one_csp, a=data(1:,1:,p), lda=br, &
                            x=xwork_bc, incx=1, beta=zero_csp, y=ywork_br, incy=1)
                        vec_y(row_i_start:row_i_start+i_max_i-1) = vec_y(row_i_start:row_i_start+i_max_i-1) + &
                            alpha_ * conjg(ywork_br(1:i_max_i))

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemv('C', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_i_start:), incx=1, beta=one_csp, &
                            y=vec_y(row_j_start:), incy=1)
                    end do
                end do
            end if
        end associate
    end subroutine

    module subroutine spmv_bsr_2d_csp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_csp_type), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:,:)
        complex(sp), intent(inout) :: vec_y(:,:)
        complex(sp), intent(in), optional :: alpha
        complex(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(sp) :: alpha_, beta_
        character(1) :: op_
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

            nbrow = nrows / br
            nrhs = size(vec_x,dim=2)

            if( storage == sparse_full .and. op_ == sparse_op_none ) then
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
                    end do
                end do

            else if( storage == sparse_full .and. op_ == sparse_op_transpose ) then
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

            else if( storage == sparse_lower .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_csp, &
                            c=vec_y(row_j_start:, :), ldc=i_max_j)
                    end do
                end do

            else if( storage == sparse_upper .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_csp, &
                            c=vec_y(row_j_start:, :), ldc=i_max_j)
                    end do
                end do

            else if( storage == sparse_full .and. op_ == sparse_op_hermitian ) then
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

            else if( storage == sparse_lower .and. op_ == sparse_op_hermitian ) then
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

                        xwork_bc(1:j_max_j, :) = conjg(vec_x(col_j_start:col_j_start+j_max_j-1, :))
                        ywork_br(1:i_max_i, :) = zero_csp
                        call gemm('N','N', m=i_max_i, n=nrhs, k=j_max_j, alpha=one_csp, a=data(1:,1:,p), lda=br, &
                            b=xwork_bc, ldb=bc, beta=zero_csp, c=ywork_br, ldc=br)
                        vec_y(row_i_start:row_i_start+i_max_i-1, :) = vec_y(row_i_start:row_i_start+i_max_i-1, :) + &
                            alpha_ * conjg(ywork_br(1:i_max_i, :))

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemm('C','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_csp, &
                            c=vec_y(row_j_start:, :), ldc=i_max_j)
                    end do
                end do

            else if( storage == sparse_upper .and. op_ == sparse_op_hermitian ) then
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

                        xwork_bc(1:j_max_j, :) = conjg(vec_x(col_j_start:col_j_start+j_max_j-1, :))
                        ywork_br(1:i_max_i, :) = zero_csp
                        call gemm('N','N', m=i_max_i, n=nrhs, k=j_max_j, alpha=one_csp, a=data(1:,1:,p), lda=br, &
                            b=xwork_bc, ldb=bc, beta=zero_csp, c=ywork_br, ldc=br)
                        vec_y(row_i_start:row_i_start+i_max_i-1, :) = vec_y(row_i_start:row_i_start+i_max_i-1, :) + &
                            alpha_ * conjg(ywork_br(1:i_max_i, :))

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemm('C','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_csp, &
                            c=vec_y(row_j_start:, :), ldc=i_max_j)
                    end do
                end do
            end if
        end associate
    end subroutine

    module subroutine spmv_bsr_1d_cdp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_cdp_type), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(dp) :: alpha_, beta_
        character(1) :: op_
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

            nbrow = nrows / br

            if( storage == sparse_full .and. op_ == sparse_op_none ) then
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
                    end do
                end do

            else if( storage == sparse_full .and. op_ == sparse_op_transpose ) then
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

            else if( storage == sparse_lower .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_i_start:), incx=1, beta=one_cdp, &
                            y=vec_y(row_j_start:), incy=1)
                    end do
                end do

            else if( storage == sparse_upper .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemv('T', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_i_start:), incx=1, beta=one_cdp, &
                            y=vec_y(row_j_start:), incy=1)
                    end do
                end do

            else if( storage == sparse_full .and. op_ == sparse_op_hermitian ) then
                do bi = 1, nbrow
                    row_i_start = (bi - 1) * br + 1
                    if (row_i_start > nrows) cycle
                    i_max_i = min(br, nrows - row_i_start + 1)
                    do p = rowptr(bi), rowptr(bi+1)-1
                        bj = col(p)
                        col_j_start = (bj - 1) * bc + 1
                        if (col_j_start > ncols) cycle
                        j_max_j = min(bc, ncols - col_j_start + 1)

                        call gemv('C', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(row_i_start:), incx=1, beta=one_cdp, &
                            y=vec_y(col_j_start:), incy=1)
                    end do
                end do

            else if( storage == sparse_lower .and. op_ == sparse_op_hermitian ) then
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

                        xwork_bc(1:j_max_j) = conjg(vec_x(col_j_start:col_j_start+j_max_j-1))
                        ywork_br(1:i_max_i) = zero_cdp
                        call gemv('N', m=i_max_i, n=j_max_j, alpha=one_cdp, a=data(1:,1:,p), lda=br, &
                            x=xwork_bc, incx=1, beta=zero_cdp, y=ywork_br, incy=1)
                        vec_y(row_i_start:row_i_start+i_max_i-1) = vec_y(row_i_start:row_i_start+i_max_i-1) + &
                            alpha_ * conjg(ywork_br(1:i_max_i))

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemv('C', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_i_start:), incx=1, beta=one_cdp, &
                            y=vec_y(row_j_start:), incy=1)
                    end do
                end do

            else if( storage == sparse_upper .and. op_ == sparse_op_hermitian ) then
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

                        xwork_bc(1:j_max_j) = conjg(vec_x(col_j_start:col_j_start+j_max_j-1))
                        ywork_br(1:i_max_i) = zero_cdp
                        call gemv('N', m=i_max_i, n=j_max_j, alpha=one_cdp, a=data(1:,1:,p), lda=br, &
                            x=xwork_bc, incx=1, beta=zero_cdp, y=ywork_br, incy=1)
                        vec_y(row_i_start:row_i_start+i_max_i-1) = vec_y(row_i_start:row_i_start+i_max_i-1) + &
                            alpha_ * conjg(ywork_br(1:i_max_i))

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemv('C', m=i_max_i, n=j_max_j, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            x=vec_x(col_i_start:), incx=1, beta=one_cdp, &
                            y=vec_y(row_j_start:), incy=1)
                    end do
                end do
            end if
        end associate
    end subroutine

    module subroutine spmv_bsr_2d_cdp(matrix,vec_x,vec_y,alpha,beta,op)
        type(BSR_cdp_type), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:,:)
        complex(dp), intent(inout) :: vec_y(:,:)
        complex(dp), intent(in), optional :: alpha
        complex(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op
        complex(dp) :: alpha_, beta_
        character(1) :: op_
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

            nbrow = nrows / br
            nrhs = size(vec_x,dim=2)

            if( storage == sparse_full .and. op_ == sparse_op_none ) then
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
                    end do
                end do

            else if( storage == sparse_full .and. op_ == sparse_op_transpose ) then
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

            else if( storage == sparse_lower .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_cdp, &
                            c=vec_y(row_j_start:, :), ldc=i_max_j)
                    end do
                end do

            else if( storage == sparse_upper .and. op_ /= sparse_op_hermitian ) then
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

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemm('T','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_cdp, &
                            c=vec_y(row_j_start:, :), ldc=i_max_j)
                    end do
                end do

            else if( storage == sparse_full .and. op_ == sparse_op_hermitian ) then
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

            else if( storage == sparse_lower .and. op_ == sparse_op_hermitian ) then
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

                        xwork_bc(1:j_max_j, :) = conjg(vec_x(col_j_start:col_j_start+j_max_j-1, :))
                        ywork_br(1:i_max_i, :) = zero_cdp
                        call gemm('N','N', m=i_max_i, n=nrhs, k=j_max_j, alpha=one_cdp, a=data(1:,1:,p), lda=br, &
                            b=xwork_bc, ldb=bc, beta=zero_cdp, c=ywork_br, ldc=br)
                        vec_y(row_i_start:row_i_start+i_max_i-1, :) = vec_y(row_i_start:row_i_start+i_max_i-1, :) + &
                            alpha_ * conjg(ywork_br(1:i_max_i, :))

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemm('C','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_cdp, &
                            c=vec_y(row_j_start:, :), ldc=i_max_j)
                    end do
                end do

            else if( storage == sparse_upper .and. op_ == sparse_op_hermitian ) then
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

                        xwork_bc(1:j_max_j, :) = conjg(vec_x(col_j_start:col_j_start+j_max_j-1, :))
                        ywork_br(1:i_max_i, :) = zero_cdp
                        call gemm('N','N', m=i_max_i, n=nrhs, k=j_max_j, alpha=one_cdp, a=data(1:,1:,p), lda=br, &
                            b=xwork_bc, ldb=bc, beta=zero_cdp, c=ywork_br, ldc=br)
                        vec_y(row_i_start:row_i_start+i_max_i-1, :) = vec_y(row_i_start:row_i_start+i_max_i-1, :) + &
                            alpha_ * conjg(ywork_br(1:i_max_i, :))

                        if (bi == bj) cycle

                        row_j_start = (bj - 1) * br + 1
                        if (row_j_start > nrows) cycle
                        i_max_j = min(br, nrows - row_j_start + 1)
                        col_i_start = (bi - 1) * bc + 1
                        if (col_i_start > ncols) cycle
                        j_max_i = min(bc, ncols - col_i_start + 1)

                        call gemm('C','N', m=j_max_j, n=nrhs, k=i_max_i, alpha=alpha_, a=data(1:,1:,p), lda=br, &
                            b=vec_x(col_i_start:, :), ldb=j_max_i, beta=one_cdp, &
                            c=vec_y(row_j_start:, :), ldc=i_max_j)
                    end do
                end do
            end if
        end associate
    end subroutine


end submodule stdlib_sparse_spmv_bsr