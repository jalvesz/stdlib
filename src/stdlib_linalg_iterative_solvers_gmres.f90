
submodule(stdlib_linalg_iterative_solvers) stdlib_linalg_iterative_gmres
    use stdlib_kinds
    use stdlib_sparse
    use stdlib_constants
    use stdlib_optval, only: optval
    implicit none

contains

    module subroutine stdlib_solve_gmres_kernel_sp(A,M,b,x,rtol,atol,maxiter,kdim,workspace)
        class(stdlib_linop_sp_type), intent(in) :: A
        class(stdlib_linop_sp_type), intent(in) :: M
        real(sp), intent(in) :: b(:), rtol, atol
        real(sp), intent(inout) :: x(:)
        integer, intent(in) :: maxiter, kdim
        type(stdlib_solver_workspace_sp_type), intent(inout) :: workspace
        integer :: i, iter, inner_iter, j
        real(sp) :: beta, denom, hnext, norm_sq, norm_sq0, temp, tolsq
        real(sp), allocatable :: cs(:), g(:), h(:,:), sn(:), x_base(:), y(:)

        allocate(h(kdim+1, kdim), cs(kdim), sn(kdim), g(kdim+1), y(kdim), x_base(size(x)))

        associate( r => workspace%tmp(:,1), &
                   w => workspace%tmp(:,2), &
                   v => workspace%tmp(:,3:kdim+3), &
                   z => workspace%tmp(:,kdim+4:2*kdim+3) )

            ! Initialize convergence targets from the right-hand side norm.
            norm_sq0 = A%inner_product(b, b)
            tolsq = max(rtol*rtol*norm_sq0, atol*atol)

            ! Form the initial residual and report the starting iterate.
            r = b
            call A%matvec(x, r, alpha=-one_sp, beta=one_sp, op='N')
            norm_sq = A%inner_product(r, r)
            if (associated(workspace%callback)) call workspace%callback(x, norm_sq, 0)

            if (norm_sq <= tolsq .or. maxiter <= 0) then
                deallocate(h, cs, sn, g, y, x_base)
                return
            end if

            iter = 0
            do while (iter < maxiter .and. norm_sq >= tolsq)
                ! Start a new GMRES cycle from the current residual.
                beta = sqrt(max(norm_sq, zero_sp))
                if (beta <= epsilon(one_sp)) exit

                inner_iter = min(kdim, maxiter - iter)
                x_base = x
                h = zero_sp
                cs = zero_sp
                sn = zero_sp
                g = zero_sp
                y = zero_sp

                ! Initialize the Krylov basis and least-squares right-hand side.
                v(:,1) = r / beta
                g(1) = beta

                do j = 1, inner_iter
                    ! Run Arnoldi with the preconditioned basis vector.
                    call M%matvec(v(:,j), z(:,j), alpha=one_sp, beta=zero_sp, op='N')
                    call A%matvec(z(:,j), w, alpha=one_sp, beta=zero_sp, op='N')

                    do i = 1, j
                        h(i,j) = A%inner_product(v(:,i), w)
                        w = w - h(i,j) * v(:,i)
                    end do

                    hnext = sqrt(max(A%inner_product(w, w), zero_sp))
                    h(j+1,j) = hnext
                    if (hnext > epsilon(one_sp)) then
                        v(:,j+1) = w / hnext
                    else
                        v(:,j+1) = zero_sp
                    end if

                    ! Apply the previously accumulated Givens rotations.
                    do i = 1, j - 1
                        temp = cs(i) * h(i,j) + sn(i) * h(i+1,j)
                        h(i+1,j) = -sn(i) * h(i,j) + cs(i) * h(i+1,j)
                        h(i,j) = temp
                    end do

                    ! Build and apply the next Givens rotation.
                    denom = sqrt(h(j,j) * h(j,j) + h(j+1,j) * h(j+1,j))
                    if (denom > epsilon(one_sp)) then
                        cs(j) = h(j,j) / denom
                        sn(j) = h(j+1,j) / denom
                    else
                        cs(j) = one_sp
                        sn(j) = zero_sp
                    end if

                    temp = cs(j) * h(j,j) + sn(j) * h(j+1,j)
                    h(j+1,j) = -sn(j) * h(j,j) + cs(j) * h(j+1,j)
                    h(j,j) = temp

                    temp = cs(j) * g(j) + sn(j) * g(j+1)
                    g(j+1) = -sn(j) * g(j) + cs(j) * g(j+1)
                    g(j) = temp

                    ! Solve the reduced system and update the current iterate.
                    call upper_triangular_solve(h, g, y, j)
                    x = x_base
                    do i = 1, j
                        x = x + y(i) * z(:,i)
                    end do

                    ! Track the residual norm estimate for convergence.
                    norm_sq = g(j+1) * g(j+1)
                    iter = iter + 1
                    if (associated(workspace%callback)) call workspace%callback(x, norm_sq, iter)

                    if (norm_sq < tolsq .or. hnext <= epsilon(one_sp) .or. iter >= maxiter) exit
                end do

                if (norm_sq < tolsq .or. iter >= maxiter) exit

                ! Refresh the residual before the next restarted cycle.
                r = b
                call A%matvec(x, r, alpha=-one_sp, beta=one_sp, op='N')
                norm_sq = A%inner_product(r, r)
            end do
        end associate

        deallocate(h, cs, sn, g, y, x_base)

    contains

        subroutine upper_triangular_solve(h, g, y, n)
            real(sp), intent(in) :: h(:,:), g(:)
            real(sp), intent(inout) :: y(:)
            integer, intent(in) :: n
            integer :: row

            y(1:n) = g(1:n)
            do row = n, 1, -1
                if (row < n) y(row) = y(row) - dot_product(h(row,row+1:n) , y(row+1:n))
                if (abs(h(row,row)) > epsilon(one_sp)) then
                    y(row) = y(row) / h(row,row)
                else
                    y(row) = zero_sp
                end if
            end do
        end subroutine
    end subroutine
    module subroutine stdlib_solve_gmres_kernel_dp(A,M,b,x,rtol,atol,maxiter,kdim,workspace)
        class(stdlib_linop_dp_type), intent(in) :: A
        class(stdlib_linop_dp_type), intent(in) :: M
        real(dp), intent(in) :: b(:), rtol, atol
        real(dp), intent(inout) :: x(:)
        integer, intent(in) :: maxiter, kdim
        type(stdlib_solver_workspace_dp_type), intent(inout) :: workspace
        integer :: i, iter, inner_iter, j
        real(dp) :: beta, denom, hnext, norm_sq, norm_sq0, temp, tolsq
        real(dp), allocatable :: cs(:), g(:), h(:,:), sn(:), x_base(:), y(:)

        allocate(h(kdim+1, kdim), cs(kdim), sn(kdim), g(kdim+1), y(kdim), x_base(size(x)))

        associate( r => workspace%tmp(:,1), &
                   w => workspace%tmp(:,2), &
                   v => workspace%tmp(:,3:kdim+3), &
                   z => workspace%tmp(:,kdim+4:2*kdim+3) )

            ! Initialize convergence targets from the right-hand side norm.
            norm_sq0 = A%inner_product(b, b)
            tolsq = max(rtol*rtol*norm_sq0, atol*atol)

            ! Form the initial residual and report the starting iterate.
            r = b
            call A%matvec(x, r, alpha=-one_dp, beta=one_dp, op='N')
            norm_sq = A%inner_product(r, r)
            if (associated(workspace%callback)) call workspace%callback(x, norm_sq, 0)

            if (norm_sq <= tolsq .or. maxiter <= 0) then
                deallocate(h, cs, sn, g, y, x_base)
                return
            end if

            iter = 0
            do while (iter < maxiter .and. norm_sq >= tolsq)
                ! Start a new GMRES cycle from the current residual.
                beta = sqrt(max(norm_sq, zero_dp))
                if (beta <= epsilon(one_dp)) exit

                inner_iter = min(kdim, maxiter - iter)
                x_base = x
                h = zero_dp
                cs = zero_dp
                sn = zero_dp
                g = zero_dp
                y = zero_dp

                ! Initialize the Krylov basis and least-squares right-hand side.
                v(:,1) = r / beta
                g(1) = beta

                do j = 1, inner_iter
                    ! Run Arnoldi with the preconditioned basis vector.
                    call M%matvec(v(:,j), z(:,j), alpha=one_dp, beta=zero_dp, op='N')
                    call A%matvec(z(:,j), w, alpha=one_dp, beta=zero_dp, op='N')

                    do i = 1, j
                        h(i,j) = A%inner_product(v(:,i), w)
                        w = w - h(i,j) * v(:,i)
                    end do

                    hnext = sqrt(max(A%inner_product(w, w), zero_dp))
                    h(j+1,j) = hnext
                    if (hnext > epsilon(one_dp)) then
                        v(:,j+1) = w / hnext
                    else
                        v(:,j+1) = zero_dp
                    end if

                    ! Apply the previously accumulated Givens rotations.
                    do i = 1, j - 1
                        temp = cs(i) * h(i,j) + sn(i) * h(i+1,j)
                        h(i+1,j) = -sn(i) * h(i,j) + cs(i) * h(i+1,j)
                        h(i,j) = temp
                    end do

                    ! Build and apply the next Givens rotation.
                    denom = sqrt(h(j,j) * h(j,j) + h(j+1,j) * h(j+1,j))
                    if (denom > epsilon(one_dp)) then
                        cs(j) = h(j,j) / denom
                        sn(j) = h(j+1,j) / denom
                    else
                        cs(j) = one_dp
                        sn(j) = zero_dp
                    end if

                    temp = cs(j) * h(j,j) + sn(j) * h(j+1,j)
                    h(j+1,j) = -sn(j) * h(j,j) + cs(j) * h(j+1,j)
                    h(j,j) = temp

                    temp = cs(j) * g(j) + sn(j) * g(j+1)
                    g(j+1) = -sn(j) * g(j) + cs(j) * g(j+1)
                    g(j) = temp

                    ! Solve the reduced system and update the current iterate.
                    call upper_triangular_solve(h, g, y, j)
                    x = x_base
                    do i = 1, j
                        x = x + y(i) * z(:,i)
                    end do

                    ! Track the residual norm estimate for convergence.
                    norm_sq = g(j+1) * g(j+1)
                    iter = iter + 1
                    if (associated(workspace%callback)) call workspace%callback(x, norm_sq, iter)

                    if (norm_sq < tolsq .or. hnext <= epsilon(one_dp) .or. iter >= maxiter) exit
                end do

                if (norm_sq < tolsq .or. iter >= maxiter) exit

                ! Refresh the residual before the next restarted cycle.
                r = b
                call A%matvec(x, r, alpha=-one_dp, beta=one_dp, op='N')
                norm_sq = A%inner_product(r, r)
            end do
        end associate

        deallocate(h, cs, sn, g, y, x_base)

    contains

        subroutine upper_triangular_solve(h, g, y, n)
            real(dp), intent(in) :: h(:,:), g(:)
            real(dp), intent(inout) :: y(:)
            integer, intent(in) :: n
            integer :: row

            y(1:n) = g(1:n)
            do row = n, 1, -1
                if (row < n) y(row) = y(row) - dot_product(h(row,row+1:n) , y(row+1:n))
                if (abs(h(row,row)) > epsilon(one_dp)) then
                    y(row) = y(row) / h(row,row)
                else
                    y(row) = zero_dp
                end if
            end do
        end subroutine
    end subroutine

    module subroutine stdlib_solve_gmres_dense_sp(A,b,x,di,rtol,atol,maxiter,restart,kdim,precond,M,workspace)
        use stdlib_linalg, only: diag
        real(sp), intent(in) :: A(:,:)
        real(sp), intent(in) :: b(:)
        real(sp), intent(inout) :: x(:)
        real(sp), intent(in), optional :: rtol, atol
        logical(int8), intent(in), optional, target :: di(:)
        integer, intent(in), optional :: maxiter, kdim
        logical, intent(in), optional :: restart
        integer, intent(in), optional :: precond
        class(stdlib_linop_sp_type), optional, intent(in), target :: M
        type(stdlib_solver_workspace_sp_type), optional, intent(inout), target :: workspace
        type(stdlib_linop_sp_type) :: op
        type(stdlib_linop_sp_type), pointer :: M_ => null()
        type(stdlib_solver_workspace_sp_type), pointer :: workspace_
        integer :: kdim_, maxiter_, n, ncols, precond_
        real(sp) :: rtol_, atol_
        logical :: restart_
        logical(int8), pointer :: di_(:)
        real(sp), allocatable :: diagonal(:)

        n = size(b)
        maxiter_ = optval(x=maxiter, default=n)
        kdim_ = max(1, min(optval(x=kdim, default=min(30, n)), n))
        restart_ = optval(x=restart, default=.true.)
        rtol_ = optval(x=rtol, default=1.e-5_sp)
        atol_ = optval(x=atol, default=epsilon(one_sp))
        precond_ = optval(x=precond, default=pc_none)
        ncols = 2 * kdim_ + stdlib_size_wksp_gmres

        if (present(M)) then
            M_ => M
        else
            allocate(M_)
            allocate(diagonal(n), source=zero_sp)

            select case(precond_)
            case(pc_jacobi)
                diagonal = diag(A)
                M_%matvec => precond_jacobi
            case default
                M_%matvec => precond_none
            end select
            where(abs(diagonal) > epsilon(zero_sp)) diagonal = one_sp / diagonal
        end if

        op%matvec => matvec

        if (present(di)) then
            di_ => di
        else
            allocate(di_(n), source=.false._int8)
        end if

        if (present(workspace)) then
            workspace_ => workspace
        else
            allocate(workspace_)
        end if
        if (.not.allocated(workspace_%tmp)) then
            allocate(workspace_%tmp(n, ncols), source=zero_sp)
        else if (size(workspace_%tmp,1) /= n .or. size(workspace_%tmp,2) < ncols) then
            deallocate(workspace_%tmp)
            allocate(workspace_%tmp(n, ncols), source=zero_sp)
        end if

        if (restart_) x = zero_sp
        x = merge(b, x, di_)
        call stdlib_solve_gmres_kernel(op, M_, b, x, rtol_, atol_, maxiter_, kdim_, workspace_)

        if (.not.present(di)) deallocate(di_)
        di_ => null()

        if (.not.present(workspace)) then
            deallocate(workspace_%tmp)
            deallocate(workspace_)
        end if
        M_ => null()
        workspace_ => null()

    contains

        subroutine matvec(x,y,alpha,beta,op)
            use stdlib_linalg_blas, only: gemv
            real(sp), intent(in) :: x(:)
            real(sp), intent(inout) :: y(:)
            real(sp), intent(in) :: alpha
            real(sp), intent(in) :: beta
            character(1), intent(in) :: op
            call gemv(op, m=size(A,1), n=size(A,2), alpha=alpha, a=A, lda=size(A,1), x=x, incx=1, beta=beta, y=y, incy=1)
            y = merge(zero_sp, y, di_)
        end subroutine

        subroutine precond_none(x,y,alpha,beta,op)
            real(sp), intent(in) :: x(:)
            real(sp), intent(inout) :: y(:)
            real(sp), intent(in) :: alpha
            real(sp), intent(in) :: beta
            character(1), intent(in) :: op
            y = merge(zero_sp, x, di_)
        end subroutine

        subroutine precond_jacobi(x,y,alpha,beta,op)
            real(sp), intent(in) :: x(:)
            real(sp), intent(inout) :: y(:)
            real(sp), intent(in) :: alpha
            real(sp), intent(in) :: beta
            character(1), intent(in) :: op
            y = merge(zero_sp, diagonal * x, di_)
        end subroutine
    end subroutine
    module subroutine stdlib_solve_gmres_dense_dp(A,b,x,di,rtol,atol,maxiter,restart,kdim,precond,M,workspace)
        use stdlib_linalg, only: diag
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in), optional :: rtol, atol
        logical(int8), intent(in), optional, target :: di(:)
        integer, intent(in), optional :: maxiter, kdim
        logical, intent(in), optional :: restart
        integer, intent(in), optional :: precond
        class(stdlib_linop_dp_type), optional, intent(in), target :: M
        type(stdlib_solver_workspace_dp_type), optional, intent(inout), target :: workspace
        type(stdlib_linop_dp_type) :: op
        type(stdlib_linop_dp_type), pointer :: M_ => null()
        type(stdlib_solver_workspace_dp_type), pointer :: workspace_
        integer :: kdim_, maxiter_, n, ncols, precond_
        real(dp) :: rtol_, atol_
        logical :: restart_
        logical(int8), pointer :: di_(:)
        real(dp), allocatable :: diagonal(:)

        n = size(b)
        maxiter_ = optval(x=maxiter, default=n)
        kdim_ = max(1, min(optval(x=kdim, default=min(30, n)), n))
        restart_ = optval(x=restart, default=.true.)
        rtol_ = optval(x=rtol, default=1.e-5_dp)
        atol_ = optval(x=atol, default=epsilon(one_dp))
        precond_ = optval(x=precond, default=pc_none)
        ncols = 2 * kdim_ + stdlib_size_wksp_gmres

        if (present(M)) then
            M_ => M
        else
            allocate(M_)
            allocate(diagonal(n), source=zero_dp)

            select case(precond_)
            case(pc_jacobi)
                diagonal = diag(A)
                M_%matvec => precond_jacobi
            case default
                M_%matvec => precond_none
            end select
            where(abs(diagonal) > epsilon(zero_dp)) diagonal = one_dp / diagonal
        end if

        op%matvec => matvec

        if (present(di)) then
            di_ => di
        else
            allocate(di_(n), source=.false._int8)
        end if

        if (present(workspace)) then
            workspace_ => workspace
        else
            allocate(workspace_)
        end if
        if (.not.allocated(workspace_%tmp)) then
            allocate(workspace_%tmp(n, ncols), source=zero_dp)
        else if (size(workspace_%tmp,1) /= n .or. size(workspace_%tmp,2) < ncols) then
            deallocate(workspace_%tmp)
            allocate(workspace_%tmp(n, ncols), source=zero_dp)
        end if

        if (restart_) x = zero_dp
        x = merge(b, x, di_)
        call stdlib_solve_gmres_kernel(op, M_, b, x, rtol_, atol_, maxiter_, kdim_, workspace_)

        if (.not.present(di)) deallocate(di_)
        di_ => null()

        if (.not.present(workspace)) then
            deallocate(workspace_%tmp)
            deallocate(workspace_)
        end if
        M_ => null()
        workspace_ => null()

    contains

        subroutine matvec(x,y,alpha,beta,op)
            use stdlib_linalg_blas, only: gemv
            real(dp), intent(in) :: x(:)
            real(dp), intent(inout) :: y(:)
            real(dp), intent(in) :: alpha
            real(dp), intent(in) :: beta
            character(1), intent(in) :: op
            call gemv(op, m=size(A,1), n=size(A,2), alpha=alpha, a=A, lda=size(A,1), x=x, incx=1, beta=beta, y=y, incy=1)
            y = merge(zero_dp, y, di_)
        end subroutine

        subroutine precond_none(x,y,alpha,beta,op)
            real(dp), intent(in) :: x(:)
            real(dp), intent(inout) :: y(:)
            real(dp), intent(in) :: alpha
            real(dp), intent(in) :: beta
            character(1), intent(in) :: op
            y = merge(zero_dp, x, di_)
        end subroutine

        subroutine precond_jacobi(x,y,alpha,beta,op)
            real(dp), intent(in) :: x(:)
            real(dp), intent(inout) :: y(:)
            real(dp), intent(in) :: alpha
            real(dp), intent(in) :: beta
            character(1), intent(in) :: op
            y = merge(zero_dp, diagonal * x, di_)
        end subroutine
    end subroutine
    module subroutine stdlib_solve_gmres_CSR_sp(A,b,x,di,rtol,atol,maxiter,restart,kdim,precond,M,workspace)
        type(CSR_sp_type), intent(in) :: A
        real(sp), intent(in) :: b(:)
        real(sp), intent(inout) :: x(:)
        real(sp), intent(in), optional :: rtol, atol
        logical(int8), intent(in), optional, target :: di(:)
        integer, intent(in), optional :: maxiter, kdim
        logical, intent(in), optional :: restart
        integer, intent(in), optional :: precond
        class(stdlib_linop_sp_type), optional, intent(in), target :: M
        type(stdlib_solver_workspace_sp_type), optional, intent(inout), target :: workspace
        type(stdlib_linop_sp_type) :: op
        type(stdlib_linop_sp_type), pointer :: M_ => null()
        type(stdlib_solver_workspace_sp_type), pointer :: workspace_
        integer :: kdim_, maxiter_, n, ncols, precond_
        real(sp) :: rtol_, atol_
        logical :: restart_
        logical(int8), pointer :: di_(:)
        real(sp), allocatable :: diagonal(:)

        n = size(b)
        maxiter_ = optval(x=maxiter, default=n)
        kdim_ = max(1, min(optval(x=kdim, default=min(30, n)), n))
        restart_ = optval(x=restart, default=.true.)
        rtol_ = optval(x=rtol, default=1.e-5_sp)
        atol_ = optval(x=atol, default=epsilon(one_sp))
        precond_ = optval(x=precond, default=pc_none)
        ncols = 2 * kdim_ + stdlib_size_wksp_gmres

        if (present(M)) then
            M_ => M
        else
            allocate(M_)
            allocate(diagonal(n), source=zero_sp)

            select case(precond_)
            case(pc_jacobi)
                call diag(A, diagonal)
                M_%matvec => precond_jacobi
            case default
                M_%matvec => precond_none
            end select
            where(abs(diagonal) > epsilon(zero_sp)) diagonal = one_sp / diagonal
        end if

        op%matvec => matvec

        if (present(di)) then
            di_ => di
        else
            allocate(di_(n), source=.false._int8)
        end if

        if (present(workspace)) then
            workspace_ => workspace
        else
            allocate(workspace_)
        end if
        if (.not.allocated(workspace_%tmp)) then
            allocate(workspace_%tmp(n, ncols), source=zero_sp)
        else if (size(workspace_%tmp,1) /= n .or. size(workspace_%tmp,2) < ncols) then
            deallocate(workspace_%tmp)
            allocate(workspace_%tmp(n, ncols), source=zero_sp)
        end if

        if (restart_) x = zero_sp
        x = merge(b, x, di_)
        call stdlib_solve_gmres_kernel(op, M_, b, x, rtol_, atol_, maxiter_, kdim_, workspace_)

        if (.not.present(di)) deallocate(di_)
        di_ => null()

        if (.not.present(workspace)) then
            deallocate(workspace_%tmp)
            deallocate(workspace_)
        end if
        M_ => null()
        workspace_ => null()

    contains

        subroutine matvec(x,y,alpha,beta,op)
            real(sp), intent(in) :: x(:)
            real(sp), intent(inout) :: y(:)
            real(sp), intent(in) :: alpha
            real(sp), intent(in) :: beta
            character(1), intent(in) :: op
            call spmv(A, x, y, alpha, beta, op)
            y = merge(zero_sp, y, di_)
        end subroutine

        subroutine precond_none(x,y,alpha,beta,op)
            real(sp), intent(in) :: x(:)
            real(sp), intent(inout) :: y(:)
            real(sp), intent(in) :: alpha
            real(sp), intent(in) :: beta
            character(1), intent(in) :: op
            y = merge(zero_sp, x, di_)
        end subroutine

        subroutine precond_jacobi(x,y,alpha,beta,op)
            real(sp), intent(in) :: x(:)
            real(sp), intent(inout) :: y(:)
            real(sp), intent(in) :: alpha
            real(sp), intent(in) :: beta
            character(1), intent(in) :: op
            y = merge(zero_sp, diagonal * x, di_)
        end subroutine
    end subroutine
    module subroutine stdlib_solve_gmres_CSR_dp(A,b,x,di,rtol,atol,maxiter,restart,kdim,precond,M,workspace)
        type(CSR_dp_type), intent(in) :: A
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in), optional :: rtol, atol
        logical(int8), intent(in), optional, target :: di(:)
        integer, intent(in), optional :: maxiter, kdim
        logical, intent(in), optional :: restart
        integer, intent(in), optional :: precond
        class(stdlib_linop_dp_type), optional, intent(in), target :: M
        type(stdlib_solver_workspace_dp_type), optional, intent(inout), target :: workspace
        type(stdlib_linop_dp_type) :: op
        type(stdlib_linop_dp_type), pointer :: M_ => null()
        type(stdlib_solver_workspace_dp_type), pointer :: workspace_
        integer :: kdim_, maxiter_, n, ncols, precond_
        real(dp) :: rtol_, atol_
        logical :: restart_
        logical(int8), pointer :: di_(:)
        real(dp), allocatable :: diagonal(:)

        n = size(b)
        maxiter_ = optval(x=maxiter, default=n)
        kdim_ = max(1, min(optval(x=kdim, default=min(30, n)), n))
        restart_ = optval(x=restart, default=.true.)
        rtol_ = optval(x=rtol, default=1.e-5_dp)
        atol_ = optval(x=atol, default=epsilon(one_dp))
        precond_ = optval(x=precond, default=pc_none)
        ncols = 2 * kdim_ + stdlib_size_wksp_gmres

        if (present(M)) then
            M_ => M
        else
            allocate(M_)
            allocate(diagonal(n), source=zero_dp)

            select case(precond_)
            case(pc_jacobi)
                call diag(A, diagonal)
                M_%matvec => precond_jacobi
            case default
                M_%matvec => precond_none
            end select
            where(abs(diagonal) > epsilon(zero_dp)) diagonal = one_dp / diagonal
        end if

        op%matvec => matvec

        if (present(di)) then
            di_ => di
        else
            allocate(di_(n), source=.false._int8)
        end if

        if (present(workspace)) then
            workspace_ => workspace
        else
            allocate(workspace_)
        end if
        if (.not.allocated(workspace_%tmp)) then
            allocate(workspace_%tmp(n, ncols), source=zero_dp)
        else if (size(workspace_%tmp,1) /= n .or. size(workspace_%tmp,2) < ncols) then
            deallocate(workspace_%tmp)
            allocate(workspace_%tmp(n, ncols), source=zero_dp)
        end if

        if (restart_) x = zero_dp
        x = merge(b, x, di_)
        call stdlib_solve_gmres_kernel(op, M_, b, x, rtol_, atol_, maxiter_, kdim_, workspace_)

        if (.not.present(di)) deallocate(di_)
        di_ => null()

        if (.not.present(workspace)) then
            deallocate(workspace_%tmp)
            deallocate(workspace_)
        end if
        M_ => null()
        workspace_ => null()

    contains

        subroutine matvec(x,y,alpha,beta,op)
            real(dp), intent(in) :: x(:)
            real(dp), intent(inout) :: y(:)
            real(dp), intent(in) :: alpha
            real(dp), intent(in) :: beta
            character(1), intent(in) :: op
            call spmv(A, x, y, alpha, beta, op)
            y = merge(zero_dp, y, di_)
        end subroutine

        subroutine precond_none(x,y,alpha,beta,op)
            real(dp), intent(in) :: x(:)
            real(dp), intent(inout) :: y(:)
            real(dp), intent(in) :: alpha
            real(dp), intent(in) :: beta
            character(1), intent(in) :: op
            y = merge(zero_dp, x, di_)
        end subroutine

        subroutine precond_jacobi(x,y,alpha,beta,op)
            real(dp), intent(in) :: x(:)
            real(dp), intent(inout) :: y(:)
            real(dp), intent(in) :: alpha
            real(dp), intent(in) :: beta
            character(1), intent(in) :: op
            y = merge(zero_dp, diagonal * x, di_)
        end subroutine
    end subroutine

end submodule stdlib_linalg_iterative_gmres