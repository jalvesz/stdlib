!! The `stdlib_linalg_iterative_solvers` module provides interfaces for iterative solvers.
!!
module stdlib_linalg_iterative_solvers
    use stdlib_kinds
    use stdlib_sparse
    implicit none
    private 

    !! workspace sizes: defined by the number of vectors used by the iterative solver.
    enum, bind(c)
        enumerator :: stdlib_size_wksp_cg = 3
        enumerator :: stdlib_size_wksp_pcg = 4
    end enum
    public :: stdlib_size_wksp_cg, stdlib_size_wksp_pcg

    !! version: experimental
    !!
    !! linop type holding the linear operator and its associated methods.
    !! The `linop` type is used to define the linear operator for the iterative solvers.
    type, public :: linop_sp_type
        procedure(vector_sub_sp), nopass, pointer    :: matvec => null()
        procedure(reduction_sub_sp), nopass, pointer :: inner_product => default_dot_sp
    end type
    type, public :: linop_dp_type
        procedure(vector_sub_dp), nopass, pointer    :: matvec => null()
        procedure(reduction_sub_dp), nopass, pointer :: inner_product => default_dot_dp
    end type

    !! version: experimental
    !!
    !! solver_workspace type holding temporal array data for the iterative solvers.
    type, public :: solver_workspace_sp_type
        real(sp), allocatable :: tmp(:,:)
        procedure(logger_sub_sp), pointer, nopass :: callback => null()
    end type 

    type, public :: solver_workspace_dp_type
        real(dp), allocatable :: tmp(:,:)
        procedure(logger_sub_dp), pointer, nopass :: callback => null()
    end type 


    abstract interface
        subroutine vector_sub_sp(x,y,alpha,beta,op)
            import :: sp
            real(sp), intent(in)  :: x(:)
            real(sp), intent(inout) :: y(:)
            real(sp), intent(in) :: alpha
            real(sp), intent(in) :: beta
            character(1), intent(in) :: op
        end subroutine
        pure real(sp) function reduction_sub_sp(x,y) result(r)
            import :: sp
            real(sp), intent(in) :: x(:)
            real(sp), intent(in) :: y(:)
        end function
        subroutine logger_sub_sp(x,norm_sq,iter)
            import :: sp
            real(sp), intent(in) :: x(:)
            real(sp), intent(in) :: norm_sq
            integer, intent(in) :: iter
        end subroutine
        subroutine vector_sub_dp(x,y,alpha,beta,op)
            import :: dp
            real(dp), intent(in)  :: x(:)
            real(dp), intent(inout) :: y(:)
            real(dp), intent(in) :: alpha
            real(dp), intent(in) :: beta
            character(1), intent(in) :: op
        end subroutine
        pure real(dp) function reduction_sub_dp(x,y) result(r)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(in) :: y(:)
        end function
        subroutine logger_sub_dp(x,norm_sq,iter)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(in) :: norm_sq
            integer, intent(in) :: iter
        end subroutine
    end interface

    !! version: experimental
    !!
    !! solve_cg_kernel interface for the conjugate gradient method.
    !! [Specifications](../page/specs/stdlib_linalg_iterative_solvers.html#solve_cg_kernel)
    interface solve_cg_kernel
        module subroutine solve_cg_kernel_sp(A,b,x,tol,maxiter,workspace)
            class(linop_sp_type), intent(in) :: A !! linear operator
            real(sp), intent(in) :: b(:) !! right-hand side vector
            real(sp), intent(inout) :: x(:) !! solution vector and initial guess
            real(sp), intent(in) :: tol !! tolerance for convergence
            integer, intent(in) :: maxiter !! maximum number of iterations
            type(solver_workspace_sp_type), intent(inout) :: workspace !! workspace for the solver
        end subroutine
        module subroutine solve_cg_kernel_dp(A,b,x,tol,maxiter,workspace)
            class(linop_dp_type), intent(in) :: A !! linear operator
            real(dp), intent(in) :: b(:) !! right-hand side vector
            real(dp), intent(inout) :: x(:) !! solution vector and initial guess
            real(dp), intent(in) :: tol !! tolerance for convergence
            integer, intent(in) :: maxiter !! maximum number of iterations
            type(solver_workspace_dp_type), intent(inout) :: workspace !! workspace for the solver
        end subroutine
    end interface
    public :: solve_cg_kernel

    interface solve_cg
        module subroutine solve_cg_dense_sp(A,b,x,di,tol,maxiter,restart,workspace)
            !! linear operator matrix
            real(sp), intent(in) :: A(:,:) 
            real(sp), intent(in) :: b(:) !! right-hand side vector
            real(sp), intent(inout) :: x(:) !! solution vector and initial guess
            real(sp), intent(in), optional :: tol !! tolerance for convergence
            logical(1), intent(in), optional, target  :: di(:) !! dirichlet conditions mask
            integer, intent(in), optional :: maxiter !! maximum number of iterations
            logical, intent(in), optional :: restart !! restart flag
            type(solver_workspace_sp_type), optional, intent(inout), target :: workspace !! workspace for the solver
        end subroutine
        module subroutine solve_cg_dense_dp(A,b,x,di,tol,maxiter,restart,workspace)
            !! linear operator matrix
            real(dp), intent(in) :: A(:,:) 
            real(dp), intent(in) :: b(:) !! right-hand side vector
            real(dp), intent(inout) :: x(:) !! solution vector and initial guess
            real(dp), intent(in), optional :: tol !! tolerance for convergence
            logical(1), intent(in), optional, target  :: di(:) !! dirichlet conditions mask
            integer, intent(in), optional :: maxiter !! maximum number of iterations
            logical, intent(in), optional :: restart !! restart flag
            type(solver_workspace_dp_type), optional, intent(inout), target :: workspace !! workspace for the solver
        end subroutine
        module subroutine solve_cg_CSR_sp(A,b,x,di,tol,maxiter,restart,workspace)
            !! linear operator matrix
            type(CSR_sp_type), intent(in) :: A
            real(sp), intent(in) :: b(:) !! right-hand side vector
            real(sp), intent(inout) :: x(:) !! solution vector and initial guess
            real(sp), intent(in), optional :: tol !! tolerance for convergence
            logical(1), intent(in), optional, target  :: di(:) !! dirichlet conditions mask
            integer, intent(in), optional :: maxiter !! maximum number of iterations
            logical, intent(in), optional :: restart !! restart flag
            type(solver_workspace_sp_type), optional, intent(inout), target :: workspace !! workspace for the solver
        end subroutine
        module subroutine solve_cg_CSR_dp(A,b,x,di,tol,maxiter,restart,workspace)
            !! linear operator matrix
            type(CSR_dp_type), intent(in) :: A
            real(dp), intent(in) :: b(:) !! right-hand side vector
            real(dp), intent(inout) :: x(:) !! solution vector and initial guess
            real(dp), intent(in), optional :: tol !! tolerance for convergence
            logical(1), intent(in), optional, target  :: di(:) !! dirichlet conditions mask
            integer, intent(in), optional :: maxiter !! maximum number of iterations
            logical, intent(in), optional :: restart !! restart flag
            type(solver_workspace_dp_type), optional, intent(inout), target :: workspace !! workspace for the solver
        end subroutine
    end interface
    public :: solve_cg

    !! version: experimental
    !!
    !! solve_pcg_kernel interface for the preconditionned conjugate gradient method.
    !! [Specifications](../page/specs/stdlib_linalg_iterative_solvers.html#solve_pcg_kernel)
    interface solve_pcg_kernel
        module subroutine solve_pcg_kernel_sp(A,M,b,x,tol,maxiter,workspace)
            class(linop_sp_type), intent(in) :: A !! linear operator
            class(linop_sp_type), intent(in) :: M !! preconditioner linear operator
            real(sp), intent(in) :: b(:) !! right-hand side vector
            real(sp), intent(inout) :: x(:) !! solution vector and initial guess
            real(sp), intent(in) :: tol !! tolerance for convergence
            integer, intent(in) :: maxiter !! maximum number of iterations
            type(solver_workspace_sp_type), intent(inout) :: workspace !! workspace for the solver
        end subroutine
        module subroutine solve_pcg_kernel_dp(A,M,b,x,tol,maxiter,workspace)
            class(linop_dp_type), intent(in) :: A !! linear operator
            class(linop_dp_type), intent(in) :: M !! preconditioner linear operator
            real(dp), intent(in) :: b(:) !! right-hand side vector
            real(dp), intent(inout) :: x(:) !! solution vector and initial guess
            real(dp), intent(in) :: tol !! tolerance for convergence
            integer, intent(in) :: maxiter !! maximum number of iterations
            type(solver_workspace_dp_type), intent(inout) :: workspace !! workspace for the solver
        end subroutine
    end interface
    public :: solve_pcg_kernel

    interface solve_pcg
        module subroutine solve_pcg_dense_sp(A,b,x,di,tol,maxiter,restart,precond,M,workspace)
            !! linear operator matrix
            real(sp), intent(in) :: A(:,:)
            real(sp), intent(in) :: b(:) !! right-hand side vector
            real(sp), intent(inout) :: x(:) !! solution vector and initial guess
            real(sp), intent(in), optional :: tol !! tolerance for convergence
            logical(1), intent(in), optional, target  :: di(:) !! dirichlet conditions mask
            integer, intent(in), optional  :: maxiter !! maximum number of iterations
            logical, intent(in), optional :: restart !! restart flag
            integer, intent(in), optional  :: precond !! preconditioner method enumerator
            class(linop_sp_type), optional , intent(in), target :: M !! preconditioner linear operator
            type(solver_workspace_sp_type), optional, intent(inout), target :: workspace !! workspace for the solver
        end subroutine
        module subroutine solve_pcg_dense_dp(A,b,x,di,tol,maxiter,restart,precond,M,workspace)
            !! linear operator matrix
            real(dp), intent(in) :: A(:,:)
            real(dp), intent(in) :: b(:) !! right-hand side vector
            real(dp), intent(inout) :: x(:) !! solution vector and initial guess
            real(dp), intent(in), optional :: tol !! tolerance for convergence
            logical(1), intent(in), optional, target  :: di(:) !! dirichlet conditions mask
            integer, intent(in), optional  :: maxiter !! maximum number of iterations
            logical, intent(in), optional :: restart !! restart flag
            integer, intent(in), optional  :: precond !! preconditioner method enumerator
            class(linop_dp_type), optional , intent(in), target :: M !! preconditioner linear operator
            type(solver_workspace_dp_type), optional, intent(inout), target :: workspace !! workspace for the solver
        end subroutine
        module subroutine solve_pcg_CSR_sp(A,b,x,di,tol,maxiter,restart,precond,M,workspace)
            !! linear operator matrix
            type(CSR_sp_type), intent(in) :: A
            real(sp), intent(in) :: b(:) !! right-hand side vector
            real(sp), intent(inout) :: x(:) !! solution vector and initial guess
            real(sp), intent(in), optional :: tol !! tolerance for convergence
            logical(1), intent(in), optional, target  :: di(:) !! dirichlet conditions mask
            integer, intent(in), optional  :: maxiter !! maximum number of iterations
            logical, intent(in), optional :: restart !! restart flag
            integer, intent(in), optional  :: precond !! preconditioner method enumerator
            class(linop_sp_type), optional , intent(in), target :: M !! preconditioner linear operator
            type(solver_workspace_sp_type), optional, intent(inout), target :: workspace !! workspace for the solver
        end subroutine
        module subroutine solve_pcg_CSR_dp(A,b,x,di,tol,maxiter,restart,precond,M,workspace)
            !! linear operator matrix
            type(CSR_dp_type), intent(in) :: A
            real(dp), intent(in) :: b(:) !! right-hand side vector
            real(dp), intent(inout) :: x(:) !! solution vector and initial guess
            real(dp), intent(in), optional :: tol !! tolerance for convergence
            logical(1), intent(in), optional, target  :: di(:) !! dirichlet conditions mask
            integer, intent(in), optional  :: maxiter !! maximum number of iterations
            logical, intent(in), optional :: restart !! restart flag
            integer, intent(in), optional  :: precond !! preconditioner method enumerator
            class(linop_dp_type), optional , intent(in), target :: M !! preconditioner linear operator
            type(solver_workspace_dp_type), optional, intent(inout), target :: workspace !! workspace for the solver
        end subroutine
    end interface
    public :: solve_pcg 

contains

    !------------------------------------------------------------------
    ! defaults
    !------------------------------------------------------------------
    pure real(sp) function default_dot_sp(x,y) result(r)
        use stdlib_intrinsics, only: stdlib_dot_product
        real(sp), intent(in) :: x(:)
        real(sp), intent(in) :: y(:)
        r = stdlib_dot_product(x,y)
    end function

    pure real(dp) function default_dot_dp(x,y) result(r)
        use stdlib_intrinsics, only: stdlib_dot_product
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: y(:)
        r = stdlib_dot_product(x,y)
    end function

    
end module stdlib_linalg_iterative_solvers
