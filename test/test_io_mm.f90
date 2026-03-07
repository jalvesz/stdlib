module test_io_mm
    use testdrive, only : new_unittest, unittest_type, error_type, check, skip_test
    use stdlib_kinds
    use stdlib_math, only: all_close
    use stdlib_io_mm
    implicit none

    integer, parameter :: MS_general = 1
    integer, parameter :: MS_symmetric = 2
    integer, parameter :: MS_skew_symmetric = 3
    integer, parameter :: MS_hermitian = 4

contains


    !> Collect all exported unit tests
    subroutine collect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest('io_mm_array', test_io_mm_array), &
            new_unittest('io_mm_coordinate', test_io_mm_coordinate) &
        ]
    end subroutine

    subroutine generate_random_sp_dense_matrix(A)
        real(sp), intent(out) :: A(:, :)

        ! Internal variables
        integer :: n

        n = size(A,dim=1)
        
        call random_number(A)
    end subroutine
    pure function compare_dense_sp(A, B) result(result)
        real(sp), intent(in) :: A(:, :), B(:, :)
        logical :: result

        result = all_close(A, B)
    end function
    subroutine generate_random_dp_dense_matrix(A)
        real(dp), intent(out) :: A(:, :)

        ! Internal variables
        integer :: n

        n = size(A,dim=1)
        
        call random_number(A)
    end subroutine
    pure function compare_dense_dp(A, B) result(result)
        real(dp), intent(in) :: A(:, :), B(:, :)
        logical :: result

        result = all_close(A, B)
    end function
    subroutine generate_random_csp_dense_matrix(A)
        complex(sp), intent(out) :: A(:, :)

        ! Internal variables
        real(sp), allocatable :: R(:, :, :)
        integer :: n

        n = size(A,dim=1)
        
        allocate(R(n,n,2))
        call random_number(R(:,:,1))
        call random_number(R(:,:,2))
        A = cmplx(R(:,:,1), R(:,:,2), kind=sp)
    end subroutine
    pure function compare_dense_csp(A, B) result(result)
        complex(sp), intent(in) :: A(:, :), B(:, :)
        logical :: result

        result = all_close(A, B)
    end function
    subroutine generate_random_cdp_dense_matrix(A)
        complex(dp), intent(out) :: A(:, :)

        ! Internal variables
        real(dp), allocatable :: R(:, :, :)
        integer :: n

        n = size(A,dim=1)
        
        allocate(R(n,n,2))
        call random_number(R(:,:,1))
        call random_number(R(:,:,2))
        A = cmplx(R(:,:,1), R(:,:,2), kind=dp)
    end subroutine
    pure function compare_dense_cdp(A, B) result(result)
        complex(dp), intent(in) :: A(:, :), B(:, :)
        logical :: result

        result = all_close(A, B)
    end function
    subroutine generate_random_int8_dense_matrix(A)
        integer(int8), intent(out) :: A(:, :)

        ! Internal variables
        real :: rnd
        integer :: i, j
        integer :: n

        n = size(A,dim=1)
        
        do j = 1, n
            do i = 1,n
                call random_number(rnd)
                A(i,j) = int(rnd * 100 - 50, kind=int8)
            end do
        end do
    end subroutine
    pure function compare_dense_int8(A, B) result(result)
        integer(int8), intent(in) :: A(:, :), B(:, :)
        logical :: result

        result = all(A==B)
    end function
    subroutine generate_random_int16_dense_matrix(A)
        integer(int16), intent(out) :: A(:, :)

        ! Internal variables
        real :: rnd
        integer :: i, j
        integer :: n

        n = size(A,dim=1)
        
        do j = 1, n
            do i = 1,n
                call random_number(rnd)
                A(i,j) = int(rnd * 100 - 50, kind=int16)
            end do
        end do
    end subroutine
    pure function compare_dense_int16(A, B) result(result)
        integer(int16), intent(in) :: A(:, :), B(:, :)
        logical :: result

        result = all(A==B)
    end function
    subroutine generate_random_int32_dense_matrix(A)
        integer(int32), intent(out) :: A(:, :)

        ! Internal variables
        real :: rnd
        integer :: i, j
        integer :: n

        n = size(A,dim=1)
        
        do j = 1, n
            do i = 1,n
                call random_number(rnd)
                A(i,j) = int(rnd * 100 - 50, kind=int32)
            end do
        end do
    end subroutine
    pure function compare_dense_int32(A, B) result(result)
        integer(int32), intent(in) :: A(:, :), B(:, :)
        logical :: result

        result = all(A==B)
    end function
    subroutine generate_random_int64_dense_matrix(A)
        integer(int64), intent(out) :: A(:, :)

        ! Internal variables
        real :: rnd
        integer :: i, j
        integer :: n

        n = size(A,dim=1)
        
        do j = 1, n
            do i = 1,n
                call random_number(rnd)
                A(i,j) = int(rnd * 100 - 50, kind=int64)
            end do
        end do
    end subroutine
    pure function compare_dense_int64(A, B) result(result)
        integer(int64), intent(in) :: A(:, :), B(:, :)
        logical :: result

        result = all(A==B)
    end function

    subroutine test_io_mm_array(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: n = 5
            real(sp), allocatable :: matrix_save(:, :), matrix_load(:, :), A(:, :)
            logical :: result
            allocate(matrix_save(n,n))
            allocate(A(n,n))

            ! General matrix
            call generate_random_sp_dense_matrix(matrix_save)
            call save_mm("test_mmio_dense.mtx", matrix_save, format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_sp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=unspecified, type=real(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_sp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=auto, type=real(sp)")
            if(allocated(error)) return

            ! Symmetric matrix
            call generate_random_sp_dense_matrix(A)
            ! Construct symmetric matrix using (A + A.T)
            matrix_save = A + transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_sp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=symmetric, type=real(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_sp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=auto, type=real(sp)")
            if(allocated(error)) return

            ! Skew-symmetric matrix
            call generate_random_sp_dense_matrix(A)
            ! Construct symmetric matrix using (A - A.T)
            matrix_save = A - transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_sp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=skew-symmetric, type=real(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_sp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=auto, type=real(sp)")
            if(allocated(error)) return

        end block
        block
            integer, parameter :: n = 5
            real(dp), allocatable :: matrix_save(:, :), matrix_load(:, :), A(:, :)
            logical :: result
            allocate(matrix_save(n,n))
            allocate(A(n,n))

            ! General matrix
            call generate_random_dp_dense_matrix(matrix_save)
            call save_mm("test_mmio_dense.mtx", matrix_save, format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_dp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=unspecified, type=real(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_dp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=auto, type=real(dp)")
            if(allocated(error)) return

            ! Symmetric matrix
            call generate_random_dp_dense_matrix(A)
            ! Construct symmetric matrix using (A + A.T)
            matrix_save = A + transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_dp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=symmetric, type=real(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_dp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=auto, type=real(dp)")
            if(allocated(error)) return

            ! Skew-symmetric matrix
            call generate_random_dp_dense_matrix(A)
            ! Construct symmetric matrix using (A - A.T)
            matrix_save = A - transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_dp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=skew-symmetric, type=real(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_dp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=auto, type=real(dp)")
            if(allocated(error)) return

        end block
        block
            integer, parameter :: n = 5
            complex(sp), allocatable :: matrix_save(:, :), matrix_load(:, :), A(:, :)
            logical :: result
            allocate(matrix_save(n,n))
            allocate(A(n,n))

            ! General matrix
            call generate_random_csp_dense_matrix(matrix_save)
            call save_mm("test_mmio_dense.mtx", matrix_save, format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_csp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=unspecified, type=complex(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_csp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=auto, type=complex(sp)")
            if(allocated(error)) return

            ! Symmetric matrix
            call generate_random_csp_dense_matrix(A)
            ! Construct symmetric matrix using (A + A.T)
            matrix_save = A + transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_csp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=symmetric, type=complex(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_csp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=auto, type=complex(sp)")
            if(allocated(error)) return

            ! Skew-symmetric matrix
            call generate_random_csp_dense_matrix(A)
            ! Construct symmetric matrix using (A - A.T)
            matrix_save = A - transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_csp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=skew-symmetric, type=complex(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_csp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=auto, type=complex(sp)")
            if(allocated(error)) return

            ! Hermitian matrix
            call generate_random_csp_dense_matrix(A)
            ! Construct symmetric matrix using (A + A.H)
            matrix_save = A + transpose(conjg(A))
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "hermitian", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_csp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=hermitian, symmetry_arg=hermitian, type=complex(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_csp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=hermitian, symmetry_arg=auto, type=complex(sp)")
            if(allocated(error)) return
        end block
        block
            integer, parameter :: n = 5
            complex(dp), allocatable :: matrix_save(:, :), matrix_load(:, :), A(:, :)
            logical :: result
            allocate(matrix_save(n,n))
            allocate(A(n,n))

            ! General matrix
            call generate_random_cdp_dense_matrix(matrix_save)
            call save_mm("test_mmio_dense.mtx", matrix_save, format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_cdp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=unspecified, type=complex(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_cdp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=auto, type=complex(dp)")
            if(allocated(error)) return

            ! Symmetric matrix
            call generate_random_cdp_dense_matrix(A)
            ! Construct symmetric matrix using (A + A.T)
            matrix_save = A + transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_cdp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=symmetric, type=complex(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_cdp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=auto, type=complex(dp)")
            if(allocated(error)) return

            ! Skew-symmetric matrix
            call generate_random_cdp_dense_matrix(A)
            ! Construct symmetric matrix using (A - A.T)
            matrix_save = A - transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_cdp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=skew-symmetric, type=complex(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_cdp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=auto, type=complex(dp)")
            if(allocated(error)) return

            ! Hermitian matrix
            call generate_random_cdp_dense_matrix(A)
            ! Construct symmetric matrix using (A + A.H)
            matrix_save = A + transpose(conjg(A))
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "hermitian", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_cdp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=hermitian, symmetry_arg=hermitian, type=complex(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_cdp(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=hermitian, symmetry_arg=auto, type=complex(dp)")
            if(allocated(error)) return
        end block
        block
            integer, parameter :: n = 5
            integer(int8), allocatable :: matrix_save(:, :), matrix_load(:, :), A(:, :)
            logical :: result
            allocate(matrix_save(n,n))
            allocate(A(n,n))

            ! General matrix
            call generate_random_int8_dense_matrix(matrix_save)
            call save_mm("test_mmio_dense.mtx", matrix_save, format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int8(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=unspecified, type=integer(int8)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int8(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=auto, type=integer(int8)")
            if(allocated(error)) return

            ! Symmetric matrix
            call generate_random_int8_dense_matrix(A)
            ! Construct symmetric matrix using (A + A.T)
            matrix_save = A + transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int8(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=symmetric, type=integer(int8)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int8(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=auto, type=integer(int8)")
            if(allocated(error)) return

            ! Skew-symmetric matrix
            call generate_random_int8_dense_matrix(A)
            ! Construct symmetric matrix using (A - A.T)
            matrix_save = A - transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int8(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=skew-symmetric, type=integer(int8)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int8(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=auto, type=integer(int8)")
            if(allocated(error)) return

        end block
        block
            integer, parameter :: n = 5
            integer(int16), allocatable :: matrix_save(:, :), matrix_load(:, :), A(:, :)
            logical :: result
            allocate(matrix_save(n,n))
            allocate(A(n,n))

            ! General matrix
            call generate_random_int16_dense_matrix(matrix_save)
            call save_mm("test_mmio_dense.mtx", matrix_save, format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int16(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=unspecified, type=integer(int16)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int16(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=auto, type=integer(int16)")
            if(allocated(error)) return

            ! Symmetric matrix
            call generate_random_int16_dense_matrix(A)
            ! Construct symmetric matrix using (A + A.T)
            matrix_save = A + transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int16(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=symmetric, type=integer(int16)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int16(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=auto, type=integer(int16)")
            if(allocated(error)) return

            ! Skew-symmetric matrix
            call generate_random_int16_dense_matrix(A)
            ! Construct symmetric matrix using (A - A.T)
            matrix_save = A - transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int16(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=skew-symmetric, type=integer(int16)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int16(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=auto, type=integer(int16)")
            if(allocated(error)) return

        end block
        block
            integer, parameter :: n = 5
            integer(int32), allocatable :: matrix_save(:, :), matrix_load(:, :), A(:, :)
            logical :: result
            allocate(matrix_save(n,n))
            allocate(A(n,n))

            ! General matrix
            call generate_random_int32_dense_matrix(matrix_save)
            call save_mm("test_mmio_dense.mtx", matrix_save, format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int32(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=unspecified, type=integer(int32)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int32(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=auto, type=integer(int32)")
            if(allocated(error)) return

            ! Symmetric matrix
            call generate_random_int32_dense_matrix(A)
            ! Construct symmetric matrix using (A + A.T)
            matrix_save = A + transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int32(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=symmetric, type=integer(int32)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int32(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=auto, type=integer(int32)")
            if(allocated(error)) return

            ! Skew-symmetric matrix
            call generate_random_int32_dense_matrix(A)
            ! Construct symmetric matrix using (A - A.T)
            matrix_save = A - transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int32(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=skew-symmetric, type=integer(int32)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int32(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=auto, type=integer(int32)")
            if(allocated(error)) return

        end block
        block
            integer, parameter :: n = 5
            integer(int64), allocatable :: matrix_save(:, :), matrix_load(:, :), A(:, :)
            logical :: result
            allocate(matrix_save(n,n))
            allocate(A(n,n))

            ! General matrix
            call generate_random_int64_dense_matrix(matrix_save)
            call save_mm("test_mmio_dense.mtx", matrix_save, format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int64(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=unspecified, type=integer(int64)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int64(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=auto, type=integer(int64)")
            if(allocated(error)) return

            ! Symmetric matrix
            call generate_random_int64_dense_matrix(A)
            ! Construct symmetric matrix using (A + A.T)
            matrix_save = A + transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int64(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=symmetric, type=integer(int64)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int64(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=auto, type=integer(int64)")
            if(allocated(error)) return

            ! Skew-symmetric matrix
            call generate_random_int64_dense_matrix(A)
            ! Construct symmetric matrix using (A - A.T)
            matrix_save = A - transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int64(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=skew-symmetric, type=integer(int64)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = compare_dense_int64(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=auto, type=integer(int64)")
            if(allocated(error)) return

        end block
    end subroutine

    subroutine generate_random_positions(pos, t_entries)
        integer, intent(inout) :: pos(:)
        integer, intent(in) :: t_entries

        ! Internal variables
        integer :: i, j, temp
        real(dp) :: rnd

        do i = 1, t_entries
            pos(i) = i
        end do
        do i = t_entries, 1, -1
            call random_number(rnd)
            j = ceiling(t_entries*rnd)
            temp = pos(i)
            pos(i) = pos(j)
            pos(j) = temp
        end do
    end subroutine

    subroutine fill_first_half_indices(index_save, pos, nnz, nrows, ncols, symmetry, j)
        integer, intent(out) :: index_save(:, :)
        integer, intent(in) :: pos(:), nnz, nrows, ncols, symmetry
        integer, intent(out) :: j

        ! Internal variables
        integer :: i, row, col

        if(symmetry == MS_symmetric .or. symmetry == MS_hermitian) then
            j = 1
            do i = 1, nnz
                row = mod(pos(i) - 1,nrows) + 1
                col = (pos(i) - 1)/ncols + 1
                if(row < col) cycle
                index_save(1,j) = row
                index_save(2,j) = col
                j = j + 1
            end do
        else
            j = 1
            do i = 1, nnz
                row = mod(pos(i) - 1,nrows) + 1
                col = (pos(i) - 1)/ncols + 1
                if(row <= col) cycle
                index_save(1,j) = row
                index_save(2,j) = col
                j = j + 1
            end do
        end if
    end subroutine

    subroutine generate_random_data_for_sp_coo(A, nnz_to_write)
        real(sp), intent(out) :: A(:)
        integer, intent(in) :: nnz_to_write

        ! Internal variables

        call random_number(A(1:nnz_to_write))
    end subroutine

    pure function compare_coo_sp(index_save, index_load, data_save, data_load) result(result)
        real(sp), intent(in) :: data_save(:), data_load(:)
        integer, intent(in) :: index_save(:, :), index_load(:,:)
        logical :: result
        
        result = all(index_save == index_load) .and. all_close(data_save, data_load)
    end function

    subroutine fill_other_half_sp(index_save, data_save, j, half_nnz, symmetry)
        integer, intent(out) :: index_save(:, :)
        real(sp), intent(out) :: data_save(:)
        integer, intent(in) :: half_nnz, symmetry
        integer, intent(inout) :: j

        ! Internal variables.
        integer :: i

        if(symmetry == MS_symmetric) then
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = data_save(i)
                j=j+1
            end do
        else
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = -data_save(i)
                j=j+1
            end do
        end if
    end subroutine
    
    subroutine generate_random_sp_coo_matrix(index_save, data_save, nrows, ncols, symmetry)
        real(sp), allocatable, intent(out) :: data_save(:)
        integer, allocatable, intent(out) :: index_save(:, :)
        integer, intent(in) :: nrows, ncols, symmetry

        ! Internal variables
        integer, allocatable :: pos(:)
        integer :: nnz, nnz_lower, nnz_diag, i, j
        real(sp) :: density

        allocate(pos(nrows * ncols))
        call generate_random_positions(pos, nrows * ncols)
        call random_number(density)
        nnz = ceiling(density*nrows*ncols)
        
        if(symmetry == MS_general) then
            allocate(index_save(2, nnz))
            allocate(data_save(nnz))
            do i = 1, nnz
                index_save(1,i) = mod(pos(i) - 1,ncols) + 1
                index_save(2,i) = (pos(i) - 1)/nrows + 1
            end do
            call generate_random_data_for_sp_coo(data_save, nnz)
        else if(symmetry == MS_symmetric) then
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            nnz_diag = count(mod(pos(1:nnz) - 1,nrows) == (pos(1:nnz) - 1)/ncols) !! diagonal
            allocate(index_save(2, 2*nnz_lower + nnz_diag))
            allocate(data_save(2*nnz_lower + nnz_diag))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_symmetric, j)
            call generate_random_data_for_sp_coo(data_save, nnz_lower + nnz_diag)
            call fill_other_half_sp(index_save, data_save, j, nnz_lower+nnz_diag, MS_symmetric)
        else
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            allocate(index_save(2, 2*nnz_lower))
            allocate(data_save(2*nnz_lower))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_skew_symmetric, j)
            call generate_random_data_for_sp_coo(data_save, nnz_lower)
            call fill_other_half_sp(index_save, data_save, j, nnz_lower, MS_skew_symmetric)
        end if
    end subroutine
    subroutine generate_random_data_for_dp_coo(A, nnz_to_write)
        real(dp), intent(out) :: A(:)
        integer, intent(in) :: nnz_to_write

        ! Internal variables

        call random_number(A(1:nnz_to_write))
    end subroutine

    pure function compare_coo_dp(index_save, index_load, data_save, data_load) result(result)
        real(dp), intent(in) :: data_save(:), data_load(:)
        integer, intent(in) :: index_save(:, :), index_load(:,:)
        logical :: result
        
        result = all(index_save == index_load) .and. all_close(data_save, data_load)
    end function

    subroutine fill_other_half_dp(index_save, data_save, j, half_nnz, symmetry)
        integer, intent(out) :: index_save(:, :)
        real(dp), intent(out) :: data_save(:)
        integer, intent(in) :: half_nnz, symmetry
        integer, intent(inout) :: j

        ! Internal variables.
        integer :: i

        if(symmetry == MS_symmetric) then
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = data_save(i)
                j=j+1
            end do
        else
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = -data_save(i)
                j=j+1
            end do
        end if
    end subroutine
    
    subroutine generate_random_dp_coo_matrix(index_save, data_save, nrows, ncols, symmetry)
        real(dp), allocatable, intent(out) :: data_save(:)
        integer, allocatable, intent(out) :: index_save(:, :)
        integer, intent(in) :: nrows, ncols, symmetry

        ! Internal variables
        integer, allocatable :: pos(:)
        integer :: nnz, nnz_lower, nnz_diag, i, j
        real(dp) :: density

        allocate(pos(nrows * ncols))
        call generate_random_positions(pos, nrows * ncols)
        call random_number(density)
        nnz = ceiling(density*nrows*ncols)
        
        if(symmetry == MS_general) then
            allocate(index_save(2, nnz))
            allocate(data_save(nnz))
            do i = 1, nnz
                index_save(1,i) = mod(pos(i) - 1,ncols) + 1
                index_save(2,i) = (pos(i) - 1)/nrows + 1
            end do
            call generate_random_data_for_dp_coo(data_save, nnz)
        else if(symmetry == MS_symmetric) then
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            nnz_diag = count(mod(pos(1:nnz) - 1,nrows) == (pos(1:nnz) - 1)/ncols) !! diagonal
            allocate(index_save(2, 2*nnz_lower + nnz_diag))
            allocate(data_save(2*nnz_lower + nnz_diag))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_symmetric, j)
            call generate_random_data_for_dp_coo(data_save, nnz_lower + nnz_diag)
            call fill_other_half_dp(index_save, data_save, j, nnz_lower+nnz_diag, MS_symmetric)
        else
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            allocate(index_save(2, 2*nnz_lower))
            allocate(data_save(2*nnz_lower))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_skew_symmetric, j)
            call generate_random_data_for_dp_coo(data_save, nnz_lower)
            call fill_other_half_dp(index_save, data_save, j, nnz_lower, MS_skew_symmetric)
        end if
    end subroutine
    subroutine generate_random_data_for_csp_coo(A, nnz_to_write)
        complex(sp), intent(out) :: A(:)
        integer, intent(in) :: nnz_to_write

        ! Internal variables
        real(sp), allocatable :: R(:, :)

        allocate(R(nnz_to_write, 2))
        call random_number(R(:,1))
        call random_number(R(:,2))
        A(1:nnz_to_write) = cmplx(R(1:nnz_to_write, 1), R(1:nnz_to_write, 2), kind=sp)
    end subroutine

    pure function compare_coo_csp(index_save, index_load, data_save, data_load) result(result)
        complex(sp), intent(in) :: data_save(:), data_load(:)
        integer, intent(in) :: index_save(:, :), index_load(:,:)
        logical :: result
        
        result = all(index_save == index_load) .and. all_close(data_save, data_load)
    end function

    subroutine fill_other_half_csp(index_save, data_save, j, half_nnz, symmetry)
        integer, intent(out) :: index_save(:, :)
        complex(sp), intent(out) :: data_save(:)
        integer, intent(in) :: half_nnz, symmetry
        integer, intent(inout) :: j

        ! Internal variables.
        integer :: i

        if(symmetry == MS_symmetric) then
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = data_save(i)
                j=j+1
            end do
        else if(symmetry == MS_hermitian) then
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = conjg(data_save(i))
                j=j+1
            end do
        else
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = -data_save(i)
                j=j+1
            end do
        end if
    end subroutine
    
    subroutine generate_random_csp_coo_matrix(index_save, data_save, nrows, ncols, symmetry)
        complex(sp), allocatable, intent(out) :: data_save(:)
        integer, allocatable, intent(out) :: index_save(:, :)
        integer, intent(in) :: nrows, ncols, symmetry

        ! Internal variables
        integer, allocatable :: pos(:)
        integer :: nnz, nnz_lower, nnz_diag, i, j
        real(sp) :: density

        allocate(pos(nrows * ncols))
        call generate_random_positions(pos, nrows * ncols)
        call random_number(density)
        nnz = ceiling(density*nrows*ncols)
        
        if(symmetry == MS_general) then
            allocate(index_save(2, nnz))
            allocate(data_save(nnz))
            do i = 1, nnz
                index_save(1,i) = mod(pos(i) - 1,ncols) + 1
                index_save(2,i) = (pos(i) - 1)/nrows + 1
            end do
            call generate_random_data_for_csp_coo(data_save, nnz)
        else if(symmetry == MS_symmetric) then
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            nnz_diag = count(mod(pos(1:nnz) - 1,nrows) == (pos(1:nnz) - 1)/ncols) !! diagonal
            allocate(index_save(2, 2*nnz_lower + nnz_diag))
            allocate(data_save(2*nnz_lower + nnz_diag))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_symmetric, j)
            call generate_random_data_for_csp_coo(data_save, nnz_lower + nnz_diag)
            call fill_other_half_csp(index_save, data_save, j, nnz_lower+nnz_diag, MS_symmetric)
        else if(symmetry == MS_hermitian) then
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            nnz_diag = count(mod(pos(1:nnz) - 1,nrows) == (pos(1:nnz) - 1)/ncols) !! diagonal
            allocate(index_save(2, 2*nnz_lower + nnz_diag))
            allocate(data_save(2*nnz_lower + nnz_diag))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_hermitian, j)
            call generate_random_data_for_csp_coo(data_save, nnz_lower+nnz_diag)
            do i = 1, nnz_lower + nnz_diag
                if(index_save(1, i) == index_save(2,i)) data_save(i) = real(data_save(i))
            end do
            call fill_other_half_csp(index_save, data_save, j, nnz_lower+nnz_diag, MS_hermitian)
        else
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            allocate(index_save(2, 2*nnz_lower))
            allocate(data_save(2*nnz_lower))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_skew_symmetric, j)
            call generate_random_data_for_csp_coo(data_save, nnz_lower)
            call fill_other_half_csp(index_save, data_save, j, nnz_lower, MS_skew_symmetric)
        end if
    end subroutine
    subroutine generate_random_data_for_cdp_coo(A, nnz_to_write)
        complex(dp), intent(out) :: A(:)
        integer, intent(in) :: nnz_to_write

        ! Internal variables
        real(dp), allocatable :: R(:, :)

        allocate(R(nnz_to_write, 2))
        call random_number(R(:,1))
        call random_number(R(:,2))
        A(1:nnz_to_write) = cmplx(R(1:nnz_to_write, 1), R(1:nnz_to_write, 2), kind=dp)
    end subroutine

    pure function compare_coo_cdp(index_save, index_load, data_save, data_load) result(result)
        complex(dp), intent(in) :: data_save(:), data_load(:)
        integer, intent(in) :: index_save(:, :), index_load(:,:)
        logical :: result
        
        result = all(index_save == index_load) .and. all_close(data_save, data_load)
    end function

    subroutine fill_other_half_cdp(index_save, data_save, j, half_nnz, symmetry)
        integer, intent(out) :: index_save(:, :)
        complex(dp), intent(out) :: data_save(:)
        integer, intent(in) :: half_nnz, symmetry
        integer, intent(inout) :: j

        ! Internal variables.
        integer :: i

        if(symmetry == MS_symmetric) then
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = data_save(i)
                j=j+1
            end do
        else if(symmetry == MS_hermitian) then
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = conjg(data_save(i))
                j=j+1
            end do
        else
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = -data_save(i)
                j=j+1
            end do
        end if
    end subroutine
    
    subroutine generate_random_cdp_coo_matrix(index_save, data_save, nrows, ncols, symmetry)
        complex(dp), allocatable, intent(out) :: data_save(:)
        integer, allocatable, intent(out) :: index_save(:, :)
        integer, intent(in) :: nrows, ncols, symmetry

        ! Internal variables
        integer, allocatable :: pos(:)
        integer :: nnz, nnz_lower, nnz_diag, i, j
        real(dp) :: density

        allocate(pos(nrows * ncols))
        call generate_random_positions(pos, nrows * ncols)
        call random_number(density)
        nnz = ceiling(density*nrows*ncols)
        
        if(symmetry == MS_general) then
            allocate(index_save(2, nnz))
            allocate(data_save(nnz))
            do i = 1, nnz
                index_save(1,i) = mod(pos(i) - 1,ncols) + 1
                index_save(2,i) = (pos(i) - 1)/nrows + 1
            end do
            call generate_random_data_for_cdp_coo(data_save, nnz)
        else if(symmetry == MS_symmetric) then
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            nnz_diag = count(mod(pos(1:nnz) - 1,nrows) == (pos(1:nnz) - 1)/ncols) !! diagonal
            allocate(index_save(2, 2*nnz_lower + nnz_diag))
            allocate(data_save(2*nnz_lower + nnz_diag))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_symmetric, j)
            call generate_random_data_for_cdp_coo(data_save, nnz_lower + nnz_diag)
            call fill_other_half_cdp(index_save, data_save, j, nnz_lower+nnz_diag, MS_symmetric)
        else if(symmetry == MS_hermitian) then
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            nnz_diag = count(mod(pos(1:nnz) - 1,nrows) == (pos(1:nnz) - 1)/ncols) !! diagonal
            allocate(index_save(2, 2*nnz_lower + nnz_diag))
            allocate(data_save(2*nnz_lower + nnz_diag))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_hermitian, j)
            call generate_random_data_for_cdp_coo(data_save, nnz_lower+nnz_diag)
            do i = 1, nnz_lower + nnz_diag
                if(index_save(1, i) == index_save(2,i)) data_save(i) = real(data_save(i))
            end do
            call fill_other_half_cdp(index_save, data_save, j, nnz_lower+nnz_diag, MS_hermitian)
        else
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            allocate(index_save(2, 2*nnz_lower))
            allocate(data_save(2*nnz_lower))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_skew_symmetric, j)
            call generate_random_data_for_cdp_coo(data_save, nnz_lower)
            call fill_other_half_cdp(index_save, data_save, j, nnz_lower, MS_skew_symmetric)
        end if
    end subroutine
    subroutine generate_random_data_for_int8_coo(A, nnz_to_write)
        integer(int8), intent(out) :: A(:)
        integer, intent(in) :: nnz_to_write

        ! Internal variables
        real :: rnd
        integer :: i

        do i = 1, nnz_to_write
            call random_number(rnd)
            A(i) = int(rnd * 100 - 50, kind=int8)
        end do
    end subroutine

    pure function compare_coo_int8(index_save, index_load, data_save, data_load) result(result)
        integer(int8), intent(in) :: data_save(:), data_load(:)
        integer, intent(in) :: index_save(:, :), index_load(:,:)
        logical :: result
        
        result = all(index_save == index_load) .and. all(data_save == data_load)
    end function

    subroutine fill_other_half_int8(index_save, data_save, j, half_nnz, symmetry)
        integer, intent(out) :: index_save(:, :)
        integer(int8), intent(out) :: data_save(:)
        integer, intent(in) :: half_nnz, symmetry
        integer, intent(inout) :: j

        ! Internal variables.
        integer :: i

        if(symmetry == MS_symmetric) then
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = data_save(i)
                j=j+1
            end do
        else
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = -data_save(i)
                j=j+1
            end do
        end if
    end subroutine
    
    subroutine generate_random_int8_coo_matrix(index_save, data_save, nrows, ncols, symmetry)
        integer(int8), allocatable, intent(out) :: data_save(:)
        integer, allocatable, intent(out) :: index_save(:, :)
        integer, intent(in) :: nrows, ncols, symmetry

        ! Internal variables
        integer, allocatable :: pos(:)
        integer :: nnz, nnz_lower, nnz_diag, i, j
        real(dp) :: density

        allocate(pos(nrows * ncols))
        call generate_random_positions(pos, nrows * ncols)
        call random_number(density)
        nnz = ceiling(density*nrows*ncols)
        
        if(symmetry == MS_general) then
            allocate(index_save(2, nnz))
            allocate(data_save(nnz))
            do i = 1, nnz
                index_save(1,i) = mod(pos(i) - 1,ncols) + 1
                index_save(2,i) = (pos(i) - 1)/nrows + 1
            end do
            call generate_random_data_for_int8_coo(data_save, nnz)
        else if(symmetry == MS_symmetric) then
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            nnz_diag = count(mod(pos(1:nnz) - 1,nrows) == (pos(1:nnz) - 1)/ncols) !! diagonal
            allocate(index_save(2, 2*nnz_lower + nnz_diag))
            allocate(data_save(2*nnz_lower + nnz_diag))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_symmetric, j)
            call generate_random_data_for_int8_coo(data_save, nnz_lower + nnz_diag)
            call fill_other_half_int8(index_save, data_save, j, nnz_lower+nnz_diag, MS_symmetric)
        else
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            allocate(index_save(2, 2*nnz_lower))
            allocate(data_save(2*nnz_lower))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_skew_symmetric, j)
            call generate_random_data_for_int8_coo(data_save, nnz_lower)
            call fill_other_half_int8(index_save, data_save, j, nnz_lower, MS_skew_symmetric)
        end if
    end subroutine
    subroutine generate_random_data_for_int16_coo(A, nnz_to_write)
        integer(int16), intent(out) :: A(:)
        integer, intent(in) :: nnz_to_write

        ! Internal variables
        real :: rnd
        integer :: i

        do i = 1, nnz_to_write
            call random_number(rnd)
            A(i) = int(rnd * 100 - 50, kind=int16)
        end do
    end subroutine

    pure function compare_coo_int16(index_save, index_load, data_save, data_load) result(result)
        integer(int16), intent(in) :: data_save(:), data_load(:)
        integer, intent(in) :: index_save(:, :), index_load(:,:)
        logical :: result
        
        result = all(index_save == index_load) .and. all(data_save == data_load)
    end function

    subroutine fill_other_half_int16(index_save, data_save, j, half_nnz, symmetry)
        integer, intent(out) :: index_save(:, :)
        integer(int16), intent(out) :: data_save(:)
        integer, intent(in) :: half_nnz, symmetry
        integer, intent(inout) :: j

        ! Internal variables.
        integer :: i

        if(symmetry == MS_symmetric) then
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = data_save(i)
                j=j+1
            end do
        else
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = -data_save(i)
                j=j+1
            end do
        end if
    end subroutine
    
    subroutine generate_random_int16_coo_matrix(index_save, data_save, nrows, ncols, symmetry)
        integer(int16), allocatable, intent(out) :: data_save(:)
        integer, allocatable, intent(out) :: index_save(:, :)
        integer, intent(in) :: nrows, ncols, symmetry

        ! Internal variables
        integer, allocatable :: pos(:)
        integer :: nnz, nnz_lower, nnz_diag, i, j
        real(dp) :: density

        allocate(pos(nrows * ncols))
        call generate_random_positions(pos, nrows * ncols)
        call random_number(density)
        nnz = ceiling(density*nrows*ncols)
        
        if(symmetry == MS_general) then
            allocate(index_save(2, nnz))
            allocate(data_save(nnz))
            do i = 1, nnz
                index_save(1,i) = mod(pos(i) - 1,ncols) + 1
                index_save(2,i) = (pos(i) - 1)/nrows + 1
            end do
            call generate_random_data_for_int16_coo(data_save, nnz)
        else if(symmetry == MS_symmetric) then
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            nnz_diag = count(mod(pos(1:nnz) - 1,nrows) == (pos(1:nnz) - 1)/ncols) !! diagonal
            allocate(index_save(2, 2*nnz_lower + nnz_diag))
            allocate(data_save(2*nnz_lower + nnz_diag))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_symmetric, j)
            call generate_random_data_for_int16_coo(data_save, nnz_lower + nnz_diag)
            call fill_other_half_int16(index_save, data_save, j, nnz_lower+nnz_diag, MS_symmetric)
        else
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            allocate(index_save(2, 2*nnz_lower))
            allocate(data_save(2*nnz_lower))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_skew_symmetric, j)
            call generate_random_data_for_int16_coo(data_save, nnz_lower)
            call fill_other_half_int16(index_save, data_save, j, nnz_lower, MS_skew_symmetric)
        end if
    end subroutine
    subroutine generate_random_data_for_int32_coo(A, nnz_to_write)
        integer(int32), intent(out) :: A(:)
        integer, intent(in) :: nnz_to_write

        ! Internal variables
        real :: rnd
        integer :: i

        do i = 1, nnz_to_write
            call random_number(rnd)
            A(i) = int(rnd * 100 - 50, kind=int32)
        end do
    end subroutine

    pure function compare_coo_int32(index_save, index_load, data_save, data_load) result(result)
        integer(int32), intent(in) :: data_save(:), data_load(:)
        integer, intent(in) :: index_save(:, :), index_load(:,:)
        logical :: result
        
        result = all(index_save == index_load) .and. all(data_save == data_load)
    end function

    subroutine fill_other_half_int32(index_save, data_save, j, half_nnz, symmetry)
        integer, intent(out) :: index_save(:, :)
        integer(int32), intent(out) :: data_save(:)
        integer, intent(in) :: half_nnz, symmetry
        integer, intent(inout) :: j

        ! Internal variables.
        integer :: i

        if(symmetry == MS_symmetric) then
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = data_save(i)
                j=j+1
            end do
        else
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = -data_save(i)
                j=j+1
            end do
        end if
    end subroutine
    
    subroutine generate_random_int32_coo_matrix(index_save, data_save, nrows, ncols, symmetry)
        integer(int32), allocatable, intent(out) :: data_save(:)
        integer, allocatable, intent(out) :: index_save(:, :)
        integer, intent(in) :: nrows, ncols, symmetry

        ! Internal variables
        integer, allocatable :: pos(:)
        integer :: nnz, nnz_lower, nnz_diag, i, j
        real(dp) :: density

        allocate(pos(nrows * ncols))
        call generate_random_positions(pos, nrows * ncols)
        call random_number(density)
        nnz = ceiling(density*nrows*ncols)
        
        if(symmetry == MS_general) then
            allocate(index_save(2, nnz))
            allocate(data_save(nnz))
            do i = 1, nnz
                index_save(1,i) = mod(pos(i) - 1,ncols) + 1
                index_save(2,i) = (pos(i) - 1)/nrows + 1
            end do
            call generate_random_data_for_int32_coo(data_save, nnz)
        else if(symmetry == MS_symmetric) then
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            nnz_diag = count(mod(pos(1:nnz) - 1,nrows) == (pos(1:nnz) - 1)/ncols) !! diagonal
            allocate(index_save(2, 2*nnz_lower + nnz_diag))
            allocate(data_save(2*nnz_lower + nnz_diag))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_symmetric, j)
            call generate_random_data_for_int32_coo(data_save, nnz_lower + nnz_diag)
            call fill_other_half_int32(index_save, data_save, j, nnz_lower+nnz_diag, MS_symmetric)
        else
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            allocate(index_save(2, 2*nnz_lower))
            allocate(data_save(2*nnz_lower))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_skew_symmetric, j)
            call generate_random_data_for_int32_coo(data_save, nnz_lower)
            call fill_other_half_int32(index_save, data_save, j, nnz_lower, MS_skew_symmetric)
        end if
    end subroutine
    subroutine generate_random_data_for_int64_coo(A, nnz_to_write)
        integer(int64), intent(out) :: A(:)
        integer, intent(in) :: nnz_to_write

        ! Internal variables
        real :: rnd
        integer :: i

        do i = 1, nnz_to_write
            call random_number(rnd)
            A(i) = int(rnd * 100 - 50, kind=int64)
        end do
    end subroutine

    pure function compare_coo_int64(index_save, index_load, data_save, data_load) result(result)
        integer(int64), intent(in) :: data_save(:), data_load(:)
        integer, intent(in) :: index_save(:, :), index_load(:,:)
        logical :: result
        
        result = all(index_save == index_load) .and. all(data_save == data_load)
    end function

    subroutine fill_other_half_int64(index_save, data_save, j, half_nnz, symmetry)
        integer, intent(out) :: index_save(:, :)
        integer(int64), intent(out) :: data_save(:)
        integer, intent(in) :: half_nnz, symmetry
        integer, intent(inout) :: j

        ! Internal variables.
        integer :: i

        if(symmetry == MS_symmetric) then
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = data_save(i)
                j=j+1
            end do
        else
            do i = 1, half_nnz
                if(index_save(1,i) == index_save(2,i)) cycle
                index_save(1,j) = index_save(2,i)
                index_save(2,j) = index_save(1,i)
                data_save(j) = -data_save(i)
                j=j+1
            end do
        end if
    end subroutine
    
    subroutine generate_random_int64_coo_matrix(index_save, data_save, nrows, ncols, symmetry)
        integer(int64), allocatable, intent(out) :: data_save(:)
        integer, allocatable, intent(out) :: index_save(:, :)
        integer, intent(in) :: nrows, ncols, symmetry

        ! Internal variables
        integer, allocatable :: pos(:)
        integer :: nnz, nnz_lower, nnz_diag, i, j
        real(dp) :: density

        allocate(pos(nrows * ncols))
        call generate_random_positions(pos, nrows * ncols)
        call random_number(density)
        nnz = ceiling(density*nrows*ncols)
        
        if(symmetry == MS_general) then
            allocate(index_save(2, nnz))
            allocate(data_save(nnz))
            do i = 1, nnz
                index_save(1,i) = mod(pos(i) - 1,ncols) + 1
                index_save(2,i) = (pos(i) - 1)/nrows + 1
            end do
            call generate_random_data_for_int64_coo(data_save, nnz)
        else if(symmetry == MS_symmetric) then
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            nnz_diag = count(mod(pos(1:nnz) - 1,nrows) == (pos(1:nnz) - 1)/ncols) !! diagonal
            allocate(index_save(2, 2*nnz_lower + nnz_diag))
            allocate(data_save(2*nnz_lower + nnz_diag))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_symmetric, j)
            call generate_random_data_for_int64_coo(data_save, nnz_lower + nnz_diag)
            call fill_other_half_int64(index_save, data_save, j, nnz_lower+nnz_diag, MS_symmetric)
        else
            nnz_lower = count(mod(pos(1:nnz) - 1,nrows) > (pos(1:nnz) - 1)/ncols) !! lower triangular part
            allocate(index_save(2, 2*nnz_lower))
            allocate(data_save(2*nnz_lower))
            call fill_first_half_indices(index_save, pos, nnz, nrows, ncols, MS_skew_symmetric, j)
            call generate_random_data_for_int64_coo(data_save, nnz_lower)
            call fill_other_half_int64(index_save, data_save, j, nnz_lower, MS_skew_symmetric)
        end if
    end subroutine

    subroutine test_io_mm_coordinate(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer :: nrows, ncols
            real(sp), allocatable :: data_save(:), data_load(:)
            integer, allocatable :: index_save(:, :), index_load(:,:)
            logical :: result

            nrows = 5
            ncols = 5

            call random_seed()
            ! General matrix
            call generate_random_sp_coo_matrix(index_save, data_save, nrows, ncols, MS_general)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_sp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=unspecified, type=real(sp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Symmetric matrix
            call generate_random_sp_coo_matrix(index_save, data_save, nrows, ncols, MS_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_sp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=symmetric, type=real(sp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Skew-symmetric matrix
            call generate_random_sp_coo_matrix(index_save, data_save, nrows, ncols, MS_skew_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_sp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=skew-symmetric, type=real(sp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)
            
        end block
        block
            integer :: nrows, ncols
            real(dp), allocatable :: data_save(:), data_load(:)
            integer, allocatable :: index_save(:, :), index_load(:,:)
            logical :: result

            nrows = 5
            ncols = 5

            call random_seed()
            ! General matrix
            call generate_random_dp_coo_matrix(index_save, data_save, nrows, ncols, MS_general)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_dp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=unspecified, type=real(dp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Symmetric matrix
            call generate_random_dp_coo_matrix(index_save, data_save, nrows, ncols, MS_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_dp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=symmetric, type=real(dp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Skew-symmetric matrix
            call generate_random_dp_coo_matrix(index_save, data_save, nrows, ncols, MS_skew_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_dp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=skew-symmetric, type=real(dp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)
            
        end block
        block
            integer :: nrows, ncols
            complex(sp), allocatable :: data_save(:), data_load(:)
            integer, allocatable :: index_save(:, :), index_load(:,:)
            logical :: result

            nrows = 5
            ncols = 5

            call random_seed()
            ! General matrix
            call generate_random_csp_coo_matrix(index_save, data_save, nrows, ncols, MS_general)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_csp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=unspecified, type=complex(sp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Symmetric matrix
            call generate_random_csp_coo_matrix(index_save, data_save, nrows, ncols, MS_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_csp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=symmetric, type=complex(sp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Skew-symmetric matrix
            call generate_random_csp_coo_matrix(index_save, data_save, nrows, ncols, MS_skew_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_csp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=skew-symmetric, type=complex(sp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)
            
            ! Hermitian matrix
            call generate_random_csp_coo_matrix(index_save, data_save, nrows, ncols, MS_hermitian)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "hermitian", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_csp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=hermitian, type=complex(sp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)
        end block
        block
            integer :: nrows, ncols
            complex(dp), allocatable :: data_save(:), data_load(:)
            integer, allocatable :: index_save(:, :), index_load(:,:)
            logical :: result

            nrows = 5
            ncols = 5

            call random_seed()
            ! General matrix
            call generate_random_cdp_coo_matrix(index_save, data_save, nrows, ncols, MS_general)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_cdp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=unspecified, type=complex(dp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Symmetric matrix
            call generate_random_cdp_coo_matrix(index_save, data_save, nrows, ncols, MS_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_cdp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=symmetric, type=complex(dp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Skew-symmetric matrix
            call generate_random_cdp_coo_matrix(index_save, data_save, nrows, ncols, MS_skew_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_cdp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=skew-symmetric, type=complex(dp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)
            
            ! Hermitian matrix
            call generate_random_cdp_coo_matrix(index_save, data_save, nrows, ncols, MS_hermitian)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "hermitian", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_cdp(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=hermitian, type=complex(dp)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)
        end block
        block
            integer :: nrows, ncols
            integer(int8), allocatable :: data_save(:), data_load(:)
            integer, allocatable :: index_save(:, :), index_load(:,:)
            logical :: result

            nrows = 5
            ncols = 5

            call random_seed()
            ! General matrix
            call generate_random_int8_coo_matrix(index_save, data_save, nrows, ncols, MS_general)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_int8(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=unspecified, type=integer(int8)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Symmetric matrix
            call generate_random_int8_coo_matrix(index_save, data_save, nrows, ncols, MS_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_int8(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=symmetric, type=integer(int8)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Skew-symmetric matrix
            call generate_random_int8_coo_matrix(index_save, data_save, nrows, ncols, MS_skew_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_int8(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=skew-symmetric, type=integer(int8)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)
            
        end block
        block
            integer :: nrows, ncols
            integer(int16), allocatable :: data_save(:), data_load(:)
            integer, allocatable :: index_save(:, :), index_load(:,:)
            logical :: result

            nrows = 5
            ncols = 5

            call random_seed()
            ! General matrix
            call generate_random_int16_coo_matrix(index_save, data_save, nrows, ncols, MS_general)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_int16(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=unspecified, type=integer(int16)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Symmetric matrix
            call generate_random_int16_coo_matrix(index_save, data_save, nrows, ncols, MS_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_int16(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=symmetric, type=integer(int16)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Skew-symmetric matrix
            call generate_random_int16_coo_matrix(index_save, data_save, nrows, ncols, MS_skew_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_int16(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=skew-symmetric, type=integer(int16)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)
            
        end block
        block
            integer :: nrows, ncols
            integer(int32), allocatable :: data_save(:), data_load(:)
            integer, allocatable :: index_save(:, :), index_load(:,:)
            logical :: result

            nrows = 5
            ncols = 5

            call random_seed()
            ! General matrix
            call generate_random_int32_coo_matrix(index_save, data_save, nrows, ncols, MS_general)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_int32(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=unspecified, type=integer(int32)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Symmetric matrix
            call generate_random_int32_coo_matrix(index_save, data_save, nrows, ncols, MS_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_int32(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=symmetric, type=integer(int32)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Skew-symmetric matrix
            call generate_random_int32_coo_matrix(index_save, data_save, nrows, ncols, MS_skew_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_int32(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=skew-symmetric, type=integer(int32)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)
            
        end block
        block
            integer :: nrows, ncols
            integer(int64), allocatable :: data_save(:), data_load(:)
            integer, allocatable :: index_save(:, :), index_load(:,:)
            logical :: result

            nrows = 5
            ncols = 5

            call random_seed()
            ! General matrix
            call generate_random_int64_coo_matrix(index_save, data_save, nrows, ncols, MS_general)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_int64(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=unspecified, type=integer(int64)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Symmetric matrix
            call generate_random_int64_coo_matrix(index_save, data_save, nrows, ncols, MS_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_int64(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=symmetric, type=integer(int64)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)

            ! Skew-symmetric matrix
            call generate_random_int64_coo_matrix(index_save, data_save, nrows, ncols, MS_skew_symmetric)
            call save_mm("test_mmio_sparse.mtx", index_save, data_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_sparse.mtx", index_load, data_load)
            result = compare_coo_int64(index_save, index_load, data_save, data_load)
            call check(error, result, .true.,&
                "MM coordinate test failed: symmetry_arg=skew-symmetric, type=integer(int64)")
            if(allocated(error)) return
            if(allocated(index_save)) deallocate(index_save)
            if(allocated(data_save)) deallocate(data_save)
            
        end block
    end subroutine

end module


program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_io_mm, only : collect_suite
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("io_mm", collect_suite) &
        ]

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if
end program
