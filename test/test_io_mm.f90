module test_io_mm
    use testdrive, only : new_unittest, unittest_type, error_type, check, skip_test
    use stdlib_kinds
    use stdlib_math, only: all_close
    use stdlib_io_mm
    implicit none

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
            call random_number(matrix_save)


            call save_mm("test_mmio_dense.mtx", matrix_save, format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=unspecified, type=real(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=auto, type=real(sp)")
            if(allocated(error)) return

            ! Symmetric matrix
            call random_number(A)


            ! Construct symmetric matrix using (A + A.T)
            matrix_save = A + transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=symmetric, type=real(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=auto, type=real(sp)")
            if(allocated(error)) return

            ! Skew-symmetric matrix
            call random_number(A)


            ! Construct symmetric matrix using (A - A.T)
            matrix_save = A - transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=skew-symmetric, type=real(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
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
            call random_number(matrix_save)


            call save_mm("test_mmio_dense.mtx", matrix_save, format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=unspecified, type=real(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=auto, type=real(dp)")
            if(allocated(error)) return

            ! Symmetric matrix
            call random_number(A)


            ! Construct symmetric matrix using (A + A.T)
            matrix_save = A + transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=symmetric, type=real(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=auto, type=real(dp)")
            if(allocated(error)) return

            ! Skew-symmetric matrix
            call random_number(A)


            ! Construct symmetric matrix using (A - A.T)
            matrix_save = A - transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=skew-symmetric, type=real(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=auto, type=real(dp)")
            if(allocated(error)) return

        end block
        block
            integer, parameter :: n = 5
            complex(sp), allocatable :: matrix_save(:, :), matrix_load(:, :), A(:, :)
            real(sp), allocatable :: R(:, :)
            real(sp), allocatable :: I(:,:)
            logical :: result
            allocate(matrix_save(n,n))
            allocate(A(n,n))
            allocate(R(n,n))
            allocate(I(n,n))

            ! General matrix

            call random_number(R)
            call random_number(I)
            matrix_save = cmplx(R, I, kind=sp)

            call save_mm("test_mmio_dense.mtx", matrix_save, format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=unspecified, type=complex(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=auto, type=complex(sp)")
            if(allocated(error)) return

            ! Symmetric matrix

            call random_number(R)
            call random_number(I)
            A = cmplx(R, I, kind=sp)

            ! Construct symmetric matrix using (A + A.T)
            matrix_save = A + transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=symmetric, type=complex(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=auto, type=complex(sp)")
            if(allocated(error)) return

            ! Skew-symmetric matrix

            call random_number(R)
            call random_number(I)
            A = cmplx(R, I, kind=sp)

            ! Construct symmetric matrix using (A - A.T)
            matrix_save = A - transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=skew-symmetric, type=complex(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=auto, type=complex(sp)")
            if(allocated(error)) return

            ! Hermitian matrix
            call random_number(R)
            call random_number(I)
            A = cmplx(R, I, kind=sp)
            ! Construct symmetric matrix using (A + A.H)
            matrix_save = A + transpose(conjg(A))
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "hermitian", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=hermitian, symmetry_arg=hermitian, type=complex(sp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=hermitian, symmetry_arg=auto, type=complex(sp)")
            if(allocated(error)) return
        end block
        block
            integer, parameter :: n = 5
            complex(dp), allocatable :: matrix_save(:, :), matrix_load(:, :), A(:, :)
            real(dp), allocatable :: R(:, :)
            real(dp), allocatable :: I(:,:)
            logical :: result
            allocate(matrix_save(n,n))
            allocate(A(n,n))
            allocate(R(n,n))
            allocate(I(n,n))

            ! General matrix

            call random_number(R)
            call random_number(I)
            matrix_save = cmplx(R, I, kind=dp)

            call save_mm("test_mmio_dense.mtx", matrix_save, format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=unspecified, type=complex(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=general, symmetry_arg=auto, type=complex(dp)")
            if(allocated(error)) return

            ! Symmetric matrix

            call random_number(R)
            call random_number(I)
            A = cmplx(R, I, kind=dp)

            ! Construct symmetric matrix using (A + A.T)
            matrix_save = A + transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=symmetric, type=complex(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=symmetric, symmetry_arg=auto, type=complex(dp)")
            if(allocated(error)) return

            ! Skew-symmetric matrix

            call random_number(R)
            call random_number(I)
            A = cmplx(R, I, kind=dp)

            ! Construct symmetric matrix using (A - A.T)
            matrix_save = A - transpose(A)
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "skew-symmetric", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=skew-symmetric, type=complex(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=skew-symmetric, symmetry_arg=auto, type=complex(dp)")
            if(allocated(error)) return

            ! Hermitian matrix
            call random_number(R)
            call random_number(I)
            A = cmplx(R, I, kind=dp)
            ! Construct symmetric matrix using (A + A.H)
            matrix_save = A + transpose(conjg(A))
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "hermitian", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=hermitian, symmetry_arg=hermitian, type=complex(dp)")
            if(allocated(error)) return
            ! Check if symmetry = auto
            call save_mm("test_mmio_dense.mtx", matrix_save, symmetry = "auto", format = "G0")
            call load_mm("test_mmio_dense.mtx", matrix_load)
            result = all_close(matrix_save, matrix_load)
            call check(error, result, .true.,&
                "MM array test failed: matrix=hermitian, symmetry_arg=auto, type=complex(dp)")
            if(allocated(error)) return
        end block
    end subroutine

    subroutine test_io_mm_coordinate(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: n = 5
            real(sp), allocatable :: index_save(:), index_load(:)
            integer, allocatable :: data_save(:, :), data_load(:,:)
        end block
        block
            integer, parameter :: n = 5
            real(dp), allocatable :: index_save(:), index_load(:)
            integer, allocatable :: data_save(:, :), data_load(:,:)
        end block
        block
            integer, parameter :: n = 5
            complex(sp), allocatable :: index_save(:), index_load(:)
            integer, allocatable :: data_save(:, :), data_load(:,:)
        end block
        block
            integer, parameter :: n = 5
            complex(dp), allocatable :: index_save(:), index_load(:)
            integer, allocatable :: data_save(:, :), data_load(:,:)
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
