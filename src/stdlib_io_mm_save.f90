! SPDX-Identifier: MIT


!> Implementation for saving multidimensional arrays to Matrix Market files
submodule (stdlib_io_mm) stdlib_io_mm_save
    use stdlib_error, only : error_stop
    use stdlib_strings, only : to_string
    use stdlib_io, only : open
    use stdlib_ascii, only : to_lower
    use stdlib_constants, only : zero_sp, zero_dp, zero_sp, zero_dp
    implicit none

    ! Matrix Market format constants
    character(len=*), parameter :: MM_BANNER = "%%MatrixMarket"
    character(len=*), parameter :: MM_COMMENT_CHAR = "%"
    
    ! Matrix Market object types
    character(len=*), parameter :: &
        MM_MATRIX = "matrix", &
        MM_VECTOR = "vector"
    
    ! Matrix Market format types  
    character(len=*), parameter :: &
        MM_COORDINATE = "coordinate", &
        MM_ARRAY = "array"
    
    ! Matrix Market data types
    character(len=*), parameter :: &
        MM_REAL = "real", &
        MM_COMPLEX = "complex", &
        MM_INTEGER = "integer", &
        MM_PATTERN = "pattern"
    
    ! Matrix Market storage schemes
    character(len=*), parameter :: &
        MM_GENERAL = "general", &
        MM_SYMMETRIC = "symmetric", &
        MM_SKEW_SYMMETRIC = "skew-symmetric", &
        MM_HERMITIAN = "hermitian"

contains

    module subroutine save_mm_dense_sp(filename, matrix, comment, symmetry, iostat, iomsg)
        !> Name of the Matrix Market file to save to
        character(len=*), intent(in) :: filename
        !> Matrix to be saved to the Matrix Market file
        real(sp), intent(in) :: matrix(:,:)
        !> Optional comment information
        character(len=*), intent(in), optional :: comment
        !> Symmetry type of the matrix (general, symmetric, skew-symmetric, hermitian)
        character(len=*), intent(in), optional :: symmetry
        !> Error status of saving, zero on success
        integer, intent(out), optional :: iostat
        !> Associated error message in case of non-zero status code
        character(len=:), allocatable, intent(out), optional :: iomsg

        integer :: io, stat, i, j, nnz
        character(len=:), allocatable :: msg
        character(len=32) :: field_type
        character(len=32) :: symmetry_

        io = open(filename, "w", iostat=stat)
        if (stat /= 0) then
            if (present(iostat)) then
                iostat = stat
                if (present(iomsg)) iomsg = "Could not create file: " // filename
                return
            else
                call error_stop("Could not create file: " // filename)
            end if
        end if

        ! Determine symmetry type
        symmetry_ = "general"
        if (present(symmetry)) then
            symmetry_ = to_lower(trim(symmetry))
        end if

        ! Determine field type based on matrix type
        field_type = MM_REAL

        catch: block
            ! Write header
            call write_mm_header(io, MM_ARRAY, field_type, symmetry_, &
                                size(matrix, 1), size(matrix, 2), nnz, comment, stat, msg)
            if (stat /= 0) exit catch

            ! Write array format (column-major order)
            if(symmetry_ == MM_GENERAL) then
                do j = 1, size(matrix, 2)
                    do i = 1, size(matrix, 1)
                        write(io, '(ES24.16E3)', iostat=stat) matrix(i, j)
                        if (stat /= 0) then
                            msg = "Error writing array element (" // &
                                to_string(i) // "," // to_string(j) // ")"
                            exit catch
                        end if
                    end do
                end do
            else
                ! For symmetric, skew-symmetric, hermitian matrices, only write the
                ! lower triangle (including diagonal)
                do j = 1, size(matrix, 2)
                    do i = j, size(matrix, 1)
                        write(io, '(ES24.16E3)', iostat=stat) matrix(i, j)
                        if (stat /= 0) then
                            msg = "Error writing array element (" // &
                                to_string(i) // "," // to_string(j) // ")"
                            exit catch
                        end if
                    end do
                end do
            end if
        end block catch
        
        close(io)

        if (present(iostat)) then
            iostat = stat
        else if (stat /= 0) then
            if (allocated(msg)) then
                call error_stop("Failed to save Matrix Market file '" // filename // "': " // msg)
            else
                call error_stop("Failed to save Matrix Market file '" // filename // "'")
            end if
        end if

        if (present(iomsg) .and. allocated(msg)) call move_alloc(msg, iomsg)
    end subroutine 
    module subroutine save_mm_dense_dp(filename, matrix, comment, symmetry, iostat, iomsg)
        !> Name of the Matrix Market file to save to
        character(len=*), intent(in) :: filename
        !> Matrix to be saved to the Matrix Market file
        real(dp), intent(in) :: matrix(:,:)
        !> Optional comment information
        character(len=*), intent(in), optional :: comment
        !> Symmetry type of the matrix (general, symmetric, skew-symmetric, hermitian)
        character(len=*), intent(in), optional :: symmetry
        !> Error status of saving, zero on success
        integer, intent(out), optional :: iostat
        !> Associated error message in case of non-zero status code
        character(len=:), allocatable, intent(out), optional :: iomsg

        integer :: io, stat, i, j, nnz
        character(len=:), allocatable :: msg
        character(len=32) :: field_type
        character(len=32) :: symmetry_

        io = open(filename, "w", iostat=stat)
        if (stat /= 0) then
            if (present(iostat)) then
                iostat = stat
                if (present(iomsg)) iomsg = "Could not create file: " // filename
                return
            else
                call error_stop("Could not create file: " // filename)
            end if
        end if

        ! Determine symmetry type
        symmetry_ = "general"
        if (present(symmetry)) then
            symmetry_ = to_lower(trim(symmetry))
        end if

        ! Determine field type based on matrix type
        field_type = MM_REAL

        catch: block
            ! Write header
            call write_mm_header(io, MM_ARRAY, field_type, symmetry_, &
                                size(matrix, 1), size(matrix, 2), nnz, comment, stat, msg)
            if (stat /= 0) exit catch

            ! Write array format (column-major order)
            if(symmetry_ == MM_GENERAL) then
                do j = 1, size(matrix, 2)
                    do i = 1, size(matrix, 1)
                        write(io, '(ES24.16E3)', iostat=stat) matrix(i, j)
                        if (stat /= 0) then
                            msg = "Error writing array element (" // &
                                to_string(i) // "," // to_string(j) // ")"
                            exit catch
                        end if
                    end do
                end do
            else
                ! For symmetric, skew-symmetric, hermitian matrices, only write the
                ! lower triangle (including diagonal)
                do j = 1, size(matrix, 2)
                    do i = j, size(matrix, 1)
                        write(io, '(ES24.16E3)', iostat=stat) matrix(i, j)
                        if (stat /= 0) then
                            msg = "Error writing array element (" // &
                                to_string(i) // "," // to_string(j) // ")"
                            exit catch
                        end if
                    end do
                end do
            end if
        end block catch
        
        close(io)

        if (present(iostat)) then
            iostat = stat
        else if (stat /= 0) then
            if (allocated(msg)) then
                call error_stop("Failed to save Matrix Market file '" // filename // "': " // msg)
            else
                call error_stop("Failed to save Matrix Market file '" // filename // "'")
            end if
        end if

        if (present(iomsg) .and. allocated(msg)) call move_alloc(msg, iomsg)
    end subroutine 
    module subroutine save_mm_dense_csp(filename, matrix, comment, symmetry, iostat, iomsg)
        !> Name of the Matrix Market file to save to
        character(len=*), intent(in) :: filename
        !> Matrix to be saved to the Matrix Market file
        complex(sp), intent(in) :: matrix(:,:)
        !> Optional comment information
        character(len=*), intent(in), optional :: comment
        !> Symmetry type of the matrix (general, symmetric, skew-symmetric, hermitian)
        character(len=*), intent(in), optional :: symmetry
        !> Error status of saving, zero on success
        integer, intent(out), optional :: iostat
        !> Associated error message in case of non-zero status code
        character(len=:), allocatable, intent(out), optional :: iomsg

        integer :: io, stat, i, j, nnz
        character(len=:), allocatable :: msg
        character(len=32) :: field_type
        character(len=32) :: symmetry_
        real(sp) :: real_part, imag_part

        io = open(filename, "w", iostat=stat)
        if (stat /= 0) then
            if (present(iostat)) then
                iostat = stat
                if (present(iomsg)) iomsg = "Could not create file: " // filename
                return
            else
                call error_stop("Could not create file: " // filename)
            end if
        end if

        ! Determine symmetry type
        symmetry_ = "general"
        if (present(symmetry)) then
            symmetry_ = to_lower(trim(symmetry))
        end if

        ! Determine field type based on matrix type
        field_type = MM_COMPLEX

        catch: block
            ! Write header
            call write_mm_header(io, MM_ARRAY, field_type, symmetry_, &
                                size(matrix, 1), size(matrix, 2), nnz, comment, stat, msg)
            if (stat /= 0) exit catch

            ! Write array format (column-major order)
            if(symmetry_ == MM_GENERAL) then
                do j = 1, size(matrix, 2)
                    do i = 1, size(matrix, 1)
                        real_part = real(matrix(i, j), kind=sp)
                        imag_part = aimag(matrix(i, j))
                        write(io, '(ES24.16E3,1X,ES24.16E3)', iostat=stat) real_part, imag_part
                        if (stat /= 0) then
                            msg = "Error writing array element (" // &
                                to_string(i) // "," // to_string(j) // ")"
                            exit catch
                        end if
                    end do
                end do
            else
                ! For symmetric, skew-symmetric, hermitian matrices, only write the
                ! lower triangle (including diagonal)
                do j = 1, size(matrix, 2)
                    do i = j, size(matrix, 1)
                        real_part = real(matrix(i, j), kind=sp)
                        imag_part = aimag(matrix(i, j))
                        write(io, '(ES24.16E3,1X,ES24.16E3)', iostat=stat) real_part, imag_part
                        if (stat /= 0) then
                            msg = "Error writing array element (" // &
                                to_string(i) // "," // to_string(j) // ")"
                            exit catch
                        end if
                    end do
                end do
            end if
        end block catch
        
        close(io)

        if (present(iostat)) then
            iostat = stat
        else if (stat /= 0) then
            if (allocated(msg)) then
                call error_stop("Failed to save Matrix Market file '" // filename // "': " // msg)
            else
                call error_stop("Failed to save Matrix Market file '" // filename // "'")
            end if
        end if

        if (present(iomsg) .and. allocated(msg)) call move_alloc(msg, iomsg)
    end subroutine 
    module subroutine save_mm_dense_cdp(filename, matrix, comment, symmetry, iostat, iomsg)
        !> Name of the Matrix Market file to save to
        character(len=*), intent(in) :: filename
        !> Matrix to be saved to the Matrix Market file
        complex(dp), intent(in) :: matrix(:,:)
        !> Optional comment information
        character(len=*), intent(in), optional :: comment
        !> Symmetry type of the matrix (general, symmetric, skew-symmetric, hermitian)
        character(len=*), intent(in), optional :: symmetry
        !> Error status of saving, zero on success
        integer, intent(out), optional :: iostat
        !> Associated error message in case of non-zero status code
        character(len=:), allocatable, intent(out), optional :: iomsg

        integer :: io, stat, i, j, nnz
        character(len=:), allocatable :: msg
        character(len=32) :: field_type
        character(len=32) :: symmetry_
        real(dp) :: real_part, imag_part

        io = open(filename, "w", iostat=stat)
        if (stat /= 0) then
            if (present(iostat)) then
                iostat = stat
                if (present(iomsg)) iomsg = "Could not create file: " // filename
                return
            else
                call error_stop("Could not create file: " // filename)
            end if
        end if

        ! Determine symmetry type
        symmetry_ = "general"
        if (present(symmetry)) then
            symmetry_ = to_lower(trim(symmetry))
        end if

        ! Determine field type based on matrix type
        field_type = MM_COMPLEX

        catch: block
            ! Write header
            call write_mm_header(io, MM_ARRAY, field_type, symmetry_, &
                                size(matrix, 1), size(matrix, 2), nnz, comment, stat, msg)
            if (stat /= 0) exit catch

            ! Write array format (column-major order)
            if(symmetry_ == MM_GENERAL) then
                do j = 1, size(matrix, 2)
                    do i = 1, size(matrix, 1)
                        real_part = real(matrix(i, j), kind=dp)
                        imag_part = aimag(matrix(i, j))
                        write(io, '(ES24.16E3,1X,ES24.16E3)', iostat=stat) real_part, imag_part
                        if (stat /= 0) then
                            msg = "Error writing array element (" // &
                                to_string(i) // "," // to_string(j) // ")"
                            exit catch
                        end if
                    end do
                end do
            else
                ! For symmetric, skew-symmetric, hermitian matrices, only write the
                ! lower triangle (including diagonal)
                do j = 1, size(matrix, 2)
                    do i = j, size(matrix, 1)
                        real_part = real(matrix(i, j), kind=dp)
                        imag_part = aimag(matrix(i, j))
                        write(io, '(ES24.16E3,1X,ES24.16E3)', iostat=stat) real_part, imag_part
                        if (stat /= 0) then
                            msg = "Error writing array element (" // &
                                to_string(i) // "," // to_string(j) // ")"
                            exit catch
                        end if
                    end do
                end do
            end if
        end block catch
        
        close(io)

        if (present(iostat)) then
            iostat = stat
        else if (stat /= 0) then
            if (allocated(msg)) then
                call error_stop("Failed to save Matrix Market file '" // filename // "': " // msg)
            else
                call error_stop("Failed to save Matrix Market file '" // filename // "'")
            end if
        end if

        if (present(iomsg) .and. allocated(msg)) call move_alloc(msg, iomsg)
    end subroutine 

    module subroutine save_mm_coo_sp(filename, matrix, comment, symmetry, iostat, iomsg)
        !> Name of the Matrix Market file to save to
        character(len=*), intent(in) :: filename
        !> Matrix to be saved to the Matrix Market file
        type(COO_sp_type), intent(in) :: matrix
        !> Optional comment information
        character(len=*), intent(in), optional :: comment
        !> Symmetry type of the matrix (general, symmetric, skew-symmetric, hermitian)
        character(len=*), intent(in), optional :: symmetry
        !> Error status of saving, zero on success
        integer, intent(out), optional :: iostat
        !> Associated error message in case of non-zero status code
        character(len=:), allocatable, intent(out), optional :: iomsg

        integer :: io, stat, i, j, nnz
        character(len=:), allocatable :: msg
        character(len=32) :: field_type
        character(len=32) :: symmetry_

        io = open(filename, "w", iostat=stat)
        if (stat /= 0) then
            if (present(iostat)) then
                iostat = stat
                if (present(iomsg)) iomsg = "Could not create file: " // filename
                return
            else
                call error_stop("Could not create file: " // filename)
            end if
        end if

        ! Determine symmetry type
        symmetry_ = "general"
        if (present(symmetry)) then
            symmetry_ = to_lower(trim(symmetry)) 
        end if

        ! Determine field type based on matrix type
        field_type = MM_REAL

        catch: block
            ! Write header
            call write_mm_header(io, MM_COORDINATE, field_type, symmetry_, &
                                matrix%nrows, matrix%ncols, matrix%nnz, comment, stat, msg)
            if (stat /= 0) exit catch

            ! Write array format (column-major order)
            if(symmetry_ == MM_GENERAL) then
                do i = 1, matrix%nnz
                    write(io, '(ES24.16E3)', iostat=stat) matrix%data(i)
                    if (stat /= 0) then
                        msg = "Error writing array element (" // to_string(i) // ")"
                        exit catch
                    end if
                end do
            else
                ! For symmetric, skew-symmetric, hermitian matrices, only write the
                ! lower triangle (including diagonal)
                do i = 1, matrix%nnz
                    if(matrix%index(1,i) > matrix%index(2,i)) cycle
                    write(io, '(ES24.16E3)', iostat=stat) matrix%data(i)
                    if (stat /= 0) then
                        msg = "Error writing array element (" // to_string(i) // ")"
                        exit catch
                    end if
                end do
            end if
        end block catch
        
        close(io)

        if (present(iostat)) then
            iostat = stat
        else if (stat /= 0) then
            if (allocated(msg)) then
                call error_stop("Failed to save Matrix Market file '" // filename // "': " // msg)
            else
                call error_stop("Failed to save Matrix Market file '" // filename // "'")
            end if
        end if

        if (present(iomsg) .and. allocated(msg)) call move_alloc(msg, iomsg)
    end subroutine 
    module subroutine save_mm_coo_dp(filename, matrix, comment, symmetry, iostat, iomsg)
        !> Name of the Matrix Market file to save to
        character(len=*), intent(in) :: filename
        !> Matrix to be saved to the Matrix Market file
        type(COO_dp_type), intent(in) :: matrix
        !> Optional comment information
        character(len=*), intent(in), optional :: comment
        !> Symmetry type of the matrix (general, symmetric, skew-symmetric, hermitian)
        character(len=*), intent(in), optional :: symmetry
        !> Error status of saving, zero on success
        integer, intent(out), optional :: iostat
        !> Associated error message in case of non-zero status code
        character(len=:), allocatable, intent(out), optional :: iomsg

        integer :: io, stat, i, j, nnz
        character(len=:), allocatable :: msg
        character(len=32) :: field_type
        character(len=32) :: symmetry_

        io = open(filename, "w", iostat=stat)
        if (stat /= 0) then
            if (present(iostat)) then
                iostat = stat
                if (present(iomsg)) iomsg = "Could not create file: " // filename
                return
            else
                call error_stop("Could not create file: " // filename)
            end if
        end if

        ! Determine symmetry type
        symmetry_ = "general"
        if (present(symmetry)) then
            symmetry_ = to_lower(trim(symmetry)) 
        end if

        ! Determine field type based on matrix type
        field_type = MM_REAL

        catch: block
            ! Write header
            call write_mm_header(io, MM_COORDINATE, field_type, symmetry_, &
                                matrix%nrows, matrix%ncols, matrix%nnz, comment, stat, msg)
            if (stat /= 0) exit catch

            ! Write array format (column-major order)
            if(symmetry_ == MM_GENERAL) then
                do i = 1, matrix%nnz
                    write(io, '(ES24.16E3)', iostat=stat) matrix%data(i)
                    if (stat /= 0) then
                        msg = "Error writing array element (" // to_string(i) // ")"
                        exit catch
                    end if
                end do
            else
                ! For symmetric, skew-symmetric, hermitian matrices, only write the
                ! lower triangle (including diagonal)
                do i = 1, matrix%nnz
                    if(matrix%index(1,i) > matrix%index(2,i)) cycle
                    write(io, '(ES24.16E3)', iostat=stat) matrix%data(i)
                    if (stat /= 0) then
                        msg = "Error writing array element (" // to_string(i) // ")"
                        exit catch
                    end if
                end do
            end if
        end block catch
        
        close(io)

        if (present(iostat)) then
            iostat = stat
        else if (stat /= 0) then
            if (allocated(msg)) then
                call error_stop("Failed to save Matrix Market file '" // filename // "': " // msg)
            else
                call error_stop("Failed to save Matrix Market file '" // filename // "'")
            end if
        end if

        if (present(iomsg) .and. allocated(msg)) call move_alloc(msg, iomsg)
    end subroutine 
    module subroutine save_mm_coo_csp(filename, matrix, comment, symmetry, iostat, iomsg)
        !> Name of the Matrix Market file to save to
        character(len=*), intent(in) :: filename
        !> Matrix to be saved to the Matrix Market file
        type(COO_csp_type), intent(in) :: matrix
        !> Optional comment information
        character(len=*), intent(in), optional :: comment
        !> Symmetry type of the matrix (general, symmetric, skew-symmetric, hermitian)
        character(len=*), intent(in), optional :: symmetry
        !> Error status of saving, zero on success
        integer, intent(out), optional :: iostat
        !> Associated error message in case of non-zero status code
        character(len=:), allocatable, intent(out), optional :: iomsg

        integer :: io, stat, i, j, nnz
        character(len=:), allocatable :: msg
        character(len=32) :: field_type
        character(len=32) :: symmetry_
        real(sp) :: real_part, imag_part

        io = open(filename, "w", iostat=stat)
        if (stat /= 0) then
            if (present(iostat)) then
                iostat = stat
                if (present(iomsg)) iomsg = "Could not create file: " // filename
                return
            else
                call error_stop("Could not create file: " // filename)
            end if
        end if

        ! Determine symmetry type
        symmetry_ = "general"
        if (present(symmetry)) then
            symmetry_ = to_lower(trim(symmetry)) 
        end if

        ! Determine field type based on matrix type
        field_type = MM_COMPLEX

        catch: block
            ! Write header
            call write_mm_header(io, MM_COORDINATE, field_type, symmetry_, &
                                matrix%nrows, matrix%ncols, matrix%nnz, comment, stat, msg)
            if (stat /= 0) exit catch

            ! Write array format (column-major order)
            if(symmetry_ == MM_GENERAL) then
                do i = 1, matrix%nnz
                    real_part = real(matrix%data(i), kind=sp)
                    imag_part = aimag(matrix%data(i))
                    write(io, '(ES24.16E3,1X,ES24.16E3)', iostat=stat) real_part, imag_part
                    if (stat /= 0) then
                        msg = "Error writing array element (" // to_string(i) // ")"
                        exit catch
                    end if
                end do
            else
                ! For symmetric, skew-symmetric, hermitian matrices, only write the
                ! lower triangle (including diagonal)
                do i = 1, matrix%nnz
                    if(matrix%index(1,i) > matrix%index(2,i)) cycle
                    real_part = real(matrix%data(i), kind=sp)
                    imag_part = aimag(matrix%data(i))
                    write(io, '(ES24.16E3,1X,ES24.16E3)', iostat=stat) real_part, imag_part
                    if (stat /= 0) then
                        msg = "Error writing array element (" // to_string(i) // ")"
                        exit catch
                    end if
                end do
            end if
        end block catch
        
        close(io)

        if (present(iostat)) then
            iostat = stat
        else if (stat /= 0) then
            if (allocated(msg)) then
                call error_stop("Failed to save Matrix Market file '" // filename // "': " // msg)
            else
                call error_stop("Failed to save Matrix Market file '" // filename // "'")
            end if
        end if

        if (present(iomsg) .and. allocated(msg)) call move_alloc(msg, iomsg)
    end subroutine 
    module subroutine save_mm_coo_cdp(filename, matrix, comment, symmetry, iostat, iomsg)
        !> Name of the Matrix Market file to save to
        character(len=*), intent(in) :: filename
        !> Matrix to be saved to the Matrix Market file
        type(COO_cdp_type), intent(in) :: matrix
        !> Optional comment information
        character(len=*), intent(in), optional :: comment
        !> Symmetry type of the matrix (general, symmetric, skew-symmetric, hermitian)
        character(len=*), intent(in), optional :: symmetry
        !> Error status of saving, zero on success
        integer, intent(out), optional :: iostat
        !> Associated error message in case of non-zero status code
        character(len=:), allocatable, intent(out), optional :: iomsg

        integer :: io, stat, i, j, nnz
        character(len=:), allocatable :: msg
        character(len=32) :: field_type
        character(len=32) :: symmetry_
        real(dp) :: real_part, imag_part

        io = open(filename, "w", iostat=stat)
        if (stat /= 0) then
            if (present(iostat)) then
                iostat = stat
                if (present(iomsg)) iomsg = "Could not create file: " // filename
                return
            else
                call error_stop("Could not create file: " // filename)
            end if
        end if

        ! Determine symmetry type
        symmetry_ = "general"
        if (present(symmetry)) then
            symmetry_ = to_lower(trim(symmetry)) 
        end if

        ! Determine field type based on matrix type
        field_type = MM_COMPLEX

        catch: block
            ! Write header
            call write_mm_header(io, MM_COORDINATE, field_type, symmetry_, &
                                matrix%nrows, matrix%ncols, matrix%nnz, comment, stat, msg)
            if (stat /= 0) exit catch

            ! Write array format (column-major order)
            if(symmetry_ == MM_GENERAL) then
                do i = 1, matrix%nnz
                    real_part = real(matrix%data(i), kind=dp)
                    imag_part = aimag(matrix%data(i))
                    write(io, '(ES24.16E3,1X,ES24.16E3)', iostat=stat) real_part, imag_part
                    if (stat /= 0) then
                        msg = "Error writing array element (" // to_string(i) // ")"
                        exit catch
                    end if
                end do
            else
                ! For symmetric, skew-symmetric, hermitian matrices, only write the
                ! lower triangle (including diagonal)
                do i = 1, matrix%nnz
                    if(matrix%index(1,i) > matrix%index(2,i)) cycle
                    real_part = real(matrix%data(i), kind=dp)
                    imag_part = aimag(matrix%data(i))
                    write(io, '(ES24.16E3,1X,ES24.16E3)', iostat=stat) real_part, imag_part
                    if (stat /= 0) then
                        msg = "Error writing array element (" // to_string(i) // ")"
                        exit catch
                    end if
                end do
            end if
        end block catch
        
        close(io)

        if (present(iostat)) then
            iostat = stat
        else if (stat /= 0) then
            if (allocated(msg)) then
                call error_stop("Failed to save Matrix Market file '" // filename // "': " // msg)
            else
                call error_stop("Failed to save Matrix Market file '" // filename // "'")
            end if
        end if

        if (present(iomsg) .and. allocated(msg)) call move_alloc(msg, iomsg)
    end subroutine 

    !> Write Matrix Market header
    subroutine write_mm_header(io, format, field, symmetry, nrows, ncols, nnz, &
                              comment, iostat, iomsg)
        integer, intent(in) :: io
        character(len=*), intent(in) :: format, field, symmetry
        integer, intent(in) :: nrows, ncols, nnz
        character(len=*), intent(in), optional :: comment
        integer, intent(out) :: iostat
        character(len=:), allocatable, intent(out) :: iomsg

        integer :: stat
        character(len=*), parameter :: iso_date_fmt = '(I4.4,"-",I2.2,"-",I2.2)'
        integer :: date_values(8)

        iostat = 0

        ! Write banner line
        write(io, '(A)', iostat=stat) MM_BANNER // " " // MM_MATRIX // " " // &
              format // " " // field // " " // symmetry
        if (stat /= 0) then
            iostat = stat
            iomsg = "Error writing Matrix Market banner"
            return
        end if

        ! Write comments (including optional header_info and generation info)
        call date_and_time(values=date_values)
        write(io, '(A)', iostat=stat) "% Generated by Fortran stdlib on " // &
              to_string(date_values(1)) // "-" // &
              to_string(date_values(2)) // "-" // &
              to_string(date_values(3))
        if (stat /= 0) then
            iostat = stat
            iomsg = "Error writing comment line"
            return
        end if

        if (present(comment)) then
            if(len_trim(comment) > 0) then
                write(io, '(A)', iostat=stat) "% " // trim(comment)
                if (stat /= 0) then
                    iostat = stat
                    iomsg = "Error writing header info"
                    return
                end if
            end if
        end if

        ! Write size line
        if (format == MM_COORDINATE) then
            write(io, '(I0,1X,I0,1X,I0)', iostat=stat) nrows, ncols, nnz
        else
            write(io, '(I0,1X,I0)', iostat=stat) nrows, ncols
        end if
        
        if (stat /= 0) then
            iostat = stat
            iomsg = "Error writing matrix dimensions"
            return
        end if
    end subroutine write_mm_header

end submodule stdlib_io_mm_save