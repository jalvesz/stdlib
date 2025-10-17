! SPDX-Identifier: MIT


!> Implementation of saving multidimensional arrays to Matrix Market files
submodule (stdlib_io_mm) stdlib_io_mm_save
    use stdlib_error, only : error_stop
    use stdlib_strings, only : to_string
    use stdlib_io, only : open
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

    module subroutine save_mm_dense_sp(filename, matrix, header_info, iostat, iomsg)
        !> Name of the Matrix Market file to save to
        character(len=*), intent(in) :: filename
        !> Matrix to be saved to the Matrix Market file
        real(sp), intent(in) :: matrix(:,:)
        !> Optional header information (comments, format preference)
        character(len=*), intent(in), optional :: header_info
        !> Error status of saving, zero on success
        integer, intent(out), optional :: iostat
        !> Associated error message in case of non-zero status code
        character(len=:), allocatable, intent(out), optional :: iomsg

        integer :: io, stat, i, j, nnz
        character(len=:), allocatable :: msg
        character(len=32) :: field_type, symmetry_type
        logical :: save_as_coordinate

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

        catch: block
            ! Determine field type based on matrix type
            field_type = MM_REAL

            ! For now, assume general symmetry (could be enhanced to detect symmetry)
            symmetry_type = MM_GENERAL

            ! Count non-zero elements to decide format
            nnz = 0
            do j = 1, size(matrix, 2)
                do i = 1, size(matrix, 1)
                    if (matrix(i, j) /= 0) nnz = nnz + 1
                end do
            end do

            ! Decide format based on sparsity (save as coordinate if < 50% non-zero)
            save_as_coordinate = (real(nnz) / real(size(matrix)) < 0.5)
            
            ! Allow override via header_info (simple implementation)
            if (present(header_info)) then
                if (index(header_info, "array") > 0) save_as_coordinate = .false.
                if (index(header_info, "coordinate") > 0) save_as_coordinate = .true.
            end if

            ! Write header
            if (save_as_coordinate) then
                call write_mm_header(io, MM_COORDINATE, field_type, symmetry_type, &
                                   size(matrix, 1), size(matrix, 2), nnz, header_info, stat, msg)
            else
                call write_mm_header(io, MM_ARRAY, field_type, symmetry_type, &
                                   size(matrix, 1), size(matrix, 2), 0, header_info, stat, msg)
            end if
            if (stat /= 0) exit catch

            ! Write data
            if (save_as_coordinate) then
                ! Write coordinate format
                do j = 1, size(matrix, 2)
                    do i = 1, size(matrix, 1)
                        if (matrix(i, j) /= 0) then
                            write(io, '(I0,1X,I0,1X,ES24.16E3)', iostat=stat) i, j, matrix(i, j)
                            if (stat /= 0) then
                                msg = "Error writing coordinate entry (" // &
                                      to_string(i) // "," // to_string(j) // ")"
                                exit catch
                            end if
                        end if
                    end do
                end do
            else
                ! Write array format (column-major order)
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
    module subroutine save_mm_dense_dp(filename, matrix, header_info, iostat, iomsg)
        !> Name of the Matrix Market file to save to
        character(len=*), intent(in) :: filename
        !> Matrix to be saved to the Matrix Market file
        real(dp), intent(in) :: matrix(:,:)
        !> Optional header information (comments, format preference)
        character(len=*), intent(in), optional :: header_info
        !> Error status of saving, zero on success
        integer, intent(out), optional :: iostat
        !> Associated error message in case of non-zero status code
        character(len=:), allocatable, intent(out), optional :: iomsg

        integer :: io, stat, i, j, nnz
        character(len=:), allocatable :: msg
        character(len=32) :: field_type, symmetry_type
        logical :: save_as_coordinate

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

        catch: block
            ! Determine field type based on matrix type
            field_type = MM_REAL

            ! For now, assume general symmetry (could be enhanced to detect symmetry)
            symmetry_type = MM_GENERAL

            ! Count non-zero elements to decide format
            nnz = 0
            do j = 1, size(matrix, 2)
                do i = 1, size(matrix, 1)
                    if (matrix(i, j) /= 0) nnz = nnz + 1
                end do
            end do

            ! Decide format based on sparsity (save as coordinate if < 50% non-zero)
            save_as_coordinate = (real(nnz) / real(size(matrix)) < 0.5)
            
            ! Allow override via header_info (simple implementation)
            if (present(header_info)) then
                if (index(header_info, "array") > 0) save_as_coordinate = .false.
                if (index(header_info, "coordinate") > 0) save_as_coordinate = .true.
            end if

            ! Write header
            if (save_as_coordinate) then
                call write_mm_header(io, MM_COORDINATE, field_type, symmetry_type, &
                                   size(matrix, 1), size(matrix, 2), nnz, header_info, stat, msg)
            else
                call write_mm_header(io, MM_ARRAY, field_type, symmetry_type, &
                                   size(matrix, 1), size(matrix, 2), 0, header_info, stat, msg)
            end if
            if (stat /= 0) exit catch

            ! Write data
            if (save_as_coordinate) then
                ! Write coordinate format
                do j = 1, size(matrix, 2)
                    do i = 1, size(matrix, 1)
                        if (matrix(i, j) /= 0) then
                            write(io, '(I0,1X,I0,1X,ES24.16E3)', iostat=stat) i, j, matrix(i, j)
                            if (stat /= 0) then
                                msg = "Error writing coordinate entry (" // &
                                      to_string(i) // "," // to_string(j) // ")"
                                exit catch
                            end if
                        end if
                    end do
                end do
            else
                ! Write array format (column-major order)
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
    module subroutine save_mm_dense_csp(filename, matrix, header_info, iostat, iomsg)
        !> Name of the Matrix Market file to save to
        character(len=*), intent(in) :: filename
        !> Matrix to be saved to the Matrix Market file
        complex(sp), intent(in) :: matrix(:,:)
        !> Optional header information (comments, format preference)
        character(len=*), intent(in), optional :: header_info
        !> Error status of saving, zero on success
        integer, intent(out), optional :: iostat
        !> Associated error message in case of non-zero status code
        character(len=:), allocatable, intent(out), optional :: iomsg

        integer :: io, stat, i, j, nnz
        character(len=:), allocatable :: msg
        character(len=32) :: field_type, symmetry_type
        logical :: save_as_coordinate
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

        catch: block
            ! Determine field type based on matrix type
            field_type = MM_COMPLEX

            ! For now, assume general symmetry (could be enhanced to detect symmetry)
            symmetry_type = MM_GENERAL

            ! Count non-zero elements to decide format
            nnz = 0
            do j = 1, size(matrix, 2)
                do i = 1, size(matrix, 1)
                    if (abs(matrix(i, j)) /= 0) nnz = nnz + 1
                end do
            end do

            ! Decide format based on sparsity (save as coordinate if < 50% non-zero)
            save_as_coordinate = (real(nnz) / real(size(matrix)) < 0.5)
            
            ! Allow override via header_info (simple implementation)
            if (present(header_info)) then
                if (index(header_info, "array") > 0) save_as_coordinate = .false.
                if (index(header_info, "coordinate") > 0) save_as_coordinate = .true.
            end if

            ! Write header
            if (save_as_coordinate) then
                call write_mm_header(io, MM_COORDINATE, field_type, symmetry_type, &
                                   size(matrix, 1), size(matrix, 2), nnz, header_info, stat, msg)
            else
                call write_mm_header(io, MM_ARRAY, field_type, symmetry_type, &
                                   size(matrix, 1), size(matrix, 2), 0, header_info, stat, msg)
            end if
            if (stat /= 0) exit catch

            ! Write data
            if (save_as_coordinate) then
                ! Write coordinate format
                do j = 1, size(matrix, 2)
                    do i = 1, size(matrix, 1)
                        if (abs(matrix(i, j)) /= 0) then
                            real_part = real(matrix(i, j), kind=sp)
                            imag_part = aimag(matrix(i, j))
                            write(io, '(I0,1X,I0,1X,ES24.16E3,1X,ES24.16E3)', iostat=stat) &
                                  i, j, real_part, imag_part
                            if (stat /= 0) then
                                msg = "Error writing coordinate entry (" // &
                                      to_string(i) // "," // to_string(j) // ")"
                                exit catch
                            end if
                        end if
                    end do
                end do
            else
                ! Write array format (column-major order)
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
    module subroutine save_mm_dense_cdp(filename, matrix, header_info, iostat, iomsg)
        !> Name of the Matrix Market file to save to
        character(len=*), intent(in) :: filename
        !> Matrix to be saved to the Matrix Market file
        complex(dp), intent(in) :: matrix(:,:)
        !> Optional header information (comments, format preference)
        character(len=*), intent(in), optional :: header_info
        !> Error status of saving, zero on success
        integer, intent(out), optional :: iostat
        !> Associated error message in case of non-zero status code
        character(len=:), allocatable, intent(out), optional :: iomsg

        integer :: io, stat, i, j, nnz
        character(len=:), allocatable :: msg
        character(len=32) :: field_type, symmetry_type
        logical :: save_as_coordinate
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

        catch: block
            ! Determine field type based on matrix type
            field_type = MM_COMPLEX

            ! For now, assume general symmetry (could be enhanced to detect symmetry)
            symmetry_type = MM_GENERAL

            ! Count non-zero elements to decide format
            nnz = 0
            do j = 1, size(matrix, 2)
                do i = 1, size(matrix, 1)
                    if (abs(matrix(i, j)) /= 0) nnz = nnz + 1
                end do
            end do

            ! Decide format based on sparsity (save as coordinate if < 50% non-zero)
            save_as_coordinate = (real(nnz) / real(size(matrix)) < 0.5)
            
            ! Allow override via header_info (simple implementation)
            if (present(header_info)) then
                if (index(header_info, "array") > 0) save_as_coordinate = .false.
                if (index(header_info, "coordinate") > 0) save_as_coordinate = .true.
            end if

            ! Write header
            if (save_as_coordinate) then
                call write_mm_header(io, MM_COORDINATE, field_type, symmetry_type, &
                                   size(matrix, 1), size(matrix, 2), nnz, header_info, stat, msg)
            else
                call write_mm_header(io, MM_ARRAY, field_type, symmetry_type, &
                                   size(matrix, 1), size(matrix, 2), 0, header_info, stat, msg)
            end if
            if (stat /= 0) exit catch

            ! Write data
            if (save_as_coordinate) then
                ! Write coordinate format
                do j = 1, size(matrix, 2)
                    do i = 1, size(matrix, 1)
                        if (abs(matrix(i, j)) /= 0) then
                            real_part = real(matrix(i, j), kind=dp)
                            imag_part = aimag(matrix(i, j))
                            write(io, '(I0,1X,I0,1X,ES24.16E3,1X,ES24.16E3)', iostat=stat) &
                                  i, j, real_part, imag_part
                            if (stat /= 0) then
                                msg = "Error writing coordinate entry (" // &
                                      to_string(i) // "," // to_string(j) // ")"
                                exit catch
                            end if
                        end if
                    end do
                end do
            else
                ! Write array format (column-major order)
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
                              header_info, iostat, iomsg)
        integer, intent(in) :: io
        character(len=*), intent(in) :: format, field, symmetry
        integer, intent(in) :: nrows, ncols, nnz
        character(len=*), intent(in), optional :: header_info
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

        if (present(header_info) .and. len_trim(header_info) > 0) then
            write(io, '(A)', iostat=stat) "% " // trim(header_info)
            if (stat /= 0) then
                iostat = stat
                iomsg = "Error writing header info"
                return
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