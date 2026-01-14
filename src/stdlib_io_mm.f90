! SPDX-Identifier: MIT


!> The Matrix Market (MM) format is a simple, human-readable, ASCII format for sparse 
!> and dense matrices. The format was developed at NIST (National Institute of Standards 
!> and Technology) for the Matrix Market, a repository of test matrices for use in 
!> comparative studies of algorithms for numerical linear algebra.
!>
!> For more information, see: https://math.nist.gov/MatrixMarket/formats.html
module stdlib_io_mm
    use stdlib_kinds, only : int8, int16, int32, int64, sp, dp, xdp, qp
    use stdlib_sparse_kinds
    implicit none
    private

    type, public :: mm_header_type
        integer :: object
        integer :: format
        integer :: qualifier
        integer :: symmetry
        character(len=1024), allocatable :: comments(:)
    end type mm_header_type

    !> Version: experimental
    !>
    !> Load a matrix from a Matrix Market file
    !> ([Specification](../page/specs/stdlib_io.html#load_mm))
    interface load_mm
        module subroutine load_mm_dense_sp(filename, matrix, iostat, iomsg)
            !> Name of the Matrix Market file to load from
            character(len=*), intent(in) :: filename
            !> Matrix to be loaded from the Matrix Market file
            real(sp), allocatable, intent(out) :: matrix(:,:)
            !> Error status of loading, zero on success
            integer, intent(out), optional :: iostat
            !> Associated error message in case of non-zero status code
            character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine
        module subroutine load_mm_dense_dp(filename, matrix, iostat, iomsg)
            !> Name of the Matrix Market file to load from
            character(len=*), intent(in) :: filename
            !> Matrix to be loaded from the Matrix Market file
            real(dp), allocatable, intent(out) :: matrix(:,:)
            !> Error status of loading, zero on success
            integer, intent(out), optional :: iostat
            !> Associated error message in case of non-zero status code
            character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine
        module subroutine load_mm_dense_csp(filename, matrix, iostat, iomsg)
            !> Name of the Matrix Market file to load from
            character(len=*), intent(in) :: filename
            !> Matrix to be loaded from the Matrix Market file
            complex(sp), allocatable, intent(out) :: matrix(:,:)
            !> Error status of loading, zero on success
            integer, intent(out), optional :: iostat
            !> Associated error message in case of non-zero status code
            character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine
        module subroutine load_mm_dense_cdp(filename, matrix, iostat, iomsg)
            !> Name of the Matrix Market file to load from
            character(len=*), intent(in) :: filename
            !> Matrix to be loaded from the Matrix Market file
            complex(dp), allocatable, intent(out) :: matrix(:,:)
            !> Error status of loading, zero on success
            integer, intent(out), optional :: iostat
            !> Associated error message in case of non-zero status code
            character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine
        module subroutine load_mm_coo_sp(filename, matrix, iostat, iomsg)
            !> Name of the Matrix Market file to load from
            character(len=*), intent(in) :: filename
            !> Matrix to be loaded from the Matrix Market file
            type(COO_sp_type), intent(out) :: matrix
            !> Error status of loading, zero on success
            integer, intent(out), optional :: iostat
            !> Associated error message in case of non-zero status code
            character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine
        module subroutine load_mm_coo_dp(filename, matrix, iostat, iomsg)
            !> Name of the Matrix Market file to load from
            character(len=*), intent(in) :: filename
            !> Matrix to be loaded from the Matrix Market file
            type(COO_dp_type), intent(out) :: matrix
            !> Error status of loading, zero on success
            integer, intent(out), optional :: iostat
            !> Associated error message in case of non-zero status code
            character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine
        module subroutine load_mm_coo_csp(filename, matrix, iostat, iomsg)
            !> Name of the Matrix Market file to load from
            character(len=*), intent(in) :: filename
            !> Matrix to be loaded from the Matrix Market file
            type(COO_csp_type), intent(out) :: matrix
            !> Error status of loading, zero on success
            integer, intent(out), optional :: iostat
            !> Associated error message in case of non-zero status code
            character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine
        module subroutine load_mm_coo_cdp(filename, matrix, iostat, iomsg)
            !> Name of the Matrix Market file to load from
            character(len=*), intent(in) :: filename
            !> Matrix to be loaded from the Matrix Market file
            type(COO_cdp_type), intent(out) :: matrix
            !> Error status of loading, zero on success
            integer, intent(out), optional :: iostat
            !> Associated error message in case of non-zero status code
            character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine
    end interface
    public :: load_mm

    !> Version: experimental
    !>
    !> Save a matrix to a Matrix Market file
    !> ([Specification](../page/specs/stdlib_io.html#save_mm))
    interface save_mm
        module subroutine save_mm_dense_sp(filename, matrix, comment, symmetry, iostat, iomsg)
        character(len=*), intent(in) :: filename
        real(sp), intent(in) :: matrix(:,:)
        character(len=*), intent(in), optional :: comment
        character(len=*), intent(in), optional :: symmetry
        integer, intent(out), optional :: iostat
        character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine 
        module subroutine save_mm_dense_dp(filename, matrix, comment, symmetry, iostat, iomsg)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: matrix(:,:)
        character(len=*), intent(in), optional :: comment
        character(len=*), intent(in), optional :: symmetry
        integer, intent(out), optional :: iostat
        character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine 
        module subroutine save_mm_dense_csp(filename, matrix, comment, symmetry, iostat, iomsg)
        character(len=*), intent(in) :: filename
        complex(sp), intent(in) :: matrix(:,:)
        character(len=*), intent(in), optional :: comment
        character(len=*), intent(in), optional :: symmetry
        integer, intent(out), optional :: iostat
        character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine 
        module subroutine save_mm_dense_cdp(filename, matrix, comment, symmetry, iostat, iomsg)
        character(len=*), intent(in) :: filename
        complex(dp), intent(in) :: matrix(:,:)
        character(len=*), intent(in), optional :: comment
        character(len=*), intent(in), optional :: symmetry
        integer, intent(out), optional :: iostat
        character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine 

        module subroutine save_coo_sp(filename, matrix, comment, symmetry, iostat, iomsg)
        character(len=*), intent(in) :: filename
        type(COO_sp_type), intent(in) :: matrix
        character(len=*), intent(in), optional :: comment
        character(len=*), intent(in), optional :: symmetry
        integer, intent(out), optional :: iostat
        character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine 
        module subroutine save_coo_dp(filename, matrix, comment, symmetry, iostat, iomsg)
        character(len=*), intent(in) :: filename
        type(COO_dp_type), intent(in) :: matrix
        character(len=*), intent(in), optional :: comment
        character(len=*), intent(in), optional :: symmetry
        integer, intent(out), optional :: iostat
        character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine 
        module subroutine save_coo_csp(filename, matrix, comment, symmetry, iostat, iomsg)
        character(len=*), intent(in) :: filename
        type(COO_csp_type), intent(in) :: matrix
        character(len=*), intent(in), optional :: comment
        character(len=*), intent(in), optional :: symmetry
        integer, intent(out), optional :: iostat
        character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine 
        module subroutine save_coo_cdp(filename, matrix, comment, symmetry, iostat, iomsg)
        character(len=*), intent(in) :: filename
        type(COO_cdp_type), intent(in) :: matrix
        character(len=*), intent(in), optional :: comment
        character(len=*), intent(in), optional :: symmetry
        integer, intent(out), optional :: iostat
        character(len=:), allocatable, intent(out), optional :: iomsg
        end subroutine 
    end interface save_mm
    public :: save_mm

end module stdlib_io_mm