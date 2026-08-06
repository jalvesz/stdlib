! Process example: Plot data with Gnuplot via stdin pipe
!
! Demonstrates passing a large stdin payload to an external process using
! `stdlib_system%run`.  The data are assembled in memory and sent to gnuplot
! through the `stdin` argument — no manual temporary file management is needed
! by the caller (the module handles the required redirection internally).
program example_gnuplot
    use stdlib_system, only: run, is_completed, process_type
    use stdlib_strings, only: to_string
    use iso_c_binding, only: c_new_line
    implicit none

    integer,          parameter :: N     = 800    ! number of data points
    real,             parameter :: XSTEP = 0.1    ! step size

    integer :: i
    real    :: x(N), y(N)

    ! Generate plotting data
    do i = 1, N
        x(i) = -30.0 + real(i - 1) * XSTEP
        y(i) = sin(x(i) * 20.0) * atan(x(i))
    end do

    ! Plot data as a line chart using gnuplot
    call plot(x, y, 'Gnuplot from Fortran')

contains

    !> Open a pipe to Gnuplot, write settings, and send x/y data via stdin.
    subroutine plot(x, y, title)
        real,             intent(in) :: x(:)   !! x values
        real,             intent(in) :: y(:)   !! y values
        character(len=*), intent(in) :: title  !! plot and window title

        type(process_type)        :: p
        character(:), allocatable :: stdin
        integer                   :: i

        ! Build gnuplot command script in memory
        stdin = 'reset session'                          // c_new_line
        stdin = stdin // 'set terminal dumb'             // c_new_line  ! terminal that needs no display
        stdin = stdin // 'set title "' // trim(title) // '"' // c_new_line
        stdin = stdin // 'set grid'                      // c_new_line
        stdin = stdin // 'unset key'                     // c_new_line
        stdin = stdin // 'plot "-" using 1:2 with lines' // c_new_line
        do i = 1, size(x)
            stdin = stdin // to_string(x(i), '(f12.6)') // ' ' &
                          // to_string(y(i), '(f12.6)') // c_new_line
        end do
        stdin = stdin // 'e' // c_new_line

        ! Send the full script to gnuplot via stdin (no manual tmp files needed)
        p = run("gnuplot", stdin=stdin, want_stdout=.true., want_stderr=.true.)

        if (is_completed(p)) then
            print *, "gnuplot finished. Output:"
            if (allocated(p%stdout)) print *, p%stdout
        else
            print *, "gnuplot did not complete as expected."
            if (allocated(p%stderr)) print *, "stderr:", p%stderr
        end if

    end subroutine plot

end program example_gnuplot
