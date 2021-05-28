program test_sorting

    use, intrinsic :: iso_fortran_env, only: &
        compiler_version
    use stdlib_kinds, only: int32, int64, dp, sp
    use stdlib_sorting
    use stdlib_string_type, only: string_type, assignment(=), operator(>), operator(<), &
        write(formatted)
    use stdlib_error, only: check

    implicit none

    integer(int32), parameter :: test_size = 2_int32**20
    integer(int32), parameter :: char_size = 26**4
    integer(int32), parameter :: string_size = 26**3
    integer(int32), parameter :: block_size = test_size/6
    integer, parameter        :: repeat = 8

    integer(int32) ::             &
        blocks(0:test_size-1),    &
        decrease(0:test_size-1),  &
        identical(0:test_size-1), &
        increase(0:test_size-1),  &
        rand0(0:test_size-1),     &
        rand1(0:test_size-1),     &
        rand2(0:test_size-1),     &
        rand3(0:test_size-1),     &
        rand10(0:test_size-1)
    character(len=4) ::               &
        char_decrease(0:char_size-1), &
        char_increase(0:char_size-1), &
        char_rand(0:char_size-1)
    type(string_type) ::                  &
        string_decrease(0:string_size-1), &
        string_increase(0:string_size-1), &
        string_rand(0:string_size-1)

    integer(int32)          :: dummy(0:test_size-1)
    character(len=4)        :: char_dummy(0:char_size-1)
    type(string_type)       :: string_dummy(0:string_size-1)
    integer(int_size)       :: index(0:test_size-1)
    integer(int32)          :: work(0:test_size/2-1)
    character(len=4)        :: char_work(0:char_size/2-1)
    type(string_type)       :: string_work(0:string_size/2-1)
    integer(int_size)       :: iwork(0:test_size/2-1)
    integer                 :: count, i, index1, index2, j, k, l, temp
    real(sp)                :: arand, brand
    character(*), parameter :: filename = 'test_sorting.txt'
    integer                 :: lun
    character(len=4)        :: char_temp
    type(string_type)       :: string_temp
    logical                 :: ltest, ldummy

! Create the test arrays
    identical(:) = 10
    do i=0, test_size-1
        increase(i) = i
        decrease(i) = test_size - 1 - i
        call random_number( arand )
        rand0(i) = int( floor( 4 * arand * test_size ), kind=int32 )
        rand1(i) = int( floor( arand * test_size / 4 ), kind=int32 )
    end do
    blocks(:) = increase(:)
    blocks(0:block_size-1) = increase(4*block_size:5*block_size-1)
    blocks(block_size:2*block_size-1) = increase(0:block_size-1)
    blocks(2*block_size:3*block_size-1) = increase(2*block_size:3*block_size-1)
    blocks(3*block_size:4*block_size-1) = increase(block_size:2*block_size-1)
    blocks(4*block_size:5*block_size-1) = increase(3*block_size:4*block_size-1)
    rand2(:) = increase(:)
    do i=0, test_size-1
        call random_number( arand )
        index1 = int( floor( arand * test_size ), kind=int32 )
        temp = rand2(i)
        rand2(i) = rand2(index1)
        rand2(index1) = temp
    end do
    rand3(:) = increase(:)
    do i=0, 2
        call random_number( arand )
        call random_number( brand )
        index1 = int( floor( arand * test_size ), kind=int32 )
        index2 = int( floor( brand * test_size ), kind=int32 )
        temp = rand3(index1)
        rand3(index1) = rand3(index2)
        rand3(index2) = temp
    end do
    rand10(:) = increase(:)
    do i=test_size-10, test_size-1
        call random_number( arand )
        rand10(i) = int( floor( arand * test_size ), kind=int32 )
    end do

    count = 0
    do i=0, 25
        do j=0, 25
            do k=0, 25
                do l=0, 25
                    char_increase(count) = achar(97+i) // achar(97+j) // &
                        achar(97+k) // achar(97+l)
                    count = count + 1
                end do
            end do
        end do
    end do

    do i=0, char_size-1
        char_decrease(char_size-1-i) = char_increase(i)
    end do

    char_rand(:) = char_increase(:)
    do i=0, char_size-1
        call random_number( arand )
        index1 = int( floor( arand * char_size ), kind=int32 )
        char_temp = char_rand(i)
        char_rand(i) = char_rand(index1)
        char_rand(index1) = char_temp
    end do

    count = 0
    do i=0, 25
        do j=0, 25
            do k=0, 25
                string_increase(count) = achar(97+i) // achar(97+j) // &
                    achar(97+k)
                count = count + 1
            end do
        end do
    end do

    do i=0, string_size-1
        string_decrease(string_size - 1 - i) = char_increase(i)
    end do

    string_rand(:) = string_increase(:)
    do i=0, string_size-1
        call random_number( arand )
        index1 = int( floor( arand * string_size ), kind=int32 )
        string_temp = string_rand(i)
        string_rand(i) = string_rand(index1)
        string_rand(index1) = string_temp
    end do

! Create and intialize file to report the results of the sortings
    open( newunit=lun, file=filename, access='sequential', action='write', &
        form='formatted', status='replace' )
    write( lun, '(a)' ) trim(compiler_version())
    write( lun, * )
    write( lun, '("|    Type     | Elements |    Array Name   |    Method ' // &
        '  |  Time (s) |")' )
    write( lun, '("|-------------|----------|-----------------|-----------' // &
        '--|-----------|")' )

! test the sorting routines on the test arrays
    ltest = .true.

    call test_int_ord_sorts( ldummy );    ltest = (ltest .and. ldummy)

    call test_char_ord_sorts(ldummy );    ltest = (ltest .and. ldummy)

    call test_string_ord_sorts( ldummy ); ltest = (ltest .and. ldummy)

    call test_int_sorts( ldummy );        ltest = (ltest .and. ldummy)

    call test_char_sorts( ldummy );       ltest = (ltest .and. ldummy)

    call test_string_sorts( ldummy );     ltest = (ltest .and. ldummy)

    call test_int_sort_indexes( )

    call test_char_sort_indexes( )

    call test_string_sort_indexes( )


    call check(ltest)

contains

    subroutine test_int_ord_sorts( ltest )
        logical, intent(out) :: ltest

        logical :: ldummy

        ltest = .true.

        call test_int_ord_sort( blocks, "Blocks", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_ord_sort( decrease, "Decreasing", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_ord_sort( identical, "Identical", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_ord_sort( increase, "Increasing", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_ord_sort( rand1, "Random dense", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_ord_sort( rand2, "Random order", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_ord_sort( rand0, "Random sparse", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_ord_sort( rand3, "Random 3", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_ord_sort( rand10, "Random 10", ldummy )
        ltest = (ltest .and. ldummy)

    end subroutine test_int_ord_sorts


    subroutine test_int_ord_sort( a, a_name, ltest )
        integer(int32), intent(in) :: a(:)
        character(*), intent(in)   :: a_name 
        logical, intent(out)       :: ltest

        integer(int64) :: t0, t1, tdiff
        real(dp)       :: rate
        integer(int64) :: i
        logical :: valid

        ltest = .true.

        tdiff = 0
        do i = 1, repeat
            dummy = a
            call system_clock( t0, rate )
            call ord_sort( dummy, work )
            call system_clock( t1, rate )
            tdiff = tdiff + t1 - t0
        end do
        tdiff = tdiff/repeat

        call verify_sort( dummy, valid, i )
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "ORD_SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a12, 2i7)') 'dummy(i-1:i) = ', dummy(i-1:i)
        end if
        write( lun, '("|     Integer |", 1x, i7, 2x, "|", 1x, a15, " |", ' // &
            'a12, " |",  F10.5, " |" )' ) &
            test_size, a_name, "Ord_Sort", tdiff/rate

        !reverse
        dummy = a
        call ord_sort( dummy, work, reverse = .true.)
        call verify_reverse_sort( dummy, valid, i )
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "reverse + work ORD_SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a12, 2i7)') 'dummy(i-1:i) = ', dummy(i-1:i)
        end if

        dummy = a
        call ord_sort( dummy, reverse = .true.)
        call verify_reverse_sort( dummy, valid, i )
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "reverse ORD_SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a12, 2i7)') 'dummy(i-1:i) = ', dummy(i-1:i)
        end if

    end subroutine test_int_ord_sort

    subroutine test_char_ord_sorts( ltest )
        logical, intent(out) :: ltest

        logical :: ldummy

        ltest = .true.

        call test_char_ord_sort( char_decrease, "Char. Decrease", ldummy )
        ltest = (ltest .and. ldummy)
        call test_char_ord_sort( char_increase, "Char. Increase", ldummy )
        ltest = (ltest .and. ldummy)
        call test_char_ord_sort( char_rand, "Char. Random", ldummy )
        ltest = (ltest .and. ldummy)

    end subroutine test_char_ord_sorts

    subroutine test_char_ord_sort( a, a_name, ltest )
        character(len=4), intent(in) :: a(0:)
        character(*), intent(in) :: a_name
        logical, intent(out) :: ltest

        integer(int64) :: t0, t1, tdiff
        real(dp)       :: rate
        integer(int64) :: i
        logical        :: valid

        tdiff = 0
        do i = 1, repeat
            char_dummy = a
            call system_clock( t0, rate )
            call ord_sort( char_dummy, char_work )
            call system_clock( t1, rate )
            tdiff = tdiff + t1 - t0
        end do
        tdiff = tdiff/repeat

        call verify_char_sort( char_dummy, valid, i )
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "ORD_SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a, 2(1x,a4))') 'char_dummy(i-1:i) = ', char_dummy(i-1:i)
        end if
        write( lun, '("|   Character |", 1x, i7, 2x, "|", 1x, a15, " |", ' // &
            'a12, " |",  F10.5, " |" )' ) &
            char_size, a_name, "Ord_Sort", tdiff/rate

        !reverse
        char_dummy = a
        call ord_sort( char_dummy, char_work, reverse = .true. )

        call verify_char_reverse_sort( char_dummy, valid, i )
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "reverse + work ORD_SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a, 2(1x,a4))') 'char_dummy(i-1:i) = ', char_dummy(i-1:i)
        end if

        char_dummy = a
        call ord_sort( char_dummy, reverse = .true. )

        call verify_char_reverse_sort( char_dummy, valid, i )
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "reverse + work ORD_SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a, 2(1x,a4))') 'char_dummy(i-1:i) = ', char_dummy(i-1:i)
        end if

    end subroutine test_char_ord_sort

    subroutine test_string_ord_sorts( ltest )
        logical, intent(out) :: ltest

        logical:: ldummy

        ltest = .true.

        call test_string_ord_sort( string_decrease, "String Decrease", ldummy )
        ltest = (ltest .and. ldummy)

        call test_string_ord_sort( string_increase, "String Increase", ldummy )
        ltest = (ltest .and. ldummy)

        call test_string_ord_sort( string_rand, "String Random" , ldummy)
        ltest = (ltest .and. ldummy)

    end subroutine test_string_ord_sorts

    subroutine test_string_ord_sort( a, a_name, ltest )
        type(string_type), intent(in) :: a(0:)
        character(*), intent(in)      :: a_name
        logical, intent(out)          :: ltest

        integer(int64) :: t0, t1, tdiff
        real(dp)       :: rate
        integer(int64) :: i
        logical        :: valid

        ltest = .true.

        tdiff = 0
        do i = 1, repeat
            string_dummy = a
            call system_clock( t0, rate )
            call ord_sort( string_dummy, string_work )
            call system_clock( t1, rate )
            tdiff = tdiff + t1 - t0
        end do
        tdiff = tdiff/repeat

        call verify_string_sort( string_dummy, valid, i )
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "ORD_SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a, 2(1x,a))') 'string_dummy(i-1:i) = ', &
                string_dummy(i-1:i)
        end if
        write( lun, '("| String_type |", 1x, i7, 2x, "|", 1x, a15, " |", ' // &
            'a12, " |",  F10.5, " |" )' ) &
            string_size, a_name, "Ord_Sort", tdiff/rate

        !reverse
        string_dummy = a
        call ord_sort( string_dummy, string_work, reverse = .true. )

        call verify_string_reverse_sort( string_dummy, valid, i )
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "reverse + work ORD_SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a, 2(1x,a))') 'string_dummy(i-1:i) = ', &
                string_dummy(i-1:i)
        end if
 
        string_dummy = a
        call ord_sort( string_dummy, reverse = .true. )

        call verify_string_reverse_sort( string_dummy, valid, i )
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "reverse ORD_SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a, 2(1x,a))') 'string_dummy(i-1:i) = ', &
                string_dummy(i-1:i)
        end if

    end subroutine test_string_ord_sort


    subroutine test_int_sorts( ltest )
        logical, intent(out) :: ltest

        logical :: ldummy

        ltest = .true.

        call test_int_sort( blocks, "Blocks", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_sort( decrease, "Decreasing", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_sort( identical, "Identical", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_sort( increase, "Increasing", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_sort( rand1, "Random dense", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_sort( rand2, "Random order", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_sort( rand0, "Random sparse", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_sort( rand3, "Random 3", ldummy )
        ltest = (ltest .and. ldummy)
        call test_int_sort( rand10, "Random 10", ldummy )
        ltest = (ltest .and. ldummy)

    end subroutine test_int_sorts

    subroutine test_int_sort( a, a_name, ltest )
        integer(int32), intent(in) :: a(:)
        character(*), intent(in)   :: a_name
        logical, intent(out) :: ltest

        integer(int64) :: t0, t1, tdiff
        real(dp)       :: rate
        integer(int64) :: i
        logical        :: valid

        ltest = .true.

        tdiff = 0
        do i = 1, repeat
            dummy = a
            call system_clock( t0, rate )
            call sort( dummy )
            call system_clock( t1, rate )
            tdiff = tdiff + t1 - t0
        end do
        tdiff = tdiff/repeat


        call verify_sort( dummy, valid, i )
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a12, 2i7)') 'dummy(i-1:i) = ', dummy(i-1:i)
        end if
        write( lun, '("|     Integer |", 1x, i7, 2x, "|", 1x, a15, " |", ' // &
            'a12, " |",  F10.5, " |" )' ) &
            test_size, a_name, "Sort", tdiff/rate


        ! reverse
        dummy = a
        call sort( dummy, .true.)
        call verify_reverse_sort(dummy, valid, i)
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "reverse SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a12, 2i7)') 'dummy(i-1:i) = ', dummy(i-1:i)
        end if

    end subroutine test_int_sort

    subroutine test_char_sorts( ltest )
        logical, intent(out) :: ltest

        logical :: ldummy

        call test_char_sort( char_decrease, "Char. Decrease", ldummy )
        ltest = (ltest .and. ldummy)

        call test_char_sort( char_increase, "Char. Increase", ldummy )
        ltest = (ltest .and. ldummy)

        call test_char_sort( char_rand, "Char. Random", ldummy )
        ltest = (ltest .and. ldummy)

    end subroutine test_char_sorts

    subroutine test_char_sort( a, a_name, ltest )
        character(len=4), intent(in) :: a(0:)
        character(*), intent(in) :: a_name
        logical, intent(out) :: ltest

        integer(int64) :: t0, t1, tdiff
        real(dp)       :: rate
        integer(int64) :: i
        logical        :: valid

        tdiff = 0
        do i = 1, repeat
            char_dummy = a
            call system_clock( t0, rate )
            call sort( char_dummy )
            call system_clock( t1, rate )
            tdiff = tdiff + t1 - t0
        end do
        tdiff = tdiff/repeat

        call verify_char_sort( char_dummy, valid, i )
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a17, 2(1x,a4))') 'char_dummy(i-1:i) = ', char_dummy(i-1:i)
        end if
        write( lun, '("|   Character |", 1x, i7, 2x, "|", 1x, a15, " |", ' // &
            'a12, " |",  F10.5, " |" )' ) &
            char_size, a_name, "Sort", tdiff/rate

        !reverse
        char_dummy = a
        call sort( char_dummy, .true.)
        call verify_char_reverse_sort( char_dummy, valid, i )
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "reverse SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a17, 2(1x,a4))') 'char_dummy(i-1:i) = ', char_dummy(i-1:i)
        end if
   
    end subroutine test_char_sort

    subroutine test_string_sorts( ltest )
        logical, intent(out) :: ltest

        logical :: ldummy

        ltest = .true.

        call test_string_sort( string_decrease, "String Decrease", ldummy )
        ltest = (ltest .and. ldummy)
        call test_string_sort( string_increase, "String Increase", ldummy )
        ltest = (ltest .and. ldummy)
        call test_string_sort( string_rand, "String Random", ldummy )
        ltest = (ltest .and. ldummy)

    end subroutine test_string_sorts

    subroutine test_string_sort( a, a_name, ltest )
        type(string_type), intent(in) :: a(0:)
        character(*), intent(in) :: a_name
        logical, intent(out) :: ltest

        integer(int64) :: t0, t1, tdiff
        real(dp)       :: rate
        integer(int64) :: i
        logical        :: valid

        ltest = .true.

        tdiff = 0
        do i = 1, repeat
            string_dummy = a
            call system_clock( t0, rate )
            call sort( string_dummy )
            call system_clock( t1, rate )
            tdiff = tdiff + t1 - t0
        end do
        tdiff = tdiff/repeat

        call verify_string_sort( string_dummy, valid, i )
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a17, 2(1x,a4))') 'string_dummy(i-1:i) = ', &
                string_dummy(i-1:i)
        end if
        write( lun, '("| String_type |", 1x, i7, 2x, "|", 1x, a15, " |", ' // &
            'a12, " |",  F10.5, " |" )' ) &
            string_size, a_name, "Sort", tdiff/rate

        ! reverse
        string_dummy = a
        call sort( string_dummy, .true.)
        call verify_string_reverse_sort(string_dummy, valid, i)
        ltest = (ltest .and. valid)
        if ( .not. valid ) then
            write( *, * ) "reverse SORT did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a17, 2(1x,a4))') 'string_dummy(i-1:i) = ', &
                string_dummy(i-1:i)
        end if


    end subroutine test_string_sort

    subroutine test_int_sort_indexes( )

        call test_int_sort_index( blocks, "Blocks" )
        call test_int_sort_index( decrease, "Decreasing" )
        call test_int_sort_index( identical, "Identical" )
        call test_int_sort_index( increase, "Increasing" )
        call test_int_sort_index( rand1, "Random dense" )
        call test_int_sort_index( rand2, "Random order" )
        call test_int_sort_index( rand0, "Random sparse" )
        call test_int_sort_index( rand3, "Random 3" )
        call test_int_sort_index( rand10, "Random 10" )

    end subroutine test_int_sort_indexes

    subroutine test_int_sort_index( a, a_name )
        integer(int32), intent(inout) :: a(:)
        character(*), intent(in)      :: a_name

        integer(int64)                 :: t0, t1, tdiff
        real(dp)                       :: rate
        integer(int64)                 :: i
        logical                        :: valid

        tdiff = 0
        do i = 1, repeat
            dummy = a
            call system_clock( t0, rate )
            call sort_index( dummy, index, work, iwork )
            call system_clock( t1, rate )
            tdiff = tdiff + t1 - t0
        end do
        tdiff = tdiff/repeat

        dummy = a(index)
        call verify_sort( dummy, valid, i )
        if ( .not. valid ) then
            write( *, * ) "SORT_INDEX did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a18, 2i7)') 'a(index(i-1:i)) = ', a(index(i-1:i))
        end if
        write( lun, '("|     Integer |", 1x, i7, 2x, "|", 1x, a15, " |", ' // &
            'a12, " |",  F10.5, " |" )' ) &
            test_size, a_name, "Sort_Index", tdiff/rate

        dummy = a
        call sort_index( dummy, index, work, iwork, reverse=.true. )
        dummy = a(index)
        call verify_reverse_sort( dummy, valid, i )
        if ( .not. valid ) then
            write( *, * ) "SORT_INDEX did not reverse sort " // &
                a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a18, 2i7)') 'a(index(i-1:i)) = ', a(index(i-1:i))
        end if

    end subroutine test_int_sort_index

    subroutine test_char_sort_indexes( )

        call test_char_sort_index( char_decrease, "Char. Decrease" )
        call test_char_sort_index( char_increase, "Char. Increase" )
        call test_char_sort_index( char_rand, "Char. Random" )

    end subroutine test_char_sort_indexes

    subroutine test_char_sort_index( a, a_name )
        character(len=4), intent(in) :: a(0:)
        character(*), intent(in) :: a_name

        integer(int64) :: t0, t1, tdiff
        real(dp)       :: rate
        integer(int64) :: i
        logical        :: valid

        tdiff = 0
        do i = 1, repeat
            char_dummy = a
            call system_clock( t0, rate )
            call sort_index( char_dummy, index, char_work, iwork )
            call system_clock( t1, rate )
            tdiff = tdiff + t1 - t0
        end do
        tdiff = tdiff/repeat

        call verify_char_sort( char_dummy, valid, i )
        if ( .not. valid ) then
            write( *, * ) "SORT_INDEX did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a17, 2(1x,a4))') 'char_dummy(i-1:i) = ', char_dummy(i-1:i)
        end if
        write( lun, '("|   Character |", 1x, i7, 2x, "|", 1x, a15, " |", ' // &
            'a12, " |",  F10.5, " |" )' ) &
            char_size, a_name, "Sort_Index", tdiff/rate

    end subroutine test_char_sort_index

    subroutine test_string_sort_indexes( )

        call test_string_sort_index( string_decrease, "String Decrease" )
        call test_string_sort_index( string_increase, "String Increase" )
        call test_string_sort_index( string_rand, "String Random" )

    end subroutine test_string_sort_indexes

    subroutine test_string_sort_index( a, a_name )
        type(string_type), intent(in) :: a(0:)
        character(*), intent(in) :: a_name

        integer(int64) :: t0, t1, tdiff
        real(dp)       :: rate
        integer(int64) :: i
        logical        :: valid

        tdiff = 0
        do i = 1, repeat
            string_dummy = a
            call system_clock( t0, rate )
            call sort_index( string_dummy, index, string_work, iwork )
            call system_clock( t1, rate )
            tdiff = tdiff + t1 - t0
        end do
        tdiff = tdiff/repeat

        call verify_string_sort( string_dummy, valid, i )
        if ( .not. valid ) then
            write( *, * ) "SORT_INDEX did not sort " // a_name // "."
            write(*,*) 'i = ', i
            write(*,'(a17, 2(1x,a4))') 'string_dummy(i-1:i) = ', &
                string_dummy(i-1:i)
        end if
        write( lun, '("| String_type |", 1x, i7, 2x, "|", 1x, a15, " |", ' // &
            'a12, " |",  F10.5, " |" )' ) &
            string_size, a_name, "Sort_Index", tdiff/rate

    end subroutine test_string_sort_index


    subroutine verify_sort( a, valid, i )
        integer(int32), intent(in) :: a(0:)
        logical, intent(out) :: valid
        integer(int64), intent(out) :: i

        integer(int64) :: n

        n = size( a, kind=int64 )
        valid = .false.
        do i=1, n-1
            if ( a(i-1) > a(i) ) return
        end do
        valid = .true.

    end subroutine verify_sort



    subroutine verify_string_sort( a, valid, i )
        type(string_type), intent(in) :: a(0:)
        logical, intent(out) :: valid
        integer(int64), intent(out) :: i

        integer(int64) :: n

        n = size( a, kind=int64 )
        valid = .false.
        do i=1, n-1
            if ( a(i-1) > a(i) ) return
        end do
        valid = .true.

    end subroutine verify_string_sort

    subroutine verify_char_sort( a, valid, i )
        character(len=4), intent(in) :: a(0:)
        logical, intent(out) :: valid
        integer(int64), intent(out) :: i

        integer(int64) :: n

        n = size( a, kind=int64 )
        valid = .false.
        do i=1, n-1
            if ( a(i-1) > a(i) ) return
        end do
        valid = .true.

    end subroutine verify_char_sort

    subroutine verify_char_reverse_sort( a, valid, i )
        character(len=4), intent(in) :: a(0:)
        logical, intent(out) :: valid
        integer(int64), intent(out) :: i

        integer(int64) :: n

        n = size( a, kind=int64 )
        valid = .false.
        do i=1, n-1
            if ( a(i-1) < a(i) ) return
        end do
        valid = .true.

    end subroutine verify_char_reverse_sort

    subroutine verify_reverse_sort( a, valid, i )
        integer(int32), intent(in) :: a(0:)
        logical, intent(out) :: valid
        integer(int64), intent(out) :: i

        integer(int64) :: n

        n = size( a, kind=int64 )
        valid = .false.
        do i=1, n-1
            if ( a(i-1) < a(i) ) return
        end do
        valid = .true.

    end subroutine verify_reverse_sort

    subroutine verify_string_reverse_sort( a, valid, i )
        type(string_type), intent(in) :: a(0:)
        logical, intent(out) :: valid
        integer(int64), intent(out) :: i

        integer(int64) :: n

        n = size( a, kind=int64 )
        valid = .false.
        do i=1, n-1
            if ( a(i-1) < a(i) ) return
        end do
        valid = .true.

    end subroutine verify_string_reverse_sort

end program test_sorting
