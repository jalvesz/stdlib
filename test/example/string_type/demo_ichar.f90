program demo_ichar
use stdlib_string_type
implicit none
type(string_type) :: string
integer :: code

string = "Fortran"
code = ichar(string)
end program demo_ichar
