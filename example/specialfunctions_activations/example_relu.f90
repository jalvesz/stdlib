program example_relu
    use stdlib_kinds, only: sp
    use stdlib_math, only: linspace
    use stdlib_specialfunctions, only: relu
    
    integer, parameter :: n = 10
    real(sp) :: x(n), y(n)
    implicit none
  
    x = linspace(-2._sp, 2._sp, n)
    y = relu( x , 1.0 )
end program example_relu
  