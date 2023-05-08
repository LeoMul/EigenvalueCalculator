program qrdecomp 
    use my_library

    implicit none 
    real * 8 :: A(5,5),G(5,5),GA(5,5)
    real * 8 :: theta 
    integer :: N 

    N = size(A,1)

    A = 0d0 
    A(1,1) = 12.0d0
    A(1,2) = -51.0d0
    A(1,3) = 4.0d0
    A(2,1) = 6.0d0
    A(2,2) = 167.0d0
    A(2,3) = -68.0d0
    A(3,1) = -4.0d0
    A(3,2) = 24.0d0
    A(3,3) = -41.0d0

    !A(1,1) = 6.0d0
    !A(1,2) = 5.0d0
    !A(1,3) = 0.0d0
    !A(2,1) = 5.0d0
    !A(2,2) = 1.0d0
    !A(2,3) = 4.0d0
    !A(3,1) = 0.0d0
    !A(3,2) = 4.0d0
    !A(3,3) = 3.0d0

    A = create_identity(5)

   
    call triangularise(A,G)
    call display_matrix(G)

end program