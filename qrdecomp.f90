program qrdecomp 
    use my_library

    implicit none 
    real * 8 :: A(3,3),G(3,3),GA(3,3),Q(3,3),AK(3,3)
    real * 8 ::eigenvalues(3)
    real * 8 :: tolerance 
    real * 8 :: work(3,3)
    integer :: N ,ii ,max_iter

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

    !A = create_identity(5)

   
    !call triangularise(A,G,Q)
    !call display_matrix(matmul(Q,G))
    !call display_matrix(A)
    !max_iter = 1
    !AK = A
    !do ii = 1,max_iter
    !    call triangularise(AK,G,Q)
    !    AK = matmul(AK,Q)
    !    AK = matmul(TRANSPOSE(Q),AK)
    !    print*, "iteration",ii
    !    call display_matrix(AK)
    !    print*,"/////"
    !end do
    tolerance = 0.000001d0
    call eigenvalue_QR(A,eigenvalues,tolerance,work)
    call display_matrix(work)
    !call display_matrix(reshape(eigenvalues,[1,N]))
    print*,(eigenvalues)
    call display_matrix(matmul(transpose(work),matmul(A,work)))

end program