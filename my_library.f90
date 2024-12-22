module my_library 
    implicit none 
    contains 

    function my_matmul_serial(x_matrix,y_matrix) result(res_matrix)
        real * 8 , intent(in) :: x_matrix(:,:), y_matrix(:,:) 
        real * 8 , allocatable :: res_matrix(:,:)
        real * 8 :: sum
        integer :: ii,jj,kk 
        integer :: n,m,o,r 
        n = size(x_matrix,1)
        m = size(x_matrix,2)
        o = size(y_matrix,1)
        r = size(y_matrix,2) 

        print*,"beginning", n,m,o,r

        if (m .ne. o ) then 
            allocate(res_matrix(1,1))
            res_matrix = 0.0d0
            print*, "These two matrices cannot be multiplied, check dimensions. Returning zero scalar."
        else 
            print*,"allocating"
            allocate(res_matrix(n,r))
            res_matrix = 0.0d0
            !call display_matrix(res_matrix)
            do  ii = 1,n
                do jj = 1,r 
                    sum = 0.0d0
                    do kk = 1,o 
                        sum = sum + x_matrix(ii,kk) * y_matrix(kk,jj)
                    end do 
                    res_matrix(ii,jj) = sum
                    !print*,ii,jj 
                end do 
            end do 
            end if 
    end function


    subroutine display_matrix(matrix) 
        real * 8 ,intent(in) :: matrix(:,:)
        integer :: n,m ,i,j

        m = size(matrix,1)
        n = size(matrix,2)
        do i = 1,m 
            do j = 1,n 
                write(*,FMT = "(10F10.5)", advance = "no") matrix(i,j)
                write (*,FMT = "(A1)",advance = "no") " "
            end do
            print*, ""
        end do 
    end subroutine

    function create_identity(n) result(identity)
        integer :: n 
        real * 8 , allocatable :: identity(:,:)
        integer :: ii
        allocate(identity(n,n))
        identity = 0.0d0
        do ii = 1,n 
            identity(ii,ii) = 1.0d0
        end do 
    end function

    function create_example(n,m) result(example)
        integer :: n,m 
        real * 8 , allocatable :: example(:,:)
        integer :: ii,jj
        allocate(example(n,m))
        do ii = 1,n 
            do jj = 1,m 
                example(ii,jj) = real((ii-1)*n + jj)
            end do 
        end do 

    end function

    subroutine givens_rotation(a, b, c, s, r)

        real*8:: a, b, c, s, r
        real*8:: h, d
        
        if (b.ne.0.0) then


            h = hypot(a, b)
            d = 1.0d0 / h
            c = (a) * d
            s = -b*d
            r = sign(1.0d0, a) * h
            !print*,a,b,c,h

        else
            c = 1.0d0
            s = 0.0d0
            r = a
        end if
        
        return
    end 

    subroutine triangularise(A,out_matrix,Q_matrix)
        real * 8 :: A(:,:),out_matrix(:,:),Q_matrix(:,:)
        real * 8,allocatable :: G(:,:)
        real * 8:: s,c,r
        integer :: N ,ii,jj
        N = size(A,1)
        allocate(G(N,N))

        !print*,"size",N
        
        out_matrix = A
        !call display_matrix(out_matrix)
        Q_matrix = create_identity(n)
        do jj = 1,N
            do ii = N,jj+1,-1
                !print*,ii,jj,out_matrix(ii,jj),out_matrix(jj,jj)

                call givens_rotation(out_matrix(jj,jj),out_matrix(ii,jj),c,s,r)

                call create_G_matrix(N,ii,jj,G,s,c)
                !print*, "G matrix:"
                !call display_matrix(G)
                !print*,"////////////"
                out_matrix = matmul(G,out_matrix)
                !print*,"Out matrix:"
                !call display_matrix(out_matrix)
                !print*,"///////"
                Q_matrix = matmul(Q_matrix,transpose(G))
            end do 

        end do 


    end 

    subroutine eigenvalue_QR(A,eigenvalues,tolerance,work) 
        use, intrinsic :: iso_fortran_env, only: REAL64, INT64
        real * 8 :: A(:,:), eigenvalues(:),tolerance,work(:,:)
        real * 8,allocatable:: AK(:,:),G(:,:),Q(:,:)
        real * 8,allocatable:: prev_eigenvalues(:),eigendiff(:)
        real * 8:: norm,prevnorm
        integer :: N ,ii, max_iter

        N = size(A,1)
        allocate(AK(N,N),G(N,N),Q(N,N))
        allocate(prev_eigenvalues(N),eigendiff(N))

        norm =  huge(0.0_REAL64)
        prevnorm = 0.0d0 
        max_iter = 1000 
        AK = A
        call triangularise(AK,G,Q)
        AK = matmul(AK,Q)
        AK = matmul(TRANSPOSE(Q),AK)
        call get_diagonals(AK,prev_eigenvalues)
        work = Q
        do ii = 1,max_iter
            call triangularise(AK,G,Q)
            AK = matmul(AK,Q)
            AK = matmul(TRANSPOSE(Q),AK)
            work = matmul(work,Q)
            call get_diagonals(AK,eigenvalues)
            eigendiff = eigenvalues - prev_eigenvalues
            norm = dot_product(eigendiff,eigendiff) 
            !print*,norm,eigenvalues,prev_eigenvalues
            if (norm < tolerance) then
                PRINT*,"exiting eigenvalue finder at iteration",ii
                EXIT
            end if 
            prev_eigenvalues = eigenvalues
        end do 

        !work = AK


    end 

    subroutine get_diagonals(A,diagonals)
        real * 8 :: A(:,:), diagonals(:)
        integer :: N ,ii 
        N = size(A,1)
        do ii = 1,N 
            diagonals(ii) = A(ii,ii)
        end do
    end

    subroutine calculate_Q_matrix(A,Q)
        real * 8 :: A(:,:),Q(:,:)
        real * 8:: s,c,r
        integer :: N ,ii,jj
        N = size(A,1)

        !print*,"size",N
        
        !Q = create_identity(N)
        !!call display_matrix(out_matrix)
!
        !do jj = 1,N
        !    do ii = N,jj+1,-1
        !        print*,ii,jj,out_matrix(ii,jj),out_matrix(jj,jj)
!
        !        call givens_rotation(out_matrix(jj,jj),out_matrix(ii,jj),c,s,r)
!
        !        call create_G_matrix(N,ii,jj,G,s,c)
        !        print*, "G matrix:"
        !        call display_matrix(G)
        !        print*,"////////////"
        !        out_matrix = matmul(G,out_matrix)
        !        print*,"Out matrix:"
        !        call display_matrix(out_matrix)
        !        print*,"///////"
        !    end do 
!
        !end do 
    end 


    subroutine create_G_matrix(N,i,j,Gmatrix,s,c)
        integer :: N,i,j, ii 
        real * 8 :: Gmatrix(N,N),s,c
        Gmatrix = 0.0d0
        do ii = 1,N 
            Gmatrix(ii,ii) = 1.0d0 
        end do 
        Gmatrix(i,i) = c 
        Gmatrix(j,j) = c 
        Gmatrix(i,j) = s 
        Gmatrix(j,i) = -s

    end subroutine



end module