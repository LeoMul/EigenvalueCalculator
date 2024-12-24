module real_symmetric
    use my_library
 
    implicit none 
    contains 

    subroutine diag_rsym_jacobi(A,v,n)
    integer :: n
    real*8 :: A(n,n)
    real*8 :: Acopy(n,n)
    real*8 :: v(n,n),b(n),d(n),z(n)
    
    real*8 :: sm ,tresh,g,h,t,theta,c,s
    real*8 :: tau
    integer :: ii,jj ,ip,iq 
    integer :: nrotations 

    v = 0.0d0 
    Acopy = A  

    nrotations = 0 
    do ii = 1,n 
        v(ii,ii) = 1.0d0
        b(ii) = Acopy(ii,ii)
        z(ii) = 0.0d0
        !d(ii) = A(ii,ii)
    end do

    d = b 

    do ii = 1,53
        sm = 0.0d0 
        do ip=1,n-1
            do iq = ip+1,n 
                sm = sm + abs(Acopy(ip,iq))
            end do
        end do
        if (sm.eq.0.0d0) then 
            call display_matrix(Acopy)
            print*,d
            return
        end if  

        if (ii .lt. 4) then 
            tresh = 0.2 * sm / (n*n)
        else 
            tresh = 0.0d0 
        end if

        do ip=1,n-1
            do iq=ip+1,n 
                g = 100.0d0 * abs(acopy(ip,iq)) 

                if ((ii .gt. 4) .and. (abs(d(ip))+g .eq. abs(d(ip))) .and. (abs(d(iq))+g .eq. abs(d(iq)))) then 
                    acopy(ip,iq) = 0.0d0 
                else if (abs(acopy(ip,iq)) .gt. tresh) then 
                    h = d(iq) - d(ip) 
                    if (abs(h) + g .eq. abs(h)) then 
                        t = acopy(ip,iq) / h 
                    else 
                        theta = 0.5 * h / acopy(ip,iq)
                        t = 1.d0 / (abs(theta) + sqrt(1.0d0+theta*theta))
                        if(theta.lt.0.0d0) t=-t 
                    end if 
                    c = 1.0 / sqrt(1.0d0 + t*t)
                    s = t*c 
                    tau = s/(1.0d0 + c)
                    h = t * acopy(ip,iq)

                    z(ip) = z(ip) - h 
                    z(iq) = z(iq) + h

                    d(ip) = d(ip) - h 
                    d(iq) = d(iq) + h     
                    
                    acopy(ip,iq) = 0.0d0 
                    
                    do jj = 1,ip-1 
                        g = acopy(jj,ip)
                        h = acopy(jj,iq)
                        acopy(jj,ip) = g - s * (h + g * tau)
                        acopy(jj,iq) = h + s * (g - h * tau)
                    end do 

                    do jj=ip+1,iq-1 
                        g = acopy(ip,jj)
                        h = acopy(jj,iq)
                        acopy(ip,jj) = g - s * (h + g * tau) 
                        acopy(jj,iq) = h + s * (g - h * tau)
                    end do 

                    do jj=iq+1,n 
                        g = acopy(ip,jj)
                        h = acopy(iq,jj)
                        acopy(ip,jj) = g - s * (h + g * tau) 
                        acopy(iq,jj) = h + s * (g - h * tau)
                    end do 

                    do jj=1,n 
                        g = v(jj,ip)
                        h = v(jj,iq)
                        v(jj,ip) = g - s * (h + g * tau) 
                        v(jj,iq) = h + s * (g - h * tau)
                    end do 
                    nrotations = nrotations + 1
                end if 

            end do 
        end do 

        do ip=1,n 
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0d0
        end do

    end do 
    PRINT*, 'too many iterations'
    return 
    end subroutine

end module 