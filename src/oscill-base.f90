! Program oblicza rozwiazania kwantowego oscylatora harmonicznego
! (c) 2016 Dominik Gronkiewicz www.gronki.pl

! Procedura do obliczania wielomianow Hermite'a:
! https://people.sc.fsu.edu/~jburkardt/f_src/hermite_polynomial/hermite_polynomial.html

program oscill
    use cgs_constants
    use qho

    implicit none


    double precision :: x_max
    double precision, allocatable :: x(:), psi(:,:)
    integer :: i, n_max, n_pts
    character (len=*), parameter :: fmt1e = "(A30,Es12.4,1X,A)"

    character (len=128) :: buf, buf2

    print '(A40,$)', 'Najwyzszy poziom wzbudzony > '
    read (*,*), n_max
    if ( n_max < 2 ) then
        print *, 'MUSI BYC CO NAJMNIEJ 2'
        stop -1
    endif

    print '(A40,$)', 'Liczba punktow obliczeniowych > '
    read (*,*), n_pts
    if ( n_pts < 10 ) then
        print *, 'MUSI BYC CO NAJMNIEJ 10'
        stop -1
    endif


    allocate( x(n_pts) )
    allocate( psi(n_pts,n_max) )

    print '(A40,$)', 'Zakres obliczen [jedn.oscyl.] > '
    read (*,*), x_max

    call qho_xrang(x,x_max)

    call qho_state_base(x,psi)

    write (buf,'(A,I0,A,I0,A)'), 'base.n1to', n_max, '.p', n_pts, '.csv'

    open(unit=33, file=buf, action='write')

    write (33,'(A8,A12,$)'), 'PRZEDZ', 'X [j.o.]'
    do i = 1,n_max
        write (buf, '(A,I0)'), 'Psi_', i
        write (33,'(A12,$)'), trim(buf)
    enddo

    write (33,*), ''
    do i = 1,n_pts
        write (33,'(I8,200Es12.4)'), i, x(i), psi(i,:)
    end do

    close(unit=33)


end program
