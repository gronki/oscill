! Program oblicza rozwiazania kwantowego oscylatora harmonicznego
! (c) 2016 Dominik Gronkiewicz www.gronki.pl

! Procedura do obliczania wielomianow Hermite'a:
! https://people.sc.fsu.edu/~jburkardt/f_src/hermite_polynomial/hermite_polynomial.html

program oscill
    use cgs_constants
    use qho

    implicit none


    double precision :: x_max, mass, omega, x_osc
    double precision, allocatable :: x(:), psi(:), pkl(:)
    integer :: i, n_max, n_pts
    character (len=*), parameter :: fmt1e = "(A30,Es12.4,1X,A)"

    character (len=128) :: buf, buf2

    print '(A40,$)', 'Poziom wzbudzony > '
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
    allocate( psi(n_pts), pkl(n_pts) )

    print '(A40,$)', 'Zakres obliczen [jedn.oscyl.] > '
    read (*,*), x_max

    print '(A40,$)', 'Masa czasteczki [g] > '
    read (*,*), mass
    if ( mass <= 0 ) then
        print *, 'MUSI BYC DODATNIA'
        stop -1
    endif

    print '(A40,$)', 'Czestosc pulapki [rad/s] > '
    read (*,*), omega
    if ( omega <= 0 ) then
        print *, 'MUSI BYC DODATNIA'
        stop -1
    endif

    x_osc = qho_x_osc(mass,omega)
    print fmt1e, 'Jednostka dlugosci', x_osc, '[cm]'

    call qho_xrang(x,x_max)
    call qho_state(n_max,x,psi)
    call qho_state_clas(n_max,x,pkl)

    write (buf,'(A,I0,A,I0,A)'), 'state.n', n_max, '.p', n_pts, '.csv'

    open(unit=33, file=buf, action='write')

    write (33,'(A8,$)'), 'PRZEDZ'
    write (buf, '(A,I0)'), 'Psi_', n_max
    write (33,'(4A12,$)'), 'X [j.o.]', trim(buf), '|Psi|^2', 'P_kl'
    write (33,'(4A12)'), 'X [cm]', trim(buf), '|Psi|^2', 'P_kl'


    do i = 1,n_pts
        write (33,'(I8,12Es12.4)'), i,      &
                x(i), psi(i), abs(psi(i))**2, pkl(i),   &
                x(i) * x_osc, psi(i) / sqrt(x_osc), abs(psi(i))**2 / x_osc, pkl(i) / x_osc
    end do

    close(unit=33)


end program
