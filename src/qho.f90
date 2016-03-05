module qho

    use cgs_constants
    implicit none

contains

    subroutine qho_state_base ( x, wave )
        real (kind=8), intent(out) :: wave(:,:)
        real (kind=8), intent(in) :: x(:)
        double precision, allocatable :: cof(:,:)
        integer :: n_max, i, j, n
        double precision ::  t

        n_max = size(wave,2)

        allocate( cof(n_max+1,n_max+1) )

        ! call h_polynomial_coefficients(n_max,cof)
        call qho_hermit(cof)

        do n = 0, n_max-1
            t = ( pi**0.25 * sqrt( 2.**n * qho_fac(n) ) )
            do j = 1, size(x)
                wave(j,n+1) = exp(-x(j)**2 / 2.) * qho_eval_polyn(x(j),cof(:,n+1)) / t
            end do
        enddo

    end subroutine

    subroutine qho_state (n,x,wave)
        real (kind=8), intent(out) :: wave(:)
        real (kind=8), intent(in) :: x(:)
        integer, intent(in) :: n
        double precision, allocatable :: cof(:,:)
        integer :: i, j
        real (kind=8) :: hermi, t


        allocate( cof(n+1,n+1) )

        call qho_hermit(cof)

        t = ( pi**0.25 * sqrt( 2.**n * qho_fac(n) ) )
        do j = 1, size(x)
            wave(j) = exp(-x(j)**2 / 2.) * qho_eval_polyn(x(j),cof(:,n+1)) / t
        end do

    end subroutine

    real (kind=8) function qho_eval_polyn(x,coeff) result(w)
        real (kind=8), intent(in) :: coeff(:)
        real (kind=8), intent(in) :: x
        integer :: i

        w = 0
        do i = size(coeff), 1, -1
            w = w * x + coeff(i)
        end do
    end function

    real (kind=8) function qho_eval_polyn0(x,coeff) result(w)
        real (kind=8), intent(in) :: coeff(:)
        real (kind=8), intent(in) :: x
        integer :: i
        w = coeff(1)
        if ( size(coeff) .ge. 2 ) then
            do i = 2,size(coeff)
                w = w + coeff(i) * x**real(i-1,kind=8)
            end do
        endif
    end function

    subroutine qho_state_clas (n,x,wave)
        real (kind=8), intent(out) :: wave(:)
        real (kind=8), intent(in) :: x(:)
        integer, intent(in) :: n
        integer :: j


        do j = 1, size(x)
            if ( abs(x(j)) < sqrt(2. * n + 1) ) then
                wave(j) = 1. / ( pi * sqrt( 2. * n + 1 - x(j)**2 ) )
            else
                wave(j) = 0.
            endif
        end do

    end subroutine


    elemental real(kind=8) function qho_x_osc (mass,omega)
        real (kind=8), intent(in) :: mass, omega
        qho_x_osc = sqrt( cgs_hbar / (mass*omega) )
    end function


    subroutine qho_xrang (x,x_max)
        real (kind=8), intent(in) :: x_max
        real (kind=8), intent(out) :: x(:)
        integer :: i,n
        n = size(x)
        do i = 1,n
            x(i) = 2. * x_max * (i-1.)/(n-1.) - x_max
        end do
    end subroutine

    real(kind=8) function qho_fac(n) result(f)
        integer, intent(in) :: n
        integer :: i
        f = 1d0
        if ( n < 2 ) return
        do i = 1,n
            f = f * i
        end do
    end function

    subroutine qho_hermit(coeff)
        real (kind=8), intent(inout) :: coeff(:,:)
        integer :: n, k, nn

        nn = size(coeff,2)

        if ( size(coeff,1) < 3 .or. size(coeff,2) < 3 .or. size(coeff,1) < size(coeff,2) ) then
            write (0,*) 'BAD ARRAY SIZE'
            return
        endif

        coeff(:,:) = 0.
        coeff(1,1) = 1.
        coeff(2,2) = 2.

        do n = 2, nn-1
            coeff(1,n+1) = -coeff(2,n)
            do k = 1, nn-2
                coeff(k+1,n+1) = 2 * coeff(k,n) - (k+1) * coeff(k+2,n)
            end do
            coeff(nn,n+1) = 2 * coeff(nn-1,n)
        end do
    end subroutine

end module
