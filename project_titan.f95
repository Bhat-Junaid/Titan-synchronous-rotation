program Titan
    implicit none
    real(8), parameter :: pi = 4*atan(1.0d0)
    real(8) :: t, p, dp, h
    integer :: i, n
    character(len=22) :: file

    ! initial values
    t = 0.0d0
    p = pi              ! Initial rotation angle rad 
    dp = 0.01643d0      ! Initial rotation rate rad/hour
    n = 38200 * 100     ! Total time steps
    h = 0.01            ! Time step in hours
    

    ! output file name
    file = 'titan_p_dp_r100_e.txt'

    ! Run RK4 integration
    call rk4(force_func, t, p, dp, h, n, file)
    

contains
    ! Runge-Kutta 4th order integration
    subroutine rk4(func, t, p, dp, h, n, file)
        implicit none
        real(8), intent(inout) :: t, p, dp, h
        integer, intent(in) :: n
        character(len=*), intent(in) :: file
        real(8) :: k1p, k2p, k3p, k4p, k1dp, k2dp, k3dp, k4dp ! Runge-Kutta coefficients
        integer :: unit

        real(8) :: t_values(n+1), p_values(n+1), dp_values(n+1)

        ! Function interface
        interface
            real(8) function func(t, p, dp)
                real(8), intent(in) :: t, p, dp
            end function func
        end interface

        ! Open the output file
        open(newunit=unit, file=file, status='replace', action='write')

        do i = 1, n+1
            !t_values(i) = t
            !p_values(i) = p
            !dp_values(i) = dp

            k1p = h * dp
            k1dp = h * func(t, p, dp)
            
            k2p = h * (dp + 0.5 * k1dp)
            k2dp = h * func(t + 0.5 * h, p + 0.5 * k1p, dp + 0.5 * k1dp)
            
            k3p = h * (dp + 0.5 * k2dp)
            k3dp = h * func(t + 0.5 * h, p + 0.5 * k2p, dp + 0.5 * k2dp)
            
            k4p = h * (dp + k3dp)
            k4dp = h * func(t + h, p + k3p, dp + k3dp)

            p = p + (k1p + 2.0 * k2p + 2.0 * k3p + k4p) / 6.0
            dp = dp + (k1dp + 2.0 * k2dp + 2.0 * k3dp + k4dp) / 6.0
            t = t + h

            ! Write values to the file
            write(unit, *) t, p - (0.01643d0)*t, dp    
        end do

        ! Close the output file
        close(unit)
    end subroutine rk4

    ! Force function
    real(8) function force_func(t, p, dp)
        real(8), intent(in) :: t, p, dp
        real(8) :: omega, e, l, I

        omega = 0.01643d0            ! rotation frequency rad/hour
        e = 0.0289                  ! eccentricity
        I = 1.35e-4                  ! moment of inertia
        l = omega * t                ! mean longitude

        force_func = -0.75 * omega**2 * I * (-e * sin(2.0 * p - l) + &
            2.0 * (1.0 - 5.0 * e**2 / 2.0) * sin(2.0 * p - 2.0 * l) + &
            7.0 * e * sin(2.0 * p - 3.0 * l) + 17.0 * e**2 * sin(2.0 * p - 4.0 * l)) 
    end function force_func

end program Titan


