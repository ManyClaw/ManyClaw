program manyclaw

    ! Test function definitions
    use test_functions_module

    ! Utility modules
    use timer_module
    use precision_module

    implicit none
    
    ! Test suite parameters
    integer, parameter :: NUM_TESTS = 4
    
    ! Local storage
    real(kind=TK) :: ave_time(NUM_TESTS)

    ! Common block parameters
    real(kind=sk) :: u,v
    real(kind=sk) :: gamma,gamma1
    common /advection_rp/ u,v
    common /euler_rp/  gamma,gamma1
    
    ! Test advection Riemann solvers
    nx = 1028
    ny = 1028
    num_states = 1
    num_aux = 0
    num_ghost = 2
    num_waves = 1
    u = 1.0_sk
    v = 1.0_sk
    ave_time(1) = time_rp_row_function(row_rp_advection,3)
    print *,"Advection row Riemann solver time = ",ave_time(1)," ms."
    ave_time(2) = time_rp_ptwise_function(ptwise_rp_advection,3)
    print *,"Advection ptwise Riemann solver time = ",ave_time(2)," ms."
    
    ! Test Euler Riemann solvers
    nx = 1028
    ny = 1028
    num_states = 4
    num_aux = 0
    num_ghost = 2
    num_waves = 4
    gamma = 1.4_sk
    gamma1 = 1.0_sk - gamma
    ave_time(3) = time_rp_row_function(row_rp_euler,3)
    print *,"Euler row Riemann solver time = ",ave_time(3)," ms."
    ave_time(4) = time_rp_ptwise_function(ptwise_rp_euler,3)
    print *,"Euler row Riemann solver time = ",ave_time(4)," ms."

end program manyclaw
