program manyclaw

    ! Test function definitions
    use test_functions_module

    ! Utility modules
    use timer_module
    use precision_module

    implicit none
    
    ! Command line arguments
    character(len=10) :: input_nx,input_ny,input_num_tests
    
    ! Local storage
    integer :: num_tests
    real(kind=TK) :: ave_time(4)

    ! Common block parameters
    real(kind=sk) :: u,v
    real(kind=sk) :: gamma,gamma1
    common /advection_rp/ u,v
    common /euler_rp/  gamma,gamma1
    
    ! Parse command line arguments
    nx = 1028
    ny = 1028
    num_tests = 3
    select case(iargc())
        case(0)
            continue
        case(1)
            call getarg(1,input_nx)
            read(input_nx,'(I10)') nx
            ny = nx
            num_tests = 3
        case(2)
            call getarg(1,input_nx)
            call getarg(2,input_ny)
            read(input_nx,'(I10)') nx
            read(input_ny,'(I10)') ny
        case(3)
            call getarg(1,input_nx)
            call getarg(2,input_ny)
            call getarg(3,input_num_tests)
            read(input_nx,'(I10)') nx
            read(input_ny,'(I10)') ny
            read(input_num_tests,'(I10)') num_tests
        case default
            print *, "***ERROR*** Too many command line arguments!"
            stop
    end select
    
    ! Test advection Riemann solvers
    num_states = 1
    num_aux = 0
    num_ghost = 2
    num_waves = 1
    u = 1.0_sk
    v = 1.0_sk
    ave_time(1) = time_rp_row1_function(row_rp_advection,num_tests)
    print *,"Advection row method 1 Riemann solver time = ",ave_time(1)," ms."
    ave_time(1) = time_rp_row2_function(row_rp_advection,num_tests)
    print *,"Advection row method 2 Riemann solver time = ",ave_time(1)," ms."
    ave_time(1) = time_rp_row3_function(row_rp_advection,num_tests)
    print *,"Advection row method 3 Riemann solver time = ",ave_time(1)," ms."
    ave_time(2) = time_rp_ptwise_function(ptwise_rp_advection,num_tests)
    print *,"Advection ptwise Riemann solver time = ",ave_time(2)," ms."
    
    ! Test Euler Riemann solvers
    num_states = 4
    num_aux = 0
    num_ghost = 2
    num_waves = 4
    gamma = 1.4_sk
    gamma1 = 1.0_sk - gamma
    ave_time(3) = time_rp_row1_function(row_rp_euler,num_tests)
    print *,"Euler row method 1 Riemann solver time = ",ave_time(3)," ms."
    ave_time(3) = time_rp_row2_function(row_rp_euler,num_tests)
    print *,"Euler row method 2 Riemann solver time = ",ave_time(3)," ms."
    ave_time(3) = time_rp_row3_function(row_rp_euler,num_tests)
    print *,"Euler row method 3 Riemann solver time = ",ave_time(3)," ms."
    ave_time(4) = time_rp_ptwise_function(ptwise_rp_euler,num_tests)
    print *,"Euler ptwise Riemann solver time = ",ave_time(4)," ms."

end program manyclaw
