module precision_module

    implicit none
    
    ! Precision of state arrays (q and aux) and Riemann solver output arrays
    integer, parameter :: sk = kind(1.d0)
    
    ! Precision of timer variables
    integer, parameter :: tk = kind(1.d0)
    
end module precision_module
