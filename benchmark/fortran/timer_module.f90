module timer_module

    use precision_module

    implicit none
    
    ! Parameter related to problem setup
    integer, save :: num_states,num_aux,nx,ny,num_ghost,num_waves
    
    integer, save, private :: clock_start,clock_finish,clock_rate
    
contains
    
    subroutine start_timer()
        implicit none
        call system_clock(clock_start,clock_rate)
    end subroutine start_timer
    
    
    real(kind=tk) function end_timer() result(time)
        use precision_module
        implicit none
        call system_clock(clock_finish,clock_rate)
        time = real(clock_finish - clock_start,kind=tk) &
                    / real(clock_rate,kind=tk)
    end function end_timer
    
    
    real(kind=tk) function time_rp_row1_function(rp_function,num_tests) &
                                result(ave_time)

        use precision_module
        
        implicit none
        
        ! Input arguments
        integer, intent(in) :: num_tests

        ! Arrays for computin'
        real(kind=sk), allocatable :: q(:,:,:)
        real(kind=sk), allocatable :: aux(:,:,:)

        ! Input function interface
        interface
            subroutine rp_function(ixy,maxm,meqn,mwaves,num_ghost,mx,ql,qr,auxl, &
                                        auxr,wave,s,amdq,apdq,num_aux)
                use precision_module
                integer :: ixy,maxm,meqn,mwaves,num_ghost,mx,num_aux
                real(kind=sk) :: wave(meqn, mwaves, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: s(mwaves, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: ql(meqn, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: qr(meqn, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: auxl(num_aux, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: auxr(num_aux, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: apdq(meqn, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: amdq(meqn, 1-num_ghost:maxm+num_ghost)
            end subroutine rp_function
        end interface
        
        ! Riemann solver input arguments
        real(kind=sk), allocatable :: ql(:,:)
        real(kind=sk), allocatable :: qr(:,:)
        real(kind=sk), allocatable :: auxl(:,:)
        real(kind=sk), allocatable :: auxr(:,:)
        
        ! Output from the Riemann solver
        real(kind=sk), allocatable :: wave(:,:,:)
        real(kind=sk), allocatable :: s(:,:)
        real(kind=sk), allocatable :: apdq(:,:)
        real(kind=sk), allocatable :: amdq(:,:)
        
        ! Local storage
        real(kind=tk) :: total_time = 0.0_tk
        integer :: m,n,i,j,maxm

        ! Allocate actual arrays
        allocate(q(num_states,1-num_ghost:nx+num_ghost,1-num_ghost:ny+num_ghost))
        allocate(aux(num_aux,1-num_ghost:nx+num_ghost,1-num_ghost:ny+num_ghost))

        ! Setup arrays and problem
        maxm = max(nx,ny)
        allocate(ql(num_states,1-num_ghost:maxm+num_ghost))
        allocate(qr(num_states,1-num_ghost:maxm+num_ghost))
        allocate(auxl(num_states,1-num_ghost:maxm+num_ghost))
        allocate(auxr(num_states,1-num_ghost:maxm+num_ghost))
        allocate(wave(num_states, num_waves, 1-num_ghost:maxm+num_ghost))
        allocate(s(num_waves, 1-num_ghost:maxm+num_ghost))
        allocate(apdq(num_states, 1-num_ghost:maxm+num_ghost))
        allocate(amdq(num_states, 1-num_ghost:maxm+num_ghost))
        
        ! Initialize random number generator
        call init_random_seed()
        
        do n=1,num_tests
            ! Setup this test with random input
            call random_number(q)
            call random_number(aux)
            
            call start_timer()
            do j=1-num_ghost,ny+num_ghost
                do m=1,num_states
                    do i=1-num_ghost,nx+num_ghost
                        ql(m,i) = q(m,i,j)
                        qr(m,i) = q(m,i,j)
                    enddo
                enddo
                do m=1,num_aux
                    do i=1-num_ghost,nx+num_ghost
                        auxl(m,i) = aux(m,i,j)
                        auxr(m,i) = aux(m,i,j)
                    enddo
                enddo
                call rp_function(1,maxm,num_states,num_waves,num_ghost,nx,ql,qr, &
                                    auxl,auxr,wave,s,amdq,apdq,num_aux)
            enddo
            do i=1-num_ghost,nx+num_ghost
                do m=1,num_states
                    do j=1-num_ghost,ny+num_ghost
                        ql(m,j) = q(m,i,j)
                        qr(m,j) = q(m,i,j)
                    enddo
                enddo
                do m=1,num_aux
                    do j=1-num_ghost,ny+num_ghost
                        auxl(m,j) = aux(m,i,j)
                        auxr(m,j) = aux(m,i,j)
                    enddo
                enddo
                call rp_function(2,maxm,num_states,num_waves,num_ghost,nx,ql,qr, &
                    auxl,auxr,wave,s,amdq,apdq,num_aux)
            enddo
            total_time = total_time + end_timer()
        enddo
        
        ! Output time (ms)
        ave_time = 1000.0_sk * (total_time / real(num_tests,kind=tk))
        
    end function time_rp_row1_function
    
    
    real(kind=tk) function time_rp_row2_function(rp_function,num_tests) &
                                result(ave_time)

        use precision_module
        
        implicit none
        
        ! Input arguments
        integer, intent(in) :: num_tests

        ! Arrays for computin'
        real(kind=sk), allocatable :: q(:,:,:)
        real(kind=sk), allocatable :: aux(:,:,:)

        ! Input function interface
        interface
            subroutine rp_function(ixy,maxm,meqn,mwaves,num_ghost,mx,ql,qr,auxl, &
                                        auxr,wave,s,amdq,apdq,num_aux)
                use precision_module
                integer :: ixy,maxm,meqn,mwaves,num_ghost,mx,num_aux
                real(kind=sk) :: wave(meqn, mwaves, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: s(mwaves, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: ql(meqn, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: qr(meqn, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: auxl(num_aux, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: auxr(num_aux, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: apdq(meqn, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: amdq(meqn, 1-num_ghost:maxm+num_ghost)
            end subroutine rp_function
        end interface
        
        ! Riemann solver input arguments
        real(kind=sk), allocatable :: ql(:,:)
        real(kind=sk), allocatable :: qr(:,:)
        real(kind=sk), allocatable :: auxl(:,:)
        real(kind=sk), allocatable :: auxr(:,:)
        
        ! Output from the Riemann solver
        real(kind=sk), allocatable :: wave(:,:,:)
        real(kind=sk), allocatable :: s(:,:)
        real(kind=sk), allocatable :: apdq(:,:)
        real(kind=sk), allocatable :: amdq(:,:)
        
        ! Local storage
        real(kind=tk) :: total_time = 0.0_tk
        integer :: m,n,i,j,maxm

        ! Allocate actual arrays
        allocate(q(num_states,1-num_ghost:nx+num_ghost,1-num_ghost:ny+num_ghost))
        allocate(aux(num_aux,1-num_ghost:nx+num_ghost,1-num_ghost:ny+num_ghost))

        ! Setup arrays and problem
        maxm = max(nx,ny)
        allocate(ql(num_states,1-num_ghost:maxm+num_ghost))
        allocate(qr(num_states,1-num_ghost:maxm+num_ghost))
        allocate(auxl(num_states,1-num_ghost:maxm+num_ghost))
        allocate(auxr(num_states,1-num_ghost:maxm+num_ghost))
        allocate(wave(num_states, num_waves, 1-num_ghost:maxm+num_ghost))
        allocate(s(num_waves, 1-num_ghost:maxm+num_ghost))
        allocate(apdq(num_states, 1-num_ghost:maxm+num_ghost))
        allocate(amdq(num_states, 1-num_ghost:maxm+num_ghost))
        
        ! Initialize random number generator
        call init_random_seed()
        
        do n=1,num_tests
            ! Setup this test with random input
            call random_number(q)
            call random_number(aux)
            
            call start_timer()
            do j=1-num_ghost,ny+num_ghost
                ql = q(:,:,j)
                qr = q(:,:,j)
                if (num_aux > 0) then
                    auxl = aux(:,:,j)
                    auxr = aux(:,:,j)
                endif
                call rp_function(1,maxm,num_states,num_waves,num_ghost,nx,ql,qr, &
                                    auxl,auxr,wave,s,amdq,apdq,num_aux)
            enddo
            do i=1-num_ghost,nx+num_ghost
                ql = q(:,i,:)
                qr = q(:,i,:)
                if (num_aux > 0) then
                    auxl = aux(:,i,:)
                    auxr = aux(:,i,:)
                endif
                call rp_function(2,maxm,num_states,num_waves,num_ghost,nx,ql,qr, &
                    auxl,auxr,wave,s,amdq,apdq,num_aux)
            enddo
            total_time = total_time + end_timer()
        enddo
        
        ! Output time (ms)
        ave_time = 1000.0_sk * (total_time / real(num_tests,kind=tk))
        
    end function time_rp_row2_function
    
    
    real(kind=tk) function time_rp_row3_function(rp_function,num_tests) &
                                result(ave_time)

        use precision_module
        
        implicit none
        
        ! Input arguments
        integer, intent(in) :: num_tests

        ! Arrays for computin'
        real(kind=sk), allocatable :: q(:,:,:)
        real(kind=sk), allocatable :: aux(:,:,:)

        ! Input function interface
        interface
            subroutine rp_function(ixy,maxm,meqn,mwaves,num_ghost,mx,ql,qr,auxl, &
                                        auxr,wave,s,amdq,apdq,num_aux)
                use precision_module
                integer :: ixy,maxm,meqn,mwaves,num_ghost,mx,num_aux
                real(kind=sk) :: wave(meqn, mwaves, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: s(mwaves, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: ql(meqn, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: qr(meqn, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: auxl(num_aux, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: auxr(num_aux, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: apdq(meqn, 1-num_ghost:maxm+num_ghost)
                real(kind=sk) :: amdq(meqn, 1-num_ghost:maxm+num_ghost)
            end subroutine rp_function
        end interface
        
        ! Riemann solver input arguments
        real(kind=sk), allocatable :: ql(:,:)
        real(kind=sk), allocatable :: qr(:,:)
        real(kind=sk), allocatable :: auxl(:,:)
        real(kind=sk), allocatable :: auxr(:,:)
        
        ! Output from the Riemann solver
        real(kind=sk), allocatable :: wave(:,:,:)
        real(kind=sk), allocatable :: s(:,:)
        real(kind=sk), allocatable :: apdq(:,:)
        real(kind=sk), allocatable :: amdq(:,:)
        
        ! Local storage
        real(kind=tk) :: total_time = 0.0_tk
        integer :: m,n,i,j,maxm

        ! Allocate actual arrays
        allocate(q(num_states,1-num_ghost:nx+num_ghost,1-num_ghost:ny+num_ghost))
        allocate(aux(num_aux,1-num_ghost:nx+num_ghost,1-num_ghost:ny+num_ghost))

        ! Setup arrays and problem
        maxm = max(nx,ny)
        allocate(ql(num_states,1-num_ghost:maxm+num_ghost))
        allocate(qr(num_states,1-num_ghost:maxm+num_ghost))
        allocate(auxl(num_states,1-num_ghost:maxm+num_ghost))
        allocate(auxr(num_states,1-num_ghost:maxm+num_ghost))
        allocate(wave(num_states, num_waves, 1-num_ghost:maxm+num_ghost))
        allocate(s(num_waves, 1-num_ghost:maxm+num_ghost))
        allocate(apdq(num_states, 1-num_ghost:maxm+num_ghost))
        allocate(amdq(num_states, 1-num_ghost:maxm+num_ghost))
        
        ! Initialize random number generator
        call init_random_seed()
        
        do n=1,num_tests
            ! Setup this test with random input
            call random_number(q)
            call random_number(aux)
            ql = q(:,:,1)
            qr = q(:,1,:)
            if (num_aux > 0) then
                auxl = aux(:,:,1)
                auxr = aux(:,1,:)
            endif
            
            call start_timer()
            do j=1-num_ghost,ny+num_ghost
                call rp_function(1,maxm,num_states,num_waves,num_ghost,nx,ql,qr, &
                                    auxl,auxr,wave,s,amdq,apdq,num_aux)
            enddo
            do i=1-num_ghost,nx+num_ghost
                call rp_function(2,maxm,num_states,num_waves,num_ghost,nx,ql,qr, &
                    auxl,auxr,wave,s,amdq,apdq,num_aux)
            enddo
            total_time = total_time + end_timer()
        enddo
        
        ! Output time (ms)
        ave_time = 1000.0_sk * (total_time / real(num_tests,kind=tk))
        
    end function time_rp_row3_function
    
    
    real(kind=tk) function time_rp_ptwise_function(rp_function,num_tests) &
                                            result(ave_time)

        use precision_module
        
        implicit none
        
        ! Input arguments
        integer, intent(in) :: num_tests

        ! Arrays for computin'
        real(kind=sk), allocatable :: q(:,:,:)
        real(kind=sk), allocatable :: aux(:,:,:)

        ! Input function interface
        interface
            subroutine rp_function(xy,num_states,num_aux,num_waves,num_ghost, &
                                        ql,qr,auxl,auxr,wave,s,amdq,apdq)
                use precision_module
                integer :: xy,num_states,num_aux,num_waves,num_ghost
                real(kind=sk) :: wave(num_states, num_waves)
                real(kind=sk) :: s(num_waves)
                real(kind=sk) :: ql(num_states)
                real(kind=sk) :: qr(num_states)
                real(kind=sk) :: auxl(num_aux)
                real(kind=sk) :: auxr(num_aux)
                real(kind=sk) :: apdq(num_states)
                real(kind=sk) :: amdq(num_states)
            end subroutine rp_function
        end interface
        
        ! Output from the Riemann solver
        real(kind=sk), allocatable :: wave(:,:)
        real(kind=sk), allocatable :: s(:)
        real(kind=sk), allocatable :: apdq(:)
        real(kind=sk), allocatable :: amdq(:)
        
        ! Local storage
        real(kind=tk) :: total_time = 0.0_tk
        integer :: n,i,j

        ! Allocate actual arrays
        allocate(q(num_states,1-num_ghost:nx+num_ghost,1-num_ghost:ny+num_ghost))
        allocate(aux(num_aux,1-num_ghost:nx+num_ghost,1-num_ghost:ny+num_ghost))

        ! Setup Riemann arrays and problem
        allocate(wave(num_states, num_waves))
        allocate(s(num_waves))
        allocate(apdq(num_states))
        allocate(amdq(num_states))
        
        ! Initialize random number generator
        call init_random_seed()
        
        do n=1,num_tests
            ! Setup this test with random input
            call random_number(q)
            call random_number(aux)
            
            call start_timer()
            do i=1-num_ghost,nx+num_ghost
                do j=1-num_ghost,ny+num_ghost
                    call rp_function(1,num_states,num_aux,num_waves,num_ghost, &
                                    q(:,i,j),q(:,i,j),aux(:,i,j),aux(:,i,j),wave,s,amdq,apdq)
                enddo
            enddo
            total_time = total_time + end_timer()
        enddo
        
        ! Output time (ms)
        ave_time = 1000.0_sk * (total_time / real(num_tests,kind=tk))
        
    end function time_rp_ptwise_function
    
    subroutine init_random_seed()
    
        implicit none
        
        integer :: i, n, clock
        integer, allocatable :: seed(:)

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(PUT = seed)

        deallocate(seed)
    end subroutine
    

end module timer_module
