module test_functions_module
    
    implicit none
    
contains

    ! =========================================================================
    !  Dummy Functions
    subroutine row_dummy_function(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr, &
                                    wave,s,amdq,apdq,num_aux)

        use precision_module
        implicit none

        integer :: ixy,maxm,meqn,mwaves,mbc,mx,num_aux

        real(kind=sk), intent(inout) :: wave(meqn, mwaves, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: s(mwaves, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: ql(meqn, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: qr(meqn, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: auxl(num_aux, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: auxr(num_aux, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: apdq(meqn, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: amdq(meqn, 1-mbc:maxm+mbc)
        
        call sleep(1)
    end subroutine row_dummy_function
    
    subroutine ptwise_dummy_function(num_states,num_aux,num_waves,num_ghost, &
                                        ql,qr,auxl,auxr,wave,s,amdq,apdq)

        use precision_module
        implicit none
        
        integer :: num_states,num_aux,num_waves,num_ghost
        real(kind=sk) :: wave(num_states, num_waves)
        real(kind=sk) :: s(num_waves)
        real(kind=sk) :: ql(num_states)
        real(kind=sk) :: qr(num_states)
        real(kind=sk) :: auxl(num_aux)
        real(kind=sk) :: auxr(num_aux)
        real(kind=sk) :: apdq(num_states)
        real(kind=sk) :: amdq(num_states)
        
        call sleep(1)
        
    end subroutine ptwise_dummy_function
    ! =========================================================================
    
    
    ! =========================================================================
    !  Constant coefficient advection Riemann solvers
    subroutine row_rp_advection(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr, &
                                    wave,s,amdq,apdq,num_aux)

        use precision_module
        implicit none

        integer :: ixy,maxm,meqn,mwaves,mbc,mx,num_aux

        real(kind=sk), intent(inout) :: wave(meqn, mwaves, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: s(mwaves, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: ql(meqn, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: qr(meqn, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: auxl(num_aux, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: auxr(num_aux, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: apdq(meqn, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: amdq(meqn, 1-mbc:maxm+mbc)

        integer :: i
        
        real(kind=sk) :: u,v
        common /advection_rp/ u,v

        do i = 2-mbc, mx+mbc
            wave(1,1,i) = ql(1,i) - qr(1,i-1)
            if (ixy == 1) then
                s(1,i) = u
            else
                s(1,i) = v
            endif
    
            ! Flux differences:
            amdq(1,i) = dmin1(s(1,i), 0.d0) * wave(1,1,i)
            apdq(1,i) = dmax1(s(1,i), 0.d0) * wave(1,1,i)
    
        enddo
        
    end subroutine row_rp_advection
    
    
    subroutine ptwise_rp_advection(xy,num_states,num_aux,num_waves,num_ghost, &
                                        ql,qr,auxl,auxr,wave,s,amdq,apdq)

        use precision_module
        implicit none
        
        integer :: xy,num_states,num_aux,num_waves,num_ghost
        real(kind=sk) :: wave(num_states, num_waves)
        real(kind=sk) :: s(num_waves)
        real(kind=sk) :: ql(num_states)
        real(kind=sk) :: qr(num_states)
        real(kind=sk) :: auxl(num_aux)
        real(kind=sk) :: auxr(num_aux)
        real(kind=sk) :: apdq(num_states)
        real(kind=sk) :: amdq(num_states)
        
        real(kind=sk) :: u,v
        common /advection_rp/ u,v
        
        wave(1,1) = ql(1) - qr(1)
        if (xy == 1) then
            s(1) = u
        else
            s(1) = v
        endif
    
        ! Flux differences:
        amdq(1) = min(s(1), 0.d0) * wave(1,1)
        apdq(1) = max(s(1), 0.d0) * wave(1,1)
        
    end subroutine ptwise_rp_advection
    ! =========================================================================
    

    ! =========================================================================
    !  Euler Riemann solvers
    subroutine row_rp_euler(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr, &
                                    wave,s,amdq,apdq,num_aux)

        use precision_module
        implicit none

        integer :: ixy,maxm,meqn,mwaves,mbc,mx,num_aux

        real(kind=sk), intent(inout) :: wave(meqn, mwaves, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: s(mwaves, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: ql(meqn, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: qr(meqn, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: auxl(num_aux, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: auxr(num_aux, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: apdq(meqn, 1-mbc:maxm+mbc)
        real(kind=sk), intent(inout) :: amdq(meqn, 1-mbc:maxm+mbc)
    
        ! Local storage
        integer :: mu,mv,i,m,mw
        double precision :: rhsqrtl,rhsqrtr,rhsq2,pl,pr,a1,a2,a3,a4
        double precision, dimension(1-mbc:mx+mbc) :: u,v,enth,u2v2,a,g1a2,euv
        double precision :: delta(4)
        logical, parameter :: entropy_fix = .false.
        
        real(kind=sk) :: gamma,gamma1
        common /euler_rp/  gamma,gamma1
    
        if (ixy == 1) then
            mu = 2
            mv = 3
        else
            mu = 3
            mv = 2
        endif
        
        ! Calculate 
        do i = 2-mbc, mx+mbc
            rhsqrtl = sqrt(qr(1,i-1))
            rhsqrtr = sqrt(ql(1,i))
            pl = gamma1*(qr(4,i-1) - 0.5d0*(qr(2,i-1)**2 + &
            qr(3,i-1)**2)/qr(1,i-1))
            pr = gamma1*(ql(4,i) - 0.5d0*(ql(2,i)**2 + &
            ql(3,i)**2)/ql(1,i))
            rhsq2 = rhsqrtl + rhsqrtr
            u(i) = (qr(mu,i-1)/rhsqrtl + ql(mu,i)/rhsqrtr) / rhsq2
            v(i) = (qr(mv,i-1)/rhsqrtl + ql(mv,i)/rhsqrtr) / rhsq2
            enth(i) = (((qr(4,i-1)+pl)/rhsqrtl &
            + (ql(4,i)+pr)/rhsqrtr)) / rhsq2
            u2v2(i) = u(i)**2 + v(i)**2
            a2 = gamma1*(enth(i) - .5d0*u2v2(i))
            a(i) = sqrt(a2)
            g1a2(i) = gamma1 / a2
            euv(i) = enth(i) - u2v2(i)
        enddo
        
        ! Primary loop 
        do i = 2-mbc, mx+mbc
            delta(1) = ql(1,i) - qr(1,i-1)
            delta(2) = ql(mu,i) - qr(mu,i-1)
            delta(3) = ql(mv,i) - qr(mv,i-1)
            delta(4) = ql(4,i) - qr(4,i-1)
            a3 = g1a2(i) * (euv(i)*delta(1) &
            + u(i)*delta(2) + v(i)*delta(3) - delta(4))
            a2 = delta(3) - v(i)*delta(1)
            a4 = (delta(2) + (a(i)-u(i))*delta(1) - a(i)*a3) / (2.d0*a(i))
            a1 = delta(1) - a3 - a4
    
            ! Computer the waves
            ! Acoustic
            wave(1,1,i) = a1
            wave(mu,1,i) = a1*(u(i)-a(i))
            wave(mv,1,i) = a1*v(i)
            wave(4,1,i) = a1*(enth(i) - u(i)*a(i))
            s(1,i) = u(i)-a(i)
    
            ! Shear
            wave(1,2,i) = 0.d0
            wave(mu,2,i) = 0.d0
            wave(mv,2,i) = a2
            wave(4,2,i) = a2*v(i)
            s(2,i) = u(i)
    
            ! Entropy
            wave(1,3,i) = a3
            wave(mu,3,i) = a3*u(i)
            wave(mv,3,i) = a3*v(i)
            wave(4,3,i) = a3*0.5d0*u2v2(i)
            s(3,i) = u(i)
    
            ! Acoustic
            wave(1,4,i) = a4
            wave(mu,4,i) = a4*(u(i)+a(i))
            wave(mv,4,i) = a4*v(i)
            wave(4,4,i) = a4*(enth(i)+u(i)*a(i))
            s(4,i) = u(i)+a(i)
        enddo

        if (.not.entropy_fix) then
            amdq = 0.d0
            apdq = 0.d0
            do m=1,meqn
                do i=2-mbc,mx+mbc
                    do mw=1,mwaves
                        if (s(mw,i) < 0.d0) then
                            amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                        else
                            apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                        endif
                    enddo
                enddo
            enddo
        else
            stop "Entropy fix not implemented!"
        endif
            
    end subroutine row_rp_euler
    
    
    subroutine ptwise_rp_euler(xy,num_states,num_aux,num_waves,num_ghost, &
                                        ql,qr,auxl,auxr,wave,s,amdq,apdq)

        use precision_module
        implicit none
        
        integer :: xy,num_states,num_aux,num_waves,num_ghost
        real(kind=sk) :: wave(num_states, num_waves)
        real(kind=sk) :: s(num_waves)
        real(kind=sk) :: ql(num_states)
        real(kind=sk) :: qr(num_states)
        real(kind=sk) :: auxl(num_aux)
        real(kind=sk) :: auxr(num_aux)
        real(kind=sk) :: apdq(num_states)
        real(kind=sk) :: amdq(num_states)
        
        ! Local storage
        integer :: mu,mv,i,m,mw
        double precision :: rhsqrtl,rhsqrtr,rhsq2,pl,pr,a1,a2,a3,a4
        double precision :: u,v,enth,u2v2,a,g1a2,euv
        double precision :: delta(4)
        logical, parameter :: entropy_fix = .false.
        
        real(kind=sk) :: gamma,gamma1
        common /euler_rp/  gamma,gamma1
    
        if (xy == 1) then
            mu = 2
            mv = 3
        else
            mu = 3
            mv = 2
        endif
        
        ! Calculate Roe variables
        rhsqrtl = sqrt(qr(1))
        rhsqrtr = sqrt(ql(1))
        pl = gamma1*(qr(4) - 0.5d0*(qr(2)**2 + &
        qr(3)**2)/qr(1))
        pr = gamma1*(ql(4) - 0.5d0*(ql(2)**2 + &
        ql(3)**2)/ql(1))
        rhsq2 = rhsqrtl + rhsqrtr
        u = (qr(mu)/rhsqrtl + ql(mu)/rhsqrtr) / rhsq2
        v = (qr(mv)/rhsqrtl + ql(mv)/rhsqrtr) / rhsq2
        enth = (((qr(4)+pl)/rhsqrtl &
        + (ql(4)+pr)/rhsqrtr)) / rhsq2
        u2v2 = u**2 + v**2
        a2 = gamma1*(enth - .5d0*u2v2)
        a = sqrt(a2)
        g1a2 = gamma1 / a2
        euv = enth - u2v2
        
        ! Calculate waves
        delta(1) = ql(1) - qr(1)
        delta(2) = ql(mu) - qr(mu)
        delta(3) = ql(mv) - qr(mv)
        delta(4) = ql(4) - qr(4)
        a3 = g1a2 * (euv*delta(1) + u*delta(2) + v*delta(3) - delta(4))
        a2 = delta(3) - v*delta(1)
        a4 = (delta(2) + (a-u)*delta(1) - a*a3) / (2.d0*a)
        a1 = delta(1) - a3 - a4
    
        ! Acoustic
        wave(1,1) = a1
        wave(mu,1) = a1*(u-a)
        wave(mv,1) = a1*v
        wave(4,1) = a1*(enth - u*a)
        s(1) = u-a
    
        ! Shear
        wave(1,2) = 0.d0
        wave(mu,2) = 0.d0
        wave(mv,2) = a2
        wave(4,2) = a2*v
        s(2) = u
    
        ! Entropy
        wave(1,3) = a3
        wave(mu,3) = a3*u
        wave(mv,3) = a3*v
        wave(4,3) = a3*0.5d0*u2v2
        s(3) = u
    
        ! Acoustic
        wave(1,4) = a4
        wave(mu,4) = a4*(u+a)
        wave(mv,4) = a4*v
        wave(4,4) = a4*(enth+u*a)
        s(4) = u+a

        if (.not.entropy_fix) then
            amdq = 0.d0
            apdq = 0.d0
            do m=1,num_states
                do mw=1,num_waves
                    if (s(mw) < 0.d0) then
                        amdq(m) = amdq(m) + s(mw)*wave(m,mw)
                    else
                        apdq(m) = apdq(m) + s(mw)*wave(m,mw)
                    endif
                enddo
            enddo
        else
            stop "Entropy fix not implemented!"
        endif
        
    
    end subroutine ptwise_rp_euler
    
end module test_functions_module