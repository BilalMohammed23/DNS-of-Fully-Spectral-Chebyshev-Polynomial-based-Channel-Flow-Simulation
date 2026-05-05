function [q_n_u,q_n_v] = find_q_0(K,dPdx,gamma,d2udx2,d2udy2,H1,u0,delt_t_initial,dPdy,d2vdx2,d2vdy2,H2,v0,Nx,theta)
    Q1_0 = - K - ((1/2).*dPdx) + (gamma/2).*(d2udx2 + d2udy2) - H1 + (u0/delt_t_initial);
    Q2_0 = - ((1/2).*dPdy) + (gamma/2).*(d2vdx2 + d2vdy2) - H2 + (v0/delt_t_initial);
    
    Q1 = Q1_0;
    Q2 = Q2_0;
    
    Q1_hatk = fft(Q1,[],2)/Nx;
    Q2_hatk = fft(Q2,[],2)/Nx;
    
    q_n_u = zeros(Nx,Nx);
    q_n_v = zeros(Nx,Nx);
    
    for ix = 1:Nx
        q_n_u(:,ix) = a_cheb_coeff(Q1_hatk(:,ix),theta);     %want this
        q_n_v(:,ix) = a_cheb_coeff(Q2_hatk(:,ix),theta);     %want this
    end
end