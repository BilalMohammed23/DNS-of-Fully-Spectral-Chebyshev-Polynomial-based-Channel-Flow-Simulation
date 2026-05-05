function [q_n_u,q_n_v] = find_q(K,dPdx,gamma,d2udx2,d2udy2,H1_new,u,delta_t_max,dPdy,d2vdx2,d2vdy2,H2_new,v,Nx,theta,H1_old,H2_old)
    Q1 = - K - ((1/2).*dPdx) + (gamma/2).*(d2udx2 + d2udy2) - (((3/2)*H1_new) - ((1/2)*H1_old)) + (u/delta_t_max);
    Q2 = - ((1/2).*dPdy) + (gamma/2).*(d2vdx2 + d2vdy2) - (((3/2)*H2_new) - ((1/2)*H2_old)) + (v/delta_t_max);
    
    Q1_hatk = fft(Q1,[],2)/Nx;
    Q2_hatk = fft(Q2,[],2)/Nx;
    
    q_n_u = zeros(Nx,Nx);
    q_n_v = zeros(Nx,Nx);
    
    for ix = 1:Nx
        q_n_u(:,ix) = a_cheb_coeff(Q1_hatk(:,ix),theta);     %want this
        q_n_v(:,ix) = a_cheb_coeff(Q2_hatk(:,ix),theta);     %want this
    end
end