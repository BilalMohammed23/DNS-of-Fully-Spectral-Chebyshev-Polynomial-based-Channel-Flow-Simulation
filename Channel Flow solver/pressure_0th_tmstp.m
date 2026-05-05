function [P, dPdx, dPdy] = pressure_0th_tmstp(Nx,k_x,H1,H2,theta,y)

    %pressure equation at 0th time-step:
    ghat = 0;
    fhat = 0;
    
    H1_hatk = fft(H1,[],2)/Nx;
    dH1dx_hatk = 1i .* k_x .* H1_hatk;
    dH1dx = real(ifft(dH1dx_hatk,[],2) * Nx);   %real is for safety
    
    % max_imag_dH1dx = max(abs(imag(dH1dx(:)))) 
    % max_real_dH1dx = max(abs(real(dH1dx(:))))
    
    a_H2 = zeros(Nx,Nx);
    b_H2 = zeros(Nx,Nx);
    dH2dy = zeros(Nx,Nx);
    
    for ix = 1: Nx
        a_H2(:,ix)  = a_cheb_coeff(H2(:,ix), theta);   
        b_H2(:,ix) = b_cheb_coeff(a_H2(:,ix), Nx-1); 
        dH2dy(:,ix) = cheb_eval_series(b_H2(:,ix), y);
    end
    
    q = -(dH1dx + dH2dy);
    
    q_hatk = fft(q,[],2)/Nx;
    
    alpha_k = k_x.^2;
    
    q_n = zeros(Nx,Nx);
    d_all = zeros(Nx,Nx);
    for ix = 1 : Nx
        q_n(:,ix) = a_cheb_coeff(q_hatk(:,ix),theta);
        [d_n] = solve_tridiag_dn(Nx,q_n(:,ix),alpha_k(ix),ghat,fhat);
        d_all(:,ix) = d_n;
    end
    
    P_hatk = zeros(Nx,Nx);
    
    for ix = 1 : Nx
        d_n = d_all(:,ix);
        P_hatk(:,ix) = cheb_eval_series(d_n,y);
    end
    
    P = ifft(P_hatk,[],2)*Nx;
    
    dPdx_hatk = 1i .* k_x .* P_hatk;
    dPdx = real(ifft(dPdx_hatk,[],2)*Nx);
    
    a_Phatk = zeros(Nx,Nx);
    b_Phatk = zeros(Nx,Nx);
    dPdy_hatk = zeros(Nx,Nx);
    for ix = 1: Nx
        a_Phatk(:,ix)  = a_cheb_coeff(P_hatk(:,ix), theta);
        b_Phatk(:,ix) = b_cheb_coeff(a_Phatk(:,ix), Nx-1);
        dPdy_hatk(:,ix) = cheb_eval_series(b_Phatk(:,ix), y);
    end
    dPdy = real(ifft(dPdy_hatk,[],2) * Nx);

end