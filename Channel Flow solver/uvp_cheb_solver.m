function [u_hatk,v_hatk,u,v,dPdx,dPdy] = uvp_cheb_solver(Nx,k_x,q_n_u,q_n_v,beta,gamma,y,theta)
    cheb_coeffs = zeros(3*Nx,Nx);
    
    for ix = 0 : Nx-1
        if ix==0
            [m] = global_matrix_general_k0(ix,Nx,k_x,q_n_u,q_n_v,beta,gamma);
            cheb_coeffs(:,ix+1) = m;
        else
            [m] = global_matrix_general_k(ix,Nx,k_x,q_n_u,q_n_v,beta,gamma);
            cheb_coeffs(:,ix+1) = m;
        end
    end
    
    a_n = zeros(Nx,Nx);
    b_n = zeros(Nx,Nx);
    d_n = zeros(Nx,Nx);
    
    for ii = 1:3*Nx
        for jj = 1:Nx
            if (ii<Nx+1)
                a_n(ii,jj) = cheb_coeffs(ii,jj);
            elseif (ii>Nx) && (ii<(2*Nx+1))
                b_n(ii-Nx,jj) = cheb_coeffs(ii,jj);
            elseif(ii>2*Nx)
                d_n(ii-2*Nx,jj) = cheb_coeffs(ii,jj);
            end
        end
    end
    
    u_hatk = zeros(Nx,Nx);
    v_hatk = zeros(Nx,Nx);
    P_hatk = zeros(Nx,Nx);
    
    for iy=1:Nx
        u_hatk(:,iy) = cheb_eval_series(a_n(:,iy),y);
        v_hatk(:,iy) = cheb_eval_series(b_n(:,iy),y);
        P_hatk(:,iy) = cheb_eval_series(d_n(:,iy),y);
    end
    
    u = real(ifft(u_hatk,[],2)*Nx);
    v = real(ifft(v_hatk,[],2)*Nx);
    
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