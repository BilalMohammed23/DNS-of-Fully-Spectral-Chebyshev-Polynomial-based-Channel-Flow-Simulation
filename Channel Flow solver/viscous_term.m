function [d2udx2, d2vdx2, d2udy2, d2vdy2] = viscous_term(k_x,uhatk,vhatk,Nx,theta,y)
    d2udx2_hatk = -(k_x.^2) .* uhatk;
    d2vdx2_hatk = -(k_x.^2) .* vhatk;
    
    d2udx2 = real(ifft(d2udx2_hatk,[],2) * Nx);
    d2vdx2 = real(ifft(d2vdx2_hatk,[],2) * Nx);
    
    d2udy2_hatk = zeros(Nx,Nx);
    d2vdy2_hatk = zeros(Nx,Nx);
    
    for ix = 1:Nx
        % uyy 
        a_u = a_cheb_coeff(uhatk(:,ix), theta);
        b_u = b_cheb_coeff(a_u, Nx-1);       
        c_u = b_cheb_coeff(b_u, Nx-1);       
        d2udy2_hatk(:,ix) = cheb_eval_series(c_u, y);
    
        % vyy
        a_v = a_cheb_coeff(vhatk(:,ix), theta);
        b_v = b_cheb_coeff(a_v, Nx-1);      
        c_v = b_cheb_coeff(b_v, Nx-1);     
        d2vdy2_hatk(:,ix) = cheb_eval_series(c_v, y);
    end
    
    d2udy2 = real(ifft(d2udy2_hatk,[],2) * Nx);
    d2vdy2 = real(ifft(d2vdy2_hatk,[],2) * Nx);
end