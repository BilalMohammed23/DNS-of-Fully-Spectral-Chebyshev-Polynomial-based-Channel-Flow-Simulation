function [H1,H2,dudx_hatk, dvdy_hatk] = H_non_linear(k_x,theta,y,Nx,uhatk,vhatk)
    
    %theory:
    %======
    %H1 = u (du/dx) + v (du/dy)
    %H2 = u (dv/dx) + v (dv/dy)

    dudx_hatk = 1i .* k_x .* uhatk;
    
    a_uhatk = zeros(Nx,Nx);
    b_uhatk = zeros(Nx,Nx);
    dudy_hatk = zeros(Nx,Nx);
    for ix = 1: Nx
        a_uhatk(:,ix)  = a_cheb_coeff(uhatk(:,ix), theta);
        b_uhatk(:,ix) = b_cheb_coeff(a_uhatk(:,ix), Nx-1);
        dudy_hatk(:,ix) = cheb_eval_series(b_uhatk(:,ix), y);
    end
    
    dvdx_hatk = 1i .* k_x .* vhatk;
    
    a_vhatk = zeros(Nx,Nx);
    b_vhatk = zeros(Nx,Nx);
    dvdy_hatk = zeros(Nx,Nx);
    for ix = 1: Nx
        a_vhatk(:,ix)  = a_cheb_coeff(vhatk(:,ix), theta);
        b_vhatk(:,ix) = b_cheb_coeff(a_vhatk(:,ix), Nx-1);
        dvdy_hatk(:,ix) = cheb_eval_series(b_vhatk(:,ix), y);
    end
    
    %%
    % psuedospectral
    Ny = size(uhatk,1); %no of rows
    Nx = size(uhatk,2); %no of columns
    Nx_pssp = (3/2)*Nx;
    
    uhatk_pssp = zeros(Ny, Nx_pssp);
    vhatk_pssp = zeros(Ny, Nx_pssp);
    
    dudx_hatk_pssp = zeros(Ny, Nx_pssp);
    dudy_hatk_pssp = zeros(Ny, Nx_pssp);
    
    dvdx_hatk_pssp = zeros(Ny, Nx_pssp);
    dvdy_hatk_pssp = zeros(Ny, Nx_pssp);
    % uhatk = rand(16); test case check if pseudospectral works fine
    
    %pad to 3N/2 grid: 
    % positive x modes
    uhatk_pssp(:,1:Nx/2) = uhatk(:,1:Nx/2);
    vhatk_pssp(:,1:Nx/2) = vhatk(:,1:Nx/2);
    
    dudx_hatk_pssp(:,1:Nx/2) = dudx_hatk(:,1:Nx/2);
    dudy_hatk_pssp(:,1:Nx/2) = dudy_hatk(:,1:Nx/2);
    
    dvdx_hatk_pssp(:,1:Nx/2) = dvdx_hatk(:,1:Nx/2);
    dvdy_hatk_pssp(:,1:Nx/2) = dvdy_hatk(:,1:Nx/2);
    
    % negative x modes
    uhatk_pssp(:,Nx+1:end) = uhatk(:,Nx/2 + 1:end);
    vhatk_pssp(:,Nx+1:end) = vhatk(:,Nx/2 + 1:end);
    
    dudx_hatk_pssp(:,Nx+1:end) = dudx_hatk(:,Nx/2 + 1:end);
    dudy_hatk_pssp(:,Nx+1:end) = dudy_hatk(:,Nx/2 + 1:end);
    
    dvdx_hatk_pssp(:,Nx+1:end) = dvdx_hatk(:,Nx/2 + 1:end);
    dvdy_hatk_pssp(:,Nx+1:end) = dvdy_hatk(:,Nx/2 + 1:end);
    
    %Inverse FFT
    u_pssp_if = real(ifft(uhatk_pssp, [], 2) * Nx_pssp);
    v_pssp_if = real(ifft(vhatk_pssp, [], 2) * Nx_pssp);  %
    dudx_pssp_if = ifft(dudx_hatk_pssp, [], 2) * Nx_pssp;
    dudy_pssp_if = ifft(dudy_hatk_pssp, [], 2) * Nx_pssp;
    dvdx_pssp_if = real(ifft(dvdx_hatk_pssp, [], 2) * Nx_pssp); %
    dvdy_pssp_if = real(ifft(dvdy_hatk_pssp, [], 2) * Nx_pssp); %

% max_imag_v      = max(abs(imag(v_pssp_if(:))))
% max_imag_dvdx   = max(abs(imag(dvdx_pssp_if(:))))
% max_imag_dvdy   = max(abs(imag(dvdy_pssp_if(:))))
% 
% max_real_v      = max(abs(real(v_pssp_if(:))))
% max_real_dvdx   = max(abs(real(dvdx_pssp_if(:))))
% max_real_dvdy   = max(abs(real(dvdy_pssp_if(:))))
    
    %non linear product
    H1_1_pssp = u_pssp_if .* dudx_pssp_if;
    H1_2_pssp = v_pssp_if .* dudy_pssp_if;
    
    H2_1_pssp = u_pssp_if .* dvdx_pssp_if;
    H2_2_pssp = v_pssp_if .* dvdy_pssp_if;
    
    %fft of non linear product
    H1_1_pssp_f = fft(H1_1_pssp,[],2)/Nx_pssp;
    H1_2_pssp_f = fft(H1_2_pssp,[],2)/Nx_pssp;
    
    H2_1_pssp_f = fft(H2_1_pssp,[],2)/Nx_pssp;
    H2_2_pssp_f = fft(H2_2_pssp,[],2)/Nx_pssp;
    
    postprod_H1_1 = zeros(Nx,Nx);
    postprod_H1_2 = zeros(Nx,Nx);
    
    postprod_H2_1 = zeros(Nx,Nx);
    postprod_H2_2 = zeros(Nx,Nx);
    
    %truncate to N grid:
    %positive quadrant
    postprod_H1_1(:,1:Nx/2) = H1_1_pssp_f(:,1:Nx/2);
    postprod_H1_2(:,1:Nx/2) = H1_2_pssp_f(:,1:Nx/2);
    
    postprod_H2_1(:,1:Nx/2) = H2_1_pssp_f(:,1:Nx/2);
    postprod_H2_2(:,1:Nx/2) = H2_2_pssp_f(:,1:Nx/2);
    %negative quadrant
    postprod_H1_1(:,Nx/2 + 1:end) = H1_1_pssp_f(:,Nx+1:end);
    postprod_H1_2(:,Nx/2 + 1:end) = H1_2_pssp_f(:,Nx+1:end);
    
    postprod_H2_1(:,Nx/2 + 1:end) = H2_1_pssp_f(:,Nx+1:end);
    postprod_H2_2(:,Nx/2 + 1:end) = H2_2_pssp_f(:,Nx+1:end);
    
    %ifft in N grid:
    H1_1 = real(ifft(postprod_H1_1, [], 2) * Nx);
    H1_2 = real(ifft(postprod_H1_2, [], 2) * Nx);
    
    H2_1 = real(ifft(postprod_H2_1, [], 2) * Nx);
    H2_2 = real(ifft(postprod_H2_2, [], 2) * Nx);
    
% max_imag_H1_1 = max(abs(imag(H1_1(:))))
% max_imag_H1_2 = max(abs(imag(H1_2(:))))
% max_imag_H2_1 = max(abs(imag(H2_1(:))))
% max_imag_H2_2 = max(abs(imag(H2_2(:))))
% 
% max_real_H1_1 = max(abs(real(H1_1(:))))
% max_real_H1_2 = max(abs(real(H1_2(:))))
% max_real_H2_1 = max(abs(real(H2_1(:))))
% max_real_H2_2 = max(abs(real(H2_2(:))))

    H1 = H1_1 + H1_2;
    H2 = H2_1 + H2_2;

end