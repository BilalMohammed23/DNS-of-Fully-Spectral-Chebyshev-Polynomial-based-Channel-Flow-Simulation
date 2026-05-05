clear
close all
clc

Nx = 32;

L = 2*pi;
h = L/Nx;
j = (0:Nx-1);
x = j*h;
dx = x(2) - x(1);
l = [0 : (Nx/2)-1 , -Nx/2 : -1];
k_x = (2*pi/L)*l;

theta = j*pi/(Nx-1);
y = cos(theta);
dy_min = y(1) - y(2);

[X , Y] = meshgrid(x,y);

Re = 350;
K = -2/Re;
gamma = 1/Re;

eps = 0.01;
A = 1;
m = 1;
alpha = 2*pi/(L);
u = (1 - Y.^2) - 4.*eps.*A.*Y.*(1 - Y.^2).*sin(alpha.*X);
v = - eps.*A.*alpha.*((1 - Y.^2).^2).*cos(alpha.*X);

%%
%0th time-step:
u0 = u;
v0 = v;

uhatk_0 = fft(u0,[],2)/Nx;
vhatk_0 = fft(v0,[],2)/Nx;

max_u = max(abs(u(:)));
max_v = max(abs(v(:)));

%non-linear psuedo-spectral method:
[H1_new , H2_new, dudx_hatk, dvdy_hatk] = H_non_linear(k_x,theta,y,Nx, uhatk_0,vhatk_0);

%viscous terms:
[d2udx2, d2vdx2, d2udy2, d2vdy2] = viscous_term(k_x,uhatk_0,vhatk_0,Nx,theta,y);

%continuity check:
dudx = real(ifft(dudx_hatk,[],2) * Nx);     %real is for safety
dvdy = real(ifft(dvdy_hatk,[],2) * Nx);     %real is for safety

COM_0tmstp = dudx + dvdy;
max_COM_0tmstp = max(abs(COM_0tmstp(:)));

%Pressure at 0th timestep:
[P_0, dP0dx, dP0dy] = pressure_0th_tmstp(Nx,k_x,H1_new,H2_new,theta,y);

%%
% time advancement
% from here directly go ahead with the coupled Velocity-Pressure solver

T = 10;
t = 0;
count = 0;

%% live 2D contour setup
figure(1);
vel_mag0 = sqrt(u0.^2 + v0.^2);
% vel_mag0 = u0;
hcont = contourf(X, Y, vel_mag0, 20, 'LineColor', 'none');
colormap(jet);
colorbar;
xlabel('x');
ylabel('y');
title(sprintf('Velocity magnitude, t = %.4f', t));
drawnow;

%%
while t < T
    
    if t==0
        delta_t_max = dt_channel(max_u, max_v, dx, dy_min, gamma);
        if t + delta_t_max > T
            delta_t_max = T - t;
        end

        beta = 2/(gamma*delta_t_max);
        
        [q_n_u_0,q_n_v_0] = find_q_0(K,dP0dx,gamma,d2udx2,d2udy2,H1_new,u0,delta_t_max,dP0dy,d2vdx2,d2vdy2,H2_new,v0,Nx,theta);
        
        H1_old = H1_new;
        H2_old = H2_new;

        [u_hatk,v_hatk,u,v,dPdx,dPdy] = uvp_cheb_solver(Nx,k_x,q_n_u_0,q_n_v_0,beta,gamma,y,theta);
        
        max_u = max(abs(u(:)));
        max_v = max(abs(v(:)));

        % Continuity check:
        dudx_hatk = 1i .* k_x .* u_hatk;
        a_vhatk = zeros(Nx,Nx);
        b_vhatk = zeros(Nx,Nx);
        dvdy_hatk = zeros(Nx,Nx);
        for ix = 1: Nx
            a_vhatk(:,ix)  = a_cheb_coeff(v_hatk(:,ix), theta);
            b_vhatk(:,ix) = b_cheb_coeff(a_vhatk(:,ix), Nx-1);
            dvdy_hatk(:,ix) = cheb_eval_series(b_vhatk(:,ix), y);
        end

        dudx = real(ifft(dudx_hatk,[],2) * Nx);     %real is for safety
        dvdy = real(ifft(dvdy_hatk,[],2) * Nx);     %real is for safety

        COM_tmstp = dudx + dvdy;
        max_COM_tmstp = max(abs(COM_tmstp(:)));

        %%
        % Live contour update
        vel_mag = sqrt(u.^2 + v.^2);
        % vel_mag = u;
        cla;
        contourf(X, Y, vel_mag, 20, 'LineColor', 'none');
        colormap(jet);
        colorbar;
        xlabel('x');
        ylabel('y');
        title(sprintf('Velocity magnitude, t = %.4f', t + delta_t_max));
        drawnow;

        %%
        t = t + delta_t_max;
        count = count + 1;
        fprintf('%4d: t = %.4f, dt = %.3e, max div = %.3e\n',count,t,delta_t_max,max_COM_tmstp)
    else
        delta_t_max = dt_channel(max_u, max_v, dx, dy_min, gamma);
        if t + delta_t_max > T
            delta_t_max = T - t;
        end
        
        beta = 2/(gamma*delta_t_max);    

        % Non-Linear terms:
        [H1_new , H2_new, dudx_hatk, dvdy_hatk] = H_non_linear(k_x,theta,y,Nx, u_hatk,v_hatk);
        
        % Viscous terms:
        [d2udx2, d2vdx2, d2udy2, d2vdy2] = viscous_term(k_x,u_hatk,v_hatk,Nx,theta,y);
        
        % Time-advanced q terms
        [q_n_u,q_n_v] = find_q(K,dPdx,gamma,d2udx2,d2udy2,H1_new,u,delta_t_max,dPdy,d2vdx2,d2vdy2,H2_new,v,Nx,theta,H1_old,H2_old);
        
        [u_hatk,v_hatk,u,v,dPdx,dPdy] = uvp_cheb_solver(Nx,k_x,q_n_u,q_n_v,beta,gamma,y,theta);

        max_u = max(abs(u(:)));
        max_v = max(abs(v(:)));
        
        % Continuity check:
        dudx_hatk = 1i .* k_x .* u_hatk;
        a_vhatk = zeros(Nx,Nx);
        b_vhatk = zeros(Nx,Nx);
        dvdy_hatk = zeros(Nx,Nx);
        for ix = 1: Nx
            a_vhatk(:,ix)  = a_cheb_coeff(v_hatk(:,ix), theta);
            b_vhatk(:,ix) = b_cheb_coeff(a_vhatk(:,ix), Nx-1);
            dvdy_hatk(:,ix) = cheb_eval_series(b_vhatk(:,ix), y);
        end

        dudx = real(ifft(dudx_hatk,[],2) * Nx);     %real is for safety
        dvdy = real(ifft(dvdy_hatk,[],2) * Nx);     %real is for safety

        COM_tmstp = dudx + dvdy;
        max_COM_tmstp = max(abs(COM_tmstp(:)));

        %%
        % Live contour update
        vel_mag = sqrt(u.^2 + v.^2);
        % vel_mag = u;
        cla;
        contourf(X, Y, vel_mag, 20, 'LineColor', 'none');
        colormap(jet);
        colorbar;
        xlabel('x');
        ylabel('y');
        title(sprintf('Velocity magnitude, t = %.4f', t + delta_t_max));
        drawnow;

        %%
        H1_old = H1_new;
        H2_old = H2_new;
        
        t = t + delta_t_max;
        count = count + 1;
        
        fprintf('%4d: t = %.4f, dt = %.3e, max div = %.3e\n',count,t,delta_t_max,max_COM_tmstp)
    end

end