clear
close all
clc

%%

Nx = 128;

L = 2*pi;
h = L/Nx;
j = (0:Nx-1);
x = j*h;
z = x;

theta = j*pi/(Nx-1);
y = cos(theta);

% [Y, X, Z] = ndgrid(y, x, z);%used for our ease of x and y dimensioning
[X, Y, Z] = ndgrid(x, y, z);

% Q = (-2 -8 .* pi.^2 .* (1 - Y.^2)).*(cos(2.*pi.*X) .* cos(2.*pi.*Z));
% Q = (-2 -8 .* pi.^2 .* (1 - Y.^2)).*(cos(2.*pi.*X));
Q = (-2 -8 .* pi.^2 .* (1 - Y.^2)).*(cos(X));

Q_fft = fft(fft(Q,[],1),[],3)/(Nx*Nx);

ghat = 0;
fhat = 0;

%%

l = [0 : (Nx/2)-1 , -Nx/2 : -1];
m = l;
n = l;

k_x = (2*pi/L)*l;
k_z = (2*pi/L)*n;

alpha = zeros(Nx,Nx);

for ik = 1:Nx
    for iz = 1:Nx
        alpha(ik,iz) = k_x(ik)^2 + k_z(iz)^2;
    end
end

a_all = zeros(Nx, Nx, Nx);
for ik = 1:Nx
    for iz = 1:Nx
        qhat_y = Q_fft(ik,:,iz);
        q_n = qn_cheb_coeff_RK4(qhat_y, theta);
        [a_n] = solve_tridiag_an(Nx,q_n,alpha(ik,iz),ghat,fhat);
        a_all(ik, :, iz) = a_n;
    end
end

%%

Phat = zeros(Nx, Nx, Nx);

for ik = 1:Nx
    for iz = 1:Nx
        a_n = a_all(ik,:,iz);   % coefficients for this Fourier mode
        Phat(ik,:,iz) = cheb_eval_series_RK4(a_n, y);
    end
end

P = ifft(ifft(Phat, [], 1), [], 3, 'symmetric');

%%
% 2D visualization

kz_mid = round(Nx/2);

P_xy = squeeze(P(:,:,kz_mid));   % size: x-by-y

figure
imagesc(x, y, P_xy.')
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('y')
title(['Pressure at z = ', num2str(z(kz_mid))])