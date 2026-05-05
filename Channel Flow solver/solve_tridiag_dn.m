function [a_n] = solve_tridiag_dn(Nx,q_n,alpha,ghat,fhat)
    
    c_n = 1;
    a_n = zeros(1,Nx);
    even_spac = 0:2:Nx-1;
    odd_spac = 1:2:Nx-1;
    
    [m2_even, m1_even] = tridiag_struc_even(Nx,q_n,alpha,ghat,fhat,even_spac,c_n);
    a_n_even = m1_even\m2_even;
    a_n_even = a_n_even.';

    [m2_odd, m1_odd] = tridiag_struc_odd(Nx,q_n,alpha,ghat,fhat,odd_spac,c_n);
    a_n_odd = m1_odd\m2_odd;
    a_n_odd = a_n_odd.';
    
    for i = even_spac
        a_n(i+1) = a_n_even(i/2 + 1);
    end

    for j = odd_spac
        a_n(j+1) = a_n_odd((j-1)/2 + 1);
    end
    
end