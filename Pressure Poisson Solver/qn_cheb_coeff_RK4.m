function q_n = qn_cheb_coeff_RK4(qhat_y, theta)
    N = length(qhat_y)-1;
    q_n = zeros(1,N+1);

    for n = 0:N
        summ = 0;
        for j = 0:N
            c_j = 1;
            if (j==0 || j==N) 
                c_j = 2; 
            end
            summ = summ + ( (1/c_j) * qhat_y(j+1) * cos(n*theta(j+1)) );
        end
        c_n = 1;
        if (n==0 || n==N) 
            c_n = 2; 
        end
        q_n(n+1) = (2/(c_n*N)) * summ;
    end
end