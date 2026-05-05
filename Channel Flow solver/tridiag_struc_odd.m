function [m2_odd_new, m1_odd] = tridiag_struc_odd(Nx,q_n,alpha,ghat,fhat,odd_spac,c_n)
    m2_odd = zeros(Nx , 1);
    m2_odd_new = zeros(length(odd_spac),1);
    m1_odd = zeros(Nx/2,Nx/2);

    for i = odd_spac
        %m2 outer loop
        c_n_1 = 1;
        c_n_2 = 1;

        if(i==1)
            m2_odd(i)= ghat + fhat;
        elseif(i==Nx-1)
            m2_odd(i)= ((c_n_1*q_n(i+1) - c_n_1*c_n_2*q_n(i-1))/(2*(i-1))) + ((c_n*q_n(i+1))/(2*(i+1)));
        else
            m2_odd(i)= ((c_n_1*q_n(i+1) - c_n_1*c_n_2*q_n(i-1))/(2*(i-1))) + ((c_n*q_n(i+1) - q_n(i+3))/(2*(i+1)));
        end

        %m1 inner loop
        i_n = ((i-1)/2) + 1;
        for j = 1 : Nx/2
            if i_n==1
                m1_odd(i_n,j) = 2*(((j-1)*2 + 1)^2);
            end
            if (i_n==j) && (i_n>1)
                m1_odd(i_n,j) = -(2*i + ((alpha*c_n_1)/(2*(i-1))) + ((alpha*c_n)/(2*(i+1))));
            elseif (i_n>1) && (j<i_n) && (j==i_n-1)
                m1_odd(i_n,j) = (alpha*c_n_2*c_n_1)/(2*(i-1));
            elseif (i_n>1) && (j>i_n) && (j==i_n+1)
                m1_odd(i_n,j) = alpha/(2*(i+1));
            end
        end
    end
    for z = 1:length(odd_spac)
        m2_odd_new(z) = m2_odd(odd_spac(z));
    end

end