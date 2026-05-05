function [m2_even_new, m1_even] = tridiag_struc_even(Nx,q_n,alpha,ghat,fhat,even_spac,c_n)
    m2_even = zeros(Nx , 1);
    m2_even_new = zeros(length(even_spac),1);
    m1_even = zeros(Nx/2,Nx/2);

    for i = even_spac
        %m2 outer loop
        c_n_1 = 1;
        c_n_2 = 1;
        
        if i-2 == 0
            c_n_2 = 2;
        end
        if alpha == 0
            if(i==0)
                m2_even(i+1)= ghat - fhat;  
            elseif(i==2)
                m2_even(i+1)= 0; %any constant
            elseif(i==Nx-2)
                m2_even(i+1)= ((c_n_1*q_n(i+1) - c_n_1*c_n_2*q_n(i-1))/(2*(i-1))) + ((c_n*q_n(i+1))/(2*(i+1)));
            else
                m2_even(i+1)= ((c_n_1*q_n(i+1) - c_n_1*c_n_2*q_n(i-1))/(2*(i-1))) + ((c_n*q_n(i+1) - q_n(i+3))/(2*(i+1)));
            end
        else
            if(i==0)
                m2_even(i+1)= ghat - fhat;
            elseif(i==Nx-2)
                m2_even(i+1)= ((c_n_1*q_n(i+1) - c_n_1*c_n_2*q_n(i-1))/(2*(i-1))) + ((c_n*q_n(i+1))/(2*(i+1)));
            else
                m2_even(i+1)= ((c_n_1*q_n(i+1) - c_n_1*c_n_2*q_n(i-1))/(2*(i-1))) + ((c_n*q_n(i+1) - q_n(i+3))/(2*(i+1)));
            end    
        end
        %m1 inner loop
        i_n = (i/2) + 1;
        for j = 1 : Nx/2
            if alpha==0
                if i_n==1
                    m1_even(i_n,j) = 2*(((j-1)*2)^2);
                elseif i_n==2 && j==1
                    m1_even(i_n,j) = 1;
                end
                if (i_n==j) && (i_n>2)
                    m1_even(i_n,j) = -(2*i + ((alpha*c_n_1)/(2*(i-1))) + ((alpha*c_n)/(2*(i+1))));
                elseif (i_n>2) && (j<i_n) && (j==i_n-1)
                    m1_even(i_n,j) = (alpha*c_n_2*c_n_1)/(2*(i-1));
                elseif (i_n>2) && (j>i_n) && (j==i_n+1)
                    m1_even(i_n,j) = alpha/(2*(i+1));
                end
            else
                if i_n==1
                    m1_even(i_n,j) = 2*(((j-1)*2)^2);
                end
                if (i_n==j) && (i_n>1)
                    m1_even(i_n,j) = -(2*i + ((alpha*c_n_1)/(2*(i-1))) + ((alpha*c_n)/(2*(i+1))));
                elseif (i_n>1) && (j<i_n) && (j==i_n-1)
                    m1_even(i_n,j) = (alpha*c_n_2*c_n_1)/(2*(i-1));
                elseif (i_n>1) && (j>i_n) && (j==i_n+1)
                    m1_even(i_n,j) = alpha/(2*(i+1));
                end
            end
        end
    end
    for z = 1:length(even_spac)
        m2_even_new(z) = m2_even(even_spac(z) + 1);
    end
end