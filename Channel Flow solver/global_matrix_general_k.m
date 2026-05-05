function [m] = global_matrix_general_k(ix,Nx,k_x,q_n_u,q_n_v,beta,gamma)
    m1 = zeros(3*Nx,3*Nx);
    m2 = zeros(3*Nx,1);

    kx_current = k_x(ix+1);
    alpha_eqn = kx_current^2 + beta;
    q_n_mu = q_n_u(:,ix+1);
    q_n_mv = q_n_v(:,ix+1);    
    c_n = 1;
    c_n_1 = 1;
    for i = 1 : 3*Nx
        % u-RHS set
        if (i>2) && i<(Nx+1)
            c_n_2 = 1;
            if i-3 == 0
                c_n_2 = 2;
            end
            RHS1_u = -((c_n_1*2)/(2*gamma*(i-2)))*(c_n_2*q_n_mu(i-2) - q_n_mu(i));
            if i==Nx || i==Nx-1
                RHS2_u = (2/(2*gamma*(i)))*(c_n*q_n_mu(i));
            else
                RHS2_u = (2/(2*gamma*(i)))*(c_n*q_n_mu(i) - q_n_mu(i+2));
            end
            m2(i) = RHS1_u + RHS2_u;
        end
        % v-RHS set
        if (i>Nx+2) && (i<(2*Nx)+1)
            c_n_2 = 1;
            if i-(Nx+1)-2 == 0
                c_n_2 = 2;
            end
            RHS1_v = (c_n_1*q_n_mv(i-Nx) - c_n_1*c_n_2*q_n_mv(i-Nx-2))/((i-Nx-2)*gamma);
            if i==2*Nx || i==2*Nx-1
                RHS2_v = -(- c_n*q_n_mv(i-Nx))/((i-Nx)*gamma);
            else
                RHS2_v = -((q_n_mv(i-Nx+2) - c_n*q_n_mv(i-Nx)))/((i-Nx)*gamma);
            end
            m2(i) = RHS1_v + RHS2_v;
        end
        for j = 1 : 3*Nx
            %% u-velcoity block
            % boundary rows
            if (i==1) && j<(Nx+1)
                m1(i,j) = (1)^(j-1);
            elseif (i==2) && j<(Nx+1)
                m1(i,j) = (-1)^(j-1);
            end
            % interior rows
            if (i>2) && (i<Nx+1) && j<(Nx+1)
                c_n_2 = 1;
                if i-3 == 0
                    c_n_2 = 2;
                end
                if  i==j 
                    m1(i,j) = (c_n_1*alpha_eqn/(2*(i-2))) + (alpha_eqn*c_n/(2*(i))) + 2*(i-1);
                end
                if j<i && j==i-2
                    m1(i,j) = -(alpha_eqn * c_n_1 * c_n_2/(2*(i-2)));
                elseif j>i && j==i+2
                    m1(i,j) = -alpha_eqn/(2*(i));
                end              
            end
            %% u-P block
            % interior rows
            if (i>2) && (i<Nx+1) && j>(2*Nx)
                c_n_2 = 1;
                if i-3 == 0
                    c_n_2 = 2;
                end
                if i==(j-2*Nx)
                    m1(i,j) = ((c_n_1*1i*kx_current)/(2*gamma*(i-2))) + ((c_n*1i*kx_current)/(2*gamma*(i)));
                end
                if j<(i+2*Nx) && j==(i+2*Nx)-2
                    m1(i,j) = -(c_n_1*1i*kx_current*c_n_2)/(2*gamma*(i-2));
                elseif j>(i+2*Nx) && j==(i+2*Nx)+2
                    m1(i,j) = -(1i*kx_current)/(2*gamma*(i));
                end
            end
            %% v-velocity block
            % boundary rows
            if (i==Nx+1) && (j>Nx) && (j<(2*Nx)+1)
                m1(i,j) = (1)^(j-Nx-1);
            elseif (i==Nx+2) && (j>Nx) && (j<(2*Nx)+1)
                m1(i,j) = (-1)^(j-Nx-1);
            end
            % interior rows
            if (i>Nx+2) && (i< 2*Nx+1) && (j>Nx) && j<(2*Nx +1)
                c_n_2 = 1;
                if i-(Nx+1)-2 == 0
                    c_n_2 = 2;
                end
                if  i==j 
                    m1(i,j) = (c_n_1*alpha_eqn/(2*(i-Nx-2))) + (alpha_eqn*c_n/(2*(i-Nx))) + 2*(i-Nx-1);
                end
                if j<i && j==i-2
                    m1(i,j) = -(alpha_eqn * c_n_1 * c_n_2/(2*(i-Nx-2)));
                elseif j>i && j==i+2
                    m1(i,j) = -alpha_eqn/(2*(i-Nx));
                end
            end
            %% v-P block
            % interior rows
            if (i>Nx+2) && (i<2*Nx+1) && j>2*Nx
                if j<(i+Nx) && j==i+Nx-1
                    m1(i,j) = -(c_n_1/gamma);
                elseif j>(i+Nx) && j==i+Nx+1
                    m1(i,j) = (1/gamma);
                end
            end
            %% continuity equation block
            % u COM
            if (i>2*Nx) && (j<Nx+1)
                if (j==i-(2*Nx))
                    m1(i,j) = 1i * kx_current;
                elseif (j==i-(2*Nx)+2)
                    m1(i,j) = -1i * kx_current;
                end
            end
            % v COM
            if (i>2*Nx) && (i<3*Nx -1) && (j>Nx) && (j<2*Nx + 1)
                if (j==i-Nx+1)
                    m1(i,j) = 2*(j-Nx-1);
                end
            elseif (i==3*Nx-1) && (j>Nx) && (j<2*Nx + 1)
                if (j==i-Nx+1)
                    m1(i,j) = (j-Nx-1);
                end
            end
        end
    end
    m = m1\m2;
end