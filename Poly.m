classdef Poly
    methods (Static)

        function V = hyp2vert(A,b)
            Vert = cddmex('extreme',struct('A',A,'B',b));
            V = Vert.V;
        end

        function [A,b] = vert2hyp(V)
            H = cddmex('hull',struct('V',V));
            A = H.A;
            b = H.B;
        end

        function V = conv(V)
            ind = convhull(V);
            V = V(ind,:);
        end

        function [H, D_nT] = getBezMatrices(order, dt)
            p = order;
            S = zeros(p,p+1);
            R = zeros(p+1,p);
            for i= 1:p
                R(i,i) = (p+1-i)/p;
                R(i+1,i) = i/p;
                S(i,i) = -p;
                S(i,i+1) = p;
            end
            H = (1/dt*S'*R');
            H_0 = H^0;
            H_1 = H^1;
            H_2 = H^2;
            z_0 = [1 zeros(1,order)]';
            z_T = [zeros(1,order) 1]';
            D = [H_0*z_0 H_1*z_0 H_0*z_T H_1*z_T];
            D_nT = inv(D');
        end

        function s = plotReachableTraj(Vert,col,x0,x1,dt,H)
            %%% Plot interpolations of the edges
            tau = linspace(0,1);
            set_x0 = false;
            set_x1 = false;
            if isempty(x0)
                set_x0 = true;
            end
            if isempty(x1)
                set_x1 = true;
            end
            s = [];
            for i = 1:size(Vert,1)-1
                for lambda = 0:0.5:1
                    X = (1-lambda)*Vert(i,1) + lambda*Vert(i+1,1);
                    Y = (1-lambda)*Vert(i,2) + lambda*Vert(i+1,2);
                    
                    if set_x0
                        x0 = [X Y];
                    end
                    if set_x1
                        x1 = [X Y];
                    end
                    
                    a = x0(1);
                    b = 1/3*(dt*x0(2)+3*x0(1));
                    c = 1/3*(-dt*x1(2)+3*x1(1));
                    d = x1(1);
                    pBez = @(t) a*(1-t).^3 + 3*b*t.*(1-t).^2 + 3*c*t.^2.*(1-t) + d*t.^3;
                    vBez = @(t) (3*d*t.^2 - 3*c*t.^2 - 3*a*(t - 1).^2 + 3*b*(t - 1).^2 - 6*c*t.*(t - 1) + 3*b*t.*(2*t - 2))*1/dt;
                    
                    s = [s; scatter([a b c d], [a b c d]*H,30,col,'filled')];
                    s = [s; plot(pBez(tau),vBez(tau),'linewidth',1,'color',col)];
                    %         p4=plot([a b c d],[a b c d]*H,'color','b','linestyle','--');
                end
            end
        end
        
        function [s, tau, pos, vel,cp]  = plotTraj(col,x0,x1,dt,H)
            %%% Plot interpolations of the edges
            tau = linspace(0,1);
            s = [];
                    a = x0(1);
                    b = 1/3*(dt*x0(2)+3*x0(1));
                    c = 1/3*(-dt*x1(2)+3*x1(1));
                    d = x1(1);
                    pBez = @(t) a*(1-t).^3 + 3*b*t.*(1-t).^2 + 3*c*t.^2.*(1-t) + d*t.^3;
                    vBez = @(t) (3*d*t.^2 - 3*c*t.^2 - 3*a*(t - 1).^2 + 3*b*(t - 1).^2 - 6*c*t.*(t - 1) + 3*b*t.*(2*t - 2))*1/dt;
                    
%                     s = [s; scatter([a b c d], [a b c d]*H,30,col,'filled')];
                    pos = pBez(tau);
                    vel = vBez(tau);
%                     s = [s; plot(pos, vel,'linewidth',1,'color',col)];
                cp = [[a b c d]; [a b c d]*H];
        end
        
        function [A_in, b_in] = backwardReachable(H, A, B,A_x, b_x,u_max,D_nT,x1)
            % x0 such that x0,x1,x2,x3 satisfy the input constraints
            H_0 = H^0;
            H_1 = H^1;
            H_2 = H^2;
            A_x_ = [];
            b_x_ = [];
            A_u = [];
            b_u = [];
            for m = 1:4
                I_m = zeros(1,4);
                I_m(m) = 1;
                A_tmp = [I_m*H_0'; I_m*H_1'];
                A_lin = ([1 0]*A*B)\([1 0]*A*A*[I_m*H_0'; I_m*H_1'] - I_m*H_2');
                A_in = A_lin(1:2)*D_nT(1:2,1:2);
                b_in = -A_lin(3:4)*D_nT(3:4,3:4)*x1';
                
                A_u = [A_u; A_in; -A_in];
                b_u = [b_u; [b_in+u_max; -b_in+u_max]];
                
                % and state constraints
                A_x_ = [A_x_; A_x*A_tmp(:,1:2)*D_nT(1:2,1:2)];
                b_x_ = [b_x_; b_x-A_x*A_tmp(:,3:4)*D_nT(3:4,3:4)*x1'];
                
            end
            A_in = [A_u; A_x_];
            b_in = [b_u; b_x_];
        end
        
        function [A_in, b_in] = dynamicTube(H, A, B,A_x, b_x,u_max,D_nT)
            % x0 such that x0,x1,x2,x3 satisfy the input constraints
            H_0 = H^0;
            H_1 = H^1;
            H_2 = H^2;
            A_x_ = [];
            b_x_ = [];
            A_u = [];
            b_u = [];
            for m = 1:4
                I_m = zeros(1,4);
                I_m(m) = 1;
                A_tmp = [I_m*H_0'; I_m*H_1'];
                A_lin = ([1 0]*A*B)\([1 0]*A*A*A_tmp - H_2(:,m)');
                A_in = A_lin*D_nT;
                
                A_u = [A_u; A_in; -A_in];
                b_u = [b_u; [u_max; u_max]];
                
                % and state constraints
                A_x_ = [A_x_; A_x*A_tmp*D_nT];
                b_x_ = [b_x_; b_x];
                
            end
            A_in = [A_u; A_x_];
            b_in = [b_u; b_x_];
        end
        
        function [A_in, b_in] = forwardReachable(H, A, B,A_x, b_x,u_max,D_nT,x0)
            % x0 such that x0,x1,x2,x3 satisfy the input constraints
            H_0 = H^0;
            H_1 = H^1;
            H_2 = H^2;
            A_x_ = [];
            b_x_ = [];
            A_u = [];
            b_u = [];
            for m = 1:4
                I_m = zeros(1,4);
                I_m(m) = 1;
                A_tmp = [I_m*H_0'; I_m*H_1'];
                A_lin = ([1 0]*A*B)\([1 0]*A*A*A_tmp - H_2(:,m)');
                A_in = A_lin(3:4)*D_nT(3:4,3:4);
                b_in = -A_lin(1:2)*D_nT(1:2,1:2)*x0';
                
                A_u = [A_u; A_in; -A_in];
                b_u = [b_u; [b_in+u_max; -b_in+u_max]];
                
                % and state constraints
                A_x_ = [A_x_; A_x*A_tmp(:,3:4)*D_nT(3:4,3:4)];
                b_x_ = [b_x_; b_x-A_x*A_tmp(:,1:2)*D_nT(1:2,1:2)*x0'];
                
            end
            A_in = [A_u; A_x_];
            b_in = [b_u; b_x_];
        end
        
        function [A_in, b_in] = backwardReachable2D(H, A, B,A_x, b_x,u_max,D_nT,x1)
            H_0 = H^0;
            H_1 = H^1;
            H_2 = H^2;
            A_u = [];
            b_u = [];
            A_x_ = [];
            b_x_ = [];
            for m = 1:4    
                I_m = zeros(2,8);
                I_m(1,(m-1)*2+1) = 1;
                I_m(2,(m-1)*2+2) = 1;
                A_tmp = [I_m*kron(H_0,eye(2))'; I_m*kron(H_1,eye(2))'];
                A_lin = ([1 0 0 0; 0 1 0 0]*A*B)\([1 0 0 0; 0 1 0 0]*A*A*A_tmp - I_m*kron(H_2,eye(2))');
                A_in = A_lin*kron(D_nT,eye(2));

                b_in = -A_in(:,5:8)*x1';
                A_in = A_in(:,1:4);
                
                A_u = [A_u; A_in; -A_in];
                b_u = [b_u; [b_in+u_max; -b_in+u_max]];
                
                A_state = A_x*A_tmp*kron(D_nT,eye(2));
                
                A_x_ = [A_x_; A_state(:,1:4)];
                b_x_ = [b_x_; b_x-A_state(:,5:8)*x1'];
            end
            A_in = [A_u; A_x_];
            b_in = [b_u; b_x_];
        end
        
        function [A_in, b_in] = forwardReachable2D(H, A, B,A_x, b_x,u_max,D_nT,x0)
            H_0 = H^0;
            H_1 = H^1;
            H_2 = H^2;
            A_u = [];
            b_u = [];
            A_x_ = [];
            b_x_ = [];
            for m = 1:4
                I_m = zeros(2,8);
                I_m(1,(m-1)*2+1) = 1;
                I_m(2,(m-1)*2+2) = 1;
                A_tmp = [I_m*kron(H_0,eye(2))'; I_m*kron(H_1,eye(2))'];
                A_lin = ([1 0 0 0; 0 1 0 0]*A*B)\([1 0 0 0; 0 1 0 0]*A*A*A_tmp - I_m*kron(H_2,eye(2))');
                A_in = A_lin*kron(D_nT,eye(2));

                b_in = -A_in(:,1:4)*x0';
                A_in = A_in(:,5:8);
                
                A_u = [A_u; A_in; -A_in];
                b_u = [b_u; [b_in+u_max; -b_in+u_max]];
                
                A_state = A_x*A_tmp*kron(D_nT,eye(2));
                
                A_x_ = [A_x_; A_state(:,5:8)];
                b_x_ = [b_x_; b_x-A_state(:,1:4)*x0'];
            end
            A_in = [A_u; A_x_];
            b_in = [b_u; b_x_];
        end
        
        function [A_in, b_in] = dynamicTube2D(H, A, B,A_x, b_x,u_max,D_nT)
            H_0 = H^0;
            H_1 = H^1;
            H_2 = H^2;
            A_u = [];
            b_u = [];
            A_x_ = [];
            b_x_ = [];
            for m = 1:4
                I_m = zeros(2,8);
                I_m(1,(m-1)*2+1) = 1;
                I_m(2,(m-1)*2+2) = 1;
                A_tmp = [I_m*kron(H_0,eye(2))'; I_m*kron(H_1,eye(2))'];
                A_lin = ([1 0 0 0; 0 1 0 0]*A*B)\([1 0 0 0; 0 1 0 0]*A*A*A_tmp - I_m*kron(H_2,eye(2))');
                A_in = A_lin*kron(D_nT,eye(2));
                
                A_u = [A_u; A_in; -A_in];
                b_u = [b_u; [u_max; u_max; u_max; u_max]]; % two control inputs
                
                A_x_ = [A_x_; A_x*A_tmp*kron(D_nT,eye(2))];
                b_x_ = [b_x_; b_x];
            end
            A_in = [A_u; A_x_];
            b_in = [b_u; b_x_];
        end
        
        function x_d = connect(x0,x1,dt)
            a = x0(1);
            b = 1/3*(dt*x0(2)+3*x0(1));
            c = 1/3*(-dt*x1(2)+3*x1(1));
            d = x1(1);
            pBez = @(t) a*(1-t).^3 + 3*b*t.*(1-t).^2 + 3*c*t.^2.*(1-t) + d*t.^3;
            vBez = @(t) (3*d*t.^2 - 3*c*t.^2 - 3*a*(t - 1).^2 + 3*b*(t - 1).^2 - 6*c*t.*(t - 1) + 3*b*t.*(2*t - 2))*1/dt;
            x_d = @(t) [pBez(t); vBez(t)];
        end
        
        function touching = adjacencyCheck(A1,A2,b1,b2)
            opts = optimset('Display','off');
            touching = false;
            [x,fval,ef] = linprog(zeros(1,size(A1,2)),[A1; A2],[b1; b2],[],[],[],[],[],opts);
            if ef==1
                touching = true;
            end
        end
    end
end
