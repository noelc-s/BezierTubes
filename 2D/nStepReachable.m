init()

x0_initial = [0.5 0.3];
x0 = x0_initial;

u_max = .1;
dt = 1;
steps =15;
% dt = 15;
% steps =1;

DI = DoubleIntegrator();

A_x = [0 1; 0 -1; 1 0; -1 0];
b_x = [1; 1; 1; 1];

order = 2*size(DI.A,1)-1;
[H, D_nT] = Poly.getBezMatrices(order, dt);

[A_dyn, b_dyn] = Poly.forwardReachable(H, DI.A, DI.B,A_x, b_x,u_max,D_nT,x0);
F = Poly.conv(Poly.hyp2vert(A_dyn, b_dyn));
[A_dyn, b_dyn] = Poly.backwardReachable(H, DI.A, DI.B,A_x, b_x,u_max,D_nT,x0);
B = Poly.conv(Poly.hyp2vert(A_dyn, b_dyn));
V_x = Poly.conv(Poly.hyp2vert(A_x, b_x));

% Poly.plotReachableTraj(Vert,'r',x0,[],dt,H)
clf;
hold on;
axis equal;
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)

patch(F(:,1),F(:,2),'g','facealpha',0.1);
patch(V_x(:,1),V_x(:,2),'k','facealpha',0.05,'linewidth',2);

%% Iterate
continue_F = true;
continue_B = true;
Vert_k_F = F;
Vert_k_B = B;
for k = 1:steps-1
    Vert_k_F_ = [];
    Vert_k_B_ = [];

    % Forward
    for i = 1:size(Vert_k_F,1)
        x0 = Vert_k_F(i,:);
        [A_dyn, b_dyn] = Poly.forwardReachable(H, DI.A, DI.B,A_x, b_x,u_max,D_nT,x0);
        Vert_2 = Poly.hyp2vert(A_dyn, b_dyn);
        if size(Vert_2,1)>2
            Vert_2 = Poly.conv(Vert_2);
            Vert_k_F_ = [Vert_k_F_; Vert_2];
        end
    end
%     Vert_k_F = Poly.conv([Vert_k_F; Vert_k_F_]);
    if isempty(Vert_k_F_)
        continue_F = false;
    end
    if continue_F
        Vert_k_F = Poly.conv(Vert_k_F_);
        patch(Vert_k_F(:,1),Vert_k_F(:,2),'g','facealpha',0.1);
    end

    % Backward
    for i = 1:size(Vert_k_B,1)
        x1 = Vert_k_B(i,:);
        [A_dyn, b_dyn] = Poly.backwardReachable(H, DI.A, DI.B,A_x, b_x,u_max,D_nT,x1);
        Vert_2 = Poly.hyp2vert(A_dyn, b_dyn);
        if size(Vert_2,1)>2
            Vert_2 = Poly.conv(Vert_2);
            Vert_k_B_ = [Vert_k_B_; Vert_2];
        end
    end
%     Vert_k_B = Poly.conv([Vert_k_B; Vert_k_B_]);
    if isempty(Vert_k_B_)
        continue_B = false;
    end
    if continue_B
        Vert_k_B = Poly.conv(Vert_k_B_);
        patch(Vert_k_B(:,1),Vert_k_B(:,2),'b','facealpha',0.1);
    end
    drawnow;
end
