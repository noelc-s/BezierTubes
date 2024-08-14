%% Connect boundary conditions of minimal curves
clf
dt = 1; % dtau
order = 7; % minimal curve
m = 1;

Xi = [0 -0.1 0.3 0.6 0.9 1.0 1.0 1.0; ...
      0 0.3 0.5 0.5 0.5 0.5 0.5 0.4];

figure(1)
scatter(Xi(1,:),Xi(2,:))
Z = Bezier.Z(order, dt);
tau = linspace(0,1);
C = Xi*Z(tau);
hold on;
% plot(C(1,:),C(2,:))
patch([C(1,:) NaN],[C(2,:) NaN],[tau 1],'facecolor','none','EdgeColor','interp','linewidth',5)
colorbar