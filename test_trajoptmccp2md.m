mccp2 = Mccpvd2DofMD();

x0 = zeros(12,1);
%x0(1)=x(5)=0.1;
%x0(2)=x(0.05;
%yf = [0.15; 0.36];

qf = [ pi/4;pi/6 ];
%xf = [ qf(1); qf(2); 0; 0; qf(1); 0.3; 0; 0; qf(2); 0.3; 0; 0 ];
%testtraj1 = TrajOptMccp2MD.reachFastIlqrOs(mccp2, x0, yf, 'tf', 1, 'dt', 0.02, 'm2min', [0,0], 'sfactor', 1);
%testtraj1.SettlingTime

testtraj1 = TrajOptMccp2MD.reachFastIlqr(mccp2, x0, qf, 'tf', 1, 'dt', 0.01, 'm2min', [0,0], 'sfactor', 1);

%%
%unat = [qf(1);0.3;0;qf(2);0.3;0 ];
%unat = repmat(unat,1,51);

%xsim = mccp2.simulate_feedforward(x0,unat,0.02);
xsim = mccp2.simulate_feedforward(x0,testtraj1.u,0.02);
%mccp2.plot(x0);

frames = mccp2.animate(xsim,0.02);
hold on
%scatter(yf(1),yf(2));
y=Arm2Dof.endpoint(qf,mccp2.L);
scatter(y(1),y(2))
%%
figure
subplot(2,1,1)
plot(testtraj1.t, testtraj1.x(5,:))
hold on
plot(testtraj1.t, testtraj1.x(1,:))
plot(testtraj1.t(1:end-1), testtraj1.u(1,:))
subplot(2,1,2)
plot(testtraj1.t, testtraj1.x(9,:))
hold on
plot(testtraj1.t, testtraj1.x(2,:))
plot(testtraj1.t(1:end-1), testtraj1.u(4,:))