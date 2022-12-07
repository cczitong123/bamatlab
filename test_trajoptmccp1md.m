mccp1md = Mccpvd1MD();
%%
% mccp1md.actuator.max_damping = 0.1;
% mccp1md.actuator.max_rege_damping = 0.05;
% mccp1md.actuator.max_rege_damping = 0.5;

%%
x0 = [0;   0;   0;   pi/24;   0;   0];
%xf = -0.3 ; 
xf = 0.7;
%xf = [0.7; 0; 0.7; 1; 0 ;0];
testtraj1 = TrajOptMccp1MD.reachFastIlqr(mccp1md, x0, xf, 'tf', 0.7, 'dt', 0.02, 'm2min', pi/24, 'sfactor', 1, 'ConstraintType', 2);
testtraj1.SettlingTime
%stepinfo(testtraj.x(1,:), testtraj.t, 0.7, 'SettlingTimeShreshold', 0.02 )
plot_mccp1md(testtraj1, 'PlotEnergy', 0, 'PlotImpedance',1)

%%
x0 = [0;   0;   0;   0.1;   0;   0];
%xf = -0.3 ; 
xf = 0.7;
testtraj2 = TrajOptMccp1MD.reachSmoothIlqr(mccp1md, x0, xf, 'tf', 0.5 , 'dt', 0.01, 'm2min', 0.1, 'sfactor', 1, 'ConstraintType', 2);
testtraj2.SettlingTime
%stepinfo(testtraj.x(1,:), testtraj.t, xf(1), 'SettlingTimeShreshold', 0.02 )
plot_mccp1md(testtraj2, 'PlotEnergy', 0)
% figure
% hold on
% plot(testtraj.t, testtraj.x(2,:))
% plot(testtraj.t, testtraj.x(5,:))
% hold off
%%
x0 = [0;   0;   0;  0.1;   0;   0];
%xf = [0.7; 0; 0.7; 0.1 ;   0;   0];
xf = 0.7;
testtraj3 = TrajOptMccp1MD.reachFastIlqr_fixu1(mccp1md, x0, xf, 'tf', 0.7, 'dt', 0.01, 'm2min', 0.5, 'sfactor', 1);
testtraj3.SettlingTime
%stepinfo(testtraj1.x(1,:), testtraj1.t, 'SettlingTimeShreshold', 0.02)
plot_mccp1md(testtraj3)

