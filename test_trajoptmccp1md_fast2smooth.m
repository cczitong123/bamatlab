mccp1md = Mccpvd1MD();

x0 = [0;   0;   0;   0.1;   0;   0];

xf = 0.7;

testtraj1 = TrajOptMccp1MD.reachFastIlqr(mccp1md, x0, xf, 'tf', 1, 'dt', 0.02, 'm2min', 0.1, 'sfactor', 1, 'ConstraintType', 2);

testtraj1.SettlingTime
%stepinfo(testtraj.x(1,:), testtraj.t, 0.7, 'SettlingTimeShreshold', 0.02 )
plot_mccp1md(testtraj1, 'PlotEnergy', 0, 'PlotImpedance',1)

%%
tf = 0.5;
length_u = (tf+0.1)/0.02; %the total time is tf+resting
u0 = testtraj1.u(:,1:length_u);
testtraj2 = TrajOptMccp1MD.reachSmoothIlqr(mccp1md, x0, xf, 'tf', tf , ...
    'dt', 0.02, 'm2min', 0.1, 'sfactor', 1, 'ConstraintType', 2, 'Resting', 0.1, 'u0', u0);
%testtraj.SettlingTime
%stepinfo(testtraj.x(1,:), testtraj.t, xf(1), 'SettlingTimeShreshold', 0.02 )
plot_mccp1md(testtraj2, 'PlotEnergy', 0)