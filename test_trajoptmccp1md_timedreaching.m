mccp1md = Mccpvd1MD();
x0 = [0;   0;   0;   0.1;   0;   0];
%xf = -0.3 ; 
xf = 0.7;
tf = 0.7;
testtraj = TrajOptMccp1MD.reachSmoothIlqr(mccp1md, x0, xf, 'tf', tf , ...
    'dt', 0.02, 'm2min', 0.1, 'sfactor', 3, 'ConstraintType', 2, 'Resting', 0.1);
%testtraj.SettlingTime
%stepinfo(testtraj.x(1,:), testtraj.t, xf(1), 'SettlingTimeShreshold', 0.02 )
plot_mccp1md(testtraj, 'PlotEnergy', 0)