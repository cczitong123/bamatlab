mccp2vc = Mccpvd2VC();
q0 = [0;0.2];
x0 = [ q0(1); q0(2); zeros(6,1) ] ;
%x0(6,:) = 0;
%x0(8,:) = 0;

qf = [0.7; 0.6];

y0 = mccp2vc.endpoint( q0 );
yf = mccp2vc.endpoint( qf );

tf = round(0.15*log2(1 + norm(yf-y0)/0.01), 2);

traj = mccp2vc.reach_ilqr(x0,qf,'tf',tf,'dt',0.01,'m2f',[0.1;0.1],'plot',1,'animate',1,'saveanimate',0);
%%
%traj = mccp2vc.reach_ilqr(x0,yf,'os',1,'tf',tf,'m2f',[0;0],'plot',1,'animate',1,'saveanimate',0);

qs = [ 0.2 -0.3 -1  0.5 0.2;
       0.2  1.5  1.2  0.3  0.2];
targets = mccp2vc.endpoint(qs);        
figure   
fig=gcf;
fig.Units='inches';
fig.Position=[2 2 3.5 3.5];
scatter(targets(1,:),targets(2,:))
text(targets(1,:),targets(2,:),{'1','2','3','4','End'});
xlim([-1 1])
ylim([-1 1])