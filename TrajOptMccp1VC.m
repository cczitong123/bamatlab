classdef TrajOptMccp1VC
    % For Maccepavd velocity controlled
    
    properties
    end
    
    methods
    end
    
    methods (Static)
        
        function [ traj ] = reachFastIlqr( robot, x0, xf, varargin )
            % x0, xf, tf, dt, m2min
            parser = inputParser();
            
            addRequired(parser, 'robot');
            addRequired(parser, 'x0');
            addRequired(parser, 'xf');
            addOptional(parser, 'tf', 1);
            addOptional(parser, 'dt', 0.02);
            addOptional(parser, 'm2min', 0);
            addOptional(parser, 'sfactor', 3);
            addOptional(parser, 'ConstraintType', 1);
            parse(parser, robot, x0, xf, varargin{:});
            % x0: scalar for q0, vector to include velocity and motor
            % positions
            % xf: scalar for qf, vector to include final position of 
            tic
            
            %T = tf;
            %if nargin < 5
            dt = parser.Results.dt;
            tf = parser.Results.tf;
            sfactor = parser.Results.sfactor;
            %end
            %if nargin < 6
            m2min = parser.Results.m2min;
            %end
            if isscalar(x0)
                x0 = [x0;0;x0;m2min;0;0];
            end
            qf = xf(1);
            if isscalar(xf)
                Hf = diag([1, 0, 0, 0, 0, 0]);
                xf = [qf;0;qf;m2min;0;0];
            else
                if parser.Results.ConstraintType == 1
                    Hf = diag([1, 1, 1, 1, 1, 1]);
                elseif parser.Results.ConstraintType == 2
                    Hf = diag([1, 1, 1, 1, 0, 0]);
                else
                    
                end
                
            end
            
            t = 0:dt:tf;
            Nt = tf/dt + 1;
            
            f = @(x,u)robot.dynamics_with_jacobian_fd( x, u);
            costpara = [];
            costpara.qf = qf;
            costpara.xf = xf;
            costpara.Hf = Hf;
            costpara.robot = robot;
            costpara.tf = tf;
            costpara.w_tf = 1e3*dt;
            costpara.w_t  = 1e3; %
            costpara.w_e  = sfactor*1e2; % 1 is good for spring force
            j = @(x,u,t)TrajOptMccp1MD.jReachFast(x, u, t, costpara);
            %%%
            opt_param = [];
            
            opt_param.lambda_init = 1;
            opt_param.lambda_max  = 1e10;
            opt_param.iter_max = 150;
            opt_param.online_plotting = 0;
            opt_param.online_printing = 2;
            opt_param.dcost_converge = 1e-3;
            opt_param.solver = 'rk4';
            opt_param.target = qf;
            
            opt_param.T = tf;
            %%%
            opt_param.umax = [pi/3; pi/2; 1];
            opt_param.umin = [-pi/3; m2min; 0];
            % u0 can be full command sequence or just initial point
            u0 = [qf; m2min; 0.1];
            %u0 = [0; 0.1; 0];
            
            %traj = ILQRController.ilqr(f, j, dt, Nt, x0, u0, opt_param);
            traj = iLQRv1(f, j, dt, Nt, x0, u0, opt_param);
            
            traj.t = t;
            %traj = val_traj_mccpvd(robot, task, traj);
            if nargin == 0
                plot_traj_mccpvd1(traj);
            end
            
            toc
            disp('ILQR done')
            %%% Compute the costs
            traj.J1 = sum((traj.x(1,:) - qf).^2, 2)*dt*costpara.w_t;
            traj.Jtf = costQuad(traj.x(:,end) - xf, Hf*costpara.w_tf);
            traj.Jp = traj.J1 + traj.Jtf;
            traj.Je = sum((traj.u(1,1:end) - traj.x(1,1:end-1)).^2 + traj.u(2,:).^2,2)*dt*costpara.w_e ;
            traj = evaluate_traj_Mccp1MD(robot, traj, 'qf', qf);
        end
        
        function [ traj ] = reachFastIlqr_fixu1( robot, x0, xf, varargin )
            % x0, xf, tf, dt, m2min
            parser = inputParser();
            
            addRequired(parser, 'robot');
            addRequired(parser, 'x0');
            addRequired(parser, 'xf');
            addOptional(parser, 'tf', 1);
            addOptional(parser, 'dt', 0.02);
            addOptional(parser, 'm2min', 0);
            addOptional(parser, 'sfactor', 1);
            parse(parser, robot, x0, xf, varargin{:});
            % x0: scalar for q0, vector to include velocity and motor
            % positions
            % xf: scalar for qf, vector to include final position of 
            tic
            
            %T = tf;
            %if nargin < 5
            dt = parser.Results.dt;
            tf = parser.Results.tf;
            sfactor = parser.Results.sfactor;
            %end
            %if nargin < 6
            m2min = parser.Results.m2min;
            %end
            if isscalar(x0)
                x0 = [x0;0;x0;m2min;0;0];
            end
            qf = xf(1);
            if isscalar(xf)
                Hf = diag([1, 0, 0, 0, 0, 0]);
                xf = [qf;0;qf;m2min;0;0];
            else
                Hf = diag([1, 1, 1, 1, 1, 1]);
                
            end
            
            t = 0:dt:tf;
            Nt = tf/dt + 1;
            
            f = @(x,u)robot.dynamics_with_jacobian_fd( x, u);
            costpara = [];
            costpara.qf = qf;
            costpara.xf = xf;
            costpara.Hf = Hf;
            costpara.robot = robot;
            costpara.tf = tf;
            costpara.w_tf = 1e2;
            costpara.w_t  = 1e3; %
            costpara.w_e  = sfactor*1e2; % 1 is good for spring force
            j = @(x,u,t)TrajOptMccp1MD.jReachFast(x, u, t, costpara);
            %%%
            opt_param = [];
            
            opt_param.lambda_init = 1;
            opt_param.lambda_max  = 1e10;
            opt_param.iter_max = 150;
            opt_param.online_plotting = 0;
            opt_param.online_printing = 2;
            opt_param.dcost_converge = 1e-3;
            opt_param.solver = 'rk4';
            opt_param.target = qf;
            
            opt_param.T = tf;
            %%%
            opt_param.umax = [qf; pi/2; 1];
            opt_param.umin = [qf; m2min; 0];
            % u0 can be full command sequence or just initial point
            u0 = [qf; m2min; 0.1];
            %u0 = [0; 0.1; 0];
            
            traj = ILQRController.ilqr(f, j, dt, Nt, x0, u0, opt_param);
            traj.t = t;
            %traj = val_traj_mccpvd(robot, task, traj);
            if nargin == 0
                plot_traj_mccpvd1(traj);
            end
            
            toc
            disp('ILQR done')
            
            traj = evaluate_traj_Mccp1MD(robot, traj, 'qf');
        end
        
        function [l, l_x, l_xx, l_u, l_uu, l_ux] = jReachFast( x, u, t, para)
            w_tf = para.w_tf;
            w_t  = para.w_t; %
            w_e  = para.w_e; % 1 is good for spring force
            %Rtf  = diag([w_tf, 0, 0, 0, 0, 0]);
            if (isnan(u))
                % final cost
                if nargout > 1
                    [l, l_x, l_xx] = costQuad(x-para.xf, para.Hf*w_tf);
                else
                    l = costQuad(x-para.xf, para.Hf*w_tf);
                end
            else
                % running cost
                %para = [];
                para.w_t = w_t;
                para.w_e = w_e;
                fl = @(x,u,t) costr_reachFast(x,u,t,para);
                l = fl(x,u,t);
                
                if nargout > 1
                    % finite difference
                    flJ=@(x,u,t)J_cost_fd ( fl, x, u, t );
                    [l_x ,l_u      ] = flJ ( x, u, t );
                    flH =@(x,u,t)H_cost_fd  ( flJ, x, u, t );
                    [l_xx,l_uu,l_ux] = flH  ( x, u, t );
                end
            end
        end
        
        function [ traj ] = reachSmoothIlqr( robot, x0, xf, varargin )
        % x0, xf, tf, dt, m2min
            parser = inputParser();
            
            addRequired(parser, 'robot');
            addRequired(parser, 'x0');
            addRequired(parser, 'xf');
            addOptional(parser, 'tf', 0.6);
            addOptional(parser, 'dt', 0.02);
            addOptional(parser, 'm2min', 0);
            addOptional(parser, 'sfactor', 1);
            addOptional(parser, 'ConstraintType', 2);
            addOptional(parser, 'Resting', 0.3);
            parse(parser, robot, x0, xf, varargin{:});
            % x0: scalar for q0, vector to include velocity and motor
            % positions
            % xf: scalar for qf, vector to include final position of 
            tic
            
            %T = tf;
            %if nargin < 5
            dt = parser.Results.dt;
            tf = parser.Results.tf;
            T = parser.Results.Resting + tf;
            
            sfactor = parser.Results.sfactor;
            %end
            %if nargin < 6
            m2min = parser.Results.m2min;
            %end
            if isscalar(x0)
                x0 = [x0;0;x0;m2min;0;0];
            end
            qf = xf(1);
            if isscalar(xf)
                Hf = diag([1, 1, 0, 0, 0, 0]);
                xf = [qf;0;qf;m2min;0;0];
            else
                if parser.Results.ConstraintType == 1
                    Hf = diag([1, 1, 1, 1, 1, 1]);
                elseif parser.Results.ConstraintType == 2
                    Hf = diag([1, 1, 1, 1, 0, 0]);
                else
                    
                end
                
            end
            
            t = 0:dt:T;
            Nt = length(t);
            
            f = @(x,u)robot.dynamics_with_jacobian_fd( x, u);
            costpara = [];
            costpara.qf = qf;
            costpara.dt = dt;
            costpara.xf = xf;
            costpara.Hf = Hf;
            costpara.robot = robot;
            costpara.tf = tf;
            costpara.w_tf = 1e3*dt;
            costpara.w_t  = 1e3; %
            costpara.w_e  = sfactor*1e3; % 1 is good for spring force
            j = @(x,u,t)TrajOptMccp1MD.jReachSmooth( x, u, t, costpara);
            %%%
            opt_param = [];
            
            opt_param.lambda_init = 1;
            opt_param.lambda_max  = 1e10;
            opt_param.iter_max = 150;
            opt_param.online_plotting = 0;
            opt_param.online_printing = 2;
            opt_param.dcost_converge = 1e-4;
            opt_param.solver = 'rk4';
            opt_param.target = qf;
            
            opt_param.T = T;
            %%%
            opt_param.umax = [pi/3; pi/2; 1];
            opt_param.umin = [-pi/3; m2min; 0];
            % u0 can be full command sequence or just initial point
            u0 = [qf; m2min; 0.1];
            %u0 = [0; 0.1; 0];
            
            traj = ILQRController.ilqr(f, j, dt, Nt, x0, u0, opt_param);
            traj.t = t;
            %traj = val_traj_mccpvd(robot, task, traj);
            if nargin == 0
                plot_traj_mccpvd1(traj);
            end
            
            toc
            disp('ILQR done')
            
            traj.J1 = sum( abs(gradient(traj.x(2,:), dt )) , 2)*dt*costpara.w_t;
            traj.Jtf = costQuad(traj.x(:,end) - xf, Hf*costpara.w_tf);
            
            for i=1:(length(t)-1)
                [~, cost1(i), coste(i)] = costr_reachSmooth(traj.x(:,i), traj.u(:,i), t(i), costpara);
                cost1(i) = cost1(i)*dt;
                coste(i) = coste(i)*dt;
            end
            traj.Jtf = costQuad(traj.x(:,end)-costpara.xf, costpara.Hf*costpara.w_tf);
            traj.J1 = sum(cost1);
            traj.Jp = traj.J1 + traj.Jtf;
            traj.Je = sum(coste) ;
            traj = evaluate_traj_Mccp1MD(robot, traj, 'qf', qf);
        end
        
        function [l, l_x, l_xx, l_u, l_uu, l_ux] = jReachSmooth( x, u, t, para)
            w_tf = para.w_tf;
            w_t  = para.w_t; %
            w_e  = para.w_e; % 1 is good for spring force
            %Rtf  = diag([w_tf, 0, 0, 0, 0, 0]);
            if (isnan(u))
                % final cost
                if nargout > 1
                    [l, l_x, l_xx] = costQuad(x-para.xf, para.Hf*w_tf);
                else
                    l = costQuad(x-para.xf, para.Hf*w_tf);
                end
            else
                % running cost
                %para = [];
                para.w_t = w_t;
                para.w_e = w_e;
                fl = @(x,u,t) costr_reachSmooth(x,u,t,para);
                l = fl(x,u,t);
                
                if nargout > 1
                    % finite difference
                    flJ=@(x,u,t)J_cost_fd ( fl, x, u, t );
                    [l_x ,l_u      ] = flJ ( x, u, t );
                    flH =@(x,u,t)H_cost_fd  ( flJ, x, u, t );
                    [l_xx,l_uu,l_ux] = flH  ( x, u, t );
                end
            end
        end
    end
    
end

function c = costr_reachFast(x, u, t, para)
            %c1 = (x(1) - para.qf).^2*t/para.tf;
           
            c1 = (x(1) - para.qf).^2;
            %ce = (para.robot.spring_force(x,u)).^2;
            ce = (u(1) - x(1))^2 + u(2)^2 + 1*(u(3)-0.5)^2;
            %ce = (u(1) - x(3))^2 + (u(2) - x(4))^2 + 1e-6*u(3)^2;
            c = c1*para.w_t + ce*para.w_e ;
end

function [c, c1_, ce_ ] = costr_reachSmooth(x, u, t, para)
tf = para.tf;
dt = para.dt;
tf = round(tf/dt)*dt;

s = t/tf;
B = min( s, 1);

c1 = (x(1) - para.qf).^2*B;
%ce = (para.robot.spring_force(x,u)).^2;
ce = (u(1) - x(1))^2 + u(2)^2 + 1e-6*u(3)^2;
%ce = (u(1) - x(3))^2 + (u(2) - x(4))^2 + 1e-6*u(3)^2;
c1_ = c1*(para.w_t*(s<1) + 1000*para.w_t*(s==1)+100*para.w_t*(s>1)  );
ce_ = ce*para.w_e ;
c =  c1_ + ce_;

end

function [c, c_x, c_xx] = costQuad(x, R)
            c = x'*R*x/2;
            
            if nargout > 1
                V = diag(R);
                c_x = V.*x;
                c_xx = R;
            end
end



function x = simulate_feedforward ( x0, f, u, p )

N = size(u,2)+1; % N of time index
x = nan(size(x0,1),N); x(:,1)=x0; % initialise x
for n=1:N-1
	x(:,n+1) = simulate_step ( f, x(:,n), u(:,n), p );
end

end

function xn = simulate_step ( f, x, u, p )

switch p.solver
	case 'euler'
		dt=p.dt;
		xn = x + dt*f(x,u); % euler step
	case 'rk4'
		dt = p.dt;
        
		g1 = dt*f(x            ,u);
		g2 = dt*f(x+.5*g1,u);
		g3 = dt*f(x+.5*g2,u);
		g4 = dt*f(x+   g3,u);
		xn = x + (1/6)*(g1 + 2*g2 + 2*g3 + g4);
    case 'ode45'
        tspan = [0 p.dt];
        odefun = @(t,y)f(y,u);
        [~,y] = ode45(odefun,tspan,x);
        xn = y(end,:)';
	otherwise
		error('Unknown solver.')
end
end

function [l_xx,l_uu,l_ux] = H_cost_fd ( J, x, u, t )

dimX = size(x,1);
dimU = size(u,1);

delta=1e-6;

l_xx = zeros(dimX,dimX);
l_ux = zeros(dimU,dimX);
for i=1:dimX
	dx=zeros(dimX,1); dx(i)=delta;

	[lxxp,luxp] = J( x+dx, u, t );
	[lxxm,luxm] = J( x-dx, u, t );

	l_xx(:,i) = (lxxp-lxxm)/(2*delta);
	l_ux(:,i) = (luxp-luxm)/(2*delta);
end

l_uu = zeros(dimU,dimU);
for i=1:dimU
	du=zeros(dimU,1); du(i)=delta;

	[dummy,luup] = J( x, u+du, t );
	[dummy,luum] = J( x, u-du, t );

	l_uu(:,i) = (luup-luum)/(2*delta);
end

end

function [l_x,l_u] = J_cost_fd ( l, x, u, t )

dimX = size(x,1);
dimU = size(u,1);

delta=1e-6;

l_x  = zeros(dimX,1);
for i=1:dimX
	dx=zeros(dimX,1); dx(i)=delta;

	lxp = l( x+dx, u, t );
	lxm = l( x-dx, u, t );

	l_x(i) = (lxp-lxm)/(2*delta);
end

l_u  = zeros(dimU,1);
for i=1:dimU
	du=zeros(dimU,1); du(i)=delta;

	lup = l( x, u+du, t );
	lum = l( x, u-du, t );

	l_u(i) = (lup-lum)/(2*delta);
end

end