classdef TrajOptMccp2MD 
    %TRAJOPTMCCP2MD Summary of this class goes here
    %  x: 
    
    properties
        
    end
    
    methods
    end
    
    methods (Static)
        % reach fast to target in operation space (yf) solved by ILQR
        function [ traj ] = throw(robot, x0, varargin)
            parser = inputParser();
            addRequired(parser, 'robot');
            addRequired(parser, 'x0');
            addOptional(parser, 'tf', 0.8);
            addOptional(parser, 'dt', 0.02);
            addOptional(parser, 'm2min', [0,0]);
            
            
        end
        
        
        
        function [traj]   = reachFastIlqrOs(robot, x0, yf, varargin)
            % reach fast to target in operation space solved by ILQR
            % x0, yf, tf, dt, m2min
            parser = inputParser();
            
            addRequired(parser, 'robot');
            addRequired(parser, 'x0');
            addRequired(parser, 'yf');
            addOptional(parser, 'tf', 0.8);
            addOptional(parser, 'dt', 0.02);
            addOptional(parser, 'm2min', [0,0]);
            addOptional(parser, 'sfactor', 3);
            %addOptional(parser, 'ConstraintType', 1);
            parse(parser, robot, x0, yf, varargin{:});
            
            dt = parser.Results.dt;
            tf = parser.Results.tf;
            sfactor = parser.Results.sfactor;
            
            m2min = parser.Results.m2min;
            % if length(x0)==2 ---- joint angles
            if length(x0) == 2
                x0 = [x0(1);x0(2);0;0; x0(1); m2min(1); 0; 0;x0(2); m2min(2); 0; 0];
            end
            
            t = 0:dt:tf;
            Nt = tf/dt + 1;
            
            f = @(x,u)robot.dynamics_with_jacobian_fd( x, u);
            %%%% Define parameters of cost function
            costpara = [];
            costpara.yf = yf;
            %costpara.xf = xf;
            %costpara.Hf = Hf;
            costpara.robot = robot;
            costpara.tf = tf;
            costpara.w_tf = 1e3;
            costpara.w_t  = 1e3; %
            costpara.w_e  = sfactor; % 1 is good for spring force
            costpara.L = robot.L ; % link lengths
            %%%%
            j = @(x,u,t)TrajOptMccp2MD.jReachFastOs(x, u, t, costpara);
            
            %%%% Define parameters of optimisation 
            opt_param = [];
            
            opt_param.lambda_init = 1;
            opt_param.lambda_max  = 1e10;
            opt_param.iter_max = 150;
            opt_param.online_plotting = 0;
            opt_param.online_printing = 1;
            opt_param.dcost_converge = 1e-3;
            opt_param.solver = 'rk4';
            
            opt_param.T = tf;
            %%%
            opt_param.umax = [pi/2; pi/2; 1; pi/2; pi/2; 1 ];
            opt_param.umin = [-pi/2; m2min(1); 0; -pi/2; m2min(2); 0];
            % u0 can be full command sequence or just initial point
            u0 = [0.3; m2min(1)+0.1; 0.1; 0.3 ;m2min(2)+0.1 ; 0.1];
            %u0 = [0; 0.1; 0];
            %%%%
            tic
            %traj = ILQRController.ilqr(f, j, dt, Nt, x0, u0, opt_param);
            
            traj = iLQRv1(f, j, dt, Nt, x0, u0, opt_param);
            traj.t = t;
            %traj = val_traj_mccpvd(robot, task, traj);
            %if nargin == 0
            %    plot_traj_mccpvd1(traj);
            %end
            
            toc
            disp('ILQR done')
            %%% Compute the costs
            %traj.J1 = sum((traj.x(1,:) - qf).^2, 2)*dt*costpara.w_t;
            %traj.Jtf = costQuad(traj.x(:,end) - xf, Hf*costpara.w_tf);
            %traj.Jp = traj.J1 + traj.Jtf;
            %traj.Je = sum((traj.u(1,1:end) - traj.x(1,1:end-1)).^2 + traj.u(2,:).^2,2)*dt*costpara.w_e ;
            %traj = evaluate_traj_Mccp1MD(robot, traj, 'qf', qf);
        end
        
        % optimise the joint space target
        function [traj]   = reachFastIlqr(robot, x0, qf, varargin)
            % reach fast to target in operation space solved by ILQR
            % x0, yf, tf, dt, m2min
            parser = inputParser();
            
            addRequired(parser, 'robot');
            addRequired(parser, 'x0');
            addRequired(parser, 'qf');
            addOptional(parser, 'tf', 0.8);
            addOptional(parser, 'dt', 0.02);
            addOptional(parser, 'm2min', [0,0]);
            addOptional(parser, 'sfactor', 3);
            %addOptional(parser, 'ConstraintType', 1);
            parse(parser, robot, x0, qf, varargin{:});
            
            dt = parser.Results.dt;
            tf = parser.Results.tf;
            sfactor = parser.Results.sfactor;
            
            m2min = parser.Results.m2min;
            % if length(x0)==2 ---- joint angles
            if length(x0) == 2
                x0 = [x0(1);x0(2);0;0;  x0(1); m2min(1); 0; 0; x0(2); m2min(2); 0; 0];
            end
            xf = [ qf(1);qf(2);0;0;...
                    qf(1);m2min(1);0;0;...
                    qf(2);m2min(2);0;0];
            Hf = diag(ones(1,12));
            t = 0:dt:tf;
            Nt = tf/dt + 1;
            
            f = @(x,u)robot.dynamics_with_jacobian_fd( x, u);
            %%%% Define parameters of cost function
            costpara = [];
            costpara.qf = qf;
            costpara.xf = xf;
            costpara.Hf = Hf;
            costpara.robot = robot;
            costpara.tf = tf;
            costpara.w_tf = 1e3;
            costpara.w_t  = 0.1; %
            costpara.w_e  = 1*sfactor; % 1 is good for spring force
            costpara.L = robot.L ; % link lengths
            %%%%
            j = @(x,u,t)TrajOptMccp2MD.jReachFast(x, u, t, costpara);
            
            %%%% Define parameters of optimisation 
            opt_param = [];
            
            opt_param.lambda_init = 0.01;
            opt_param.lambda_max  = 1e10;
            opt_param.iter_max = 150;
            opt_param.online_plotting = 0;
            opt_param.online_printing = 2;
            opt_param.dcost_converge = 1e-3;
            opt_param.solver = 'rk4';
            
            opt_param.T = tf;
            %%%
            opt_param.umax = [pi/2; pi/2; 1; pi/2; pi/2; 1 ] ;
            opt_param.umin = [-pi/2; m2min(1); 0; -pi/2; m2min(2); 0] ;
            % u0 can be full command sequence or just initial point
            %u0 = [qf(1); m2min(1)+0.1; 0; qf(2) ;m2min(2)+0.1 ; 0];
            u0 = [2*pi*rand(1)/2-0.5; rand(1)*pi/2 + m2min(1) ; rand(1) ; 2*pi*rand(1)/2-0.5 ; rand(1)*pi/2 ; rand(1) ];
            %u0 = zeros(6,1);
            %u0 = [0; 0.1; 0];
            %%%%
            tic
            %traj = ILQRController.ilqr(f, j, dt, Nt, x0, u0, opt_param);
            
            traj = iLQRv1(f, j, dt, Nt, x0, u0, opt_param);
            traj.t = t;
            %traj = val_traj_mccpvd(robot, task, traj);
            %if nargin == 0
            %    plot_traj_mccpvd1(traj);
            %end
            
            toc
            disp('ILQR done')
            %%% Compute the costs
            %traj.J1 = sum((traj.x(1,:) - qf).^2, 2)*dt*costpara.w_t;
            %traj.Jtf = costQuad(traj.x(:,end) - xf, Hf*costpara.w_tf);
            %traj.Jp = traj.J1 + traj.Jtf;
            %traj.Je = sum((traj.u(1,1:end) - traj.x(1,1:end-1)).^2 + traj.u(2,:).^2,2)*dt*costpara.w_e ;
            %traj = evaluate_traj_Mccp1MD(robot, traj, 'qf', qf);
        end
        
        function [l, l_x, l_xx, l_u, l_uu, l_ux] = jReachFastOs( x, u, t, para)
            w_tf = para.w_tf;
            w_t  = para.w_t; %
            w_e  = para.w_e; % 1 is good for spring force
            %Rtf  = diag([w_tf, 0, 0, 0, 0, 0]);
            
            if (isnan(u))
                % final cost
                fh = @(x)costf_reachFast_os(x, para);
                l = fh(x);
                if nargout > 1
                    l_x = get_jacobian_fd(fh, x);
                    l_xx = get_hessian_fd(fh, x);
                end
            else
                % running cost
                %para = [];
                para.w_t = w_t;
                para.w_e = w_e;
                fl = @(x,u,t) costr_reachFast_os(x,u,t,para);
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
        
        function [l, l_x, l_xx, l_u, l_uu, l_ux] = jReachFast( x, u, t, para)
            w_tf = para.w_tf;
            w_t  = para.w_t; %
            w_e  = para.w_e; % 1 is good for spring force
            %Rtf  = diag([w_tf, 0, 0, 0, 0, 0]);
            
            if (isnan(u))
                % final cost
                fh = @(x)costf_reachFast(x, para);
                l = fh(x);
                if nargout > 1
                    l_x = get_jacobian_fd(fh, x);
                    l_xx = get_hessian_fd(fh, x);
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
        
    end
end

function c = costr_reachFast_os(x, u, t, para)
            %c1 = (x(1) - para.qf).^2*t/para.tf;
            y = Arm2Dof.endpoint(x, para.L);
            c1 = sum((y - para.yf).^2); %squared distance
            
            %ce1 = (para.robot.actuator1.spring_force(x(1),x(5),x(6))).^2;
            %ce2 = (para.robot.actuator2.spring_force(x(2),x(9),x(10))).^2;
            %ce1 = (para.robot.actuator1.spring_force(x(1),u(1),u(2))).^2;
            %ce2 = (para.robot.actuator2.spring_force(x(2),u(4),u(5))).^2;
            
            
            ce1 = (u(1) - x(1))^2 + u(2)^2 + 0.5*(u(3)-0.5)^2;
            ce2 = (u(4) - x(2))^2 + u(5)^2 + 0.5*(u(6)-0.5)^2;
            
            %ce = (u(1) - x(3))^2 + (u(2) - x(4))^2 + 1e-6*u(3)^2;
            ce = ce1+ce2;
            
            c = c1*para.w_t + ce*para.w_e + norm(u)*1e-6 ;
            
end

function c = costf_reachFast_os(x, para)
    y = Arm2Dof.endpoint(x(1:2,:), para.L);
    c = sum((y - para.yf).^2)*para.w_tf ;
end

function c = costr_reachFast(x, u, t, para)
            %c1 = (x(1) - para.qf).^2*t/para.tf;
            %y = Arm2Dof.endpoint(x, para.L);
            %c1 = sum((y - para.yf).^2); %squared distance
            c1 = sum((x(1:2,:) - para.qf).^2); %squared distance
            %ce1 = (para.robot.actuator1.spring_force(x(1),x(5),x(6))).^2;
            %ce2 = (para.robot.actuator2.spring_force(x(2),x(9),x(10))).^2;
            %ce1 = (para.robot.actuator1.spring_force(x(1),u(1),u(2))).^2;
            %ce2 = (para.robot.actuator2.spring_force(x(2),u(4),u(5))).^2;
            
            
            ce1 = (u(1) - para.qf(1))^2 + u(2)^2 + 0.5*(u(3)-0.5)^2;
            ce2 = (u(4) - para.qf(2))^2 + u(5)^2 + 0.5*(u(6)-0.5)^2;
            
            %ce = (u(1) - x(3))^2 + (u(2) - x(4))^2 + 1e-6*u(3)^2;
            ce = ce1+ce2;
            %ce = sum(u.^2);
            %c1 = 0;
            c = c1*para.w_t + ce*para.w_e +1e-6*norm(u) ;
            
end

function c = costf_reachFast(x, para)
    %y = Arm2Dof.endpoint(x(1:2,:), para.L);
    c = 0.5*para.w_tf*x'*para.Hf*x ;
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
