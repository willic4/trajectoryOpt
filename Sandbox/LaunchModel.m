
% THE GOAL OF THIS CODE IS TO GENERATE AN OPTIMIZED LAUNCH ORBIT FROM
% WALLOPS ISLAND TO AN ORBIT WITH A SPECIFIC SET OF KEPLERIAN ORBITAL
% PARAMETERS
clear;
constants
global gridN maxacc initialstate initialvelocity w impactpos legendlist idx im fig generateplots mu Thrust a_desired



%% user inputs - mission
%initialpositionlla = [-75.45,37.91,0]; % NASA Wallops Island
initialstate = [0, r_earth, 0, 0, 0]; 
impactpos = [100,100,10];
Thrust = 10000;
a_desired = r_earth + 10e3; % target semi-major axis


%% user inputs - code
gridN = 50;
generateplots = 1;
TolFun = 0.0001;
TolX = 0.0000001;
MaxIter = 100;
MaxFunEvals = 100000;
DiffMinChange = 0.001;
Algorithm = 'sqp';

%% intialization
% initialpositionecef = lla2ecef(initialpositionlla)
% initialvelocity = [0, 0, 0];
legendlist = [];
idx=1;


%% set up optimization
tic
% Minimize the trajectory time
time_min = @(x) x(1);
% The initial parameter guess; 1 second, gridN x, gridN y, gridN xdot,
% gridN ydot, gridN alpha

% alpha is the angle between body and position x vectors
x0 = [1; [0:gridN:gridN^2-1]'; [r_earth:(a_desired/2-r_earth)/gridN:a_desired/2 - (a_desired/2-r_earth)/gridN]'; linspace(1,1,gridN)'; linspace(1,1,gridN)'; linspace(1,1,gridN)'];
% No linear inequality or equality constraints
A = [];
b = [];
% initial state is determined here
Aeq = [0,1,0,0,0,0,zeros(1,(gridN-1)*5);... 
       0,0,1,0,0,0,zeros(1,(gridN-1)*5);... 
       0,0,0,1,0,0,zeros(1,(gridN-1)*5);... 
       0,0,1,0,1,0,zeros(1,(gridN-1)*5);... 
       0,0,0,0,0,1,zeros(1,(gridN-1)*5)];   
Beq = [initialstate(1);initialstate(2);initialstate(3);initialstate(4);initialstate(5)];
% Lower bound the simulation time at zero seconds, and bound the
% angles between 0 and 90
lb = [0;    ones(gridN * 4, 1) * -Inf; ones(gridN, 1) * 0];
ub = [Inf;  ones(gridN * 4, 1) * Inf ; ones(gridN, 1) * 90];
% Options for fmincon
optionslist = strcat("@fmincon, 'TolFun', ",num2str(TolFun),", 'TolX',",num2str(TolX),", 'MaxIter',",num2str(MaxIter),",'MaxFunEvals', ",num2str(MaxFunEvals),...
    ", 'Display', 'iter','DiffMinChange',",num2str(DiffMinChange),", 'Algorithm', '", Algorithm,"'");
if generateplots
    optionslist = strcat(optionslist,",'OutputFcn', @out");
    fig=figure();
end
optionscmd = strcat('options = optimoptions(',optionslist,');');
eval(optionscmd);
% options = optimoptions(@fmincon, 'TolFun', 0.001, 'TolX',0.001, 'MaxIter', 100, ...
%                        'MaxFunEvals', 100000, 'Display', 'iter', ...
%                        'DiffMinChange', 0.001, 'Algorithm', 'sqp','OutputFcn', @out);
% Solve for the best simulation time + control input


optimal = fmincon(time_min, x0, A, b, Aeq, Beq, lb, ub, ...
              @double_integrator_constraints2, options);

close();
% Discretize the times
sim_time = optimal(1);
delta_time = sim_time / gridN;
times = 0 : delta_time : sim_time - delta_time;
% Get the state + accelerations (control inputs) out of the vector
xs          = optimal(2             : 1 + gridN);
ys          = optimal(2 + gridN     : 1 + gridN * 2);
xds         = optimal(2 + gridN * 2 : 1 + gridN * 3);
yds         = optimal(2 + gridN * 3 : 1 + gridN * 4);
alphas      = optimal(2 + gridN * 4 : 1 + gridN * 5);

%% plotting

if generateplots
    filename = 'testanimated.gif'; % specify file name
    for idxx = 2:idx
        [A,map] = rgb2ind(im{idxx},256);
        if idxx - 1 == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.2);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.2);
        end
    end

    figure()
    plot(xs,ys)
    hold on
    circle(0,0,r_earth)
    axis equal
    % Make the plots
    % figure();
    % plot(times, accs);
    % title('Control Input (Acceleration) vs Time');
    % xlabel('Time (s)');
    % ylabel('Acceleration (m/s^2)');
    % figure();
    % plot(times, vels);
    % title('Velocity vs Time');
    % xlabel('Time (s)');
    % ylabel('Velocity (m/s)');
%     figure();
%     plot(times, xs)
%     hold on
%     plot(times, ys)
%     % hold on
%     % plot(times, zs)
%     title('Position vs Time');
%     legend('x', 'y', 'z');
    % xlabel('Time (s)');
    % ylabel('Position (m or whatever I guess)');
    % disp(sprintf('Finished in %f seconds', toc));
    % 
    % figure();
    % plot(times, xds)
    % hold on
    % plot(times, yds)
    % hold on
    % plot(times, zds)
    % title('Relative Velocity vs Time');
    % legend('x', 'y', 'z');
    % xlabel('Time (s)');
    % ylabel('Velocity (m/s or whatever I guess)');
    % 
    % 
    % figure();
    % plot(times, alphas);
    % hold on
    % plot(times, phis)
    % title('angles vs Time');
    % xlabel('Time (s)');
    % ylabel('angles (deg)');
    % legend('alpha','phi')
    % 
    % figure();
    % plot(times, Ts);
    % title('Thrust vs Time');
    % xlabel('Time (s)');
    % ylabel('Thrust (N)');
    % 
    % figure()
    % plot3(xs,ys,zs)
    % title('3D relative motion plot');
    % xlabel('x');
    % ylabel('y');
    % zlabel('z');
end

function [ c, ceq ] = double_integrator_constraints2( x )
global gridN maxacc initialposition initialvelocity w impactpos Thrust mu a_desired
    c=[];
    % Calculate the timestep
    sim_time = x(1);
    delta_time = sim_time / gridN;
    % Get the states / inputs out of the vector
    xs = x(2             : 1 + gridN);
    ys      = x(2 + gridN     : 1 + gridN * 2);
    xds      = x(2 + gridN * 2 : 1 + gridN * 3);
    yds     = x(2 + gridN * 3 : 1 + gridN * 4);
    alphas     = x(2 + gridN * 4 : 1 + gridN * 5);

%     c=[alphas(1);phis(1)]
%     for i = 1:length(xs)-1
%         c=[c; (alphas(i+1)-alphas(i))/delta_time - 5; (phis(i+1)-phis(i))/delta_time - 5]; %set max rate of change of angles to 5 degrees/sec
%     end
    % Constrain initial position and velocity
    
%     ceq = [initialposition(1) - xs(1); initialposition(2) - ys(1);  initialposition(3) - zs(1); initialvelocity(1) - xds(1); initialvelocity(2) - yds(1); initialvelocity(3) - zds(1); alphas(1); phis(1)];
    ceq=[];
    for i = 1 : length(xs) - 1
        % The state at the beginning of the time interval
        x_i = [xs(i); ys(i); xds(i); yds(i)];
        r_i = [xs(i);ys(i)];
        gravity_i = -mu/norm(r_i)^3 * r_i; % constrained 2BP
        theta_i = atan2(xs(i),ys(i));
        Thrustvecbody_i = [Thrust*cos(theta_i-alphas(i)),Thrust*sin(theta_i-alphas(i))];
        xdot_i = [xds(i); yds(i); gravity_i(1) + Thrustvecbody_i(1); gravity_i(2) + Thrustvecbody_i(2)];

        % What the state should be at the start of the next time interval
        x_n = [xs(i+1); ys(i+1); xds(i+1); yds(i+1)];
        r_n = [xs(i);ys(i)];
        gravity_n = -mu/norm(r_n)^3 * r_n; % constrained 2BP
        theta_n = atan2(xs(i),ys(i));
        Thrustvecbody_n = [Thrust*cos(theta_n-alphas(i)),Thrust*sin(theta_n-alphas(i))];
        xdot_n = [xds(i+1); yds(i+1); gravity_n(1) + Thrustvecbody_n(1); gravity_n(2) + Thrustvecbody_n(2)];
        x_end = x_i + delta_time*(xdot_i+xdot_n)/2;
        % The time derivative of the state at the beginning of the time
        % interval
%         xdot_i = [u(i); acc * sin(alphas(i)); acc * cos(alphas(i))];
%         % The time derivative of the state at the end of the time interval
%         xdot_n = [u(i+1); acc * sin(alphas(i+1)); acc * cos(alphas(i+1))];
%         
%         A11 = zeros(3,3); A12 = eye(3);
%         A21 = [3*w^2 0 0; 0 0 0; 0 0 w^2];
%         A22 = [0 2*w 0; -2*w 0 0; 0 0 0];
%         Atot = [A11, A12; A21, A22];
%         STM = expm(Atot*delta_time);
        % The end state of the time interval calculated using quadrature
%         xend = STM * x_i;
        % Constrain the end state of the current time interval to be
        % equal to the starting state of the next time interval
        %CHANGE HERE TO INCORPERATE THRUST WITH TRAP
        ceq = [ceq ; x_n - x_end];
    end
    % Constrain end position and end velocity to 0
    ceq = [ceq ; xs(end)^2 + ys(end)^2 - a_desired];
end



%Output Function
function stop = out(x,omtimValue,state)
global gridN legendlist im idx fig 
if strcmp(state,'iter')
    idx=idx+1;
    xs = x(2             : 1 + gridN);
    ys      = x(2 + gridN     : 1 + gridN * 2);
    plot(xs,ys,'.','MarkerSize',8);
    grid on
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);
end
stop = false;
end

function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end