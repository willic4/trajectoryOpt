global gridN maxacc initialposition initialvelocity w impactpos legendlist idx im fig
gridN = 50;
initialposition = [0, 0, 0];
initialvelocity = [0, 500, 0];
impactpos = [-40,60,20];
legendlist = [];
idx=1;
% spacecraft parameters
% mass0 = 4535.9; % initial mass (kilograms)
% mdot = 5.85;    % propellant flow rate (kilograms/day)
% beta = (mdot / mass0) * TU; % normalized propellant flow rate
% thrmag = 3.781; % thrust (newtons)
% maxacc = (thrmag / mass0); % normalized thrust acceleration
% rc = 7000;
% w = sqrt(398600/rc^3);
tic
% Minimize the simulation time
time_min = @(x) x(1);
speed_max = @(x) -(x(1 + gridN * 4)^2 + x(1 + gridN * 5)^2 + x(1 + gridN * 6)^2);
J = speed_max;
% The initial parameter guess; 1 second, gridN x, gridN y, gridN z, gridN
% xdot, gridN ydot, gridN zdot, gridN alphas, gridN phis, gridN Ts
x0 = [1; linspace(500,1,gridN)'; linspace(1,1,gridN)'; linspace(1,1,gridN)'; linspace(1,1,gridN)'; linspace(50,1,gridN)'; linspace(1,1,gridN)'];
% No linear inequality or equality constraints
A = [];
b = [];
Aeq = [];
Beq = [];
% Aeq = [0,1,0,0,0,0,0,zeros(1,(gridN-1)*6);...
%        0,0,1,0,0,0,0,zeros(1,(gridN-1)*6);...
%        0,0,0,1,0,0,0,zeros(1,(gridN-1)*6)];
% Beq = [initialposition(1);initialposition(2);initialposition(3)];
% Lower bound the simulation time at zero seconds, and bound the
% accelerations between -10 and 30
lb = [1;    ones(gridN * 3, 1) * -Inf; -50 ; ones(gridN-1, 1) * -Inf; -50 ; ones(gridN-1, 1) * -Inf; -50 ; ones(gridN-1, 1) * -Inf];
ub = [Inf;  ones(gridN * 3, 1) * Inf ; 50; ones(gridN-1, 1) *  Inf; 50; ones(gridN-1, 1) *  Inf; 50; ones(gridN-1, 1) *  Inf];
% Options for fmincon
options = optimoptions(@fmincon, 'TolFun', 0.001, 'TolX',0.001, 'MaxIter', 100, ...
                       'MaxFunEvals', 100000, 'Display', 'iter', ...
                       'DiffMinChange', 0.001, 'Algorithm', 'sqp','OutputFcn', @out);
% Solve for the best simulation time + control input
fig=figure();
optimal = fmincon(J, x0, A, b, Aeq, Beq, lb, ub, ...
              @double_integrator_constraints2, options);
close();
% Discretize the times
sim_time = optimal(1);
delta_time = sim_time / gridN;
times = 0 : delta_time : sim_time - delta_time;
% Get the state + accelerations (control inputs) out of the vector
xs = optimal(2             : 1 + gridN);
ys      = optimal(2 + gridN     : 1 + gridN * 2);
zs      = optimal(2 + gridN * 2 : 1 + gridN * 3);
xds     = optimal(2 + gridN * 3 : 1 + gridN * 4);
yds     = optimal(2 + gridN * 4 : 1 + gridN * 5);
zds     = optimal(2 + gridN * 5 : 1 + gridN * 6);

%% plotting


filename = 'testanimated.gif'; % specify file name
for idxx = 2:idx
    [A,map] = rgb2ind(im{idxx},256);
    if idxx - 1 == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.2);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
end


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
figure();
plot(times, xs)
hold on
plot(times, ys)
hold on
plot(times, zs)
title('Position vs Time');
legend('x', 'y', 'z');
figure();
plot(times, xds)
hold on
plot(times, yds)
hold on
plot(times, zds)
title('velocity vs Time');
legend('x', 'y', 'z');
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
figure()
plot3(xs,ys,zs,'.')
axis equal
title('3D Relative Motion Plot - TOF Minimized');
xlabel('x');
ylabel('y');
zlabel('z');
grid on
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
disp('final speed is')
disp(sqrt(xds(end)^2 + yds(end)^2 + zds(end)^2))
disp('TOF is')
disp(sim_time)
function [ c, ceq ] = double_integrator_constraints2( x )
global gridN maxacc initialposition initialvelocity w impactpos
    c=[];
    % Calculate the timestep
    sim_time = x(1);
    delta_time = sim_time / gridN;
    % Get the states / inputs out of the vector
    xs = x(2             : 1 + gridN);
    ys      = x(2 + gridN     : 1 + gridN * 2);
    zs      = x(2 + gridN * 2 : 1 + gridN * 3);
    xds     = x(2 + gridN * 3 : 1 + gridN * 4);
    yds     = x(2 + gridN * 4 : 1 + gridN * 5);
    zds     = x(2 + gridN * 5 : 1 + gridN * 6);

%     c=[alphas(1);phis(1)]
%     for i = 1:length(xs)-1
%         c=[c; (alphas(i+1)-alphas(i))/delta_time - 5; (phis(i+1)-phis(i))/delta_time - 5]; %set max rate of change of angles to 5 degrees/sec
%     end
    % Constrain initial position and velocity
    
%     ceq = [initialposition(1) - xs(1); initialposition(2) - ys(1);  initialposition(3) - zs(1); initialvelocity(1) - xds(1); initialvelocity(2) - yds(1); initialvelocity(3) - zds(1); alphas(1); phis(1)];
    ceq=[initialposition(1) - xs(1); initialposition(2) - ys(1);  initialposition(3) - zs(1)];
    for i = 1 : length(xs) - 1
        % The state at the beginning of the time interval
        x_i = [xs(i); ys(i); zs(i); xds(i); yds(i); zds(i)];
        xdot_i = [xds(i); yds(i); zds(i);0;0;-9.81];
        % What the state should be at the start of the next time interval
        x_n = [xs(i+1); ys(i+1); zs(i+1); xds(i+1); yds(i+1); zds(i+1)];
        xdot_n = [xds(i+1); yds(i+1); zds(i+1);0;0;-9.81];
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
    ceq = [ceq ; xs(end) - impactpos(1); ys(end) - impactpos(2); zs(end) - impactpos(3)];
end



%Output Function
function stop = out(x,omtimValue,state)
global gridN legendlist im idx fig 
if strcmp(state,'iter')
    idx=idx+1;
    xs = x(2             : 1 + gridN);
    ys      = x(2 + gridN     : 1 + gridN * 2);
    zs      = x(2 + gridN * 2 : 1 + gridN * 3);
    plot3(xs,ys,zs,'.','MarkerSize',8);
    axis equal
    grid on
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);
end
stop = false;
end