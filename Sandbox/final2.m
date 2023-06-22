



% 5 state problem solver - x pos, y pos, x vel, y vel, lateral acceleration
% command

%% initialization
close all;
clear;
constants
global initialstate gridN im idx fig finalr finalstate mu
idx = 0;

filename = 'newcs.gif'; % specify file name
outFolderName = "C:\Users\billy\Documents\GitHub\trajectoryOpt\Sandbox\Presentationworthy\excursion_dalphalimit_2_minthrustandtime\";   % Your destination folder
filename = strcat(outFolderName,filename)
mkdir(outFolderName)

%% trajectory options
finalr = 1500+r_earth; %km
initialstate = [0,r_earth,0,0];
finalstate = [500,r_earth+1500,0,0];
% maxacc = 0.5; %max lateral acceleration
vehiclemass = 30000; %mass in kg
maxthrust = 15000; %max thrust in kN
maxacc = maxthrust/vehiclemass; %max acceleration in km/s^2


%% optimizer options
gridN = 25;
generateplots = 1;
generategif = 1;
TolFun = 0.0001;
TolX = 0.0000000001;
TolCon = 1e-3;

MaxIter = 1e4;
MaxFunEvals = 1e7;
DiffMinChange = 0.001;
Algorithm = 'sqp';
%weights for cost function
w1 = 1;
w2 = 0;
w3 = 1;





% J = @(x) w1*(abs(finalr - sqrt((finalstate(1)-x(gridN+1))^2 + (finalstate(2)-x(2*gridN+1))^2)) + w2*(finalstate(3)-x(3*gridN+1))^2 + (finalstate(4)-x(4*gridN+1))^2)+w3*(x(1)); %cost function - designed to minimize insertion error
% J = @(x) w1*(abs(finalr - sqrt(x(gridN+1)^2 + x(2*gridN+1)^2))) + w3*(x(1)); %cost function - designed to minimize insertion error

J = @(x) w1*(sum(abs(x(2 + gridN * 4 : 1 + gridN * 5))))+ w2*(sum(abs(x(2 + gridN * 5 : 1 + gridN * 6))))... minimize total accelerations
    +w3*(x(1)); %cost function - designed to minimize insertion error

% initial state guess

x_guess = [1; linspace(initialstate(1),finalstate(1),gridN)'; linspace(initialstate(2),finalstate(2),gridN)'; linspace(1,1,gridN)'; linspace(1,1,gridN)'; linspace(maxacc,maxacc,gridN)'; linspace(1,1,gridN)'];

A=[];
b=[];
Aeq=[];
beq=[];
lb = [1;    ones(gridN * 1, 1) * 0; ones(gridN * 1, 1) * r_earth; 0 ; ones(gridN-1, 1) * -Inf; 0 ; ones(gridN-1, 1) * -Inf; ones(gridN, 1) * 0; ones(gridN, 1) * -pi];
ub = [Inf;  ones(gridN * 1, 1) * Inf ; ones(gridN * 1, 1) * finalr; 50; ones(gridN-1, 1) *  Inf; 50; ones(gridN-1, 1) *  Inf; ones(gridN, 1) *  maxacc; ones(gridN, 1) *  pi];


%% setting up and calling optimizer (fmincon)
% Options for fmincon
optionslist = strcat("@fmincon, 'TolFun', ",num2str(TolFun),", 'TolX',",num2str(TolX),", 'TolCon',",num2str(TolCon),", 'MaxIter',",num2str(MaxIter),",'MaxFunEvals', ",num2str(MaxFunEvals),...
    ", 'Display', 'iter','DiffMinChange',",num2str(DiffMinChange),", 'Algorithm', '", Algorithm,"'");
if generategif
    optionslist = strcat(optionslist,",'OutputFcn', @out");
    fig=figure('Position',[10 10 1200 800]);
end
optionscmd = strcat('options = optimoptions(',optionslist,');');
eval(optionscmd);

optimal = fmincon(J, x_guess, A, b, Aeq, beq, lb, ub, ...
              @double_integrator_constraints2, options);

close();


%% process result
% Discretize the times
sim_time = optimal(1);
delta_time = sim_time / gridN;
times = 0 : delta_time : sim_time - delta_time;
% Get the state + accelerations (control inputs) out of the vector
xs          = optimal(2             : 1 + gridN);
ys          = optimal(2 + gridN     : 1 + gridN * 2);
xds         = optimal(2 + gridN * 2 : 1 + gridN * 3);
yds         = optimal(2 + gridN * 3 : 1 + gridN * 4);
xcmd        = optimal(2 + gridN * 4 : 1 + gridN * 5);
ycmd        = optimal(2 + gridN * 5 : 1 + gridN * 6);

% fprintf("x pos Error is: " + num2str(xs(end) - finalstate(1)) + "\n");
% fprintf("y pos Error is: " + num2str(ys(end) - finalstate(2)) + "\n");
% fprintf("x vel Error is: " + num2str(xds(end) - finalstate(3)) + "\n");
% fprintf("y vel Error is: " + num2str(yds(end) - finalstate(4)) + "\n");
fprintf("position Error is: " + num2str(finalr - sqrt(ys(end)^2 + xs(end)^2)) + "\n");

fprintf("Time of flight is: " + num2str(sim_time) + "\n");
format long
final_state = [xs(end);ys(end);xds(end);yds(end)]
%% plotting 
if generategif
    for idxx = 2:idx
        [A,map] = rgb2ind(im{idxx},256);
        if idxx - 1 == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        elseif idxx == idx
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',2);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.01);
        end
    end
end
if generateplots
    % positions
    figure()
    plot(xs,ys, '.')
    axis equal
    title('2D relative motion plot');
    xlabel('x');
    ylabel('y');

    t0 = 0;
    N = 10000;
    t = linspace(t0,N,N/10);
    options=odeset('RelTol',1e-13,'AbsTol',1e-15);
    y = [xs(end), ys(end), xds(end),yds(end)];
    [T, X] = ode113(@CR2BP, t, y);
    
    figure()
    subplot(1,3,1)
    plot(xs,ys,'.','MarkerSize',8);
    hold on
    plot(X(:,1), X(:,2), 'r')
    circle(0,0,r_earth);
    hold off
    xlabel('X (m)')
    ylabel('Y (m)')
    legend('launch trajectory','resulting orbit','earth')
    axis equal
    grid on
    subplot(1,3,2)
    plot(xs,ys,'.','MarkerSize',8);
    title('Trajectory')
    axis equal
    subplot(1,3,3)
    plot(times,xcmd);
    hold on
    plot(times,ycmd);
    title('Acceleration commands (km/s^2)')
    legend('x','y')
    ylim([-0.5,0.5])
    hold off

    % variables vs time
    figure();
    plot(times, xs)
    hold on
    plot(times, ys)
    title('Position vs Time');
    legend('x', 'y', 'z');

    figure();
    plot(times, xds)
    hold on
    plot(times, yds)
    title('velocity vs Time');
    legend('x', 'y', 'z');

    figure();
    plot(times, xcmd)
    hold on
    plot(times, ycmd)
    title('acceleration command vs Time');
    legend('x', 'y');

    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = num2str(get(FigHandle, 'Number'));
      set(0, 'CurrentFigure', FigHandle);
      savefig(fullfile(outFolderName, [FigName '.fig']));
    end
end



%% nonlinear constraint 
function [ c, ceq ] = double_integrator_constraints2( x )
global gridN initialstate finalstate finalr mu r_earth
    c=[];
    % Calculate the timestep
    sim_time = x(1);
    delta_time = sim_time / gridN;
    % Get the states / inputs out of the vector
    xs = x(2             : 1 + gridN);
    ys      = x(2 + gridN     : 1 + gridN * 2);
    xds      = x(2 + gridN * 2 : 1 + gridN * 3);
    yds     = x(2 + gridN * 3 : 1 + gridN * 4);
    xcmd      = x(2 + gridN * 4 : 1 + gridN * 5); %axial acceleration
    ycmd     = x(2 + gridN * 5 : 1 + gridN * 6); %lateral acceleration
    % constrain initial position
    ceq=[initialstate(1) - xs(1); initialstate(2) - ys(1); initialstate(3) - xds(1); initialstate(4) - yds(1)];
    for i = 1 : length(xs) - 1
        c = [c; abs(ycmd(i+1) - ycmd(i)) - 10*pi/180];
        %forces at the beginning of the time interval
        r_i = [xs(i);ys(i)];
        grav_i = -mu/norm(r_i)^3 * r_i; % constrained 2BP
        % convert thrust into x-y frame
%         v_i = [xds(i);yds(i)];
%         thrust_body_i = [xcmd(i);ycmd(i)];
%         theta_i = atan2(v_i(2),v_i(1)); %angle between ecef frame and body velocity frame
%         R_ecef_from_body = [cos(theta_i),-sin(theta_i); sin(theta_i),cos(theta_i)];
%         thrust_ecef_i = R_ecef_from_body * thrust_body_i;


        thrust_ecef_i = [xcmd(i)*cos(ycmd(i));xcmd(i)*sin(ycmd(i))];
        % The state at the beginning of the time interval\

        x_i = [xs(i); ys(i); xds(i); yds(i)];
        xdot_i = [xds(i); yds(i);thrust_ecef_i(1)+grav_i(1);thrust_ecef_i(2)+grav_i(2)];
        %forces at the end of the time interval
        r_n = [xs(i+1);ys(i+1)];
        grav_n = -mu/norm(r_n)^3 * r_n; % constrained 2BP
        % convert thrust into x-y frame
%         v_n = [xds(i+1);yds(i+1)];
%         thrust_body_n = [xcmd(i+1);ycmd(i+1)];
%         theta_n = atan2(v_n(2),v_n(1)); %angle between ecef frame and body velocity frame
%         R_ecef_from_body = [cos(theta_n),sin(theta_n); -sin(theta_n),cos(theta_n)];
%         thrust_ecef_n = R_ecef_from_body * thrust_body_n;

        thrust_ecef_n = [xcmd(i)*cos(ycmd(i));xcmd(i)*sin(ycmd(i))];

        % What the state should be at the start of the next time interval
        x_n = [xs(i+1); ys(i+1); xds(i+1); yds(i+1)];
        xdot_n = [xds(i+1); yds(i+1);thrust_ecef_n(1)+grav_n(1);thrust_ecef_n(2)+grav_n(1)];
        x_end = x_i + delta_time*(xdot_i+xdot_n)/2;

        ceq = [ceq ; x_n - x_end];
    end
    % Constrain end position and end velocity to 0
    r_end = [xs(end);ys(end)];
    vmag = sqrt(mu/norm(r_end));
    v_end = [r_end(2);-r_end(1)]/norm(r_end)*vmag;
    ceq = [ceq ; sqrt(xs(end)^2 + ys(end)^2) - finalr; xds(end) - v_end(1); yds(end) - v_end(2)];
end


function stop = out(x,omtimValue,state)
global gridN legendlist im idx fig
if strcmp(state,'iter')
    idx=idx+1;
    sim_time = x(1);
    delta_time = sim_time / gridN;
    times = 0 : delta_time : sim_time - delta_time;
    r_earth = 6371;
    xs = x(2             : 1 + gridN);
    ys      = x(2 + gridN     : 1 + gridN * 2);
    xds      = x(2 + gridN*2     : 1 + gridN * 3);
    yds      = x(2 + gridN*3     : 1 + gridN * 4);
    xcmd      = x(2 + gridN*4     : 1 + gridN * 5);
    ycmd      = x(2 + gridN*5     : 1 + gridN * 6);

    t0 = 0;
    N = 10000;
    t = linspace(t0,N,N/10);
    options=odeset('RelTol',1e-13,'AbsTol',1e-15);
    y = [xs(end), ys(end), xds(end),yds(end)];
    [T, X] = ode113(@CR2BP, t, y);
    
    subplot(2,2,1)
    plot(xs,ys,'.','MarkerSize',8);
    hold on
    plot(X(:,1), X(:,2), 'r')
    circle(0,0,r_earth);
    hold off
    title(strcat('Iteration: ',num2str(idx),'; TOF: ',num2str(x(1))))
    xlabel('X (km)')
    ylabel('Y (km)')
    legend('Launch Trajectory','Resulting Orbit','Earth')
    axis equal
    grid on
    subplot(2,2,2)
    plot(xs,ys,'.','MarkerSize',10);
    title('Trajectory')
    grid on
    axis equal
    xlabel('X (km)')
    ylabel('Y (km)')
    subplot(2,2,3:4)
    plot(times,xcmd);
    grid on
    hold on
    plot(times,ycmd);
    title('Commands')
    legend('Axial Acceleration (km/s^2)','alpha (rad)')
    ylim([-pi,pi])
    xlabel('time (s)')
    hold off
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);
end
stop = false;
end


function ydot=CR2BP(t,y)
global mu

r1=sqrt((y(1))^2+(y(2))^2);

ydot=[y(3); 
    y(4);  
    -mu/norm(r1)^3 * y(1); 
    -mu/norm(r1)^3 * y(2)];
end

function h = circle(xvar,yvar,rvar)
hold on
th = 0:pi/50:2*pi;
xunit = rvar * cos(th) + xvar;
yunit = rvar * sin(th) + yvar;
h = plot(xunit, yunit,'color', [0 0.5 0]);
hold off
end