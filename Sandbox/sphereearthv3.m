



% 5 state problem solver - x pos, y pos, x vel, y vel, lateral acceleration
% command

%% initialization
close all;
clear;
constants
global initialstate gridN im idx fig finalr finalstate mu
idx = 0;

filename = 'newcs.gif'; % specify file name
%% trajectory options
finalr = 1500+r_earth; %km
initialstate = [0,r_earth,0,0];
finalstate = [500,r_earth+1500,0,0];
maxacc = 0.5; %max lateral acceleration

%% optimizer options
gridN = 50;
generateplots = 1;
TolFun = 0.0001;
TolX = 0.0000000001;
TolCon = 1e-3;

MaxIter = 1e4;
MaxFunEvals = 1e7;
DiffMinChange = 0.001;
Algorithm = 'sqp';
%weights for cost function
w1 = 0.4;
w2 = 0.4;
w3 = 0.2;





% J = @(x) w1*(abs(finalr - sqrt((finalstate(1)-x(gridN+1))^2 + (finalstate(2)-x(2*gridN+1))^2)) + w2*(finalstate(3)-x(3*gridN+1))^2 + (finalstate(4)-x(4*gridN+1))^2)+w3*(x(1)); %cost function - designed to minimize insertion error
% J = @(x) w1*(abs(finalr - sqrt(x(gridN+1)^2 + x(2*gridN+1)^2))) + w3*(x(1)); %cost function - designed to minimize insertion error

J = @(x) w1*(sum(abs(x(2 + gridN * 4 : 1 + gridN * 5))))+ w2*(sum(abs(x(2 + gridN * 5 : 1 + gridN * 6))))... minimize total accelerations
    +w3*(x(1)); %cost function - designed to minimize insertion error

% initial state guess

x_guess = [1; linspace(initialstate(1),finalstate(1),gridN)'; linspace(initialstate(2),finalstate(2),gridN)'; linspace(1,1,gridN)'; linspace(1,1,gridN)'; linspace(1,1,gridN)'; linspace(1,1,gridN)'];

A=[];
b=[];
Aeq=[];
beq=[];
lb = [1;    ones(gridN * 1, 1) * -Inf; ones(gridN * 1, 1) * r_earth; 0 ; ones(gridN-1, 1) * -Inf; 0 ; ones(gridN-1, 1) * -Inf; ones(2*gridN, 1) * -maxacc];
ub = [Inf;  ones(gridN * 1, 1) * Inf ; ones(gridN * 1, 1) * finalr; 50; ones(gridN-1, 1) *  Inf; 50; ones(gridN-1, 1) *  Inf; ones(2*gridN, 1) *  maxacc];


%% setting up and calling optimizer (fmincon)
% Options for fmincon
optionslist = strcat("@fmincon, 'TolFun', ",num2str(TolFun),", 'TolX',",num2str(TolX),", 'TolCon',",num2str(TolCon),", 'MaxIter',",num2str(MaxIter),",'MaxFunEvals', ",num2str(MaxFunEvals),...
    ", 'Display', 'iter','DiffMinChange',",num2str(DiffMinChange),", 'Algorithm', '", Algorithm,"'");
if generateplots
    optionslist = strcat(optionslist,",'OutputFcn', @out");
    fig=figure();
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
if generateplots

    
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

    % positions
    figure()
    plot(xs,ys)
    axis equal
    title('3D relative motion plot');
    xlabel('x');
    ylabel('y');


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

end


%% nonlinear constraint 
function [ c, ceq ] = double_integrator_constraints2( x )
global gridN initialstate finalstate finalr mu
    c=[];
    % Calculate the timestep
    sim_time = x(1);
    delta_time = sim_time / gridN;
    % Get the states / inputs out of the vector
    xs = x(2             : 1 + gridN);
    ys      = x(2 + gridN     : 1 + gridN * 2);
    xds      = x(2 + gridN * 2 : 1 + gridN * 3);
    yds     = x(2 + gridN * 3 : 1 + gridN * 4);
    xcmd      = x(2 + gridN * 4 : 1 + gridN * 5);
    ycmd     = x(2 + gridN * 5 : 1 + gridN * 6);
    % constrain initial position
    ceq=[initialstate(1) - xs(1); initialstate(2) - ys(1); initialstate(3) - xds(1); initialstate(4) - yds(1)];
    for i = 1 : length(xs) - 1
        % The state at the beginning of the time interval
        x_i = [xs(i); ys(i); xds(i); yds(i)];
        xdot_i = [xds(i); yds(i);xcmd(i);ycmd(i)-9.81e-3];
        % What the state should be at the start of the next time interval
        x_n = [xs(i+1); ys(i+1); xds(i+1); yds(i+1)];
        xdot_n = [xds(i+1); yds(i+1);xcmd(i+1);ycmd(i)-9.81e-3];
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
    xs = x(2             : 1 + gridN);
    ys      = x(2 + gridN     : 1 + gridN * 2);
    plot(xs,ys,'.','MarkerSize',8);
    title(strcat('Iteration: ',num2str(idx),'; TOF: ',num2str(x(1))))
    xlabel('X (m)')
    ylabel('Y (m)')
    axis equal
    grid on
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);
end
stop = false;
end