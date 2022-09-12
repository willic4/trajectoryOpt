global G m1 m2 m3 M
%% known values
G = 6.673e-20; 
rho = 1e6; % km 
rOV = 4.60e6; %sec
m1 = 5.967e23;
% m2 = 5.967e23;
% m3 = 5.967e23;

m2 = 7.35e22;
m3 = 3.675e22; %kg
fpa = 40; % degrees


M = m1 + m2 + m3; %kg

%% Known constraints

% r0(1,1) = rho/2;
% r0(3,1) = -rho/2;
% r0(4,1) = m3*rho*sind(60)/(m1+m2+m3);
% r0(2,1) = r0(4);
% 
% r0(5,1) = 0;
% r0(6,1) = (m1 + m2)*r0(2)/m3;

% r0(1,1) = -rho*(m2 + m3*cosd(60))/M;
% r0(3,1) = r0(1) + rho;
% r0(5,1) = r0(1) + rho * cosd(60);
% r0(2,1) = rho*(-m3)*cosd(60)/M;
% r0(4,1) = -rho*sind(60) + r0(2);
% r0(6,1) = -rho*sind(60) + r0(2);


r0(1,1) = -rho*(m2 + m3*cosd(60))/M;
r0(3,1) = r0(1) + rho;
r0(5,1) = r0(1) + rho * cosd(60);
r0(6,1) = rho*(m1 + m2)*sind(60)/M;
r0(4,1) = -rho*sind(60) + r0(6);
r0(2,1) = -rho*sind(60) + r0(6);

R3 = [cosd(fpa),   sind(fpa);...
      -sind(fpa), cosd(fpa)];
v0(1:2,1) = R3*-r0(1:2)/rOV;
v0(3:4,1) = R3*-r0(3:4)/rOV;
v0(5:6,1) = R3*-r0(5:6)/rOV;


%% Integrate and plot
Y0 = [r0;v0];
tspan = [0 5e7];
options = odeset('RelTol',1e-11,'AbsTol',1e-13);

[t,Y] = ode113(@eqTriangleSolution, tspan, Y0, options);

figure()
hold on
plot(Y(:,1), Y(:,2))
plot(Y(:,3), Y(:,4))
plot(Y(:,5), Y(:,6))

grid on
scatter(r0(1,1),r0(2,1))
scatter(r0(3,1),r0(4,1))
scatter(r0(5,1),r0(6,1))

legend('m1','m2','m2')
axis equal


function [dYdt] = eqTriangleSolution(t,Y)
global G m1 m2 m3 M
% Y is a [12x1] state vector consisting of positions and velocities of the
% 3 bodies in this system 
% x = [ x1, y1, x2, y2, x3, y3, dx1, dy1, dx2, dy2, dx3, dy3]';
% t is time in seconds

%% Known values
eta1 = G/M^2*(m2^2+m3^2+m2*m3)^(3/2);
eta2 = G/M^2*(m1^2+m3^2+m1*m3)^(3/2);
eta3 = G/M^2*(m2^2+m1^2+m2*m1)^(3/2);

dYdt = [Y(7:12);...
        -eta1/norm(Y(1:2,1))^3*Y(1);...
        -eta1/norm(Y(1:2,1))^3*Y(2);...
        -eta2/norm(Y(3:4,1))^3*Y(3);...
        -eta2/norm(Y(3:4,1))^3*Y(4);...
        -eta3/norm(Y(5:6,1))^3*Y(5);...
        -eta3/norm(Y(5:6,1))^3*Y(6)];
end


