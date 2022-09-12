r_earth = 6371e3; % radius of earth in m
M = 5.972e24; % kg % mass of the Earth
G = 6.67e-11; % gravitational constant
mu = M*G; % earth gravitational parameter

t0 = 0;
N = 10000;
t = linspace(t0,N,N);


options=odeset('events',@fcn,'RelTol',1e-13,'AbsTol',1e-15);
y = [0.5819596353742e6, 7.849456285818026e6, 0.007094430247474e6,-0.000525981913864e6];
[T, X] = ode113(@CR2BP, t, y);


error = abs(X(end,1) - X(1,1));

figure()

circle(0,0,r_earth)
hold on 
plot(X(:,1), X(:,2), 'r')
axis equal



function ydot=CR2BP(t,y)
M = 5.972e24; % kg % mass of the Earth
G = 6.67e-11; % gravitational constant
mu = M*G; % earth gravitational parameter

r1=sqrt((y(1))^2+(y(2))^2);

ydot=[y(3); 
    y(4);  
    -mu/norm(r1)^3 * y(1); 
    -mu/norm(r1)^3 * y(2)];
end

function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end

function [position,isterminal,direction] = fcn(t,y)
    position = y(2); % We want position_y = 0 (crossing the x axis)
    isterminal = 1;  % Halt integration 
    direction = +1;   % The zero can be approached from either direction
end