npoints=500;
dt = 0.002; % time step in years
x=1; % initialise position of planet in AU
y=0;
v_x=0; % initialise velocity of planet in AU/yr
v_y=2*pi;
% Plot the Sun at the origin
plot(0,0,'oy','MarkerSize',30, 'MarkerFaceColor','yellow');
axis([-1 1 -1 1]);
xlabel('x(AU)');
ylabel('y(AU)');
hold on 
num = xlsread('code.xlsx');
a = num(:,1);
b = num(:,2);
plot(a,b,'g')
title('Orbit Simulation')
hold off
hold on;
for step = 1:npoints-1
% loop over the timesteps
radius=sqrt(x^2+y^2);
% Compute new velocities in the x and y directions
v_x_new=v_x - (4*pi^2*x*dt)/(radius^2);
v_y_new=v_y - (4*pi^2*y*dt)/(radius^2);
% Euler Cromer Step - update positions using newly calculated velocities
x_new=x+v_x_new*dt;
y_new=y+v_y_new*dt;
% Plot planet position immediately
plot(x_new,y_new,'ob','MarkerSize',3, 'MarkerFaceColor','yellow');
drawnow;
% Update x and y velocities with new velocities
v_x=v_x_new;
v_y=v_y_new;
% Update x and y with new positions
x=x_new;
y=y_new;
end
legend('Sun','data( F ∝ 1/r) )','analytical solution( F ∝ 1/r'))