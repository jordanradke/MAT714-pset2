%% simple implementation of 2d wave equation
function [w] = wave_equation(N,T,f,g,bdy)
% parameters:
% N: grid resolution
% T: number of time steps
% f: u_0
% g: d/dt(u_0)

dx = 1/N;
dy = 1/N;
dt = 1/N^2;

x = 0:dx:1;
y = 0:dy:1;

[X,Y] = meshgrid(x,y);

% initial conditions
%f  = zeros(size(X));
%g  = exp(-400*(X-0.5).^2).*exp(-400*(Y-0.5).^2);  % change back when done testing

% boundary conditions
%bdy = zeros(size(X));

u = zeros([size(X),T]);

% initialization to generate U^0 and U^1. Just via taylor's theorem.
u(:,:,2) = f + dt*g;

for t=3:T
    % ensure boundary conditions are satisified
    u(:,1,t)   = bdy(:,1);
    u(:,N+1,t) = bdy(:,N+1);
    u(1,:,t)   = bdy(1,:);
    u(N+1,:,t) = bdy(N+1,:);
    
    for i = 2:N
        for j =2:N
            u(i,j,t) = dt^2/(dx)^2*(u(i,j-1,t-1)+u(i-1,j,t-1)+u(i+1,j,t-1)+u(i,j+1,t-1)) ...
                       + (2 - 4*dt^2/(dx)^2)*u(i,j,t-1) - u(i,j,t-2);
        end
    end
end

w = u;

end

%{
% animating surface to test
figure
s = surf(x,y,w(:,:,1));

light               % add a light
lighting gouraud    % preferred lighting for a curved surface
axis equal off      % set axis equal and remove axis
view(40,30)         % set viewpoint
camzoom(1.5)        % zoom into scene

scale = [linspace(0,1,20) linspace(1,-1,40)];    % surface scaling (0 to 1 to -1)

for t = 1:T
    z = w(:,:,t);

    s.XData = x;    % replace surface x values
    s.YData = y;    % replace surface y values
    s.ZData = z;    % replace surface z values

    pause(0.5)     % pause to control animation speed
end
%}
