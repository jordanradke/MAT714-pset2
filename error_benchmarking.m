%% error testing against analytic fourier series

pp = 8;
err = zeros(pp,1);
h   = zeros(pp,1);

for p = 1:pp
    p
    % grid refinement + time steps
    N = 2^p;
    T = 50;
    dx = 1/N;
    dy = 1/N;
    dt = 1/N^2;

    x = 0:dx:1;
    y = 0:dy:1;

    [X,Y] = meshgrid(x,y);

    % initial conditions
    f  = zeros(size(X));
    g  = exp(-4*(X-0.5).^2).*exp(-4*(Y-0.5).^2);  
    
    % boundary conditions
    bdy = zeros(size(X));
    
    % numerical solution
    fprintf('numerical solution...')
    u = wave_equation(N,T,f,g,bdy);
    
    % analytic solution sampled on grid (only for f=0 init. data)
    fprintf('analytic solution...')
    u_sol = wave_fourier(X,Y,N,T,dt,g);
    
    err(p) = max(max(max(abs(u-u_sol))));
    h(p) = dx;
    
end

% plot loglog error 
figure(1); clf();
loglog(h,err,'o-', 'LineWidth', 2)
hold on;
loglog(h, h.^2, 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

title('Error','FontSize',24);
xlabel('$h$','Interpreter','latex','FontSize',24)
ylabel('relative $\ell^\infty$ error', 'Interpreter','latex','FontSize',24)
    
    

    
    