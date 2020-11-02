%% analytic solution of wave eq. via fourier on [0,1]x[0,1]
function [w] = wave_fourier(X,Y,N,T,dt,g)

dx = 1/N;
dy = 1/N;

% fourier coefficients of C^infty function decay O(nmodes^{-p}) for any
% p we like.
nmodes = [10,10];
b = zeros(nmodes);

% numerically finding fourier coefficients
for m = 1:nmodes
    for n = 1:nmodes
        xfun = @(x) exp(-4*(x-.5).^2).*sin(m*x.*pi);
        yfun = @(y) exp(-4*(y-.5).^2).*sin(n*y.*pi);
        xint = integral(xfun,0,1);
        yint = integral(yfun,0,1);
        b(n,m) = 4/(pi*sqrt(n^2+m^2))*xint*yint;
    end
end

% series solution for u
u = zeros([size(X),T]);

u(:,:,1) = g;

% is there a way to write in this via matrix operations?
for t = 1:T
    for i = 1:N+1
        for j = 1:N+1
            sum = 0;
            for m = 1:nmodes
                for n = 1:nmodes
                    sum = sum + b(m,n)*sin(pi*sqrt(m^2+n^2)*(dt*(t-1)))...
                        *sin(m*pi*(dx*(i-1)))*sin(n*pi*(dy*(j-1)));
                end
            end
            u(i,j,t) = sum;
        end
    end
end

w = u;

end

                    
 
