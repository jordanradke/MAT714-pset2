%% interpolation of f(x)
 
Nq  = 1000000;   
hq  = 1/Nq;
xq  = 0:hq:1; % query grid

fq  = exp(-400*(xq-0.5).^2); % the exact function, numerically

% initialize error vector. run loop until you find N that gets max_err 
% under required tolerance
err = 10^6;
max_err = 10^(-2);
errors = zeros(1);
N = 0;

while err > max_err
   
    N = N+1;
    h   = 1/N;
    x   = 0:h:1;  % spatial grid

    f   = exp(-400*(x-0.5).^2);  % the sampled function
    fint = interp1(x,f,xq);      % f interpolated from the sample function
    
    err = max(abs(fq-fint));  % compute l_infty norm of error vector
    errors(N) = err;
    
end

plot(1:N, errors)
