function [xk, k, XK] = lsBFGS( f, g, x0, maxiter, tol )
% Purpose: approximate a local min of f using the linesearch algorithm
% and the (iBGFS) update formula (to avoid the solution of linear systems)
% in a modified version suggested by Xie, Byrd , and Nocedal in "Analysis
% of the BFGS method with errors". The lengthening parameter l is
% introduced.
%
% In :  f       : noisy objectve function (handle)
%       g       : noisy gradient from objective function (handle)
%       x0      : starting point
%       maxiter : maximum number of iterations to perform
%       tol     : tolerance for the norm of the noisy gradient
%
% Out:  xk    : approximation to the minimum
%       k     : number of iterations done
%       XK    : extras 

    k = 0;
    n = length( x0 );
    I = speye( n );
    H = speye( n );
    gk = g(x0);
    gnew = gk;
    xk = x0;

    % extras
    xstar = zeros(4,1); %analytic minimum
    fk(1) = f(xk); % f(xk)
    norm_g(1) = norm(gnew,2);
    abs_e(1) = norm(xstar - xk,2);
    alphas(1) = 1;


    while norm( gnew ) > tol
        dk = -H*gk;
        
        [alpha] = backtracking (f, g, xk, dk);
        
        s = alpha*dk;
        xk = xk + s;
        gnew = g(xk + s);
        gamma = gnew - gk;
        rhoinv = dot( s, gamma );
        
        H = ( I - s*gamma'/rhoinv )*H*( I - gamma*s'/rhoinv ) + ( s*s' )/rhoinv;  % O(n^3)
        gk = gnew;
        k = k+1;   % STOP 
        
        % extras
        fk(k + 1) = f(xk); % f(xk)
        norm_g(k + 1) = norm(gnew,2);
        abs_e(k + 1) = norm(xstar - xk,2);
        alphas( k+1 ) = alpha;
        
        if  k >= maxiter
            break
        end
        

        
        
    end
%   extras
    XK = [(0:k)', norm_g', fk', abs_e', alphas'];
           
end