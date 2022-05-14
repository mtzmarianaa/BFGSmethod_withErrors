function [xk, k, XK] = lsBFGSwithErrors( f, g, x0, l, maxiter, tol, fexact, gexact )
% Purpose: approximate a local min of f using the linesearch algorithm
% and the (iBGFS) update formula (to avoid the solution of linear systems)
% in a modified version suggested by Xie, Byrd , and Nocedal in "Analysis
% of the BFGS method with errors". The lengthening parameter l is
% introduced.
%
% In :  f       : noisy objectve function (handle)
%       g       : noisy gradient from objective function (handle)
%       l       : lengthening parameter
%       x0      : starting point
%       maxiter : maximum number of iterations to perform
%       tol     : tolerance for the norm of the noisy gradient
%       fexact  : exact function
%       gexact  : exact gradient
%
% Out:  xk    : approximation to the minimum
%       k     : number of iterations done
%       XK    : extras

    k = 0;
    n = length( x0 );
    I = speye( n );
    H = speye( n );
    gk = g( x0 );
    gnew = gk;
    xk = x0;
    n_alpha0 = 0;

    % extras
    N = length(xk);
    xstar = zeros(N,1); %analytic minimum
    fk_noisy(1) = f(xk); % f(xk)
    fk(1) = fexact(xk);
    norm_g(1) = norm(gnew,2);
    abs_e(1) = norm(xstar - xk,2);
    alphas(1) = 1;
    lengthening(1) = 1;

    while norm( gk ) > tol
        dk = -H*gk;

        [alpha] = backtracking (f, g, xk, dk);


        if norm(alpha*dk) >= l
            s = alpha*dk;
            alp = alpha;
            lengthening_n = 0;
        else
            s = l*(dk/norm(dk));
            alp = 0;
            lengthening_n = 1;
        end

        if alpha ~= 0
            n_alpha0 = 0;
        end

        if alpha == 0
            n_alpha0 = n_alpha0 + 1;
        end

        if n_alpha0 >= 30
            break
        end
        
        gamma = g(xk + s) - g(xk);
        rho = 1/dot( s, gamma );
        
        H = (I - rho*(s*gamma'))*H*(I - rho*(gamma*s')) + rho*(s*s');
        gk = g(xk + alpha*dk);
        xk = xk + alpha*dk;
        k = k+1; 
        
        % extras
        alphas(k + 1) = alp;
        lengthening(k+1) = lengthening_n;
        fk_noisy(k + 1) = f(xk); % f(xk)
        fk( k+1 ) = fexact(xk);
        norm_g(k + 1) = norm(gexact(xk + alpha*dk),2);
        abs_e(k + 1) = norm(xstar - xk,2);
        
        if  k >= maxiter  % forces to stop if we've reached a lot of iterations
            break
        end
        
        
    end
%   extras
    XK = [(0:k)', norm_g', fk_noisy', fk', abs_e', alphas', lengthening'];
           
end