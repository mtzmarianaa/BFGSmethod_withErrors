function [alpha] = backtracking (f, g, xk, dk)
alpha = 1; % start at 1
c1 = 0.01; % condition set in HW
c2 = 0.5; % condition set in HW
ft = f(xk + alpha*dk); % we need this value but not this gradient
fk = f(xk); % we need this value and the gradient
gk = g(xk);
n_iter = 1;
% backtracking part
while( ft > fk + c1*alpha*gk'*dk && n_iter < 64 )
    alpha = c2*alpha;
    ft = f(xk + alpha*dk);
    n_iter = n_iter + 1;
end

if n_iter >= 64
    alpha = 0;
end
gt = g(xk + alpha*dk);
if (ft > fk + c1*alpha*gk'*dk || gt'*dk < c2*gk'*dk)
    alpha = 0;
end