function[g] = gTest_error(x)
% Test function's gradient as in Xie, Byrd, and Nocedal's paper with error
% sampled uniformly from the closed unit ball
%rng(16)
T = diag([10e-2, 1, 10e2, 10e4]);
error_g = unifrnd(-1,1, [4, 1]);
error_g = error_g./norm(error_g);
R = nthroot(unifrnd(0,1),4);
error_g = error_g.*R;
g = T*x + error_g;

end