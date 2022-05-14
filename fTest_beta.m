function[y] = fTest_beta(x)
% Test function as in Xie, Byrd, and Nocedal's paper with error from a
% beta distribution alpha = 0.5 = beta
%rng(16)
T = diag([10e-2, 1, 10e2, 10e4]);
error = betarnd(0.5,0.5);
y = 0.5*x'*T*x + error;

end