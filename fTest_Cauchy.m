function[y] = fTest_Cauchy(x)
% Test function as in Xie, Byrd, and Nocedal's paper with error from a
% Cauchy distribution
%rng(16)
T = diag([10e-2, 1, 10e2, 10e4]);
error = trnd(1,1,1);
y = 0.5*x'*T*x + error;

end