function[y] = fTest_Normal(x)
% Test function as in Xie, Byrd, and Nocedal's paper with error from a
% standard normal distribution
%rng(16)
T = diag([10e-2, 1, 10e2, 10e4]);
error = normrnd(0,1);
y = 0.5*x'*T*x + error;

end