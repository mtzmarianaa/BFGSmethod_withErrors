function[y] = fTest_Pareto(x)
% Test function as in Xie, Byrd, and Nocedal's paper with error from a
% Pareto distribution with alpha = 0.5
%rng(16)
T = diag([10e-2, 1, 10e2, 10e4]);
error = gprnd(100,1,0.01);
y = 0.5*x'*T*x + error;

end