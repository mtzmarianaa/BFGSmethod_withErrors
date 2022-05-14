function[y] = fTest(x)
% Test function as in Xie, Byrd, and Nocedal's paper

T = diag([10e-2, 1, 10e2, 10e4]);

y = 0.5*x'*T*x;

end