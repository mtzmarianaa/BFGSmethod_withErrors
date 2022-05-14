function[g] = gTest(x)
% Test function's gradient as in Xie, Byrd, and Nocedal's paper

T = diag([10e-2, 1, 10e2, 10e4]);

g = T*x;

end