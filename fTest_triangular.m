function[y] = fTest_triangular(x)
% Test function as in Xie, Byrd, and Nocedal's paper with error from a
% triangular distribution on [-1,1] with peak at c = 0.5
%rng(16)
T = diag([10e-2, 1, 10e2, 10e4]);
pd = makedist('Triangular','A',-1,'B',0.5,'C',1);
error = random(pd);
y = 0.5*x'*T*x + error;

end