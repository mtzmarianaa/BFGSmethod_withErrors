function[y] = fTest_Bernoulli(x)
% Test function as in Xie, Byrd, and Nocedal's paper with error from a
% Bernoulli distribution, p = 0.8
%rng(16)
T = diag([10e-2, 1, 10e2, 10e4]);
pd = makedist('Binomial','N',1,'p',0.8);
error = random(pd);
y = 0.5*x'*T*x + error;

end