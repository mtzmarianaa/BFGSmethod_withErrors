%% Generation of the plots of the PDF of the probabilities used

% For the continuous uniform distribution on [-1,1]

figure(1)
x11 = linspace(-1, 1, 200);
x01 = linspace(0, 1, 100);

subplot(2, 2, 1)
p_cdf1 = plot(x11, unifpdf(x11,-1,1));
p_cdf1.Color = [0, 236/255, 255/255];
title('Uniform distribution')
xlabel('x')
ylabel('PDF')

subplot(2, 2, 2)
p_cdf2 = plot(x01, betapdf(x01,0.5,0.5));
p_cdf2.Color = [0, 236/255, 255/255];
title('CDF beta alpha=beta=0.5')
xlabel('x')
ylabel('PDF')


subplot(2, 2, 3)
x = 0:1;
p_cdf4 = bar( x, binopdf(0:1,1,0.8) , 1 );
p_cdf4.FaceColor = [0, 236/255, 255/255];
title('Bernoulli distribution, p = 0.8')
xlabel('x')
ylabel('PDF')

subplot(2, 2, 4)
xplus = linspace(0,10,1000);
pdf_par = gppdf(xplus, 100,1,0.01);
p_cdf3 = plot(xplus, pdf_par);
p_cdf3.Color = [0, 236/255, 255/255];
title('Pareto distribution alpha = 0.01')
xlabel('x')
ylabel('PDF')

sgtitle('Probability distribution functions')



%% Optimization part

n_runs = 20;
fexact = @fTest;
gexact = @gTest;
g = @gTest_error;
x0 = 10e5*[1,1,1,1]';
tol = 10e-5;
l = 4*10e2;
maxiter = 60;

%% Uniform error

f_unif = @fTest_unif;



f2 = figure(2);
f2.Position = [10 10 1200 600];

subplot(1,2, 1)
[xk, k, XK] = lsBFGSwithErrors( f_unif, g, x0, l, maxiter, tol, fexact, gexact );
n_iter = size(XK);
diff_phi = [];
bound = [];
for i = 1:( n_iter(1) - 1)
    diff_phi(i) = XK(i+1, 4) - XK(i, 4);
    bound(i) = -3.125*10e-8 * XK(i+1, 2).^2;
end
pltBound = plot(1:(n_iter - 1), diff_phi);
pltBound.Color = [0, 236/255, 155/255];
hold on
pltBound2 = plot(1:(n_iter - 1), bound);
pltBound2.Color = [124/255 0 200/255];
title('Bound on (16)')
xlabel('Iteration')
ylabel('Bound on (16) ')
hold off
subplot(1,2, 2)

for i = 1:n_runs
    [xk, k, XK] = lsBFGSwithErrors( f_unif, g, x0, l, maxiter, tol, fexact, gexact );
    iter = 1:length(XK(:, 6));
    p_unif_er = plot(iter, XK(:, 6) );
    p_unif_er.Color = [1, 11/255, 172/255];
    hold on
end
title('Alpha found')
xlabel('Iteration')
ylabel('alpha')
hold off

sgtitle('Bound on (16) and choice of alpha')


f3 = figure(3);
f3.Position = [10 10 800 1700];

%%% With lsBFGSwithErrors
subplot(4,2, 1)

for i = 1:n_runs
    [xk, k, XK] = lsBFGSwithErrors( f_unif, g, x0, l, maxiter, tol, fexact, gexact );
    iter = 1:length(XK(:, 4));
    p_unif_er = plot(iter, log10(XK(:, 4) ) );
    p_unif_er.Color = [0, (i*5)/255, (i*10)/255];
    hold on
end
title('Decrease in exact function value, uniform')
xlabel('Iteration')
ylabel('Log10( phi(xk) - phi(x*) )')
hold off


subplot(4,2, 2)

for i = 1:n_runs
    [xk, k, XK] = lsBFGSwithErrors( f_unif, g, x0, l, maxiter, tol, fexact, gexact );
    iter = 1:length(XK(:, 2));
    p_unif_er = plot(iter, log10(XK(:, 2) ) );
    p_unif_er.Color = [(i*5)/255, 0, (i*10)/255];
    hold on
end
title('Decrease in norm of exact gradient, uniform')
xlabel('Iteration')
ylabel('Log10( norm(exact gradient)  )')
hold off




%% Beta distribution

f_beta = @fTest_beta;


subplot(4,2, 3)

for i = 1:n_runs
    [xk, k, XK] = lsBFGSwithErrors( f_beta, g, x0, l, maxiter, tol, fexact, gexact );
    iter = 1:length(XK(:, 4));
    p_unif_er = plot(iter, log10(XK(:, 4) ) );
    p_unif_er.Color = [0, (i*5)/255, (i*10)/255];
    hold on
end
title('Decrease in exact function value, beta')
xlabel('Iteration')
ylabel('Log10( phi(xk) - phi(x*) )')
hold off


subplot(4,2, 4)

for i = 1:n_runs
    [xk, k, XK] = lsBFGSwithErrors( f_beta, g, x0, l, maxiter, tol, fexact, gexact );
    iter = 1:length(XK(:, 2));
    p_unif_er = plot(iter, log10(XK(:, 2) ) );
    p_unif_er.Color = [(i*5)/255, 0, (i*10)/255];
    hold on
end
title('Decrease in norm of exact gradient, beta')
xlabel('Iteration')
ylabel('Log10( norm(exact gradient)  )')
hold off




%% Bernoulli distribution

f_Bernoulli = @fTest_Bernoulli;

subplot(4,2, 5)

for i = 1:n_runs
    [xk, k, XK] = lsBFGSwithErrors( f_Bernoulli, g, x0, l, maxiter, tol, fexact, gexact );
    iter = 1:length(XK(:, 4));
    p_unif_er = plot(iter, log10(XK(:, 4) ) );
    p_unif_er.Color = [0, (i*5)/255, (i*10)/255];
    hold on
end
title('Decrease in exact function value, Bernoulli')
xlabel('Iteration')
ylabel('Log10( phi(xk) - phi(x*) )')
hold off


subplot(4,2, 6)

for i = 1:n_runs
    [xk, k, XK] = lsBFGSwithErrors( f_Bernoulli, g, x0, l, maxiter, tol, fexact, gexact );
    iter = 1:length(XK(:, 2));
    p_unif_er = plot(iter, log10(XK(:, 2) ) );
    p_unif_er.Color = [(i*5)/255, 0, (i*10)/255];
    hold on
end
title('Decrease in norm of exact gradient, Bernoulli')
xlabel('Iteration')
ylabel('Log10( norm(exact gradient)  )')
hold off




%% Standard normal distribution

f_Normal = @fTest_Normal;

subplot(4,2, 7)

for i = 1:n_runs
    [xk, k, XK] = lsBFGSwithErrors( f_Normal, g, x0, l, maxiter, tol, fexact, gexact );
    iter = 1:length(XK(:, 4));
    p_unif_er = plot(iter, log10(XK(:, 4) ) );
    p_unif_er.Color = [0, (i*5)/255, (i*10)/255];
    hold on
end
title('Decrease in exact function value, Standard Normal')
xlabel('Iteration')
ylabel('Log10( phi(xk) - phi(x*) )')
hold off


subplot(4,2, 8)

for i = 1:n_runs
    [xk, k, XK] = lsBFGSwithErrors( f_Normal, g, x0, l, maxiter, tol, fexact, gexact );
    iter = 1:length(XK(:, 2));
    p_unif_er = plot(iter, log10(XK(:, 2) ) );
    p_unif_er.Color = [(i*5)/255, 0, (i*10)/255];
    hold on
end
title('Decrease in norm of exact gradient, Standard Normal')
xlabel('Iteration')
ylabel('Log10( norm(exact gradient)  )')
hold off

sgtitle('20 runs, different distributions')

%% Pareto distribution

n_runs = 20;
fexact = @fTest;
gexact = @gTest;
g = @gTest_error;
x0 = 10e5*[1,1,1,1]';
tol = 10e-5;
l = 4*10e2;
maxiter = 60;


f_Pareto = @fTest_Pareto;

figure(4)
subplot(2,2, 1)

for i = 1:n_runs
    [xk, k, XK] = lsBFGSwithErrors( f_Pareto, g, x0, l, maxiter, tol, fexact, gexact );
    iter = 1:length(XK(:, 4));
    p_unif_er = plot(iter, log10(XK(:, 4) ) );
    p_unif_er.Color = [0, (i*5)/255, (i*10)/255];
    hold on
end
title('Decrease in exact function value')
xlabel('Iteration')
ylabel('Log10( phi(xk) - phi(x*) )')
hold off


subplot(2,2, 2)

for i = 1:n_runs
    [xk, k, XK] = lsBFGSwithErrors( f_Pareto, g, x0, l, maxiter, tol, fexact, gexact );
    iter = 1:length(XK(:, 2));
    p_unif_er = plot(iter, log10(XK(:, 2) ) );
    p_unif_er.Color = [(i*5)/255, 0, (i*10)/255];
    hold on
end
title('Decrease in norm of exact gradient')
xlabel('Iteration')
ylabel('Log10( norm(exact gradient)  )')
hold off

subplot(2,2, 3)
[xk, k, XK] = lsBFGSwithErrors( f_Pareto, g, x0, l, maxiter, tol, fexact, gexact );
n_iter = size(XK);
diff_phi = [];
bound = [];
for i = 1:( n_iter(1) - 1)
    diff_phi(i) = XK(i+1, 4) - XK(i, 4);
    bound(i) = -3.125*10e-8 * XK(i+1, 2).^2;
end
pltBound = plot(1:(n_iter - 1), diff_phi);
pltBound.Color = [0 200/255 200/255];
hold on
pltBound2 = plot(1:(n_iter - 1), bound);
pltBound2.Color = [124/255 0 200/255];
title('Bound on (16)')
xlabel('Iteration')
ylabel('Bound on (16) ')
hold off
subplot(2,2, 4)

for i = 1:n_runs
    [xk, k, XK] = lsBFGSwithErrors( f_Pareto, g, x0, l, maxiter, tol, fexact, gexact );
    iter = 1:length(XK(:, 6));
    p_unif_er = plot(iter, XK(:, 6)  );
    p_unif_er.Color = [1, 11/255, 172/255];
    hold on
end
title('Alpha found')
xlabel('Iteration')
ylabel('alpha')
hold off

sgtitle('Errors from Pareto (alpha = 0.01) distribution')


