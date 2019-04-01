inlier = 100;
outlier = 100;
sig = 0.02;

X = randn(inlier,2);
LX = inlier;
LY = inlier + outlier;
Yout = randn(outlier,2);
Y = [X+ sig*randn(inlier,2);Yout];
max_size = max(max(abs(X(:))),max(abs(Y(:))));

order = randperm(LX);
X = X(order,:);

XX = X/max_size;
YY = Y/max_size;

%% ATGM
option.Mexist = 0.01;
option.maxiter_nonvex = 200;
option.maxiter_convex = 100;
option.lam_nonvex = 1-0.01;
option.lam_sparse = 10;
option.lam_convex = 0.01;
option.geofunc = '1.21';
option.unary = 1;
option.full = 1;

asgATGM = FRGM_E(XX,YY,order,option);

