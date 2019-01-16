%% experiments on pascal dataset without FGM

remove = 1;
removal_time = 10;
dataset = ['Carss';'Motor'];
dataL = [30,20];

datasets = 3;%3=Carss,4=Motors
inds = 18;%index of graph pairs

[Pts,Ims,nF] = fileload(datasets,inds);
X = Pts{1,1};
Y = Pts{1,2};
I1 = Ims{1,1};
I2 = Ims{1,2};

outliers = 25;
LX = nF;
LY_all = length(Y(:,1));
LY = min(LX + outliers,LY_all);
order = randperm(LX);
%order = 1:LX;

max_size1 = max(size(rgb2gray(I1)));
max_size2 = max(size(rgb2gray(I2)));
max_size = max(max_size1,max_size2);

XX = X(order,:)/max_size;
if outliers > 0
    F2_rest = Y(LX+1:end,:);
    out_inds = randperm(LY_all-LX);
    out_inds = out_inds(1:LY-LX);
    Y_out = F2_rest(out_inds,:);% add outliers ramdonly
    %Y_out = F2_rest(1:LY-LX,:);
    YY = [Y(1:LX,:);Y_out]/max_size;
else
    YY = Y(1:LY,:)/max_size;
end

if outliers > 0 && remove > 0
    [YY,Y_rest,removerate, rest_out] = car_remove(XX,YY,order,removal_time,[]);    
    LY = length(Y_rest);
else
    Y_rest = 1:LY;
end
removerate

%% ATGM
option.Mexist = 0.1;
option.maxiter_nonvex = 200;
option.maxiter_convex = 100;
option.lam_nonvex = 10;
option.lam_sparse = 0*1.0e2;
option.lam_convex = 0.1;
option.geofunc = '1.21';
option.unary = 1;
option.full = 1;



asgATGM = ATGM_testing(XX,YY,'del',order,option);
asgATGM.acc = ratemap(asgHun(asgATGM.X),order,Y_rest,0);
