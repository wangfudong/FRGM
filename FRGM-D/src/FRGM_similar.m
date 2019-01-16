function [ASSIGN,para] = FRGM_similar(X,Y,option,order)
% calculate the similarity transformation parameters and correspondence between two graphs
%Input:
%    X,Y: point sets
%    option:options of the algorithm
%    order: order of points in X
%    
%Output:
%    ASSIGN:the correspondence maps
%    para: the calculated parameters of transformation


[m,dim1] = size(X);
[n,~] = size(Y);

if nargin < 4 || isempty(order)
    order = 1:m;
end

%-------------------initialization-----------%
if option.regist_normalize > 0
    [X,sx0] = normalize_point(X,1);
    [Y,sy0] = normalize_point(Y,1);
else
    sx0 = 1;
    sy0 = 1;
end

cnt = 1;
regist_tol_diff = 2;
regist_tol_old = 3;
regist_tol_out = zeros(1,1,1);
regist_err_out = zeros(1,1,1);
max_regist = option.regist_it;

ASSIGN = zeros(m,n,1);
R = zeros(dim1,dim1,1);
t = zeros(1,2,1);
s = zeros(1,1,1);
alter = option.GM_convex_or_non;

option.tol_cg = option.GM_tol;
option.Mexist = option.GM_lambda1;
option.maxiter_nonvex = option.GM_nonvex_iter;
option.maxiter_convex = option.GM_convex_iter;
option.lam_nonvex = option.GM_lam_nonvex;
option.lam_convex = option.GM_lam_convex;
option.lam_sparse = option.GM_lam_sparse;
option.geofunc = option.GM_geofunc;
option.unary = option.GM_unary;

%---------------------display and write-----------%
LX = m;
if option.regist_display > 0
    filename = '.\save\similar.gif';
    hh = figure;subplot(1,2,1);
    set(hh,'position',[500, 450 ,1000 ,300]);
    plot(X(:,1),X(:,2),'r.',Y(:,1),Y(:,2),'b.','markersize',20);axis off;
    drawnow;
    subplot(1,2,2);
    plot(X(:,1),X(:,2),'r.',Y(:,1),Y(:,2),'b.','markersize',20);axis off;
    drawnow;
    frame = getframe(hh);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','overwrite', 'Loopcount',inf);
end

%-----------------alternative calculation of ASSIGN and para-------------% 
while cnt <= max_regist && regist_tol_diff >= 1.0e-7
    
    if alter(cnt) == 0
        rota = option.regist_rota;
        M = M_shape(X,Y,1/8,2,rota);
        option.GM_convex = 0;
        option.GM_geofunc = '1.21';
    else
        M = M_points(X,Y);
        option.GM_convex = 1;
        option.GM_geofunc = '0.11';
    end
    switch lower(option.GM_initial)
        case 'sinkhorn'
            lambda = 1/50;
            Map_ini = sinkhorn_OT(ones(m,1),ones(n,1)*(m/n),M,lambda,1.0e-6,1000,0);
        case 'lap'
            Map_ini = asgHun(-M);
        case 'rand'
            Map_ini = double_stochastic(rand(m,n),1000);
    end
    DXX = M_points(X,X);
    switch lower(option.GM_connected)
        case 'full'
            SX = exp(-DXX/1)-eye(size(DXX));
        case 'del'
            SX = connect_del(X);
            SX = SX.*exp(-DXX/1)-eye(size(SX));
        case 'knn'
            knei = option.GM_knn;
            SX = connect_knn(X,knei);
            SX = SX.*exp(-DXX/1)-eye(size(SX));
    end
    option.convex = alter(cnt);
    
    %-----------calculate the correspondence between two graphs---------%
    [assign] = FW_Min(X,Y,M,SX,Map_ini,option,1);
    ASSIGN(:,:,cnt) = assign;
    
    %----------calculate the geometry parameters------------------%
    if cnt <= 5
        [R1,t1,s1] = aff_sim_paremeters(X,Y,assign,SX,1,1,'rigid2d');
    else
        [R1,t1,s1] = aff_sim_paremeters(X,Y,asgHun(assign),SX,1,1,'rigid2d');
    end
    
    R(:,:,cnt) = eye(2,2)/R1;
    t(:,:,cnt) = -t1*(eye(2,2)/R1)/s1;
    s(:,:,cnt) = 1/s1;
    
    Xtmp = X;
    X = rigid_affine_transform2D(X,R(:,:,cnt),t(:,:,cnt),s(:,:,cnt),[]);% updating X by geometry para
    
    %     XX = asg_assign*Y;% 
    %     regist_err = measurement(XX,Y,order,[]);
    
    regist_err = measurement(X,Y,order,[]);
    
    regist_tol = measurement(X,Xtmp,1:LX,[]);
    regist_tol_diff = abs(regist_tol_old-regist_tol);
    regist_tol_old = regist_tol;
    regist_tol_out(cnt) = regist_tol;
    regist_err_out(cnt) = regist_err;
    
    disp(['regist_err at iter ' num2str(cnt) ': ' num2str(regist_err)]);
    disp(['regist_tol_diff at iter ' num2str(cnt) ': ' num2str(regist_tol_diff)]);
    
    cnt = cnt + 1;
    
    if option.regist_display > 0
        plot(X(:,1),X(:,2),'ro','linewidth',1.5,'markersize',6);axis off;
        hold on,plot(Y(:,1),Y(:,2),'b.','markersize',10);
        hold off;
        drawnow;
        frame = getframe(hh);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
end

%-------------------output---------------%
s(:,:,end) = s(:,:,end)*(sy0/sx0);
para.R = R;
para.t = t;
para.s = s;
para.regist_err = regist_err_out;
para.regist_tol = regist_tol_out;
para.X = X;
para.Y = Y;



