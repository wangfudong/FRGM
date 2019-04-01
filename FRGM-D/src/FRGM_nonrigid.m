function [ASSIGN,para] = FRGM_nonrigid(X,Y,GX,option,order)
% calculate the nonrigid transformation parameters and correspondence between two graphs
%Input:
%    X,Y: point sets
%    GX: raw RBF kernel
%    option:options of the algorithm
%    order: order of points in X
%    
%Output:
%    ASSIGN:the correspondence maps
%    para: the calculated parameters of transformation

[m,dim1] = size(X);
[n,~] = size(Y);

if nargin < 5 || isempty(order)
    order = 1:m;
end

if option.regist_normalize > 0
    [X,~] = normalize_point(X,1);
    [Y,~] = normalize_point(Y,1);
end
LX = length(X(:,1));
LY = length(Y(:,1));

%-------------------initialization-----------%
cnt = 1;
regist_tol_diff = 2;
regist_tol_old = 3;
regist_tol_out = zeros(1,1,1);
regist_err_out = zeros(1,1,1);
max_regist = option.regist_it;

ASSIGN = zeros(m,n,1);
W = zeros(m,dim1,1);
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
if option.regist_display > 0
    filename = '.\save\nonrigid.gif';
    hh = figure;subplot(1,2,1);
    set(hh,'position',[500 450 1000 300]);
    plot(X(:,1),X(:,2),'r.',Y(:,1),Y(:,2),'b.','markersize',20);axis off;
    drawnow;
    subplot(1,2,2);
    plot(X(:,1),X(:,2),'r.',Y(:,1),Y(:,2),'b.','markersize',20);axis off;
    drawnow;
    if option.regist_save > 0
        frame = getframe(hh);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','overwrite', 'Loopcount',inf);
    end
end

%-----------------alternative calculation of ASSIGN and para-------------% 
while cnt <= max_regist && regist_tol_diff >= 1.0e-6
    
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
            %SX = 1./(DXX+eye(size(DXX)))-eye(size(DXX));
        case 'del'
            SX = connect_del(X);
            SX = SX.*exp(-DXX/1)-eye(size(SX));
        case 'knn'
            knei = option.GM_knn;
            SX = connect_knn(X,knei);
            SX = SX.*exp(-DXX/1)-eye(size(SX));
    end
    %%%SX = (SX>0);
    option.convex = alter(cnt);
    
    %-----------calculate the correspondence between two graphs---------%
    [assign] = FW_Min(X,Y,M,SX,Map_ini,option,1);
    ASSIGN(:,:,cnt) = assign;
    asg_assign = asgHun(assign);
    
    %----------calculate the geometry parameters------------------%
    if cnt <= 5
        sigma2 = max(sum(sum(assign.*(bsxfun(@minus,X(:,1),Y(:,1)').^2 + bsxfun(@minus,X(:,2),Y(:,2)').^2)))/(LX*2*LY),1.0e-10);
        [W1] = nonrighd_weight(X,Y,assign,SX,GX,1,0.1,sigma2);% W1*GX + X_maped = X1£»
    else
        sigma2 = max(sum(sum(assign.*(bsxfun(@minus,X(:,1),Y(:,1)').^2 + bsxfun(@minus,X(:,2),Y(:,2)').^2)))/(LX*2*LY),1.0e-10);
        [W1] = nonrighd_weight(X,Y,asg_assign,SX,GX,1,0.1,sigma2);% W1*GX + X_maped = X1£»
    end
    W(:,:,cnt) = -W1;
    
    Xtmp = X;
    X = nonrigid_kernel_trans(X,-W1,GX,eye(2,2),[]);% XX = X + W1*GX% updating X by geometry para
    
%     XX = asg_assign*Y;X_i is matched to (asg_assign*Y)_i
%     regist_err = measurement(XX,Y,order,[]);

    regist_err = measurement(X,Y,order,[]);
    regist_tol = measurement(X,Xtmp,1:LX,[]);
    regist_tol_diff = abs(regist_tol_old-regist_tol);
    regist_tol_old = regist_tol;
    regist_tol_out(cnt) = regist_tol;
    regist_err_out(cnt) = regist_err;
    
    disp(['regist_err at iter ' num2str(cnt) ': ' num2str(regist_err)]);
    disp(['regist_tol_diff at iter ' num2str(cnt) ': ' num2str(regist_tol_diff)]);
    
    cnt = cnt +1;
    
    if option.regist_display > 0
        %X = XX;
        plot(X(:,1),X(:,2),'ro','linewidth',1.5,'markersize',6);axis off;
        hold on,plot(Y(:,1),Y(:,2),'b.','markersize',10);
        hold off;
        drawnow;
        if option.regist_save > 0
            frame = getframe(hh);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
        end
    end
end

%-------------------output---------------%
para.W = W;
para.regist_err = regist_err_out;
para.regist_tol = regist_tol_out;
para.X = X;
para.Y = Y;





