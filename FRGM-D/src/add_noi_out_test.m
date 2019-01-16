function [XT] = add_noi_out_test(X,opt,type)

switch type
    case 'affine'
%         opt.noise = 1;
%         opt.noise_sigma = noi;
%         opt.noise_type = 'uniform';%'guassian'
%         opt.outlier = out;
%         opt.out_sigma = 0.25;
%         opt.outlier_type = 'guassian';%'uniform1''guassian''guassian1'
        
        theta = opt.theta*pi;
        v = [cos(theta),sin(theta);
            -sin(theta),cos(theta)];
        t = [0,0];
        %                 s = [10*rand(1),0.1*rand(1)];
        %                 sss = max(s*(rand(2,1)>=0.5),0.1);
        s = 1;
        
        XT = rigid_affine_transform2D(X,v,t,s,opt);
        
    case 'nonrigid'
        sigma = opt.sigma;
        ker_type = opt.kernel;
        GX = nonrigid_ker(X,sigma,ker_type);
        theta = opt.theta*pi;
        R = [cos(theta),sin(theta);
            -sin(theta),cos(theta)];
        W = 0.5*randn(length(X(:,1)),2);
        %t = 1*rand(1,2) + 0.5;
        %t=[0,0];
        scale = 1;
        
%         opt.noise = 1;
%         opt.noise_sigma = 0.00;
%         opt.noise_type = 'uniform';%'guassian'
%         opt.outlier = 00;
%         opt.out_sigma = 0.25;
%         opt.outlier_type = 'guassian';%'uniform1''guassian''guassian1'
        
        
        XT = scale*nonrigid_kernel_trans(X,W,GX,R,opt);
        
end