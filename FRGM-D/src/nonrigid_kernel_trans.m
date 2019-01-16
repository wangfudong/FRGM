function [XT] = nonrigid_kernel_trans(X,W,ker,v,option)

if nargin < 4 || isempty(v)
    XT = X + ker*W;
else
    XT = X*v + ker*W;
end
if ~isempty(option)
    if option.noise > 0
        switch option.noise_type
            
            case 'uniform'
                sig_noise = option.noise_sigma;
                XT = XT + sig_noise*(rand(size(X))-0.5);
                
            case 'gaussian'
                sig_noise = option.noise_sigma;
                XT = XT + sig_noise*randn(size(X));
                
        end
        if option.outlier > 0
            switch option.outlier_type
                case 'uniform'
                    lx = length(X(:,1));
                    sig_out = option.out_sigma;
                    out = option.outlier;
                    index = floor(lx*rand(1,out))+1;
                    index = index(randperm(out));
                    X_out = XT(index,:) + sig_out*(rand(out,2)-0.5);
                    
                    XT = [XT;X_out];
                case 'gaussian'
                    lx = length(X(:,1));
                    sig_out = option.out_sigma;
                    out = option.outlier;
                    index = floor(lx*rand(1,out));
                    index = index(randperm(out))+1;
                    X_out = XT(index,:) + sig_out*randn(out,2);
                    XT = [XT;X_out];
                case 'uniform1'
                    % lx = length(X(:,1));
                    sig_out = option.out_sigma;
                    out = option.outlier;
                    XT = [XT;repmat(mean(XT),out,1)+sig_out*(rand(out,2)-0.5)];
                case 'gaussian1'
                    % lx = length(X(:,1));
                    sig_out = option.out_sigma;
                    out = option.outlier;
                    XT = [XT;repmat(mean(XT),out,1)+sig_out*randn(out,2)];
            end
        end
    end
    
    
    
end