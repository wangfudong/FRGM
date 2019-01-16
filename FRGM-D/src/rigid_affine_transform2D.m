function XT = rigid_affine_transform2D(X,v,t,s,option)
% transformed graph XT = s*X*v+t + noise or outliers

[m,n] = size(v);
if m~=n
    error('v should be a d*d square matrix');
end
if size(X,2)~=m
    X = X';
end
if size(t,1)~=1 || size(t,2) ~= m
    error('t should be a 1*d row vector');
end

XT = s*X*v + repmat(t,length(X(:,1)),1);

if nargin >= 5 && ~isempty(option)
    if option.noise > 0
        switch option.noise_type
            
            case 'uniform'
                sig_noise = option.noise_sigma;
                XT = XT + s*sig_noise*(rand(size(X))-0.5);
                
            case 'gaussian'
                sig_noise = option.noise_sigma;
                XT = XT +s* sig_noise*randn(size(X));
                
        end
        if option.outlier > 0
            switch option.outlier_type
                case 'uniform'
                    lx = length(X(:,1));
                    sig_out = option.out_sigma;
                    out = option.outlier;
                    index = floor(lx*rand(1,out))+1;
                    index = index(randperm(out));
                    X_out = XT(index,:) + s*sig_out*(rand(out,2)-0.5);
                  
                    XT = [XT;X_out];
                case 'gaussian'
                    lx = length(X(:,1));
                    sig_out = option.out_sigma;
                    out = option.outlier;
                    index = floor(lx*rand(1,out));
                    index = index(randperm(out))+1;
                    X_out = XT(index,:) + s*sig_out*randn(out,2);
                    XT = [XT;X_out];
                case 'uniform1'
                   % lx = length(X(:,1));
                    sig_out = option.out_sigma;
                    out = option.outlier;
                    XT = [XT;repmat(mean(XT),out,1)+s*sig_out*(rand(out,2)-0.5)];
                case 'gaussian1'
                   % lx = length(X(:,1));
                    sig_out = option.out_sigma;
                    out = option.outlier;
                    XT = [XT;repmat(mean(XT),out,1)+s*sig_out*randn(out,2)];
            end
        end
    end
end




