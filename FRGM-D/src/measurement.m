function [mse,acc] = measurement(X,Y,cor_GT,cor_re)
% compute the average error between the transformed nodes and groundtruth
[~,max_bound] = normalize_point(X,1);
X = X/max_bound;
Y = Y/max_bound;
mse = mean(sqrt(sum((X - Y(cor_GT,:)).^2,2)),1);
if nargin > 3 && ~isempty(cor_re)
    acc = sum((cor_GT == cor_re))/length(cor_GT);
else
    acc = [];
end
