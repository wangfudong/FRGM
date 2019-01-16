function ker_map = nonrigid_ker(X,sigma1,type)
% compute the kernel matrix
if nargin < 2 || isempty(sigma1)
    sigma1 = 1;
end
if nargin < 3 || isempty(type)
    type = 'rbf';
end
DX2 = (bsxfun(@minus,X(:,1),X(:,1)').^2 + bsxfun(@minus,X(:,2),X(:,2)').^2)/(2*sigma1.^2);

switch type 
    case 'rbf'
        ker_map = exp(-DX2);
    case 'tps'
        ker_map = DX2.*log(DX2.^0.5+eye(size(DX2)));
end

