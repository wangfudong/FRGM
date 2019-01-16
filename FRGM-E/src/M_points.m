function MP = M_points(X,Y,q)
% L_q distance matrix between X and Y
if nargin < 3 || isempty(q)
    q = 0.5;
end
MP =  (bsxfun(@minus,X(:,1),Y(:,1)').^(1/q) + bsxfun(@minus,X(:,2),Y(:,2)').^(1/q)).^q;
%MP = MP/max(MP(:));