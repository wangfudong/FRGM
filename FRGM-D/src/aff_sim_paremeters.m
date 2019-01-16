function [V,t,s] = aff_sim_paremeters(X,Y,P,S,lambda1,lambda2,type)
% X = s*(P*Y)*R + t or X = (P*Y)*V + t
%Input:
%     X,Y: two point sets
%     P: m*n correspondence between X and Y
%     S: m*m weighted adjacency matrix
%     lambda1,lambda2: weight of the unary or pairwise term
%     type: geometric deformation type
%Output:
%     V: d*d transformation matrix
%     t: translation vector
%     s:scale

m = length(X(:,1));
n = length(Y(:,1));
if size(P,1)~=m || size(P,2)~=n
    error('check the input');
end
if size(P,1)~=m || size(P,2)~=n
    error('check the input');
end
if lambda1 < 0 || lambda2 < 0
    error('lambda1 and lambda2 should be >= 0');
end

switch type
    case 'affine2d'
        ux = sum(X)/m;
        PY = P*Y;
        upy = sum(PY)/m;
        XX = X-repmat(ux,m,1);
        PYY = PY - repmat(upy,m,1);
        LS = diag(sum(S,2))-S;
        V = (lambda1*XX'*PYY + lambda2*X'*LS*PYY)/(lambda1*(PYY'*PYY) + lambda2*PYY'*LS*PYY);
        V = V';
        t = ux - upy*V;
        
        s = [];
    case 'rigid2d'
        ux = sum(X)/m;
        PY = P*Y;
        upy = sum(PY)/m;
        XX = X-repmat(ux,m,1);
        PYY = PY - repmat(upy,m,1);
        LS = diag(sum(S,2))-S;
        
        AA = lambda1*XX'*PYY + lambda2*XX'*LS*PYY;
        [U,~,VV] = svd(AA');
        dim = size(X,2);
        
        C = eye(dim,dim);
        C(end,end) = det(U*VV');
        
        V = U*C*VV';
        s = trace(lambda1*XX'*PYY*V + lambda2*XX'*LS*PYY*V)/trace(PYY'*(lambda1*eye(size(LS)) + lambda2*LS)*PYY);
        t = ux - s*upy*V;
end


