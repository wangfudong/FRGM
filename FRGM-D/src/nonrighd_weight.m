function [W] = nonrighd_weight(X,Y,P,S,kerX,lambda1,lambda2,lambda3)
% calculate the weight of nonrigid deformation
LS = diag(sum(S,2))-S;
XX = X;
PYY = P*Y;
W = ((XX-PYY)'*(lambda1*eye(size(S)) + lambda2*LS))/(kerX'*(lambda1*eye(size(S)) + lambda2*LS) + lambda3*eye(size(S)));
W = W';
