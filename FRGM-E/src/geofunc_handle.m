function [Fgeo,grad_geo] = geofunc_handle(X,Y,SX,option)
% handle for computing the gradient and function value of the
% objective function
% X,Y: points of two graphs
% SX: neighborhood matrix of GX and GY

geofunc = option.geofunc;
LengthX = length(X(:,1));
%SX = double(( SX >0));
sma_eps = 1.0e-10;

Norm2 = @(loc)(bsxfun(@minus,loc(:,1),loc(:,1)').^2 + bsxfun(@minus,loc(:,2),loc(:,2)').^2);
switch geofunc
    
    case '1.21'
        Fgeo = @(P)(sum(sum(SX.*(Norm2(X).^0.5 - Norm2(P*Y).^0.5).^2)));
        LSX = @(P)(SX.*(Norm2(P*Y).^0.5 - Norm2(X).^0.5).*max(Norm2(P*Y),sma_eps).^(-0.5));
        LSX = @(P)(diag(LSX(P)*ones(LengthX,1))-LSX(P));
        grad_geo = @(P)(2*(LSX(P)*P*Y)*Y');
        
    case '0.11'
        Fgeo = @(P)(sum(sum(SX.*Norm2(X-P*Y))));
        LSX = diag(sum(SX,2))-SX;
        grad_geo = @(P)(2*(LSX*(P*Y-X))*Y');
        
end

