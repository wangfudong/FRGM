function [F,gF] = func_handle(M,P,SX,DX,DY,option)
% handle for computing the gradient and function value of the
% objective function
M_exist = option.M_exist;
alpha2 = option.alpha2;
LX = length(P(:,1));
LY = length(P(1,:));

q = option.q_norm;

F = sum(sum(P.*M*M_exist)) + alpha2*sum(sum((SX.*(DX-P.^q*DY*(P').^q).^2)));
if option.gra > 0
    gF = zeros(size(P));
    if q >= 1
        for t = 1:LX
            gF(t,:) = 4*sum(repmat(SX(t,:).*(DX(t,:)-P(t,:).^q*DY*(P').^q),LY,1).*...
                repmat(-q*(P(t,:).^(q-1))',1,LX).*(DY*(P').^q),2)';
        end
    else
        for t = 1:LX
            gF(t,:) = 4*sum(repmat(SX(t,:).*(DX(t,:)-P(t,:).^q*DY*(P').^q),LY,1).*...
                repmat(-q*(max(P(t,:),1.0e-10).^(q-1))',1,LX).*(DY*(P').^q),2)';
        end
    end
    
    gF = alpha2*gF + M*M_exist;
else
    gF = zeros(size(P));
end
