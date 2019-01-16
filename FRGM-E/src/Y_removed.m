function [Y_new,rest] = Y_removed(X,Y_new,Y,Map)
% remove the outliers in Y
LX = length(X(:,1));
LY = length(Y(:,1));
X_maped = Map*Y_new;
dis = sqrt(bsxfun(@minus,X_maped(:,1),Y(:,1)').^2 + bsxfun(@minus,X_maped(:,2),Y(:,2)').^2);
rest = zeros(1,LY);
for i = 1:LX
    [d,ind] = sort(dis(i,:));
    if d(2) > 3*d(1) 
        rest(ind(1)) = 1;
    else 
        rest(ind(1:2)) = 1;
    end
end

if sum(rest) < LX
    tmp = LX - sum(rest);
    index = find(rest == 0);    
    value = zeros(length(index),1);
    for i = 1:length(index)
        mindis = min(dis(:,index));
        value(i) = mindis(1);
    end
    [~,ind] = sort(value);
    ind = ind(1:tmp);
    ind = index(ind);
    index = find(rest == 1);
    rest = [index,ind];
else
    rest = find(rest == 1);
end

Y_new = Y(rest,:);




