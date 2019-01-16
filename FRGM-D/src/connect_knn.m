function SX = connect_knn(X,knei)
% calculate the adjacency matrix of X via knn

LX = length(X(:,1));
[Ind,D] = knnsearch(X,X,'K',knei+1);
Ind = Ind(:,2:end);
D = D(:,2:end);

SX = zeros(LX,LX);
for i=1:LX
    SX(i,Ind(i,:)) = D(i,:);
end
SX = SX + SX';