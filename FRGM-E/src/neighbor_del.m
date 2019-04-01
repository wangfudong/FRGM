function [SX2] = neighbor_del(X)
LX = length(X(:,1));
[ntriX] = delaunay(X(:,1),X(:,2));
SX2 = zeros(LX,LX);
for i = 1:length(ntriX(:,1))
    tri3 = sort(ntriX(i,:),'ascend');
    SX2(tri3(1),tri3(2)) = 1;
    SX2(tri3(1),tri3(3)) = 1;
    SX2(tri3(2),tri3(3)) = 1;
end
disX = sqrt(bsxfun(@minus,X(:,1),X(:,1)').^2 + bsxfun(@minus,X(:,2),X(:,2)').^2);
index = find(SX2>0);
disX = disX(index);
[AX,CX] = kmeans(disX,2);
if CX(1) > CX(2)
    AX = (AX == 1);
    AX = AX + 1;
end
index1 = index(AX == 1);
SX2 = zeros(LX,LX);
SX2(index1) = 1;
