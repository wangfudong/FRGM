function [SX] = connect_del(X)
% calculate the adjacency matrix of X via delaunay triangulation
LX = length(X(:,1));
[ntriX] = delaunay(X(:,1),X(:,2));
SX = zeros(LX,LX);
for i = 1:length(ntriX(:,1))
    tri3 = sort(ntriX(i,:),'ascend');
    SX(tri3(1),tri3(2)) = 1;
    SX(tri3(1),tri3(3)) = 1;
    SX(tri3(2),tri3(3)) = 1;
end
SX = SX + SX';