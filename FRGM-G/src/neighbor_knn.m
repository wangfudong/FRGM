function SX = neighbor_knn(D,knei)
m = length(D(1,:));
knei = floor(knei);
knei = min(max(knei,1),m-1);

SX = zeros(size(D));

for i = 1:m
    
    [val,~] = sort(D(i,:));
    basek = val(knei+1);
    ind = (D(i,:) <= basek);
    SX(i,ind) = 1;   
    
end
SX = SX + SX';
SX = (SX>0)-eye(size(D));

    