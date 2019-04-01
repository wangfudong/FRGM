function [MS,BH1,BH2] = M_shape(X,Y,r0,r1,rota_inv)
% compute the dissimilarity between X and Y using Shape Context
LX = length(X(:,1));
LY = length(Y(:,1));

if nargin > 4 && rota_inv > 0% using rotation invariant Shape Context
    [BH1] = sc_compute_rota(X',zeros(1,LX),[],12,5,r0,r1,zeros(1,LX));%BH1(i,:)' is formed form a 5*12 matrix.
    [BH2] = sc_compute_rota(Y',zeros(1,LY),[],12,5,r0,r1,zeros(1,LY));
else
    [BH1] = sc_compute(X',zeros(1,LX),[],12,5,r0,r1,zeros(1,LX));%BH1(i,:)' is formed form a 5*12 matrix.
    [BH2] = sc_compute(Y',zeros(1,LY),[],12,5,r0,r1,zeros(1,LY));
end
BH1 = full(BH1);
BH2 = full(BH2);

%%% Euclidean distance
%     MS = zeros(LX,LY);
%     for i = 1:60
%         MS = MS + abs(bsxfun(@minus,BH1(:,i),BH2(:,i)')).^2;
%     end
%     MS = MS.^(1/2);
%     MS = MS/max(MS(:));

%%% KL divergence
MS = hist_cost_2(BH1,BH2);
MS = MS/max(MS(:));


