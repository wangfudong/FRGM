function rate = rateremove(order,rest)
% compute the rate of inliers
LX = length(order);
rate = 0;
for i = 1:length(rest)
    rate = rate + ismember(rest(i),order);
end
rate = rate/LX;