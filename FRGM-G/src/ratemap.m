function [rate,matches] = ratemap(Map,order,rest,play)
Map = asgHun(Map);
matches = zeros(1,length(order));

for i = 1:length(order)
    if max(Map(i,:))>= 0.1
        in = find(Map(i,:) == max(Map(i,:)));
        matches(i) = in(1);
    end
end
if isempty(rest)
    rate = sum(abs(matches - order)==0)/sum(matches>0);
else
    
    rest_index = rest(matches);
    rate = sum(abs(rest_index - order)==0)/sum(matches>0);
    
end
if play > 0
    figure,imagesc(Map);hold on,
    for i = 1:length(order)
        ind = find(rest == order(i));
        if isempty(ind) == 0
            plot(ind,i,'r.','markersize',15);hold on,
        end
    end
end
