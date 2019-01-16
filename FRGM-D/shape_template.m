%% plot 5shape
h=figure;
xshift = [0,0.85,2.1,3,3.9];
yshift = [0.12,0,0,-0.1,0];
for i=1:5
data_name = i;%1=beijing;2=whale;3=chinese;4=cpd_fish;5=fish_2d;
X = load_testdata(data_name);
[X,~] = normalize_point(X,1);
plot(X(:,1)+xshift(i),X(:,2)+yshift(i),'r.','markersize',20);hold on;
end
set(h,'position',[300,500,1200,300]);
set(h,'color','w');
set(gca,'position',[-0.02,-0.05,1.05,1.05])
axis off