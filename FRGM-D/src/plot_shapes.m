function plot_shapes(X_raw,Y_raw,X_tran)
hh = figure;subplot(1,2,1);
set(hh,'position',[500 450 1000 300]);
plot(X_raw(:,1),X_raw(:,2),'r.',Y_raw(:,1),Y_raw(:,2),'b.','markersize',20);
subplot(1,2,2),plot(X_tran(:,1),X_tran(:,2),'ro','linewidth',1.5,'markersize',6);
hold on,plot(Y_raw(:,1),Y_raw(:,2),'b.','markersize',10);
clear hh