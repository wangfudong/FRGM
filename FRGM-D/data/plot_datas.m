%% polt datas of all methods
%% cdfHC
clear
clc
load('cdfhc_data2D_beijing.mat');

figure,
for i = 1:7
    subplot(2,4,i),plot(beijing{i}(:,1),beijing{i}(:,2),'r.');hold on;
end

%% cdfHC
clear
clc
load('cdfhc_data2D_CC.mat');
figure,
for i = 1:7
    subplot(2,4,i),plot(CC7{i}(:,1),CC7{i}(:,2),'r.');hold on;
end
%% cdfHC
clear
clc
load('cdfhc_data2D_fish.mat');
figure,
for i = 1:7
    subplot(2,4,i),plot(fish2d{i}(:,1),fish2d{i}(:,2),'r.');hold on;
end
% whale = fish2d;
% save('cdfhc_data2D_whale.mat','whale');
%% gemreg 3dface

load('gmmreg_3dface.mat');
figure,scatter3(face_3d(:,1),face_3d(:,2),face_3d(:,3),'ro');
%% gmmreg 2dface
load('gmmreg_fish2d.mat');
figure,scatter(gmmreg_fish2d(:,1),gmmreg_fish2d(:,2),'ro');
%% majiayi L2E
load('save_chinese_def.mat');
figure,plot(x1(:,1),x1(:,2),'r.',y2(:,1),y2(:,2),'bo');

