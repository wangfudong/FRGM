function [X] = load_testdata(filename)
switch filename
    case 1
    beijing=load('cdfHC_beijing.mat');
    X = beijing.test;
    case 2
        load('cdfhc_data2D_whale.mat')
        X = whale{1};
    case 3
        load('cpd_chinese.mat')
        X = chinese;
    case 4
        load('cpd_fish.mat');
        X = X;
    case 5
        load('gmmreg_fish2d.mat');
        X = gmmreg_fish2d;
        
    case 6
        load('gmmreg_3dface.mat');
        X = face_3d;
        
end
