function [Map,para] = FRGM_pr(X,Y,option)

if strcmp(option.regist_trans,'nonrigid')
    [X,~] = normalize_point(X,1);
    sigma = option.sigma;
    GX = nonrigid_ker(X,sigma,'rbf');
else
    GX = [];
end
[Map,para] = FRGM_registration(X,Y,GX,option);
