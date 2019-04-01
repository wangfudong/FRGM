function Map = FRGM_nore(X,Y,SX1,SX2,option)
% main function to compute the correspondence map
geofunc = option.geofunc;
LX = length(X(:,1));
LY = length(Y(:,1));
M = M_shape(X,Y,1/8,2);
Map_ini = asgHun(-M);
option.convex = 0;
option.geofunc = geofunc;

%% non-convex objective function
[Map_cg1] = FW_ROT(X,Y,M,SX1,Map_ini,option);
if LX == LY
    %% convex objective function
    option.convex = 1;
    option.geofunc = '0.11';
    X_maped = Map_cg1*Y;
    M_new = M_points(X_maped,Y);
    M_new = M_new/max(M_new(:));
    Map = FW_ROT(X_maped,Y,M_new,SX2,Map_cg1,option);
else
    %% non-convex objective function
    X_maped = Map_cg1*Y;
    M_new = M_points(X_maped,Y);% adjust the unary term for unequal-sized cases
    M_new = M_new/max(M_new(:));
    
    Map_ini = asgHun(-M_new);%initialization
    option.convex = 0;
    option.geofunc = geofunc;
    
    [Map_cg2] = FW_ROT(X,Y,M_new,SX1,Map_ini,option);
    %% convex objective function
    option.convex = 1;
    option.geofunc = '0.11';
    X_maped = Map_cg2*Y;
    
    M_new = M_points(X_maped,Y);
    M_new = M_new/max(M_new(:));
    
    Map = FW_ROT(X_maped,Y,M_new,SX2,Map_cg2,option);
end



