function [Map,Para] = FRGM_registration(X,Y,GX,option)

if nargin < 2, error('No enough inputs!'); end
if nargin < 4, option.trans = 'similar'; end

if ~isfield(option,'regist_trans') || isempty(option.regist_trans), option.regist_trans = 'similar'; end;
if ~isfield(option,'regist_it') || isempty(option.regist_it), option.regist_it = 10; end;
if ~isfield(option,'regist_tol') || isempty(option.regist_tol), option.regist_tol = 1.0e-6; end;
if ~isfield(option,'regist_normalize') || isempty(option.regist_normalize), option.regist_normalize = 1; end;
if ~isfield(option,'regist_rota') || isempty(option.regist_rota), option.regist_rota = 1; end;
if ~isfield(option,'regist_display') || isempty(option.regist_display), option.regist_display = 0; end;


if ~isfield(option,'GM_nonvex_iter') || isempty(option.GM_nonvex_iter), option.GM_nonvex_iter = 100; end;
if ~isfield(option,'GM_convex_iter') || isempty(option.GM_convex_iter), option.GM_convex_iter = 100; end;
if ~isfield(option,'GM_tol') || isempty(option.GM_tol), option.GM_tol = 1.0e-6; end;
if ~isfield(option,'GM_lambda1') || isempty(option.GM_lambda1), option.GM_lambda1 = 1; end;
if ~isfield(option,'GM_lam_nonvex') || isempty(option.GM_lam_nonvex), option.GM_lam_nonvex = 100; end;
if ~isfield(option,'GM_lam_convex') || isempty(option.GM_lam_convex), option.GM_lam_convex = 1; end;
if ~isfield(option,'GM_lam_sparse') || isempty(option.GM_lam_sparse), option.GM_lam_sparse = 0; end;
if ~isfield(option,'GM_initial') || isempty(option.GM_initial), option.GM_initial = 'lap'; end;
if ~isfield(option,'GM_convex') || isempty(option.GM_convex), option.GM_convex = 0; end;
if ~isfield(option,'GM_geofunc') || isempty(option.GM_geofunc), option.GM_geofunc = '1.21'; end;
if ~isfield(option,'GM_unary') || isempty(option.GM_unary), option.GM_unary = 1; end;
if ~isfield(option,'GM_connected') || isempty(option.GM_connected), option.GM_connected = 'full'; end;
if ~isfield(option,'GM_knn') || isempty(option.GM_knn), option.GM_knn = 10; end;
if ~isfield(option,'GM_remove') || isempty(option.GM_remove), option.GM_remove = 0; end;

if ~isfield(option,'GM_convex_or_non') || isempty(option.GM_convex_or_non)
    alter = zeros(1,option.regist_it);
    alter(1:2) = 1;
    option.GM_convex_or_non = alter;
end;

order = option.order;

switch option.regist_trans 
    case 'similar'
        [Map,Para] = FRGM_similar(X,Y,option,order);
    case 'affine' 
        [Map,Para] = FRGM_affine(X,Y,option,order);       
    case 'nonrigid'
        [Map,Para] = FRGM_nonrigid(X,Y,GX,option,order);
end









