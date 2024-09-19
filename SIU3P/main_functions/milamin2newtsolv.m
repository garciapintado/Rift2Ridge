% Script to make MILAMIN variables compatible with Newton's solver
% mechanical2d_N

% TODO change MILAMIN mechanical2d and mechanical2d_N so that they have
% compatible variables and this script and newt2milamin can be removed

YIELD.Phi       = Phi; % TODO CHANGE inside the solver so that sind is sin and cosd cos
YIELD.Cohesion  = Cohesion;
YIELD.Phases    = ismember(Phases,PLASTICITY.Phases);
YIELD.yc        = PLASTICITY.yc;
YIELD.scheme    = PLASTICITY.t;

g       = G;
dP_dz   = Rho(1)*G(2);

MESH.EL2NOD = ELEM2NODE([1 2 3 6 4 5 7],:);
MESH.GCOORD = GCOORD;
MESH.nnod   = nnod;
MESH.nel    = nel;
MESH.Point_id   = Point_id;
MESH.Phases = Phases;
MESH.Corner_id  = Corner_id;

nnodel      = size(ELEM2NODE,1);

ifix = Bc_ind_fs;
vfix = Bc_val_fs;

% INITIAL GUESS FOR VERTICAL VELOCITY PROFILE
EL2NODP = reshape(1:(3*MESH.nel),3,[]);
nUg     = 2*size(GCOORD,2);
nPg     = length(EL2NODP(:));
U       = zeros(nUg,1);
P       = zeros(nPg,1);

FSSA.alpha  = alpha;
FSSA.beta   = beta;

Load_el = zeros(1,sum(Point_id==max(Point_id)-1)+2);

ELASTICITY.Shearm   = Shearm;
ELASTICITY.s        = elasticity_s;
ELASTICITY.dt       = dt;

TAU.xx_old          = TAU_xx_old;
TAU.zz_old          = TAU_yy_old;
TAU.xz_old          = TAU_xy_old;

SS.F_xx     = F_xx;
SS.F_xz     = F_xy;
SS.F_zx     = F_yx;
SS.F_zz     = F_yy;

MU.eff      = zeros(MESH.nel,nip);
MU.plast    = zeros(MESH.nel,nip);
MU.creep    = zeros(MESH.nel,nip);