% SCALE VARIABLES

% Load the scaling factors
L0 = SCALE.h;
V0 = SCALE.vel;
t0 = L0/V0;
mu0 = SCALE.mu;
P0 = mu0/t0;

G = G/V0*t0;
dt = dt/t0;
time_int = time_int/t0;

BREAKUP.min_thick = BREAKUP.min_thick/L0;
BREAKUP.break_time = BREAKUP.break_time/t0;

% TODO Look at the scaling of the thermal part

% TODO Look at the scaling of the sediments

% TODO Look at scaling for non-newtonian rheologies
RHEOL.Adis = RHEOL.Adis/mu0;
RHEOL.Adif = RHEOL.Adif/mu0;

Shearm = Shearm/P0;

Cohesion = Cohesion/P0;
SS.C = SS.C/P0;

Rho = Rho*L0^2/(mu0*t0);

ext_rate = ext_rate/V0;

x_max = x_max/L0;
x_min = x_min/L0;
y_min1 = y_min1/L0;
y_max1 = y_max1/L0;
y_max2 = y_max2/L0;
base_lithos = base_lithos/L0;
Ylim = Ylim/L0;
UCR = UCR/L0^2;
LCR = LCR/L0^2;
Elsizes = Elsizes/L0^2;
shift = shift/L0;