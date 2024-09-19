function Temp = load_Tini_Tsolver(Tini_dir,GCOORD,temp_bc)
% TEMP = LOAD_TINI_TSOLVER(TINI_DIR,GCOORD,X_MIN,TEMP_BC) calculates the
% values of temperatures for a scattered mesh defined by GCOORD node
% coordinates from temperatures of a regular mesh calculated with
% 'Temperature Solver' developed by Albert de Monserrat. TINI_DIR is the
% file where the temperatures and the regular mesh are stored

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 29-10-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

disp(['Loading initial temperatures from ',Tini_dir])

% Load files
load(Tini_dir,'T','X','Y');

% Transforming coordinates
X = X+min(GCOORD(1,:));
Y = Y+max(GCOORD(2,:));

% Interpolation function
F = TriScatteredInterp(X(:),Y(:),T(:)-273.15);
% Interpolate temperatures
Temp = F(GCOORD(1,:),GCOORD(2,:));

% Fix not interpolated temperatures in the boundaries
Temp(isnan(Temp)) = temp_bc;