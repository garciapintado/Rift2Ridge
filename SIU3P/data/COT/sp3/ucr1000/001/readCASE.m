% list of available components
%
% minimal components: one of the three following has to be selected
% "nchem" : indicates no chemistry
% "schem" : simple chemistry
% "xchem" : extended chemistry
%
% optional components

% "hydr" : get hydrothermal domain and allows to activate hydrothermal model or parameterised hydrothermal enhanced thermal conductivity
% "melt"
% "serp"  : activates serpentinization. As it requires the calculation of the hydrothermal domain, it sets
%           SETTINGS.HYDR.make = true;
%           PHY.HYDR.PARAM.cfmax = 1.0;
% "spro"  : activates surface processes
%
%
% the component vector of strings activates the corresponding ones
%
% The minimal COMPSET is "nchem", which does not activate any specific components
% Note that COMPSET "hydro_nchem" is equivalent to "hydro".

components = ["nchem","hydr","melt","serp","spro"]
COMPSET = join(components,'_');
