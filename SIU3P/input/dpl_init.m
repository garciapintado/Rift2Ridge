function Dpl = dpl_init(GCOORD, ELEM2NODE, Phases, base_lithos)
   %
   % +++ purpose +++
   % obtain initial conditions of depletion profile
   % Mantle rheology changes from wet olivine with a water content of 125 ppm H/Si (Hirth & Kohlstedt, 1996; Rüpke at al., 2006).  
   % As melting occurs in the model, the increasing depletion affects the mantle rheology and its density. 
   %
   % GCOORD      :: REAL [2,nnod]
   % base_lithos :: REAL, Lithosphere-asthenosphere boundary (LAB).
   %                      Level at which depletion is 4% [dry olivine, assuming all water has been extracted at this stage (Morgan, 1997)].
   %
   % Javier GP: reorganized previous stepwise interpolation as one-liner

   km = 1000;

   % Transition points from bottom to top
   zmin     = min(GCOORD(2,:));                            % base of box : no depletion
   top_wet  = -130*km;                                     % top level at which we assume no melted mantle [wet olivine rheology] so initial depletion is still 0 [Elena indicated -160 km in her paper]
   % base_lithos                                           % initiation of melt
   zDpl_01  =  -60*km;                                     % level at which initial depletion is 10%
   zmax     = max(GCOORD(2,:));

   if zmin > top_wet || top_wet > base_lithos || base_lithos > zDpl_01 || zDpl_01 > zmax
      error("Dpl_init :: non-monotonic interpolation support levels");
   end

   %==========================================================================
   % CALCULATE RHEOLOGIC VARIATION FACTORS
   %==========================================================================
   % Wet asthenosphere
   Dpl = interp1 ([zmin, top_wet, base_lithos, zDpl_01], ... 
                  [0.0,  0.0,     0.04,        0.04   ], ...
                  GCOORD(2,:), 'linear', 'extrap');
   notmantle = ~ismember(1:size(GCOORD,2), unique(ELEM2NODE(:, Phases == 1))); 
   Dpl(notmantle) = 0.0;
end

