function [S] = get_sed_thickness(Basement, Topography)
  % calculate sediment thickness st topography coordinates
  Topo_base = interp1(Basement(1,:), Basement(2,:), Topography(1,:));
  S         = Topography(2,:) - Topo_base;
  S(S < 0.)    = 0.;
  S = [Topography(1,:); S];
end
