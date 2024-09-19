function [ipx, ipw] = ip_line(nip)
% function [ipx, ipw] = ip_line(nip)
% integration points & weights for 1D elements
%
% canonical element in [-1,1]
%
% Author: Javier Garcia-Pintado, MARUM, 2020

  switch nip
    case 1
        ipx = 0.;
        ipw = 2.;
    case 2
        ipx(1) = -sqrt(1/3.);
        ipx(2) =  sqrt(1/3.);

        ipw(1) = 1.;
        ipw(2) = 1.;
    case 3
        ipx(1) = -sqrt(3/5.);
        ipx(2) = 0.;
        ipx(3) =  sqrt(3/5.);

        ipw(1) = 5/9.;
        ipw(2) = 8/9.;
        ipw(3) = 5/9.;
    otherwise
        error('integration rule not prepared')
  end
end % function
