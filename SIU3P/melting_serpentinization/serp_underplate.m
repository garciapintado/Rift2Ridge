function [dFserp, Dserp] = serp_underplate(GCOORD, ELEM2NODE, Temp, Dserp, dt, ...
                                  TRACKP_melt, srfac, PHY, istep)
  %
  % this function uses the underplating criterion by Ros_al2017 (originally coded by Miguel Andres Martinez),
  % the serpentinization equation by Emmanuel and Berkowitz (2006)
  % and a simple serpentinization possibility criteria based on the
  % existence of hydrothermal domain and temperature range
  % 
  % TRACKP_melt :: melt tracking points to evaluate underplatting. Empty implies no melting.
  %                [1:2,:] :: x,y coordinates
  %                [3,:]   :: step index for which the melt trackpoints were generated 
  % srfac :: REAL [1,nnod], \in (0,1). Multiplier of the kinetics of the serpentinization reaction. 
  %          This factor considers the plastic strain rate as proxy for porosity, so exercing a control on the kinetics of serpentinization.should already
  %          It is assumed here that this factor already considers the effect of hydrothermal masking and phase (srfac=0.0 for non-mantle nodes)
  % PHY   :: physical parameter struct for serpentinization, with the elements
  %    .Tmin
  %    .Tmax
  %--------------------------------------------------------------------------
  
  % Javier Garcia_Pintado, 09.09.2020
  
 
  Tk = Temp+273.15;
  
  %Serpentinization for low temperature 
  switch PHY.equation
      case 'E2006'                              % Emmanuel and Berkowitz (2006). For example, used in Elena Ros's Thesis
          b_s = 2.5e-4;                                 % [Kelvin^-2]
          c_s = 543.;                                   % [Kelvin] (= 270 ºC) temperature for maximum kinetics
          Krl = exp(-b_s*(Tk-c_s).^2);                  % [1,nnod7] [s-1]
      case 'M2012'                              % Malvoisin et al. (2012) after Lasaga (1995). For example, used in Rupke and Hasenclever (2017)
          a  = 808.3;                                   % [-]
          b  = 3640.;                                   % [K]
          c  = 8759.;                                   % [K]
          T0 = 623.6;                                   % [K] [350.45 + 273.15]
          Krl = a.*exp(-b./Tk) .* (1 - exp(-c.*(1./Tk - 1/T0)));
          Krl = max(Krl,0.0);
      case 'SKGAU'                              % Javier GP. MARUM, 2021 fits 'M2012' but with a smooth decay toward high temperatures [approaches 0.0 around 400 ºC]          
          T0 = 340 + 273.15;
          w = 94;
          alpha = -4;
          x = (Tk - T0)/w;
          phi = exp(-0.5*x.^2);
          Phi = 1/2 * (1 + erf(alpha*x/sqrt(2)));                          % \in [0,1]
          Krl = 2 * phi .* Phi;
          Krl = Krl/max(Krl);
  end
  % Temp = 0:1:500;
  % figure(); plot(Temp,Krl_1,'color','blue'); hold on; plot(Temp,Krl_2,'color','red'); hold on; plot(Temp,kn,'color','green'); hold on; plot(Temp,ksk,'color',[0 .7 .0])
  Krl = PHY.k0 * Krl .* srfac;                                             % apply strainrate factor [proxy for porosity/permeability decay]
 
                  

  % Serpentinisation for high temperature [currently not used]
  % T_max = c_s-273;
  % T_lim = 520;
  % Aht = 1/(T_max - T_lim);
  % Bht = -Aht*520;
  % Krh = A_s*(Aht*Temp + Bht);

  nnod7 = size(GCOORD,2);

  % Localize triangular elements where there is melting as well as nodes in these elements
  nounderp = true(1,nnod7);                                                % LOGICAL, not underplatted
  if ~isempty(TRACKP_melt)
      trackpmeltm = TRACKP_melt(1:2,:);                        % [2,:] matrix
      nnewmeltp = sum(ismember(TRACKP_melt(4,:), istep));      % INTEGER, number of points TRACKP_melt has from the current time step
      if nnewmeltp > 0
          trackpmeltm = trackpmeltm(:,1:end-nnewmeltp);        % remove current time step tracked points
      end
      trackpmeltm = trackpmeltm(:,all(~isnan(trackpmeltm)));   % disregard columns with NaN values
      if ~isempty(trackpmeltm)
          Tris = tsearch2(GCOORD,uint32(ELEM2NODE(1:3,:)),trackpmeltm);
          Ind22 = find(Tris==0);                             % indices for which envelope elements were not found
          if ~isempty(Ind22)
              for i=1:length(Ind22)
                  [~, Tris(Ind22(i))] = min(sqrt((GCOORD(1,ELEM2NODE(7,:)) - trackpmeltm(1,Ind22(i))).^2 + (GCOORD(2,ELEM2NODE(7,:)) - trackpmeltm(2,Ind22(i))).^2)); % nearest element centroid for this trackpoint
              end
          end
          els_underp   = unique(Tris);                      % index of elements surrounding tracked melt points
          nodes_underp = unique(ELEM2NODE(1:7,els_underp)); % index of nodes touching elements around emplaced  magma 
          nounderp  = ~ismember(1:nnod7, nodes_underp);     % [1,nnod7] true for nodes not in an element with magma emplacement
          % figure(4);plot_mesh; hold on; plot(trackpmeltm(1,:)/1000,trackpmeltm(2,:)/1000,'*m')
          % hold on; plot(GCOORD(1,nounderp)/1000,GCOORD(2,nounderp)/1000,'*g')
      end
  end % ~isempty(TRACKP_melt)

  nodboo = nounderp & Temp >= PHY.Tmin & Temp <= PHY.Tmax;                 % underplatting & hydrothermal domain & temperature criteria

  Kr = Krl;
  % Kr(logical_temp) = Krl(logical_temp);
  % Kr(~logical_temp) = Krh(~logical_temp);
  Kr(~nodboo) = 0.;
  dFserp = Kr .* (1-Dserp)*dt;                                             % [1,nnod7] % serpentinization ammount [0/1]
  dFserp(dFserp + Dserp > 1) = 1 - Dserp(dFserp + Dserp > 1);              % [1,nnod7] such that the increment wont allow for  Dserp > 1.
  Dserp = dFserp + Dserp;                                                  % [1,nnod7] \in [0,1]
  
end % function
