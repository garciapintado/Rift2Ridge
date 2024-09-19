function GEO = resampleGEO(GEO, GEOres, layshift, check_monotony)
  % +++ purpose +++
  % Add nodes to GEO according to a given soft resolution [GEOres].
  % resulting GEO segments will a length in [0.75,1.5]*GEOres(i) for each interface i 
  %
  % GEOres   :: REAL (ngeo) [m] nominal resolution for each interface
  % layshift :: REAL (scalar) [m] lateral node collapse distance toward junctions 
  
  % Author: Javier García-Pintado, MARUM, 2020-03
  
  ngeo = length(GEO);
  if (length(GEOres) ~= ngeo)
      error("resampleGEO:: ---ERR001---")
  end
  
  for i=ngeo:-1:1
      Lib = GEO(i).coo;                                                    % [2,nb]
      Ldiff = diff2D(Lib);                                                 % [nb-1]
      nnn  = floor(Ldiff / (GEOres(i)*1.5));                               % [nb-1] number of new nodes to be added to each segment: new segments will range in [0.75,1.5]*GEOres(i)
      idsb = 1:size(Lib,2);                                                % [nb] 
      idsbexp = repelem(idsb, [nnn+1 1]);                                  % [na] expanded indices
      idsbnn = find(nnn ~= 0);                                             % indices in Lib previous to insertion of new nodes
      Lia = Lib(:,idsbexp);                                                % [2,na] expanded coordinates
      for k=1:length(idsbnn)                                               % fill block of new coordinates
          xyb = Lib(:,idsbnn(k)+(0:1));                                    % bounds for new segment
          xya = [linspace(xyb(1,1),xyb(1,2),nnn(idsbnn(k))+2);
                 linspace(xyb(2,1),xyb(2,2),nnn(idsbnn(k))+2)];
          idsa = find(idsbexp==idsbnn(k));
          Lia(:,idsa(2:end)) = xya(:,2:end-1);                             % plotGEO(GEO)
      end                                                                  % hold on; plot(Lia(1,:)/1000,Lia(2,:)/1000,'x','color','red','markersize',15)
      %if ~isequal(Lib,Lia)
          GEO(i).coo = Lia;
          GEO(i).n = size(Lia,2);
          if GEO(i).horizontal
              GEO = addMissingGEO(GEO, layshift, false, check_monotony, i); % TODO: add imax/i to addMissingGEO to make it faster, so it does not need to reiterate over all layers
          end
      %end
  end
  GEO = linkGEO(GEO);                                                      % remake GEO.to
end % function resampleGEO()
