function GEO = addMissingGEO(GEO, layshift, addto, check_monotony, hidmax)
  % add missing nodes to isolated edges detected within collapsed horizontal interfaces in GEO
  %
  % GEO       :: STRUCT with interface definitions
  % addto     :: whether to add the "to" field to GEO, where GEO.to contains internal GEO(i) -> GEO(j) mapping information, with j>i
  % monotonic :: 
  
  % Author: Javier Garcia-Pintado, 2020-03., MARUM
  
  
  if nargin < 3
      addto = true;                                                        
  end
  if nargin < 4
      check_monotony = true;                                                    % check for IO monotonicity
  end
  if nargin < 5
      hidmax = length(GEO);
  end
  
  for i=1:length(GEO)
    if GEO(i).horizontal
        if any(diff(GEO(i).coo(1,:)) <= 0.) && check_monotony
            error("addMissingGEO: non-monotonic x input in horizontal layers")      
        end
    end
  end
  
  hids = find([GEO.horizontal]);                                           % bottom-top horizontal interface indices within GEO
  
  for ih = length(hids):-1:2                                               % from top to bottom phases    
      i = hids(ih-1);                                                      % indices in GEO - bottom of this subdomain
      j = hids(ih);                                                        % indices in GEO - top     "   "
      if j > hidmax
          continue;
      end
      L1 = GEO(i).coo;                                                     % mutable [lower] interface
      L2 = GEO(j).coo;                                                     % fixed [upper] interface
      L = addMissingNodes({L1,L2}, false);                                 % false => L1 -> L2 mapping
      if ~isequal(L{2},L2)
          %error("addMissingGEO: addMissingNodes modified upper interface for assymetric upward mapping")
          GEO(j).coo = L{2};
          GEO(j).n   = size(GEO(j).coo,2);
      end
      if ~isequal(L{1},L1)                                                 % hold on; plot(L{k}(1,:)/1000, L{k}(2,:)/1000,'--','Color','blue')
          GEO(i).coo = L{1};
          GEO(i).n   = size(GEO(i).coo,2);
      end
  end
  if addto
      GEO = linkGEO(GEO);
  end
  if check_monotony                                                        % check and enforce x-monotony
      xmono = checkMonotonyGEO(GEO, layshift/3.);
      if ~xmono
          disp("addMissingGEO: process changed x-monotony in horizontal layers. Correcting...")
          GEO = makeMonotonicGEO(GEO, 1.0); 
      end
  end
end % function addMissingGEO
