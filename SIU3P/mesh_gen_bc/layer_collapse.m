function GEO = layer_collapse(GEO, Base, SETTINGS)
% GEO = layer_collapse(GEO,BASE,SHIFT) 
% Takes the horizontal interfaces defined by GEOMETRY 
% and GEO_ID, find where they interfere, or get vertically closer than 
% SHIFT and merge these interfaces together. Collapsing of interfaces makes
% Triangle library to understand that it should not generate elements for 
% the collapsed layer in this area, so that discontinuous layers are 
% created. In this way a reduction of the number of excessively
% small elements is achieved together with the reduction of the amount of
% times remesh is needed. When BASE is the top basement coordinates (and
% not empty) updates the upper layer (i.e. the crust) to be of the same
% thickness of the sediments when sediment thickness is >= SHIFT.
%
% Note I: It is recommended that SHIFT is bigger than at least 1/5 of the
% maximum resolution of the interfaces to avoid elements with high angles
% and constant remeshings (i.e. where resolution at interfaces is 1 km
% SHIFT should be at least 200 m).
%
% Note II: This algorithm can only work with layers that have at least 4
% interfaces, because of the numbering assumed for the identification of
% interfaces 1 and 2.
%
% For this, the function finds intersections of the upper interfaces (2)
% with the lower interfaces (1) + SHIFT. This intersections are then used 
% to define segments where interface 1 + SHIFT is above interface 2, and
% poitns surrounding the intersections (i11,i12,i21,i22). Then, points
% between i11 and i12 are substituted by points between i21 and i22
% (including them). 
%
% Initial state:
%___.____.                           .____ Interphase 2
%         \                         /
%          \    .______.____.______/_.____ Interphase 1 + SHIFT
%           \  /                  /  i12
%            \/                  /
%            /\.______.____.___./
%____._____./  i21            i22
%         i11
%
% Final state:
%___.____.                           .____ Interphase 2
%         \                         /
%          \                       / .____ Interphase 1
%           \                     / /i12
%            \                   /_/
%             \.______.____.___./
%____._____._/ i21            i22
%         i11

  %--------------------------------------------------------------------------
  % Authors
  % Miguel Andrés _Martinez; 2018 first version
  % Javier GP:               2020-02-20: 
  %      Overall new function for layer collapse.
  %      Sediment addition and angle  corrections at the end as Miguel's. 
  %        process:
  %      - check monotonocity - this maybe be alleviated in a future
  %      - sediment-based update of surface layers (based on Miguel's)
  %      - clip2GEOM: upward horizontal layer collapse  
  %      - addMissingGEO(): upward horizontal missing node inclussion in collapsed interfaces
  %      - small angle identification/correction (Miguel's original stays)
  %      - addMissingGEO()
  % Notes:
  %      - simplifyGEO(), a within layer clipping, is conducted here and within remesh()
  %      - resampleGEO(), a GEO node augmentation to keep resolution, is conducted within remesh()
  %      - The code tries to be as much unforgiving as possible. If input
  %        is monotonic, monotonicity is required for addMissingGEO(). If input is not monotonic but min(dx) < - shift/3,  monotonicity violation is accepted
  %        in the first (out of two) instance of addMissingGEO(). As this temporary violation may be solved the angle correction. If eventually the
  %        second call to addMissingGEO() gives an error, also
  %        monotonicity allevation should be likely be applied there...
  %
  %--------------------------------------------------------------------------
  sedtocrust = SETTINGS.sedtocrust; 
  shift      = SETTINGS.layshift; 
  
  minan = 15;                                                                % threshold for small angle correction

  hboo = [GEO.horizontal];
  hids = find(hboo);                                                         % bottom-top horizontal interface indices within GEO
  ngeo = length(GEO);

  monotonic = checkMonotonyGEO(GEO, SETTINGS.layshift/3.);
  %GEO0 = GEO;
  if ~monotonic
      disp("layer collapse:: forcing monotony in GEO input")
      GEO = makeMonotonicGEO(GEO, 1.0);                                        % force 1 meter  >= dx monotonic increase   hold on; plotGEO(GEOm,[GEO.gid],'o') 
  end
  %GEO1 = GEO;
  % sediment-based update of surface layer (Note: this operation breaks links in GEO.to)
  if sedtocrust
      sed_shift = 4./3;   % Factor to avoid constant remeshing when sediments are included in the upper layer
      if ~isempty(Base)
    
          if any(diff(Base(1,[1,end])) < 0.) 
              error("layer collapse: x non-monotonic in Basement points")
          end
          Top = GEO(hids(end)).coo;
          BelowTop = GEO(hids(end-1)).coo;                                      % [2,nbt] REAL
            
          top2botx = interp1(Top(1,:), Top(2,:),  BelowTop(1,:));               % [1,nbt] REAL, Interpolate topography into x of interface below
          %isthin = (top2botx-BelowTop(2,:)) < shift;                           % [1,nbt] LOGICAL: where these layers may need updating
          base2botx = interp1(Base(1,:),Base(2,:),BelowTop(1,:));               % [1,nbt] REAL, Interpolate basement [given by track points] into x of the interface below
          isthin = BelowTop(2,:) - base2botx > 0.1 * shift;                     % [1,nbt] LOGICAL: where these layers may need updating
          isthickSed = (top2botx-base2botx) >= (shift*sed_shift);               % [1,nbt] LOGICAL: for sediment thickness height
          boo = isthin & isthickSed;
          
          % hold on; 
          % hold on; plot(BelowTop(1,isthickSed)/1000,BelowTop(2,isthickSed)/1000,'o','markersize',20)
          % hold on; plot(BelowTop(1,:)/1000,base2botx/1000,'x-','markersize',20)
          % hold on; plot(BelowTop(1,boo)/1000,base2botx(boo)/1000,'.','markersize',20)
          
          % Where both conditions apply, update GEO y-coordinates to those of the tracked basement:        
          YcorrSed = base2botx(boo);                            % [ncor <= nbt] REAL
          tol = SETTINGS.tolerance(hids(end-1));
          if ~isempty(YcorrSed)
              for i=1:length(GEO)                                          
                  if i == hids(end)
                      continue;
                  end
                  [tobase, ycorrid] = ismember(GEO(i).coo', ...
                                          BelowTop(:,boo)','rows');        % tobase :: LOGICAL [ncoo], ycorrid:: INTEGER 
                  if any(tobase)
                      GEO(i).coo(2,tobase) = YcorrSed(ycorrid(tobase));
                      % simplify: remove nodes within tol
                      transects = [[1; find(tobase)]+1 [find(tobase); length(tobase)]-1];    % INTEGER [ntr,2] transects of potentially removable nodes     
                      transects = [transects [Inf diff(GEO(i).coo(1,tobase)) Inf]' < tol];       % LOGICAL [ntr] true for transects to be removed
                      transects = transects(transects(:,2) >= transects(:,1), :);               % remove empty transects
                      transects = transects(transects(:,3)==1,1:2);                             % leave only removable transects, one per row
                      rmids = [];
                      for itr=1:size(transects,1)
                          rmids = [rmids transects(itr,1):transects(itr,2)];
                      end
                      GEO(i).coo(:,rmids) = [];
                  end
              end
                  
              for i=1:length(GEO)
                  if GEO(i).horizontal
                      if any(diff(GEO(i).coo(1,:)) < 0.)
                          error("layer collapse: basement modification changed monotonicity")      
                      end
                  end
              end
          end % ~isempty(YcorrSed)
      end % ~isempty(Base)
  end % sedtocrust
  %GEO2 = GEO;
  %monotonic = checkMonotonyGEO(GEO, 0);
  GEO = rmfield(GEO,'pids');
  for ih=length(hids):-1:2                                             % top to bottom pairwise horizontal interfaces
      i = hids(ih-1);                                                  % indices in GEO - bottom of this subdomain
      j = hids(ih);                                                    % indices in GEO - top     "   "
      L1 = GEO(i).coo;                                                 % mutable [lower] interface
      L2 = GEO(j).coo;                                                 % fixed [upper] interface
      L1(2,:) = min(L1(2,:),interp1(L2(1,:),L2(2,:),L1(1,:)));         % truncate (in the accidental case the lower was on top)
      L1a = clip2GEOM(L1, L2, SETTINGS.layshift);                            % clip: return surviving coordinates and corresponding surviving indexes in L1
      if ~isequal(L1, L1a)
          GEO(i).coo = L1a;                                                  % hold on; plot(L1a(1,:)/1000, L1a(2,:)/1000,'x','markersize',20)
          GEO(i).n   = length(GEO(i).coo);
      end
  end
  %GEO3 = GEO;
  monotonic = checkMonotonyGEO(GEO, SETTINGS.layshift);
  if ~monotonic
      GEO = makeMonotonicGEO(GEO, 1.0);
      monotonic = checkMonotonyGEO(GEO, 0.);
      if ~monotonic
          error("layer_collapse: clipping & simplify() made x-non monotonic horizontal GEO interfaces unsolved by makeMonotonicGEO()")
      end
  end
  GEO = simplifyGEO(GEO, SETTINGS.tolerance, SETTINGS.layshift);           % hold on; plotGEO(GEO,[3 6],'+')
  %GEO4 = GEO;
  monotonic = checkMonotonyGEO(GEO, 0.);                                   % simplyfyGEO() should never break x-monotony
  if ~monotonic
      GEO = makeMonotonicGEO(GEO, 1.0);
      monotonic = checkMonotonyGEO(GEO, 0.);
      if ~monotonic
          error("layer_collapse: clipping & simplify() made x-non monotonic horizontal GEO interfaces unsolved by makeMonotonicGEO()")
      end
  end
  % plotGEO(GEO0)
  % plotGEO(GEO)
  % hold on; plotGEO(GEO0,3,'+',20)
  % hold on; plotGEO(GEO,[3,6,9],'o')
  % hold on; plot(GEO(i).coo(1,backids)/1000, GEO(i).coo(2,backids)/1000,'x','markersize',15) 
  % xy = GEO(i).coo;
  % anglesi = getAnglesPolyline(xy); text(xy(1,2:end-1)/1000, xy(2,2:end-1)/1000+0.02, string(anglesi),'color','red') 
  GEO = addMissingGEO(GEO, SETTINGS.layshift, true, false);                % include GEO.to. Do not enforce x-monotony
  if ~checkMonotonyGEO(GEO, SETTINGS.layshift/3.)
      error("layer_collapse: addMissingGEO() made x-non monotonic horizontal GEO interface")
  end 
  %GEO5 = GEO;
  Gid = ngeo-1:-3:3; % Miguel set this order. TODO: change from bottom to top and make GEO-based the following angle-correction code
  [GEOMETRY, Geo_id] = GEO2geometry(GEO); % TODO: use directly GEO in the functions below and get rid of GEOMETRY completely
   
  % small angle identification
  %disp("layer collapse: conducting angle corrections")
  if 3 > 2 % do not do small angle- it does not respect original coordinates nor seems to do great
      for n = 2:length(Gid)
          cid = Gid(n);
          int1x = GEOMETRY(1,Geo_id==cid);
          int1y = GEOMETRY(2,Geo_id==cid);
          
          % Calculate angles
          vect = [diff(int1y); diff(int1x)];
          v1 = -vect(:,1:end-1);
          v2 = vect(:,2:end);
          an = acosd(sum(v1.*v2)./(sqrt(sum(v1.^2)).*sqrt(sum(v2.^2))));
          an = [Inf an Inf];
          smallan = find(an<minan);
          midp = [sum(int1x(smallan+[1; -1]))/2; sum(int1y(smallan+[1; -1]))/2];
          [Gp,~] = ismember(GEOMETRY',[int1x(smallan);int1y(smallan)]','rows');
          GEOMETRY(:,Gp) = [];
          Geo_id(Gp) = [];
      end % small angle identification
      
      % small angle correction
      for n = 2:length(Gid)
          cid = Gid(n);
          int1 = GEOMETRY(:,Geo_id==cid);
          for m = 1:n-1
              int2 = [GEOMETRY(:,Geo_id==Gid(n-m))];
              consecutive_juntions = true;
              countwhile = 0;
              while consecutive_juntions && countwhile<=3
                  countwhile = countwhile+1;
                  repin1 = ismember(int1',int2','rows')';
                  right = [false repin1(1:end-2) false];
                  left = [false repin1(3:end) false];
                  junctions1 = find(right~=left & repin1);
                  repin2 = ismember(int2',int1','rows')';
                  right = [false repin2(1:end-2) false];
                  left = [false repin2(3:end) false];
                  junctions2 = find(right~=left & repin2);
                  junctions2 = unique([junctions2 ...
                      find(ismember(int2',int1(:,junctions1)','rows'))']);
                  junctions1 = unique([junctions1 ...
                      find(ismember(int1',int2(:,junctions2)','rows'))']);
                  indepint2 = ~ismember(int2(:,junctions2+1)',int1','rows')'.*(junctions2+1);
                  indepint2(indepint2==0) = junctions2(indepint2==0)-1;
                  indepint1 = ~ismember(int1(:,junctions1+1)',int2','rows')'.*(junctions1+1);
                  indepint1(indepint1==0) = junctions1(indepint1==0)-1;
                  
                  junc_int = find(diff(junctions1)==1);
                  if ~isempty(junc_int)
                      int1_fj = int1(:,1:junctions1(junc_int(1)));
                      for fj = 1:length(junc_int)
                          segment_above = int2(:,junctions2(junc_int(fj))+1: ...
                              junctions2(junc_int(fj)+1)-1);
                          if fj==length(junc_int)
                              next_int1 = int1(:,junctions1(junc_int(fj))+1:end);
                          else
                              next_int1 = int1(:,junctions1(junc_int(fj))+1: ...
                                  junctions1(junc_int(fj+1)));
                          end
                          int1_fj = [int1_fj segment_above next_int1];
                      end
                      int1 = int1_fj;
                  else
                      consecutive_juntions = false;
                  end
              end
              
              vect = [diff(int1y); diff(int1x)];
              v1 = int1(:,indepint1)-int1(:,junctions1);
              v2 = int2(:,indepint2)-int1(:,junctions1);
              an = acosd(sum(v1.*v2)./(sqrt(sum(v1.^2)).*sqrt(sum(v2.^2))));
              
              smallan = find(an<minan);
              % Small angle junctions in other interfaces
              Bel_ind = Gid(Gid<cid);
              clear IntBelow JLI JO
              for ib = 1:length(Bel_ind)
                  IntBelow{ib} = GEOMETRY(:,Geo_id==Bel_ind(ib));
                  [JuncLowerInt,JuntOrd] = ...
                      ismember(IntBelow{ib}',int1(:,junctions1(smallan))','rows');
                  JLI{ib} = find(JuncLowerInt');
                  JO{ib} = JuntOrd(JuncLowerInt)';
              end
              for o = 1:length(smallan)
                  jun1 = junctions1(smallan(o));
                  jun2 = junctions2(smallan(o));
                  ind1 = indepint1(smallan(o));
                  ind2 = indepint2(smallan(o));
                  juncleft = jun1<ind1 & jun2<ind2;
                  juncright = jun1>ind1 & jun2>ind2;
                  
                  midp = (int1(:,jun1)+int2(:,ind2))/2;
                  if juncleft
                      int1 = [int1(:,1:jun1) midp int1(:,jun1+1:end)];
                      int2 = [int2(:,1:jun2) midp int2(:,jun2+1:end)];
                      for ib = 1:length(Bel_ind)
                          O = find(JO{ib}==o);
                          if ~isempty(O)
                              IntBelow{ib} = [IntBelow{ib}(:,1:JLI{ib}(O)) midp IntBelow{ib}(:,JLI{ib}(O)+1:end)];
                              JLI{ib}(JLI{ib}>JLI{ib}(O)) = JLI{ib}(JLI{ib}>JLI{ib}(O))+1;
                          end
                      end
                      junctions1(junctions1>jun1) = ...
                          junctions1(junctions1>jun1)+1;
                      junctions2(junctions2>jun2) = ...
                          junctions2(junctions2>jun2)+1;
                      indepint1(indepint1>ind1) = ...
                          indepint1(indepint1>ind1)+1;
                      indepint2(indepint2>ind2) = ...
                          indepint2(indepint2>ind2)+1;
                  elseif juncright
                      int1 = [int1(:,1:jun1-1) midp int1(:,jun1:end)];
                      int2 = [int2(:,1:jun2-1) midp int2(:,jun2:end)];
                      for ib = 1:length(Bel_ind)
                          O = find(JO{ib}==o);
                          if ~isempty(O)
                              IntBelow{ib} = [IntBelow{ib}(:,1:JLI{ib}(O)-1) midp IntBelow{ib}(:,JLI{ib}(O):end)];
                              JLI{ib}(JLI{ib}>JLI{ib}(O)) = JLI{ib}(JLI{ib}>JLI{ib}(O))+1;
                          end
                      end
                      junctions1(junctions1>jun1) = ...
                          junctions1(junctions1>jun1)+1;
                      junctions2(junctions2>jun2) = ...
                          junctions2(junctions2>jun2)+1;
                      indepint1(indepint1>ind1) = ...
                          indepint1(indepint1>ind1)+1;
                      indepint2(indepint2>ind2) = ...
                          indepint2(indepint2>ind2)+1;
                  end
              end % for o
              GEOMETRY = [GEOMETRY(:,Geo_id<cid) int1 GEOMETRY(:,Geo_id>cid)];
              Geo_id = [Geo_id(Geo_id<cid) cid*ones(1,size(int1,2)) ...
                  Geo_id(Geo_id>cid)];
              GEOMETRY = [GEOMETRY(:,Geo_id<Gid(n-m)) int2 ...
                  GEOMETRY(:,Geo_id>Gid(n-m))];
              Geo_id = [Geo_id(Geo_id<Gid(n-m)) Gid(n-m)*ones(1,size(int2,2)) ...
                  Geo_id(Geo_id>Gid(n-m))];
              for ib = 1:length(Bel_ind)
                  GEOMETRY = [GEOMETRY(:,Geo_id<Bel_ind(ib)) IntBelow{ib} ...
                      GEOMETRY(:,Geo_id>Bel_ind(ib))];
                  Geo_id = [Geo_id(Geo_id<Bel_ind(ib)) ...
                      Bel_ind(ib)*ones(1,size(IntBelow{ib},2)) ...
                      Geo_id(Geo_id>Bel_ind(ib))];
              end
          end % for m
      end % for small angle correction
      
      for n = 1:length(Gid)
          GEOuni = unique(GEOMETRY(:,Geo_id==Gid(n))','rows','stable')';
          GEOMETRY = [GEOMETRY(:,Geo_id<Gid(n)) GEOuni ...
              GEOMETRY(:,Geo_id>Gid(n))];
          Geo_id = [Geo_id(Geo_id<Gid(n)) ...
              Gid(n)*ones(1,size(GEOuni,2)) ...
              Geo_id(Geo_id>Gid(n))];
      end
      
      % plot_geomF(GEOMETRY, Geo_id)
      % plot_geomF(GEOMETRY, Geo_id, 6, 'x')
      % plot_geomF(GEOMETRY, Geo_id, 3, 'o')
      
      monotonic = true;
      for i=1:length(Gid)
           if any(diff(GEOMETRY(1,Geo_id==Gid(i))) < 0.)
               monotonic = false;
           end
      end
      if monotonic
          GEO = geometry2GEO(GEOMETRY,Geo_id);
      else
          disp('layer_collapse:: small angle correction would change monotonicity and is not applied here')
      end
  end % un-switch Miguel's small angle correction

  GEO = addMissingGEO(GEO, SETTINGS.layshift, true, false);                % include GEO.to and force monotonic horizontal layers - maybe needs to be modified in a future
  if ~checkMonotonyGEO(GEO, SETTINGS.layshift/3.)
      GEO = makeMonotonicGEO(GEO, 1.0);
      monotonic = checkMonotonyGEO(GEO, SETTINGS.layshift/3.);
      if ~monotonic
          error("layer_collapse: small angle correciton & addMissingGEO() made x-non monotonic horizontal GEO interfaces unsolved by makeMonotonicGEO()")
      end
  end
end % function



