function [GCOORD, ELEM2NODE, Point_id, Phase_id] = generate_meshGEO(GEO, Elsizes, triangle_mode, triangle_call, meshname)
  %  function [GCOORD, ELEM2NODE, Point_id, Phase_id] = generate_meshGEO(GEO, Elsizes, triangle_mode, triangle_call, meshname)
  %
  %  Create a triangular mesh to be used by rift2ridge2D. This functions differs from MILAMIN versions and mesh management as: 
  %  - The input has to be structured as as GEO struct [see getGEO()]. 
  %  - The flags to triangle are such that they are coordinates to work with a GEO struct. In GEO, overlapping interfaces a repeated for each of the touching regions.
  %  - The output Point_id follow the rules in GEO, such that higher order polylines in shared GEO interfaces have strict priority in Point_ID. This eases hierarchical mesh operations.
  %
  %  INPUT
  %  ---
  %  GEO           :: STRUCT defining mesh
  %  Elsizes       :: REAL, DIM(nboxes), array defining el size [m2] for each box (i.e. phase)
  %  triangle_mode :: STRING in {"binary","text"}
  %  triangle_call :: STRING indicating the triangle executable file, including path
  %  meshname      :: STRING indicating the triangle output basename in system files [$meshname.ele, $meshname.poly, $meshname.node]
  %   
  %  OUTPUT
  %  ---
  %  GCOORD    :: REAL, DIM(2,nnod) 2D coordinates of the mesh nodes
  %  ELEM2NODE :: INTEGER, DIM(6,nel) node indices at each element
  %  Point_id  :: INTEGER, DIM(1,nnod) node class
  %  Phase_id  :: INTEGER, DIM(1,nel) element class
  %
  %  where 'nnod' is the number of nodes, and 'nel' the number of elements 
  % 
  %  The convention here is to follow the output as given by triangle. For comparison:
  %
  %                                                                                    in any of the two conventions: 
  %                                                                                    |  \ 3 IPs     |     \ 6 IPs
  %  3                   triangle convention [e.g. rift2ridge2D].                      |   \          |      \
  %  | \                 start by local node 1. Counterclockwise triangles and nodes   |    \         | 5.    \
  %  5  4                                                                              | .3  \        |        \
  %  |    \                                                                            |      \       | 1.  3.  \
  %  1--6--2                                                                           |       \      |          \
  %  3                   Burges convention [e.g. Kinedyn]                              | .1  .2 \     | 6.  2.  4.\
  %  | \                                                                               |_________\    |____________\
  %  6  5
  %  |   \
  %  1--4--2

  %--------------------------------------------------------------------------
  %
  % Author: Javier Garcia-Pintado: 2020-02-20
  %
  %--------------------------------------------------------------------------



  nodelist    = [];
  nodemarker  = [];
  segmentlist = [];
  segmmarker  = [];
  regionlist  = [];
  holelist    = [];
                                
  ngeo = length(GEO);
  % intbs = max(1, ((1:nboxes)-1)*3);
  % intrs = (1:nboxes)*3 - 1; 
  % intts = intrs + 1
  % intls = intts + 1

  nboxes = floor(ngeo/3);
  ptsu = 0;
  
  for ib=1:nboxes
    intb = max(1, (ib-1)*3);                                              % bottom interface
    intr = ib*3 - 1;                                                      % right    " "
    intt = intr + 1;                                                      % top      " "
    intl = intt + 1;                                                      % left     " "
    BOX_n = [GEO(intb).coo GEO(intr).coo fliplr(GEO(intt).coo) fliplr(GEO(intl).coo)]; % BOX nodes
     
    Id_nodb = max(GEO(intb).to(2,:), GEO(intb).gid);                      % max(...) as .to can only point to higher order interfaces in GEO [and is 0 when no pointer exists for the node]
    Id_nodb(1)   = GEO(intl).gid;                                         % bottom to left interface
    Id_nodb(end) = GEO(intr).gid;                                         %   "     " right  " "
    Id_nodr = max(GEO(intr).to(2,:), GEO(intr).gid);
    Id_nodt = max(flip(GEO(intt).to(2,:)), GEO(intt).gid);
    Id_nodt(end) = GEO(intl).gid;
    Id_nodl = max(flip(GEO(intl).to(2,:)), GEO(intl).gid);
    Id_nod = [Id_nodb Id_nodr Id_nodt Id_nodl];
    Id_seg = min(Id_nod, Id_nod([2:end,1]));
     
    Id_nodb([1,end]) = Id_nodb(2); % set corners back to their gid - GEO compliant
    Id_nodt([1,end]) = Id_nodt(2); %  "     "     "        "
    Id_nod = [Id_nodb Id_nodr Id_nodt Id_nodl];

    npts = length(BOX_n);
    ptsl = ptsu + 1;
    ptsu = ptsl + npts - 1;
    BOX_s = [ptsl:ptsu; 1+(ptsl:ptsu)];                                   % BOX segments
    BOX_s(2,end) = ptsl;                                                  % close box
    
    % BOX regional attributes. Here I take 2 samples per box to improve assignation for broken layers
    if isfield(GEO,'region_samplers')
        ns = size(GEO(intb).region_samplers,2);
        BOX_r = [GEO(intb).region_samplers; repmat([ib; Elsizes(ib)],1,ns)];
    else
        BOX_r = [[GEO(intb).coo(:,1)+1.e-02;   ib; Elsizes(ib)] ...      % [4,nsamples=2] sample 1: lower left 
            [GEO(intt).coo(:,end)-1.e-02; ib; Elsizes(ib)]];             %                sample 2: upper right
    end

    nodelist    = [nodelist    BOX_n];
    nodemarker  = [nodemarker  Id_nod];
    segmentlist = [segmentlist BOX_s];
    segmmarker  = [segmmarker  Id_seg];
    regionlist  = [regionlist  BOX_r];
  end
  %  segm_xc = 0.5 * (nodelist(1,segmentlist(1,:)) + nodelist(1,segmentlist(2,:)));
  %  segm_yc = 0.5 * (nodelist(2,segmentlist(1,:)) + nodelist(2,segmentlist(2,:)));
  %  hold on; text(segm_xc/1000, segm_yc/1000 + .05, sprintfc('%d',segmmarker),'Fontsize',20);
  %  hold on; text(nodelist(1,:)/1000, nodelist(2,:)/1000+.05, sprintfc('%d',nodemarker),'Fontsize',20,'color',[0.8,0.0,0.0])
  
  %TOTAL
  n_nod = size(nodelist,2);
  n_seg = size(segmentlist,2);
  n_reg = size(regionlist,2);
  n_hol = size(holelist,2);

  nodelist    = [1:n_nod; nodelist; double(nodemarker)];
  segmentlist  = [1:n_seg; segmentlist; segmmarker];
  regionlist   = [1:n_reg; regionlist];
  holelist     = [1:n_hol; holelist];

  if exist(meshname + ".ele") == 2
      delete(meshname + ".*")                                                % avoid matlab continuation in case triangle fails and previous files exist
  end
  
  fid     = fopen(meshname + ".poly_input", 'w');                            % write triangle input file [a .poly file: Planar Straight Line Graph]
  fprintf(fid,'%d 2 0 1\n', n_nod);
  fprintf(fid,'%d %15.12f %15.12f %d\n', nodelist);
  fprintf(fid,'%d 1\n', n_seg);
  fprintf(fid,'%d %d %d %d\n', segmentlist);
  fprintf(fid,'%d\n',n_hol);
  if n_hol > 0
    fprintf(fid,'%d %15.12f %15.12f\n', holelist);
  end
  fprintf(fid,'%d 0\n', n_reg);
  fprintf(fid,'%d %15.12f %15.12f %d %15.12f\n', regionlist);
  fclose(fid);

  % area_glob = 1e7; %(2*pi*radius/no_pts_incl)^2;

  if triangle_mode == "binary"
      binary_flag = "-b";
  else
      binary_flag = "";
  end

  % flag -j ::  delete repeated vertexes
  % system(['triangle ',binary_flag,' -pQIo2q33Aa',num2str(area_glob,'%12.12f'),' ',modelname,'.poly']);
  %system(join(["cp",meshname + ".poly_input",meshname + ".poly"]));        % preserve input for debugging
  copyfile(meshname + ".poly_input", meshname + ".poly");
  system(join([triangle_call, binary_flag, "-pQIo2q33FAaj", meshname + ".poly"])); % (-p: input is a .poly file) writes *.ele, *.node, *.poly
                                                                                  % output: a conforming constrained Delaunay triangulation, with holes and concavities removed
										       
  if triangle_mode == "binary"
      %NODES READING
      fid=fopen(meshname + ".node", 'r');
      fseek(fid, 0, 1);
      file_size	= ftell(fid);
      fseek(fid, 0, -1);
      dummy	= fread(fid,file_size/8,'double');
      fclose(fid);

      GCOORD		= [dummy(6:4:end)';dummy(7:4:end)'];
      Point_id	= dummy(8:4:end)';

      %ELEMS READING
      fid=fopen(meshname + ".ele", 'r');
      fseek(fid, 0, 1);
      file_size	= ftell(fid);
      fseek(fid, 0, -1);
      dummy	= fread(fid,file_size/8,'double');
      fclose(fid);

      ELEM2NODE	= [dummy(5:8:end)';dummy(6:8:end)';dummy(7:8:end)';dummy(8:8:end)';dummy(9:8:end)';dummy(10:8:end)'];
      ELEM2NODE	= int32(ELEM2NODE);
      Phase_id		= dummy(11:8:end)';

  else % text triangle output

      %NODES READING
      fid = fopen(meshname + ".node", 'r');
      tmp = fscanf(fid, '%d',4);
      nnod = tmp(1);
      GCOORD = fscanf(fid, '%e', [4, nnod]);
      fclose(fid);

      Point_id = GCOORD(end,:);          % node idenfication
      GCOORD(1,:)   = [];                % just leave 2D coordinates
      GCOORD(end,:) = [];

      %ELEMS READING
      fid =fopen(meshname + ".ele", 'r');
      tmp = fscanf(fid, '%d',3);
      nel = tmp(1);
      ELEM2NODE = fscanf(fid, '%d',[8, nel]);
      fclose(fid);
    
      Phase_id  = ELEM2NODE(end,:);
      ELEM2NODE(1,:)   = [];
      ELEM2NODE(end,:) = [];
      ELEM2NODE	= int32(ELEM2NODE);
  end
end % function
