% PLOT PALEO-WATER DEPTH
[ISOCHRONSp,~] = remesh_isoc(GCOORD,Point_id,ELEM2NODE,ISOCHRONS, ...
    Basement,tp_isoc);
[ISOCx,ISOCy_ini,ISOCyy,ISOCpwd,ero_bool]=find_ero(ISOCHRONSp,Topography);
[ISOCpwd_new] = modify_pwd(ISOCHRONSp,Topography);

np = size(ISOCx,2);
nt = size(ISOCx,1);
lines = 1:length(ISOCx(:));
lines = [lines(1:end-1)' lines(2:end)'];
del = np*(1:nt-1);
lines(del,:) = [];
Ix = ISOCx';
Iy = ISOCyy';

e2n_pwd = lines(1:find(lines(:,2)==np*(nt-1)),:);
up_line = np+1:np*nt;
up_line(np:np:end) = [];
e2n_pwd = [e2n_pwd up_line'];

e2n_pwd_u = lines(find(lines(:,1)==np+1):end,:);
down_line = 2:np*(nt-1);
down_line(np:np:end) = [];
e2n_pwd_u = [e2n_pwd_u down_line'];

e2n_pwd = [e2n_pwd; e2n_pwd_u];

P = ISOCpwd_new';


patch('faces',e2n_pwd,'vertices',[Ix(:) Iy(:)]/1000,'facevertexcdata', ...
P(:),'FaceColor','flat')
hold on
shading interp
isoc_type = 'isolines';
plot_isoc
colorbar






