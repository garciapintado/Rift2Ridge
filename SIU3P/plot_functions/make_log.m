ID = 654:729;

file_name = 'main';

for n = 1:length(ID)
    newFileData = fileread([file_name,num2str(ID(n)),'.m']);
    ii = strfind(newFileData,'directory = ');
    jj = strfind(newFileData,';'); jj  = jj(jj>ii); jj = jj(1);
    row_s = newFileData(ii:jj);
    k = strfind(row_s,'/'); k = k(end)+1;
    l = jj-ii-2;
    
    fields{n,1} = row_s(k:l);
    fields{n,3} = ID(n);
end