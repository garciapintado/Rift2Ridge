ID = 882:995 ;

clear not_run not_finished

file_name = 'main';

[data_dir,~,~] = load_dir(4);
nr = 0;
nf = 0;
for nnn = 1:length(ID)
    disp([num2str(nnn),'/',num2str(length(ID))])
    newFileData = fileread([file_name,num2str(ID(nnn)),'.m']);
    ii = strfind(newFileData,'directory = ');
    jj = strfind(newFileData,';'); jj  = jj(jj>ii); jj = jj(1);
    row_s = newFileData(ii:jj);
    k = strfind(row_s,','); k = k(1)+2;
    l = jj-ii-2;
    
    dir_c = [data_dir,'/',row_s(k:l)];

    lt = lastest(dir_c);
    if isempty(lt)
        nr = nr+1;
        not_run(nr) = ID(nnn);
    else
        load([dir_c,'/_',num2str(lt)],'BREAKUP')
        if ~BREAKUP.bool
            nf = nf+1;
            not_finished(nf) = ID(nnn);
        end
    end
end