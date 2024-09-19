function ltest = lastest(directory)

dir_w = 'directory = ';
ldir = length(dir_w)+1;

if ~strcmp(directory(1),'/')
    % Read main
    main = fileread([directory,'.m']);
    % Find directory = '
    dir_b = strfind(main, dir_w);
    dir_c = strfind(main,';');
    dir_s = dir_c>dir_b;
    dir_e = dir_c(find(dir_s,1,'first'))-2;
    
    % Read whole path
    directory = main(dir_b+ldir:dir_e);
end

files = dir(directory);

numbers =  zeros(length(files)-3,1);
tsteps = zeros(length(files)-3,1);

for n = 4:length(files)
    nmat = strfind(files(n).name,'.mat');
    name = files(n).name;
    tsteps(n-3) = str2double(name(2:nmat-1));
end

ltest = max(tsteps);

