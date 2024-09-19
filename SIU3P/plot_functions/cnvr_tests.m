function cnvr_tests(directory)

tests   = dir(directory);
tests   = tests(3:end);

bars_c = [];

for n = 1:length(tests)
    clc
    disp(['Reading data: ',num2str(n),'/',num2str(length(tests))])
    
    dir_steps   = dir([directory,tests(n).name]);
    dir_steps   = dir_steps(3:end);
    iter_max    = lastest([directory,tests(n).name]);
    
    tPic        = 0;
    tNew        = 0;
    c           = 0;
    for m = 1:length(dir_steps)
        point   = strfind(dir_steps(m).name,'.');
        if point>2
            try
                load([directory,tests(n).name,'/',dir_steps(m).name], ...
                    'residue')
                tPic    = tPic + sum(residue(2,:)==1);
                tNew    = tNew + sum(residue(2,:)==2);
                c       =  c+1;
            end
        end
    end
    aPic    = tPic/c;
    aNew    = tNew/c;
    
    bars_c      = [bars_c; aPic aNew];
    lowerb      = strfind(tests(n).name,'_');
    name_c      = tests(n).name;
    name_c(lowerb)  = ' ';
    xlabel_c{n} = name_c;
end

% Plot bars
bar(bars_c,'stacked')
set(gca,'XtickLabel',xlabel_c)
set(gca,'XTickLabelRotation',45)