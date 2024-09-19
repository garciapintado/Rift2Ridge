function tecplot_all(dirp,nproc)
% TRANSFORM ALL THE TESTS IN A GIVEN FOLDER IN TECPLOT FILES

addpath(genpath('../main_functions'))

[data_dir,~,~] = load_dir(4);

tests = dir(dirp);

test_index = [];

for n = 1:length(tests)
    if ~strcmp(tests(n).name(1),'.')
        test_index = [test_index n];
    end
end

h = waitbar(0,'Converting files into tecplot format');

for n = 1:length(test_index)
    waitbar((n-1)/length(test_index))
    % Check worker availability
    cluster = parcluster;
    jobs = findJob(cluster);
    running = 0;
    for m = 1:length(jobs)
        if strcmp(jobs(m).State,'running')
            running = running+1;
        end
    end

    while running>=nproc
        pause(10)
        
        % Check worker availability
        cluster = parcluster;
        jobs = findJob(cluster);
        running = 0;
        for m = 1:length(jobs)
            if strcmp(jobs(m).State,'running')
                running = running+1;
            end
        end
    end
    batch('var2tec',0,{[dirp,'/',tests(test_index(n)).name], ...
        [data_dir,'/tecplot/'], ...
        {'Temp','E2all','Mu_all','Vel','I2.f','I2.p','I2.c'}})
end

delete(h)

h = msgbox('CONVERSION COMPLETED');