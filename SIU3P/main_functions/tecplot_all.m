function tecplot_all(dirp,nproc)
% TRANSFORM ALL THE TESTS IN A GIVEN FOLDER IN TECPLOT FILES

addpath(genpath('main_functions'))
addpath(genpath('plot_functions'))

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
    if n>1
        finished = sum(strcmp({j(:).State},'finished'));
        waitbar(finished/length(test_index))
    end
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
    j(n) = batch('var2tec',0,{[dirp,'/',tests(test_index(n)).name], ...
        [data_dir,'/tecplot/'], ...
        {'Temp','E2all','ErP','ErC','Mu_all','Vel','I2.f','I2.p','I2.c','RHO'}});
end

while 1
    % Check worker availability
    cluster = parcluster;
    jobs = findJob(cluster);
    running = 0;
    for m = 1:length(jobs)
        if strcmp(jobs(m).State,'running')
            running = running+1;
        end
    end
    
    if running==0
        break
    end
    finished = sum(strcmp({j(:).State},'finished'));
    waitbar(finished/length(test_index))
    
    pause(10)
end

delete(h)
delete(j)

h = msgbox('CONVERSION COMPLETED');