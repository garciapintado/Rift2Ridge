% Loads plots and put them together (i.e. for power point slides)

%==========================================================================
% INPUT
%==========================================================================
% Directory
directory_p = '/home/miguel/Desktop/';
% Saved directory
savedirp = '/home/miguel/Desktop/ShFB_GR_C35nC_WSTemp';

% Folders or names of files to open
columns_p = {'ShFB_GR_C35nC_WSTemp'};
rows_p = {'eri','mu'};

pos_size = [0 0 1500 1000];

visiblep = 'visible';

%==========================================================================
% PLOT
%==========================================================================
length_test = [];
for l = 1:length(rows_p)
    for m = 1:length(columns_p)
        files_p{l,m} = dir([directory_p,columns_p{m},'/',rows_p{l},'*.fig']);
        length_test = max([length_test length(files_p{l,m})]);
    end
end

c = length(columns_p);
r = length(rows_p);

% Time loop
for l=1:length_test
    h = figure('visible','on');
    set(h,'Position',pos_size);
    saveas(h,[savedirp,'/dummy.fig'],'fig')
    % Loop for columns
    for m = 1:c
        % Loop for rows
        for n = 1:r
           % try
                % Plot non-scaled SS
                h1 = openfig([directory_p,columns_p{m},'/',files_p{m,n}(l+3).name],'new',visiblep);
                ax1 = gca;
                fig1 = get(ax1,'children');
                h = openfig([savedirp,'/dummy.fig']);
                s1 = subplot(c,r,c*(m-1)+n);
                copyobj(fig1,s1)
                close(h1)
                saveas(h,[savedirp,'/dummy.fig'],'fig')
           % end
        end
    end
    h = openfig([savedirp,'/dummy.fig']);
    %suptitle(['Strain rate at ',num2str(int_s*(k-1)),' Myr'])
    saveas(h,[savedirp,'/',num2str(l),'.png'],'png')
    close(h)
end