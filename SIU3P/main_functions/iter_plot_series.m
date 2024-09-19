% Iterate plot series
plot_opt

SETp.basement = 1;
SETp.line_widthb = 2;
SETp.color_base = [0 0 1];

inp = {'GR_C35_RHOT_3L_c1e-2','GR_C35_RHOT_3L_c1' ...
    'GR_C40_RHOT_3L_c1e-2','GR_C40_RHOT_3L_c1' ...
    'AN_C35_RHOT_3L_c1e-2','AN_C35_RHOT_3L_c1' ...
    'AN_C40_RHOT_3L_c1e-2','AN_C40_RHOT_3L_c1'};

for n = 1:length(inp)
    try
    plot_series(['/data/miguel/SEDIMENTS/',inp{n}],'_', ...
        'eri',['Figures/Review/',inp{n}],'png',1:100:10001,SETp)
    end
end

SETp.basement = 0;

inp = {'GR_C35_RHOT_3L','GR_C40_RHOT_3L', ...
    'AN_C35_RHOT_3L','AN_C40_RHOT_3L'};

for n = 1:length(inp)
    try
    plot_series(['/data/miguel/SEDIMENTS/',inp{n}],'_', ...
        'eri',['Figures/Review/',inp{n}],'png',1:100:10001,SETp)
    end
end