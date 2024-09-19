function matlab2pdf(function_name)

% MATLAB2PDF(FUNCTION_NAME) exports the function FUCNTION_NAME to a pdf
% file located in '~/Windows_Documents/Test_figures/pdf'

opts.outputDir = 'pdf';

opts.evalCode = false;
opts.format = 'pdf';

publish(function_name,opts)

system('mv pdf ~/Windows_Documents/Test_figures/')