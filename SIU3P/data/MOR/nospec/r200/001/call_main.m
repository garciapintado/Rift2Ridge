% matlab -nosplash -nodisplay                                                      # command line - interactive
% matlab -nosplash -nodisplay -batch "pid=0;call_rift2ridge2D" -logfile siu3p.MOR.nospec.r200.log   # batch mode

HOME    = getenv('HOME');
MODEL   = "rift2ridge2D";
VERSION = "SIU3P";

region  = "MOR";
event   = "nospec";
meshnam = "r200";
nnn     = "001";               % unique test identifier: [aaa] format

dsnmod = fullfile(HOME,"docs",MODEL,VERSION);
addpath(dsnmod);

loadsave = true;
fname = "last";                                                            % TRUE  : load a previous save and run
% fname = "MOR.r200.000.00403.mat";

main;
