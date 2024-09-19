% matlab -nosplash -nodisplay                                                      # command line - interactive
% matlab -nosplash -nodisplay -batch "pid=0;call_main" -logfile rift2ridge2D.log   # batch mode

HOME    = getenv('HOME');
MODEL   = "rift2ridge2D";
VERSION = "SIU3P";

region  = "COT";
event   = "sp3";
meshnam = "ucr1000";
nnn     = "001";               % unique test identifier: [aaa] format

dsnmod = fullfile(HOME,"docs",MODEL,VERSION);
addpath(dsnmod);

loadsave = false;
fname = "last";

main;
