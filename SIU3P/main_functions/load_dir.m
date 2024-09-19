function [MUTILS_PATH, TRIANGLE_PATH, SSPARSE_PATH] = load_dir(n)
  % FIND WHICH COMPUTER THE TEST IS RUNNING AND SET UP THE DIRECTORIES
  % 
  % Javier GP: this has been modified from locally defined paths (i.e. within this function) to
  % environment variable-defined paths.
  % The following environment variables need now to be externally defined in the system 
  % (e.g in .bash_profile for linux/mac users using bash shell) and accessible to matlab:
  %
  % variable         example (Javier's iMAC)                               
  % MUTILS_PATH   :: /Users/jgp/Library/MatlabLib/mutils-0.4-2
  % TRIANGLE_PATH :: /Users/jgp/Library/MatlabLib/triangle
  % SSPARSE_PATH  :: /Users/jgp/Library/MatlabLib/SuiteSparse
  %
  % If for some undesired reason, the user is not able to provide matlab with externally defined
  % environment variables, the specific case should be added here. As an example,
  % a specific case is defined for the cluster cluster computer "GLNXA64"

  HOME    = getenv('HOME');
  LOGNAME = getenv('LOGNAME');

  switch computer
    case "GLNXA64"           % cluster cluster [9 non-interconnected nodes, ppn:40]
        MUTILS_PATH   = "/llocal1/matlabr2019a/toolbox/mutils-0.4-2";
        TRIANGLE_PATH = [];                 % triangle binary already in path
        SSPARSE_PATH  = "/llocal1/SuiteSparse-5.6.0";
        return
    otherwise                % general case; e.g. "MACI64" => new iMAC desktops [Javier & Leila]
        MUTILS_PATH   = getenv('MUTILS_PATH');
        TRIANGLE_PATH = getenv('TRIANGLE_PATH');
        SSPARSE_PATH  = getenv('SSPARSE_PATH');
        return
    end % switch
end % function
