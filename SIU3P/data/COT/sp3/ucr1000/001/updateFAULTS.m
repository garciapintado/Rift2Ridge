function FAULTS = updateFAULTS(FAULTS)
    % this file can overwrite faults given by readFAULTS() and add new faults

    global year
    myr = year*1.0E06;

    FAULTS(1).id    = 101;               % INTEGER 100 as minimum fault identifier to guarantee parent priority for layer collapse operations
    FAULTS(1).scoo   = [16.5 -4.4]*1000; % [m] 'soft' initial surface coordinates.
    FAULTS(1).dip    = 65.;              % [\deg] initial dip angle, from horizontal
    FAULTS(1).strike = 180.;             % [\deg] strike angle, clockwise from North, in {0,180}. In 2D fault 0: deepening to the right 180: deepening to the left
    FAULTS(1).atime  = 5.80*myr;         % [s]
    FAULTS(1).dtime  = 7.00*myr;         % [s] deactivation time for kinetic [via accumulated strain] fault effect
    FAULTS(1).init   = false;            % whether to apply initFaultsLines() to this fault 
    FAULTS(1).ai2pmaxb = 0.5;            % if not empty, takes over PHY.FAUL.i2pmaxb at .atime
    FAULTS(1).di2pmaxb = [];             % if not empty, takes over PHY.FAUL.i2pmaxb at .dtime
    
    FAULTS(2).id     = 102;              % INTEGER 100 as minimum fault identifier to guarantee parent priority for layer collapse operations
    FAULTS(2).scoo   = [4.5 -5.6]*1000;  % [m] 'soft' x initial surface coordinates. This is projected towards the surface by the given angle
    FAULTS(2).dip    = 65.;              % [\deg] initial dip angle, from horizontal
    FAULTS(2).strike = 0.;               % [\deg] strike angle, clockwise from North, in {0,180}. In 2D fault 0: deepening to the right 180: deepening to the left
    FAULTS(2).atime  = 6.20*myr;         % [s] activation time
    FAULTS(2).dtime  = 7.10*myr;         % [s] deactivation time
    FAULTS(2).init   = true;             % 
    FAULTS(2).ai2pmaxb = 0.5;            % if not empty, takes over PHY.FAUL.i2pmaxb at .atime
    FAULTS(2).di2pmaxb = 0.5;            % if not empty, takes over PHY.FAUL.i2pmaxb at .dtime
    
    FAULTS(3).id     = 103;              % INTEGER 100 as minimum fault identifier to guarantee parent priority for layer collapse operations
    FAULTS(3).scoo   = [-1.6 -5.87]*1000;% [m] 'soft' x initial surface coordinates. This is projected towards the surface by the given angle
    FAULTS(3).dip    = 60.;              % [\deg] initial dip angle, from horizontal
    FAULTS(3).strike = 180.;             % [\deg] strike angle, clockwise from North, in {0,180}. In 2D fault 0: deepening to the right 180: deepening to the left
    FAULTS(3).atime  = 7.00*myr;         % [s] activation time
    FAULTS(3).dtime  = 8.00*myr;         % [s] deactivation time
    FAULTS(3).init   = true;             % empty to apply initFaultsLines() to this fault
    FAULTS(3).ai2pmaxb = 0.5;            % if not empty, takes over PHY.FAUL.i2pmaxb at .atime
    FAULTS(3).di2pmaxb = 0.8;            % if not empty, takes over PHY.FAUL.i2pmaxb at .dtime
end
