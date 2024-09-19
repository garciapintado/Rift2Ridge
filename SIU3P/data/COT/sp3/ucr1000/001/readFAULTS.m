function FAULTS = readFAULTS()
    global year
    myr = year*1.0E06;

    FAULTS(1).id     = 101;               % INTEGER 100 as minimum fault identifier to guarantee parent priority for layer collapse operations
    FAULTS(1).scoo   = [16.5 -4.4]*1000;  % [m] 'soft' x initial surface coordinates.
    FAULTS(1).dip    = 65.;               % [\deg] initial dip angle, from horizontal
    FAULTS(1).strike = 180.;              % [\deg] strike angle, clockwise from North, in {0,180}. In 2D fault 0: deepening to the right 180: deepening to the left
    FAULTS(1).atime  = 5.80*myr;          % [s] activation time for kinetic [via accumulated strain] fault effect
    FAULTS(1).dtime  = 6.80*myr;          % [s] deactivation time for kinetic [via accumulated strain] fault effect
    FAULTS(1).init   = true;              % empty to apply initFaultsLines() to this fault
    FAULTS(1).ai2pmaxb = [];              % if not empty, takes over PHY.FAUL.i2pmaxb at .atime
    FAULTS(1).di2pmaxb = 0.8;             % if not empty, takes over PHY.FAUL.i2pmaxb at .dtime
end
