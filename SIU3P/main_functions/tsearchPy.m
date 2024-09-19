function xyel = tsearchPy(GCOORD, ELEM2NODE, xy, verify)
    % python-based tsearch. way faster than mutils::tsearch()
    % this function uses system files for testing.
    % If works fine; TODO could be faster by shared memory adaptation
    
    % GCOORD   : [2,ncoo]
    % ELEM2NODE: [3,nel] uint32
    % xy :       [2,:]
    % verify: whether scipy should verify input matlab file
    
    global pypath              % path to python script tsearchPy.py [defined in the global environment]

    if nargin < 4
        verify = 1;
    else
        verify = verify * 1; % cast into integer
    end
    sav_var = {'GCOORD', 'ELEM2NODE', 'xy'};
    fname = "tsearchPy_" + sprintf('%05d',floor(rand(1)*10000)); 
    fnamei = fname + "_i.mat";
    
    save(fnamei, sav_var{:});
    syscmd = join(["python3", fullfile(pypath,"tsearchPy.py"), fnamei, verify]," ");     % python tsearchPy()
    status = system(syscmd);
    
    fnameo =  fname + "_o.mat";
    load(fnameo);
    xyel = uint32(xyel + 1);                                                   % to 1-based matlab indexing
    
    delete(fnamei);
    delete(fnameo);

end % function