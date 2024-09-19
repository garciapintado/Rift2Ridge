function Pc = call_el2nod_pressure(GCOORD, ELEM2NODE, nel, Pd, method)
  % JGP: make kinedyn' el2nod_pressure() callable from rift2ridge2D

  if nargin < 5
    method = 'opt'; 
  end
    
  MODEL = "rift2ridge2D";
  
  MESH.GCOORD = GCOORD;
  MESH.EL2NOD = ELEM2NODE;
  MESH.nel    = nel;

  Pc = el2nod_pressure(MESH, Pd, MODEL, method);

  %  3                rift2ridge2D node convention
  %  | \
  %  5   4
  %  |     \
  %  1 - 6 - 2
  
end
