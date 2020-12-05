function  [nq,nv] = jinfo( jtyp)

if ischar( jtyp )
  code = jtyp;
else
  code = jtyp.code;
end

switch code
  case {'Rx','Ry','Rz','R','Px','Py','Pz','P','px','py','r','H'}
    nq = 1;
    nv = 1;
  case {'S'} % spherical
    nq = 4;
    nv = 3;
  case 'Fb'
    nq = 7;
    nv = 6;
  otherwise
    error( 'unrecognised joint code ''%s''', code );
end
