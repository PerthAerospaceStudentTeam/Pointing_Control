% Rotation matrix from body frame to inertial frame using queler angles
function rotationmatrix = TIB(phi, theta, psi)

cph = cos(phi);
sph = sin(phi);

cth = cos(theta);
sth = sin(theta);

cps = cos(psi);
sps = sin(psi);


rt = [ cth*cps, sph*sth*cps - cph*sps, cph*sth*cps + sph*sps;
       cth*sps, sph*sth*sps + cph*cps, cph*sth*sps - sph*cps;
       -sth,    sph*cth,               cph*cth               ];


rotationmatrix = rt;