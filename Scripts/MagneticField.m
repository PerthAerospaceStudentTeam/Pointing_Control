% This function returns the magnetic field components in the inertial frame

function B = MagneticField(x, y, z)
    current_date = '01-Jan-2020';   
    coord = 'geocentric';
    rho = norm([x, y, z]);
    theta_E = acos(z/rho);
    psi_E = atan2(y, x);
    phi_E = 0;
    
    latitude = 90 - theta_E*(180/pi);
    longitude = psi_E*(180/pi);
    altitude = rho/1000;

    % get magnetic field vectors %
    [Bn, Be, Bd] = igrf(current_date, latitude, longitude, altitude, coord);

    BNED = [Bn; Be; -Bd];
    BI = TIB(phi_E, theta_E+pi, psi_E)*BNED;
    Bx = BI(1)*1e-9;
    By = BI(2)*1e-9;
    Bz = BI(3)*1e-9;

B = [Bx By Bz];