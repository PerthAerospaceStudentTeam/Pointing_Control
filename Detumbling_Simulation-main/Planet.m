% Earth Attributes %

R = earthRadius; % metres
M = 5.972e24; % kg
G = 6.67e-11;
earth_mu = 3.986e14;

% for Earth plotting %
[Xe, Ye, Ze] = sphere;
Xe = Xe*R;
Ye = Ye*R;
Ze = Ze*R;
