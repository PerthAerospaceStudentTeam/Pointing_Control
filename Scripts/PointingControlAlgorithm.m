%% Creating Control Law as per "Optimal Fixed-time Sliding Mode Control for 
% Spacecraft Constrained Reorientation"
%% Setting up Values
% Setting Variables for FIX-APF Algorithm
varpho11 = 1;
varpho2 = 0.02;
delta = 50;
dbar = 0.5;
k21 = 0.01;
k22 = 0.3;
k31 = 0.001;
k32 = 0.6;
a2 = 1.2;
b2 = 0.8;
a3 = 1.2;
b3 = 0.8;
theta_sun = 10*(pi/180); % In radians. Note n=1
theta_ant = 70*(pi/180); % In radians
e0 = 0.0001;

InertialParams
J = Inertia;

x1 = [1; 1; 1];
y1 = [1; 1; 1];
x2 = [-1; 1; 1];
y2 = [-1; 1; 1];
q = quaternion(1, 0, 0, 0);
qd = quaternion(1, 0, 0, 0);
qdot = quaternion(1, 0, 0, 0);
qv = quaternion(1, 0, 0, 0);
w = [1; 1; 1];
%% Setting Quaternions

% Note d represents gradient operator
% dd represents gradient squared
% dot  represents time derivative

func = Functions;
control_law_func = ControlLawFunction;

% Leave in only for testing
K21 = k21;
K22 = k22;
u = control_law_func.control_input(q, qdot, qv, qd, w, x1, y1, theta_sun, x2, y2, ...
            theta_ant, J, a2, a3, b2, b3, e0, K21, K22, k31, k32, delta, dbar ...
            , varpho11, varpho2)