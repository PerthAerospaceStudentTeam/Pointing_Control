%% Terminal Initialisation %%
clear
clc
close all

disp('Simulation Started')

% Init Planet attributes
Planet
%

% Initial conditions for Position & Velocity
alt = 400e3; % 400km alt [metres]
inc = deg2rad(51.6);
sma = alt + earthRadius;
% 

% Start cubesat on the x-axis
x0 = alt + earthRadius;
y0 = 0;
z0 = 0;
position = [x0 y0 z0];
%

% Define the initial velocity of the cubesat
% Assume orbit is circular, so it is moving 
% around with a given speed at the given distance
% to earth
r = norm([x0; y0; z0]);
circ_vel = sqrt(earth_mu/r);
%

% Not moving in the X direction, 
% but is moving in the y and z
xdot0 = 0; 
ydot0 = circ_vel*cos(inc); 
zdot0 = circ_vel*sin(inc);
%

% Angular Velocity of Cubesat
phi0 = 0;
theta0 = 0;
psi0 = 0;
euler0 = [phi0 theta0 psi0];
quart0 = transpose(eul2quat(euler0));
%

% Initial angular velocity params in rads/s
p0 = 0*((2*pi)/360);
q0 = 0*((2*pi)/360);
r0 = 0*((2*pi)/360);
%

% Time Window
orbit_period = 2*pi*sqrt((sma^3)/(earth_mu));
num_orbits = 0.05;
tfinal = orbit_period*num_orbits;
tstep = 0.1; % seconds
tout = 0:tstep:tfinal;
%

% Magnetic Field
Bxout = 0*tout;
Byout = 0*tout;
Bzout = 0*tout;
%

% Initial Conditions for magnetorquer
moment = [0; 0; 0];
%

% Main State Initialisation
% Format: Column vector : position xyz (1:3);
%                         velocity in xyz directions (4:6);
%                         quarternion orientation (7:10);
%                         angular velocity (11:13);
state = [x0; y0; z0; xdot0; ydot0; zdot0; quart0; p0; q0; r0];
stateinitial = [x0; y0; z0; xdot0; ydot0; zdot0; quart0; p0; q0; r0];
stateout = zeros(length(tout), length(state));
%

% current_angular_position = [pi/2 0 0];
target_angular_position = [pi pi pi];
% ang_vel = [0 0 0];
total_error = [0 0 0];
max_ang_vel = 0.01;

% Numerical Integration
for idx = 1:length(tout)
    % Keep track of sim
    stateout=state;
    message = sprintf("%.1f / %.1f \n", tout(idx), max(tout));
    fprintf(message)

    % Get Magnetic field components
    magfieldcurrent = MagneticField(x0, y0, z0);
    magfieldbodycurrent = TBIquat(state(7:10))*(magfieldcurrent');
    
    % Plotting
    Bxout(idx) = magfieldcurrent(1);
    Byout(idx) = magfieldcurrent(2);
    Bzout(idx) = magfieldcurrent(3);
    momentout = zeros(length(tout), length(moment));

    % Get Moment
    moment = Magnetorquer(state(11:13), magfieldbodycurrent);
    momentout(idx, :) = transpose(moment);
    
    % Calculating Required Angular Acceleration
    current_angular_position=quat2eul(transpose(stateout(7:10)));
    angular_position_change = target_angular_position - current_angular_position;
    total_error = angular_position_change + total_error;
    angular_acceleration = angular_position_change;
    ang_vel=state(11:13);
    for i = 1:length(ang_vel)
        ang_vel(i) = ang_vel(i) + angular_acceleration(i)*tstep;
        if (ang_vel(i)) > max_ang_vel
            ang_vel(i) = max_ang_vel;
        end
        stateout(11:13)=ang_vel;
        current_angular_position = current_angular_position+transpose(ang_vel)*tstep;
        disp(current_angular_position)
        stateout(7:10)=eul2quat(current_angular_position);
        disp("Error is "+angular_position_change)
    end
end                                                                                                                                                 