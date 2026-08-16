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

%% Setting Quaternions

%% Functions
% Note d represents gradient operator
% dd represents gradient squared
% dot  represents time derivative

% Creating sig^gamma function
function result = sig(x,gamma)
    result = norm(x).^gamma.*sign(x);
end

% Creating cross product function
function y_cross = cross_prod(y)
    y1 = y(1);
    y2 = y(2);
    y3 = y(3);
    y_cross = [0, -y3, y2; y3, 0, -y1; -y2, y1, 0];
end

% Creating quartenion transpose outer product function
function result = quat_trans_outer_prod(q, p)
    [p0, pv1, pv2 ,pv3] = parts(p);
    pv = [pv1; pv2; pv3];
    [q0, qv1, qv2 ,qv3] = parts(q);
    qv = [qv1; qv2; qv3];
    result = [q0*p0+transpose(qv)*pv; q0*pv-p0*qv-cross(qv, pv)];
end
   
% Creating control input function
function u = control_input(q, qdot, qv, w, x1, y1, theta_sun, x2, y2, ...
    theta_ant, J, a2, a3, b2, b3, e0, K21, K22, k31, k32)
    % These control the repulsion and attraction to forbidden and mandatory
    % zones respectfully (11 forbidden, 2 mandatory)
    M11 = Calc_M(x1, y1, theta_sun);
    M2 = Calc_M(x2, y2, theta_ant);
    [psi1, psi2] = Calc_psi(M11, M2, q);

    % Potentials
    Va = Calc_Va(q, qd);
    Vr = Calc_Vr(varpho11, varpho2, psi1, psi2, delta);
    Vp = Calc_Vp(Va, Vr);

    % Gradient of Potentials
    dVr = Calc_dVr(q, varpho11, varpho2, psi1, psi2, delta, M11, M2);
    dVp = Calc_dVp(Vr, Va, dVr, q, qd);

    % Laplacian of Potential
    ddVp = Calc_ddVp(q, qd, varpho11, varpho2, psi1, psi2, delta, M11, M2);
    

    P = Calc_P(Vr, dVr, q, qd);
    G = Calc_G(P, w, qv);
    E = Calc_E(P, qv);

    % Time Derivative of Potentials
    Vpdot = Calc_Vpdot(dVp, E, w);
    Vrdot = Calc_Vrdot(dVr, E, w);

    wtilda = Calc_wtilda(dVp);
    wtildadot = Calc_wtildadot(ddVp, E, w);
    fS = Calc_fS(S);
    dtilda = Calc_dtilda(E, J, dbar);
    
    alpha = Calc_alpha(k31, k32, a3, b3, Vp);
    beta = Calc_beta(wtilda, e0, alpha, qdot);
    S = Calc_S(qdot, alpha, beta);
    alphadot = Calc_alphadot(a3, b3, k31, k32, Vp, Vpdot);
    betadot = Calc_betadot(wtilda, wtildadot, e0, S);

    % E^+=Einv which is the generalised inverse of E
    Einv = pinv(E); %Checked

    u = cross(w, J*w) + J*Einv*(0.5*G*w - (Vrdot/Vr)*S +2*alphadot*beta + ...
        2*alpha*betadot-2*( Vr^((a2-1)/2)*K21*sig(S, a2) +  ...
        Vr^((b2-1)/2)*K22*sig(S, b2) ) + dtilda*fS); %Checked
    
end

% Creating parameters for control input function
function M = Calc_M(x, y, theta)
    A = x*transpose(y)+y*transpose(x)-transpose(x)*y1*eye(3); %Checked
    M = [transpose(x)*y transpose(cross_prod(y)*x); (cross_prod(y)*x) (A)] - cos(theta)*eye(4); %Checked
end

% psi1 and psi2
function [psi1, psi2] = Calc_psi(M11, M2, q)
    psi1 = -(transpose(q)*M11*q)/2; %Checked
    psi2 = (transpose(q)*M2*q)/2; %Checked
end

% Va
function Va = Calc_Va(q, qd)
    Va = norm(q-qd)^2; %Checked
end

% Vr
function Vr = Calc_Vr(varpho11, varpho2, psi1, psi2, delta)
    Vr = varpho11*exp((delta*psi1)^(-1)) + varpho2*exp((delta*psi2)^(-1)); %Checked
end
% Vp
function Vp = Calc_Vp(Va, Vr)
    Vp = Va+Va*Vr; %Checked
end

% dVr
function dVr = Calc_dVr(q, varpho11, varpho2, psi1, psi2, delta, M11, M2)
    dVr = ( ((delta*varpho11)/(delta*psi1)^2)*exp((delta*psi1)^(-1))*M11*q + ...
           ((-delta*varpho2)/(delta*psi2)^2)*exp((delta*psi2)^(-1))*M2*q ); %Checked
end

% dVp
function dVp = Calc_dVp(Vr, Va, dVr, q, qd)
    dVp = 2*(Vr+1)*(q-qd)+Va*dVr; %Checked 
end

% ddVp
function ddVp = Calc_ddVp(q, qd, varpho11, varpho2, psi1, psi2, delta, M11, M2)
    %l1-5 are the lines from top to bottom of equation (34)
    % being line 1 - l1, lines 2-3 - l2, line 4 -l3, lines 5-6 - l4,
    % lines 7-8 - l5
    l1 = 2*(varpho11*exp((delta*psi1)^(-1)) + varpho2*exp((delta*psi2)^(-1))+1)*eye(4); %Checked

    l2 = (2*delta*varpho11)/((delta*psi1)^2)*exp((delta*psi1)^-1)*( ...
        (q-qd)*transpose(q)*M11+M11*q*transpose(q-qd)); %Checked

    l3 = -(2*delta*varpho2)/((delta*psi2)^2)*exp((delta*psi2)^(-1))*( ...
        (q-qd)*transpose(q)*M2+M2*q*transpose(q-qd) ); %Checked

    l4 = norm(q-qd)^2*((delta*varpho11*exp((delta*psi1)^(-1)))/((delta*psi1)^4) * ( ...
        M11*(delta*psi1)^2+2*delta*M11*q*transpose(q)*M11*(delta*psi1) + ...
        delta*M11*q*transpose(q)*M11) ); %Checked

    l5 =  norm(q-qd)^2*((-delta*varpho2*exp((delta*psi2)^(-1)))/((delta*psi2)^4) * ( ...
        M2*(delta*psi2)^2-2*delta*M2*q*transpose(q)*M2*(delta*psi2) - ...
        delta*M2*q*transpose(q)*M2) ); %Checked
 
    ddVp = l1+l2+l3+l4+l5; %Checked
end

% Vpdot
function Vpdot = Calc_Vpdot(dVp, E, w)
    Vpdot = 1/2*transpose(dVp)*E*w; %Checked
    % Check if transpose(gradient(A))=gradient(transpose(A))
end

% Vrdot
function Vrdot = Calc_Vrdot(dVr, E, w)
    Vrdot = 1/2*transpose(dVr)*E*w; % Checked
    % Check same as for Vpdot function
end

% wtilda
function wtilda = Calc_wtilda(dVp)
    wtilda = -dVp; %Checked
end

% wtildadot
function wtildadot = Calc_wtildadot(ddVp, E, w)
    wtildadot = -1/2*ddVp*E*w; %Checked
end

% fS
function fS = Calc_fS(S)
    if norm(S) == 0
        fS = 0;
    else
        fS = S/norm(S);
    end
    %Checked
end

% dtilda
function dtilda = Calc_dtilda(E, J, dbar)
    dtilda = norm(E)*norm(inv(J))*dbar/2; %Checked
end

% P
function P = Calc_P(Vr, dVr, q, qd)
    P = 2*(Vr+1)*eye(4)+dVr*transpose(q-qd); %Checked
end

% G
function G = Calc_G(P, w, qv)
    G = [transpose(P*w) ; transpose(qv)*w*eye(3)-cross_prod(P*w) ]; %Checked
end

% E
function E = Calc_E(P, qv)
    E = [-transpose(qv); P]; %Checked
end

% alpha
function alpha = Calc_alpha(k31, k32, a3, b3, Vp)
    alpha = k31*Vp^(a3) + k32*Vp^(b3); %Checked
end

% beta
function beta = Calc_beta(wtilda, e0, alpha, qdot)
    Sbar = qdot-alpha*wtilda/norm(wtilda)^2; %Checked
    if (Sbar ~= 0) & (norm(wtilda) < e0)
        beta = ((2*e0^2-transpose(wtilda)*wtilda)/e0^4)*wtilda;
    else
        beta = wtilda/norm(wtilda)^2;
    end
    %Checked
end

% S
function S = Calc_S(qdot, alpha, beta)
    S = qdot -alpha*beta; %Checked
end

% alphadot
function alphadot = Calc_alphadot(a3, b3, k31, k32, Vp, Vpdot)
    alphadot = (a3*k31*Vp^(a3-1)+b3*k32*Vp^(b3-1))*Vpdot; %Checked
end

% betadot
function betadot = Calc_betadot(wtilda, wtildadot,e0, S)
    if (S ~= 0) & (norm(wtilda) < e0)
        betadot = (((2*e0^2-transpose(wtilda)*wtilda)*eye(4)-(2*wtilda*transpose(wtilda)))/(e0^4))*wtildadot;
    else
        betadot = ((norm(wtilda)^2-2*wtilda*transpose(wtilda))/(norm(wtilda)^4))*wtildadot;
    end
    %Checked
end