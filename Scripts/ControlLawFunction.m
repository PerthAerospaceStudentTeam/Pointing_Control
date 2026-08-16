% Control Law
% Creating control input function
classdef ControlLawFunction
    methods(Static)
        function u = control_input(q, qdot, qv, qd, w, x1, y1, theta_sun, x2, y2, ...
            theta_ant, J, a2, a3, b2, b3, e0, K21, K22, k31, k32, delta, dbar ...
            , varpho11, varpho2)
            mathfunc = MathFunctions;
            func = Functions;    
            
            % These control the repulsion and attraction to forbidden and mandatory
            % zones respectfully (11 forbidden, 2 mandatory)
            M11 = func.Calc_M(x1, y1, theta_sun);
            M2 = func.Calc_M(x2, y2, theta_ant);
            [psi1, psi2] = func.Calc_psi(M11, M2, q);
       
            % Potentials
            Va = func.Calc_Va(q, qd);
            Vr = func.Calc_Vr(varpho11, varpho2, psi1, psi2, delta);
            Vp = func.Calc_Vp(Va, Vr);
        
            % Gradient of Potentials
            dVr = func.Calc_dVr(q, varpho11, varpho2, psi1, psi2, delta, M11, M2);
            dVp = func.Calc_dVp(Vr, Va, dVr, q, qd);
        
            % Laplacian of Potential
            ddVp = func.Calc_ddVp(q, qd, varpho11, varpho2, psi1, psi2, delta, M11, M2);
            
        
            P = func.Calc_P(Vr, dVr, q, qd);
            G = func.Calc_G(P, w, qv);
            E = func.Calc_E(P, qv);
        
            % Time Derivative of Potentials
            Vpdot = func.Calc_Vpdot(dVp, E, w);
            Vrdot = func.Calc_Vrdot(dVr, E, w);
        
            wtilda = func.Calc_wtilda(dVp);
            wtildadot = func.Calc_wtildadot(ddVp, E, w);
            fS = func.Calc_fS(S);
            dtilda = func.Calc_dtilda(E, J, dbar);
            
            alpha = func.Calc_alpha(k31, k32, a3, b3, Vp);
            beta = func.Calc_beta(wtilda, e0, alpha, qdot);
            S = func.Calc_S(qdot, alpha, beta);
            alphadot = func.Calc_alphadot(a3, b3, k31, k32, Vp, Vpdot);
            betadot = func.Calc_betadot(wtilda, wtildadot, e0, S);
        
            % E^+=Einv which is the generalised inverse of E
            Einv = pinv(E); %Checked
        
            u = cross(w, J*w) + J*Einv*(0.5*G*w - (Vrdot/Vr)*S +2*alphadot*beta + ...
                2*alpha*betadot-2*( Vr^((a2-1)/2)*K21*mathfunc.sig(S, a2) +  ...
                Vr^((b2-1)/2)*K22*mathfunc.sig(S, b2) ) + dtilda*fS); %Checked
            
        end
    end
end