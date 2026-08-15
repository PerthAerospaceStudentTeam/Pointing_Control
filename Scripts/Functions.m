% Functions

classdef Functions
    
    methods(Static)
        % Creating parameters for control input function
        function M = Calc_M(x, y, theta)
            mathfunc = MathFunctions;
            A = x*transpose(y)+y*transpose(x)-transpose(x)*y*eye(3); %Checked
            M = [transpose(x)*y transpose(mathfunc.cross_prod(y)*x); (mathfunc.cross_prod(y)*x) (A)] - cos(theta)*eye(4); %Checked
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
            disp(psi1)
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
            mathfunc = MathFunctions;
            G = [transpose(P*w) ; transpose(qv)*w*eye(3)-mathfunc.cross_prod(P*w) ]; %Checked
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
    end
end