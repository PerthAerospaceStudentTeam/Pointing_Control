% Math Functions
classdef MathFunctions
    
    methods(Static)
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
    end
end