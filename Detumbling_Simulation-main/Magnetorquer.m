%%% Monitors the current being supplied to the magnetorquers and
%%% limits it to a specified maximum

function moment = Magnetorquer(angvel, magfieldbody)
    
    % get k value
    GainValue

    moment = k*cross(angvel, magfieldbody)/(norm(magfieldbody)^2);
    

