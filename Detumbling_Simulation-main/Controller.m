%%% System controller that produces a torque with magnetorquers to
%%% detumble a cubesat using the equation:

%%% torque = moment x magnetic_field(body frame)
%%% where: moment = k*(ang_velocity x magnetic_field(body frame))
%%% therefore:
%%% torque = (k*(ang_velocity x magnetic_field(body frame))) x magnetic_field(body frame)

function torque = Controller(magfieldbody, angvel)


    % optimal k value found using k = 2n(1+sin(i))*Imin
    GainValue

    moment = k*cross(angvel, magfieldbody)/(norm(magfieldbody)^2);


torque = cross(moment, magfieldbody);

