%%% Returns a angular velocity reading that mimics the jumpy
%%% and discrete nature of actual sensor readings

%%% first functionality: discrete readings, in this function, a reading is
%%% taken every n iterations (n*4 because the satellite function gets 
%%% called 4 times each iteration because of RK4 integration), 
%%% mimicking time when the cubesat is moving
%%% but not being able to read, so values may jump and not be as smooth
function reading = GryoSensorModule(angvel)

    %message = sprintf("%d, %d, %d", magfield_x_ideal, magfield_y_ideal, magfield_z_ideal);

    %%% this section takes a new value when the limiter has elapsed
    %%% but holds onto that value and keeps using it until a new one
    %%% is requested
    persistent limiter angular_velocity

    if isempty(angular_velocity)
        angular_velocity = angvel;
    end

    if isempty(limiter)
        limiter = 0;
    end

    if limiter >= 0
        limiter = limiter + 1;
    end

    %%% mimicking the dicrete function call by limiting how often the
    %%% sensor can get a new reading
    frequency = 10*(4);
    if mod(limiter, frequency) == 0
        angular_velocity = angvel; % columnn
        
        %%% adds in random error, it could be plus or minus
        %%% the error added below at random
        n = 20;
        % 10 percent error where n = 1
        % total error = n*10
        randomiser_x = randi([-1000*n 1000*n])/10000;
        randomiser_y = randi([-1000*n 1000*n])/10000;
        randomiser_z = randi([-1000*n 1000*n])/10000;

        bias = [2; 2; 2];

        randomiser_matrix = [randomiser_x; randomiser_y; randomiser_z] + bias; % column

        modifier_matrix = randomiser_matrix.*angular_velocity; % column

        angular_velocity = (angular_velocity + modifier_matrix); % column

    end 

    %magfieldinertial = MagneticField(x, y, z)';

reading = angular_velocity;



