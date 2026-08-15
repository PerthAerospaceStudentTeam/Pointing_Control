% rotation matrix from inertial frame to body frame using quatenions
function rotationmatrix = TBIquat(quaternion)

q0 = quaternion(1);
q1 = quaternion(2);
q2 = quaternion(3);
q3 = quaternion(4);

q0s = q0^2;
q1s = q1^2;
q2s = q2^2;
q3s = q3^2;

rt = [ q0s + q1s - q2s - q3s,  2*(q1*q2 + q0*q3),      2*(q1*q3 - q0*q2);
       2*(q1*q2 - q0*q3),      q0s - q1s + q2s - q3s,  2*(q0*q1 + q2*q3);
       2*(q0*q2 + q1*q3),      2*(q2*q3 - q0*q1),      q0s - q1s - q2s + q3s];

rotationmatrix = rt;