% Quaternion_Test

% Functions
function noise = randmon_noise
    noise=0.9+rand*(1.1-0.9);
end

% Definining Target and Current Orientations
target_orientation = quaternion(0,0.67,0,0.742);
current_orientation = quaternion(0,0.58,0.58,0.57);

% Defining Desired Manoeuvre
manoevre = (target_orientation)*(conj(current_orientation));
manoevre_noise = quaternion(randmon_noise(),randmon_noise(),randmon_noise(),randmon_noise());

% Extracting Individual Parts of Quaternion for Element Wise Multiplication
[mp1,mp2,mp3,mp4] = parts(manoevre);
[mpn1,mpn2,mpn3,mpn4] = parts(manoevre_noise);
manoevre_updated = normalize(quaternion(mp1*mpn1, mp2*mpn2, mp3*mpn3, mp4*mpn4));

% Not Relevant Anymore
% angle_rotation=pi/4;
% rotation_axis=[1 0 0];
% desired_change = quaternion(cos(angle_rotation),sin(angle_rotation)*rotation_axis(1),sin(angle_rotation)*rotation_axis(2),sin(angle_rotation)*rotation_axis(3));
% conj_desired_change = conj(desired_change);

% Calculating New Orientation
new_orientation = (manoevre_updated)*current_orientation
