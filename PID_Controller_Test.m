% PID Controller Test
s = tf('s');
P = 1/(s^2 + 10*s + 20);
step(P);

% Creating the PID
C = pid(300,70,10);
T=feedback(C*P,1);
t=0:0.01:2;
step(T,t);