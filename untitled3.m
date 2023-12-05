s = tf('s');
tau = 1e-4;
K = 3;
B = s*tau/(1+K*s*tau+(s*tau)^2);
T = -K*B;
nyquistplot(T)