%{
% Define symbolic variables for s and z
syms s z

% Continuous transfer function
H_s = 10 / (s + 10);
T = 0.05;

% a) Trapezoidal rule approximation
H_z = c2d(tf([10], [1 10]), T, 'tustin');

% b) Manual Forward difference approximation (s ≈ (z-1)/T)
s_forward = (z-1)/T;
H_z_forward = subs(H_s, s, s_forward);  % Substitute s in H_s

% c) Manual Backward difference approximation (s ≈ (z-1)/(T*z))
s_backward = (z-1)/(T*z);
H_z_backward = subs(H_s, s, s_backward);  % Substitute s in H_s

% Simplify the transfer functions
H_z_forward = simplify(H_z_forward);
H_z_backward = simplify(H_z_backward);

% Display results
disp('Trapezoidal Rule Approximation:')
disp(H_z)

disp('Forward Rule Approximation:')
disp(H_z_forward)

disp('Backward Rule Approximation:')
disp(H_z_backward)

%}


%{

% a) Define a discrete transfer function
% Coefficients of the numerator and denominator (ensure proper and causal)
num = [1, 0.5];  % Coefficients of the numerator
den = [1, 0];    % Coefficients of the denominator

% Create the discrete transfer function
H_z2 = tf(num, den);  % Define the discrete transfer function H(z)

% Display the discrete transfer function
disp('Discrete Transfer Function H(z):');
disp(H_z2);

% b) Calculate the pulse response of the system
[y2, t2] = step(H_z2);  % Step response of the discrete system

% c) Plotting the discrete transfer function and its pulse response
figure;

% Plot discrete transfer function response
subplot(2, 1, 1);
step(H_z2);  % Step response for the discrete transfer function
title('Pulse Response of the Discrete Transfer Function');
xlabel('Time (samples)');
ylabel('Output');

% Overlay the discrete transfer function on the same plot
hold on;
t = 0:T:5; % Time vector for plotting (adjust as needed)
u = ones(size(t));  % Unit step input
y = lsim(H_z2, u, t);  % Simulate the response
plot(t, y, 'r--', 'DisplayName', 'Discrete Transfer Function');
hold off;

legend('Pulse Response', 'Discrete Transfer Function');
grid on;

% Final adjustments
sgtitle('Question 2: Discrete Transfer Function and Pulse Response');
%}



% a) Define the discrete damped sinusoid
r = 0.5;
theta = pi/2;
k = 0:100;  % Time indices
ek = r.^k .* sin(k * theta);  % Damped sinusoid

% b) Plot the signal
figure;
subplot(2, 1, 1);
stem(k, ek, 'filled');
xlabel('k');
ylabel('e(k)');
title('Discrete Damped Sinusoid');
grid on;

% c) Find the Z-transform using symbolic variables
syms z;  % Define symbolic variable z
E_z = sum(ek .* z.^(-k));  % Compute Z-transform manually

% d) Find the pole/zero locations
[num, den] = numden(E_z);  % Separate into numerator and denominator
num = sym2poly(num);  % Convert symbolic to polynomial coefficients
den = sym2poly(den);

% Create the transfer function (no 'Variable' option)
H_z = tf(num, den);  % Create the transfer function

% Find pole/zero locations
[z, p, k] = tf2zp(num, den);

% e) Find the number of samples per cycle
omega = 2 * pi * theta;  % Angular frequency
T = 2 * pi / omega;  % Period
samples_per_cycle = T;  % Samples per cycle

% Print the results
disp('Z-transform:');
disp(E_z);

disp('Pole/Zero Locations:');
disp('Zeros:');
disp(z);
disp('Poles:');
disp(p);

disp('Number of Samples per Cycle:');
disp(samples_per_cycle);

% c) Plot the Z-transform response (discrete)
subplot(2, 1, 2);
step(H_z);  % Step response of the discrete transfer function
title('Response of the Z-Transform');
xlabel('Time (samples)');
ylabel('Output');
grid on;

% Add legend for clarity
legend('Z-Transform Response');


