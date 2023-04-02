% Cleanup
clear all;
close all;
clc;

% Print the help
help spine
help fidelity

% Set the Larmor frequencies of the qubits (3 qubits)
f0 = [0.99e9, 1e9, 1.01e9];

% Set the Rabi frequencies (bit high, some RWA artifacts visible)
fR = [10e6, 10e6, 10e6];

% Generate the signal, driving the 2nd qubit (pi-rotation)
dt = 0.4e-10;
tpi = 0.5 / fR(2);
Npi = tpi / dt;
t = (1:Npi)*dt;
signal = cos(2*pi*f0(2)*t);

% Simulate and plot
[U, probabilities] = spine(fR, f0, dt, signal, 10);

% Print the probability of ending up in the spin up state
disp('Spin up probability for the 3 qubits:');
disp(squeeze(probabilities(3, end, :)));

% Print the fidelity of an identity operation/X-rotation
F = zeros(3, 1);
for n=1:3
    if (n == 2)
        % X-rotation
        Uideal = [0 1;
                  1 0];
        F(n) = fidelity(U(:, :, n), Uideal);
    else
        % Identity
        Uideal = eye(2);
        F(n) = fidelity(U(:, :, n), Uideal);
    end
end
disp('Infidelity (I/X/I-gates):');
disp(1-F);
