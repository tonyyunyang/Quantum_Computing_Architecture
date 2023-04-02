% Set the Larmor frequencies of the qubits (3 qubits)
f0 = [0.99e9, 1e9, 1.01e9];
% Set the Rabi frequencies
fR = [10e6, 10e6, 10e6];
% Generate the signal, driving the 2nd qubit (pi-rotation)
dt = 0.4e-10;
tpi = 0.5 / fR(2);
Npi = tpi / dt;
t = (1:Npi)*dt;
signal = cos(2*pi*f0(2)*t);
% Simulate and plot
[U, probabilities] = spine(fR, f0, dt, signal, 10);