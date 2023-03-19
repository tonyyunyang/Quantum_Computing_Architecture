% Cleanup
clear all;
close all;
clc;

clear;
close all;


% define input signal
% input amplitude 
A=[10e-3 0 0]; %[V]
% input frequency 
fin=[1 2 3]*1e9; %[Hz]

%define microwave frequency
fosc=10e9; %[Hz]

% Set the Larmor frequencies of the qubits (3 qubits)
f0 = [11e9, 12e9, 13e9];

% Set the Rabi frequencies (bit high, some RWA artifacts visible)
fR = [1e6, 1e6, 1e6];

% total simulation time. It is comnputed as the time needed for a pi
% rotation for Q1
T=1/fR(1)/2; %[s]
% simulation time step
dt=1e-12; %[s]
% define time vector
t=dt:dt:T; %[s]

% input signal
Vin=A(1)*cos(2*pi*fin(1)*t)+...
    A(2)*cos(2*pi*fin(2)*t)+...
    A(3)*cos(2*pi*fin(3)*t); %[V]

G1=200; %[-]

% these are the two values determined from the previous question
Gtemp1 = 805.7;
Gtemp2 = 161162;

Vout1=G1*Vin + Gtemp1 * Vin.^2 + Gtemp2 * Vin.^3;
 
Vmix=Vout1.*cos(2*pi*t*fosc);

G2=1; %[-]

% we modify this value as well as the equation for Vout2 below to see
% whether 2nd and 3rd order distortion of AMPL2 has effect
Gtemp3 = 0;
Gtemp4 = 0.5;

Vout2=G2*Vmix + Gtemp3 * Vmix.^2 + Gtemp4 * Vmix.^3;


%compute fft
X=abs(fft(Vout2));
XdB=20*log10(X); %[dB]
Y=abs(fft(Vin));
YdB=20*log10(Y); %[dB]
% FFT resolution
fres=1/T; %[Hz]
% maximum FFT frequency
fmax=1/dt; %[Hz]
% frequency vector
f=0:fres:(fmax-fres); %[Hz]

%plot FFT
figure
semilogx(f,XdB,f,YdB)
hold off
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
legend('V_{out2}','V_{in}')
xlim([0.1 14]*1e9)
grid


% Simulate and plot
[U, probabilities] = spine(fR, f0, dt, Vout2, 10);


% Compute fidelity of the operation on qubit 1
% Ideal transformation representing a pi-rotation around X
 Uideal = [0 1;
           1 0];
% Compute fidelity by comparing the actual transformation to the ideal one
F1 = fidelity(U(:, :, 1), Uideal);

% For the other 2 qubits, we would like to compare the actual transormation
% to the identity operation (I). However, any leakage form the driving
% signal is going to induce Z-rotations on the qubits. Such rotations can
% be easily cancelled by software, i.e. by not applying any extra
% electrical pulse. This correction is currently not implemented in this
% code to make it easier for the students. However, for the simple case in
% which the qubits ar einitialized in |0>, the fidelity on idel qubits can
% be simply computed by measuring the qubit along Z. This measurmenet is
% insensitive with respect to rotazion around Z.
F2=probabilities(3,end,2);
F3=probabilities(3,end,3);


disp(['Fidelity of qubit 1: ' num2str(F1*100) '%']);
disp(['Fidelity of qubit 2: ' num2str(F2*100) '%']);
disp(['Fidelity of qubit 3: ' num2str(F3*100) '%']);