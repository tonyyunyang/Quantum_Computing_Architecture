% Cleanup
clear all;
close all;
clc;

clear;
close all;

% Set the Larmor frequencies of the qubits (3 qubits)
f0 = [5.1e9, 4.9e9, 0];

% Set the Rabi frequencies
fR = [1e6, 1e6, 0];

%define microwave frequency
fosc=5e9; %[Hz]


%number of bits of the DAC
N=5;

%sample frequency of the DAC
fsample=1e9;

% define input signal
% input amplitude 
A=[10e-3 0 0]; %[V]
% input frequency 
fin=[-0.1 0 0]*1e9; %[Hz]
% inout phases
ph_in=[-0.32 0 0]; %[rad]

% total simulation time. It is computed as the time needed for a pi
% rotation for Q1
T=1/fR(1)/2; %[s]
% simulation time step
dt=1e-12; %[s]
% define time vector
t=dt:dt:T; %[s]


% input signal
% Am = 1;          % Amplitude
% mu = 0;         % Mean
% sigma = 0.1;    % Standard deviation
% in = Am * gaussmf(t, [sigma*sqrt(2*pi) mu]);

% in=A(1)*cos(2*pi*fin(1)*t+ph_in(1))+...
%     A(2)*cos(2*pi*fin(2)*t+ph_in(2))+...
%     A(3)*cos(2*pi*fin(2)*t+ph_in(3)); %[V]

% Use the Hilbert transform to extract the I and Q components
% s_hilbert = hilbert(in);
% I = real(s_hilbert);
% Q = imag(s_hilbert);

% s_hilbert = hilbert(in);
I = A(1)*cos(2*pi*fin(1)*t+ph_in(1)) + A(2)*cos(2*pi*fin(2)*t+ph_in(2));
Q = A(1)*sin(2*pi*fin(1)*t+ph_in(1)) + A(2)*sin(2*pi*fin(2)*t+ph_in(2));

%sampling the input signal I
I_sampled=I(1:1/(fsample*dt):end);

%sampling the input signal Q
Q_sampled=Q(1:1/(fsample*dt):end);

%quantize the signal on N-bit
%first, scale the signal to fit in [0:2^N-1] range, then round it
D_I=round((I_sampled-min(I))/(max(I)-min(I))*(2^N-1));
D_Q=round((Q_sampled-min(Q))/(max(Q)-min(Q))*(2^N-1));

%DAC output range
DAC_or=1; %[V]
%DAC step
LSB=DAC_or/(2^N-1);
%define DAC array
DAC_I=LSB*[0 ones(1,2^N-1)];
output_levels_I=cumsum(DAC_I);

DAC_Q=LSB*[0 ones(1,2^N-1)];
output_levels_Q=cumsum(DAC_Q);

VDAC_I=output_levels_I(D_I+1);
VDAC_Q=output_levels_Q(D_Q+1);

%add S/H
tmp=ones(1/(fsample*dt),1)*VDAC_I;
VoutDAC_I=tmp(:)'-DAC_or/2;
clear tmp

tmp=ones(1/(fsample*dt),1)*VDAC_Q;
VoutDAC_Q=tmp(:)'-DAC_or/2;
clear tmp

figure
plot(t(1:1000),VoutDAC_I(1:1000),t(1:1000),(I(1:1000)-min(I))/(max(I)-min(I))*DAC_or-DAC_or/2);
xlabel('Time [s]')
ylabel('Amplitude [-]')
legend('DAC_I input','ideal voltage (scaled for comparison)')

figure
plot(t(1:1000),VoutDAC_Q(1:1000),t(1:1000),(Q(1:1000)-min(Q))/(max(Q)-min(Q))*DAC_or-DAC_or/2);
xlabel('Time [s]')
ylabel('Amplitude [-]')
legend('DAC_Q input','ideal voltage (scaled for comparison)')

G1=1/sinc(fin(1)/fsample); %[-]

I_mix = 2 * VoutDAC_I .* cos(2 * pi * fosc * t);
Q_mix = 2 * VoutDAC_Q .* sin(2 * pi * fosc * t);

Mix_out = I_mix + Q_mix;

G2=1; %[-]

Vout2=G1 * G2 * Mix_out;


%compute fft
V2=abs(fft(Vout2));
V2dB=20*log10(V2); %[dB]
V1=abs(fft(I_mix)+1e-10);
V1dB=20*log10(V1); %[dB]
% FFT resolution
fres=1/T; %[Hz]
% maximum FFT frequency
fmax=1/dt; %[Hz]
% frequency vector
f=0:fres:(fmax-fres); %[Hz]

%% plot FFT
figure
subplot(2,1,1)
plot(f,V1dB)
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
legend('V_{I}')
grid

subplot(2,1,2)
plot(f,V2dB)
hold off
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
legend('V_{out2}')
grid

%compute fft
V2=abs(fft(Vout2));
V2dB=20*log10(V2); %[dB]
V1=abs(fft(Q_mix)+1e-10);
V1dB=20*log10(V1); %[dB]
% FFT resolution
fres=1/T; %[Hz]
% maximum FFT frequency
fmax=1/dt; %[Hz]
% frequency vector
f=0:fres:(fmax-fres); %[Hz]

%% plot FFT
figure
subplot(2,1,1)
plot(f,V1dB)
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
legend('V_{Q}')
grid

subplot(2,1,2)
plot(f,V2dB)
hold off
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
legend('V_{out2}')
grid

%%
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