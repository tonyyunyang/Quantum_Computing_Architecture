% Cleanup
clear all;
close all;
clc;

clear;
close all;

% Set the Larmor frequencies of the qubits (3 qubits)
f0 = [11e9, 12e9, 13e9];

% Set the Rabi frequencies
fR = [1e6, 1e6, 1e6];

%define microwave frequency
fosc=10e9; %[Hz]


%number of bits of the DAC
N=10;

%sample frequency of the DAC
fsample=10e9;

% define input signal
% input amplitude 
A=[10e-3 0 0]; %[V]
% input frequency 
fin=[1 2 3]*1e9; %[Hz]
% inout phases
ph_in=[0 0 0]; %[rad]

% total simulation time. It is computed as the time needed for a pi
% rotation for Q1
T=1/fR(1)/2; %[s]
% simulation time step
dt=1e-12; %[s]
% define time vector
t=dt:dt:T; %[s]

%define the mismatch
delta = 0:0.001:10;
F1 = zeros(length(delta), 0);
F2 = zeros(length(delta), 0);
F3 = zeros(length(delta), 0);

for i = 1:length(delta)
    % input signal
    in=A(1)*cos(2*pi*fin(1)*t+ph_in(1))+...
        A(2)*cos(2*pi*fin(2)*t+ph_in(2))+...
        A(3)*cos(2*pi*fin(2)*t+ph_in(3)); %[V]
    
    %sampling the input signal
    in_sampled=in(1:1/(fsample*dt):end);
    
    %quantize the signal on N-bit
    %first, scale the signal to fit in [0:2^N-1] range, then round it
    Din=round((in_sampled-min(in))/(max(in)-min(in))*(2^N-1));
    
    %DAC output range
    DAC_or=1; %[V]
    %DAC step
    LSB=DAC_or/(2^N-1);
    %define DAC array
    %introduce a 10% mismatch
    DAC = (1 + delta(i)) * LSB * [0 ones(1, 2^N-1)];
    output_levels=cumsum(DAC);
    
    VDAC=output_levels(Din+1);
    
    %add S/H
    tmp=ones(1/(fsample*dt),1)*VDAC;
    VoutDAC=tmp(:)'-DAC_or/2;
    clear tmp
    
    % figure
    % plot(t(1:1000),VoutDAC(1:1000),t(1:1000),(in(1:1000)-min(in))/(max(in)-min(in))*DAC_or-DAC_or/2);
    % xlabel('Time [s]')
    % ylabel('Amplitude [-]')
    % legend('DAC input','ideal voltage (scaled for comparison)')
    
    G1=4/sinc(fin(1)/fsample); %[-]
    
    Vout1=G1*VoutDAC;
    
    Vmix=Vout1.*cos(2*pi*t*fosc);
    
    G2=1; %[-]
    
    Vout2=G2*Vmix;
    
    
    %compute fft
    V2=abs(fft(Vout2));
    V2dB=20*log10(V2); %[dB]
    V1=abs(fft(Vout1)+1e-10);
    V1dB=20*log10(V1); %[dB]
    % FFT resolution
    fres=1/T; %[Hz]
    % maximum FFT frequency
    fmax=1/dt; %[Hz]
    % frequency vector
    f=0:fres:(fmax-fres); %[Hz]
    
    %% plot FFT
    % figure
    % subplot(2,1,1)
    % plot(f,V1dB)
    % xlabel('Frequency [Hz]')
    % ylabel('Amplitude [dB]')
    % legend('V_{out1}')
    % xlim([0 20]*1e9)
    % grid
    % 
    % subplot(2,1,2)
    % plot(f,V2dB)
    % hold off
    % xlabel('Frequency [Hz]')
    % ylabel('Amplitude [dB]')
    % legend('V_{out2}')
    % xlim([0 20]*1e9)
    % grid
    
    %%
    % Simulate and plot
    [U, probabilities] = spine_no_plot(fR, f0, dt, Vout2, 10);
    
    
    % Compute fidelity of the operation on qubit 1
    % Ideal transformation representing a pi-rotation around X
     Uideal = [0 1;
               1 0];
    % Compute fidelity by comparing the actual transformation to the ideal one
    F1(i) = fidelity(U(:, :, 1), Uideal);
    
    % For the other 2 qubits, we would like to compare the actual transormation
    % to the identity operation (I). However, any leakage form the driving
    % signal is going to induce Z-rotations on the qubits. Such rotations can
    % be easily cancelled by software, i.e. by not applying any extra
    % electrical pulse. This correction is currently not implemented in this
    % code to make it easier for the students. However, for the simple case in
    % which the qubits ar einitialized in |0>, the fidelity on idel qubits can
    % be simply computed by measuring the qubit along Z. This measurmenet is
    % insensitive with respect to rotazion around Z.
    F2(i)=probabilities(3,end,2);
    F3(i)=probabilities(3,end,3);
    
    % disp(['Fidelity of qubit 1: ' num2str(F1*100) '%']);
    % disp(['Fidelity of qubit 2: ' num2str(F2*100) '%']);
    % disp(['Fidelity of qubit 3: ' num2str(F3*100) '%']);
end

figure
subplot(3,1,1)
plot(delta,F1)
xlabel('Mismatch delta')
ylabel('Q1 Fidelity')
legend('Q1 VS Mismatch delta')
xlim([-5 5])
grid

subplot(3,1,2)
plot(delta,F2)
xlabel('Mismatch delta')
ylabel('Q2 Fidelity')
legend('Q2 VS Mismatch delta')
xlim([-5 5])
grid

subplot(3,1,3)
plot(delta,F3)
hold off
xlabel('Mismatch delta')
ylabel('Q3 Fidelity')
legend('Q3 VS Mismatch delta')
xlim([-5 5])
grid