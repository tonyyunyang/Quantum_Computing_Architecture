% Cleanup
clear all;
close all;
clc;

clear;
close all;


% define input signal
% input amplitude 
Question_2_A=[10e-3 0 0]; %[V]
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
Question_2_Vin=Question_2_A(1)*cos(2*pi*fin(1)*t)+...
    Question_2_A(2)*cos(2*pi*fin(2)*t)+...
    Question_2_A(3)*cos(2*pi*fin(3)*t); %[V]

G1=200; %[-]

% define the range of possible coeffients for both 2nd and 3rd order term
Gtemp3_range = 0:20:4000;
Gtemp4_range = 0:.025:5;
% define an array filled with zeros to store Fidelity of q2 and q3
% the size of these two arrays has to be the same as the number of possible
% coeffients for both 2nd and 3rd order term
Question_2_F2_G3 = zeros(length(Gtemp3_range), 0);
Question_2_F2_G4 = zeros(length(Gtemp4_range), 0);
Question_2_F3_G3 = zeros(length(Gtemp3_range), 0);
Question_2_F3_G4 = zeros(length(Gtemp4_range), 0);

Question_3_F2_G3 = zeros(length(Gtemp3_range), 0);
Question_3_F2_G4 = zeros(length(Gtemp4_range), 0);
Question_3_F3_G3 = zeros(length(Gtemp3_range), 0);
Question_3_F3_G4 = zeros(length(Gtemp4_range), 0);

% start the first loop to collect all possible F2 and F3 (This corresponds to only having the second order term)
for i = 1:length(Gtemp3_range)
    Gtemp1 = 805.7;
    Gtemp2 = 161162;
    Vout1=G1*Question_2_Vin + Gtemp1 * Question_2_Vin.^2 + Gtemp2 * Question_2_Vin.^3;
    Vmix=Vout1.*cos(2*pi*t*fosc);
    G2=1; %[-]
    Gtemp3 = Gtemp3_range(i);
    Gtemp4 = 0;
    Vout2=G2*Vmix + Gtemp3 * Vmix.^2 + Gtemp4 * Vmix.^3;

    %compute fft
    X=abs(fft(Vout2));
    XdB=20*log10(X); %[dB]
    Y=abs(fft(Question_2_Vin));
    YdB=20*log10(Y); %[dB]
    % FFT resolution
    fres=1/T; %[Hz]
    % maximum FFT frequency
    fmax=1/dt; %[Hz]
    % frequency vector
    f=0:fres:(fmax-fres); %[Hz]

    % Simulate and plot
    [U, probabilities] = spine_no_plot(fR, f0, dt, Vout2, 10);

    % Compute fidelity of the operation on qubit 1
    % Ideal transformation representing a pi-rotation around X
     Uideal = [0 1;
               1 0];
    % Compute fidelity by comparing the actual transformation to the ideal one
    %F1 = fidelity(U(:, :, 1), Uideal);
    
    % For the other 2 qubits, we would like to compare the actual transormation
    % to the identity operation (I). However, any leakage form the driving
    % signal is going to induce Z-rotations on the qubits. Such rotations can
    % be easily cancelled by software, i.e. by not applying any extra
    % electrical pulse. This correction is currently not implemented in this
    % code to make it easier for the students. However, for the simple case in
    % which the qubits ar einitialized in |0>, the fidelity on idel qubits can
    % be simply computed by measuring the qubit along Z. This measurmenet is
    % insensitive with respect to rotazion around Z.
    Question_2_F2_G3(i)=probabilities(3,end,2);
    Question_2_F3_G3(i)=probabilities(3,end,3);
end

% start the second loop to collect all possible F2 and F3 (This corresponds to only having the third order term)
for j = 1:length(Gtemp4_range)
    Gtemp1 = 805.7;
    Gtemp2 = 161162;
    Vout1=G1*Question_2_Vin + Gtemp1 * Question_2_Vin.^2 + Gtemp2 * Question_2_Vin.^3;
    Vmix=Vout1.*cos(2*pi*t*fosc);
    G2=1; %[-]
    Gtemp3 = 0;
    Gtemp4 = Gtemp4_range(j);
    Vout2=G2*Vmix + Gtemp3 * Vmix.^2 + Gtemp4 * Vmix.^3;

    %compute fft
    X=abs(fft(Vout2));
    XdB=20*log10(X); %[dB]
    Y=abs(fft(Question_2_Vin));
    YdB=20*log10(Y); %[dB]
    % FFT resolution
    fres=1/T; %[Hz]
    % maximum FFT frequency
    fmax=1/dt; %[Hz]
    % frequency vector
    f=0:fres:(fmax-fres); %[Hz]

    % Simulate and plot
    [U, probabilities] = spine_no_plot(fR, f0, dt, Vout2, 10);

    % Compute fidelity of the operation on qubit 1
    % Ideal transformation representing a pi-rotation around X
     Uideal = [0 1;
               1 0];
    % Compute fidelity by comparing the actual transformation to the ideal one
    %F1 = fidelity(U(:, :, 1), Uideal);
    
    % For the other 2 qubits, we would like to compare the actual transormation
    % to the identity operation (I). However, any leakage form the driving
    % signal is going to induce Z-rotations on the qubits. Such rotations can
    % be easily cancelled by software, i.e. by not applying any extra
    % electrical pulse. This correction is currently not implemented in this
    % code to make it easier for the students. However, for the simple case in
    % which the qubits ar einitialized in |0>, the fidelity on idel qubits can
    % be simply computed by measuring the qubit along Z. This measurmenet is
    % insensitive with respect to rotazion around Z.
    Question_2_F2_G4(j)=probabilities(3,end,2);
    Question_2_F3_G4(j)=probabilities(3,end,3);
end


% input amplitude 
Question_3_A=[10e-3 10e-3 0]; %[V]
% input signal
Question_3_Vin=Question_3_A(1)*cos(2*pi*fin(1)*t)+...
    Question_3_A(2)*cos(2*pi*fin(2)*t)+...
    Question_3_A(3)*cos(2*pi*fin(3)*t); %[V]

% start the first loop to collect all possible F2 and F3 (This corresponds to only having the second order term)
for i = 1:length(Gtemp3_range)
    Gtemp1 = 94.9;
    Gtemp2 = 18891.25;
    Vout1=G1*Question_3_Vin + Gtemp1 * Question_3_Vin.^2 + Gtemp2 * Question_3_Vin.^3;
    Vmix=Vout1.*cos(2*pi*t*fosc);
    G2=1; %[-]
    Gtemp3 = Gtemp3_range(i);
    Gtemp4 = 0;
    Vout2=G2*Vmix + Gtemp3 * Vmix.^2 + Gtemp4 * Vmix.^3;

    %compute fft
    X=abs(fft(Vout2));
    XdB=20*log10(X); %[dB]
    Y=abs(fft(Question_3_Vin));
    YdB=20*log10(Y); %[dB]
    % FFT resolution
    fres=1/T; %[Hz]
    % maximum FFT frequency
    fmax=1/dt; %[Hz]
    % frequency vector
    f=0:fres:(fmax-fres); %[Hz]

    % Simulate and plot
    [U, probabilities] = spine_no_plot(fR, f0, dt, Vout2, 10);

    % Compute fidelity of the operation on qubit 1
    % Ideal transformation representing a pi-rotation around X
     Uideal = [0 1;
               1 0];
    % Compute fidelity by comparing the actual transformation to the ideal one
    %F1 = fidelity(U(:, :, 1), Uideal);
    
    % For the other 2 qubits, we would like to compare the actual transormation
    % to the identity operation (I). However, any leakage form the driving
    % signal is going to induce Z-rotations on the qubits. Such rotations can
    % be easily cancelled by software, i.e. by not applying any extra
    % electrical pulse. This correction is currently not implemented in this
    % code to make it easier for the students. However, for the simple case in
    % which the qubits ar einitialized in |0>, the fidelity on idel qubits can
    % be simply computed by measuring the qubit along Z. This measurmenet is
    % insensitive with respect to rotazion around Z.
    Question_3_F2_G3(i)=probabilities(3,end,2);
    Question_3_F3_G3(i)=probabilities(3,end,3);
end

% start the second loop to collect all possible F2 and F3 (This corresponds to only having the third order term)
for j = 1:length(Gtemp4_range)
    Gtemp1 = 94.9;
    Gtemp2 = 18891.25;
    Vout1=G1*Question_3_Vin + Gtemp1 * Question_3_Vin.^2 + Gtemp2 * Question_3_Vin.^3;
    Vmix=Vout1.*cos(2*pi*t*fosc);
    G2=1; %[-]
    Gtemp3 = 0;
    Gtemp4 = Gtemp4_range(j);
    Vout2=G2*Vmix + Gtemp3 * Vmix.^2 + Gtemp4 * Vmix.^3;

    %compute fft
    X=abs(fft(Vout2));
    XdB=20*log10(X); %[dB]
    Y=abs(fft(Question_3_Vin));
    YdB=20*log10(Y); %[dB]
    % FFT resolution
    fres=1/T; %[Hz]
    % maximum FFT frequency
    fmax=1/dt; %[Hz]
    % frequency vector
    f=0:fres:(fmax-fres); %[Hz]

    % Simulate and plot
    [U, probabilities] = spine_no_plot(fR, f0, dt, Vout2, 10);

    % Compute fidelity of the operation on qubit 1
    % Ideal transformation representing a pi-rotation around X
     Uideal = [0 1;
               1 0];
    % Compute fidelity by comparing the actual transformation to the ideal one
    %F1 = fidelity(U(:, :, 1), Uideal);
    
    % For the other 2 qubits, we would like to compare the actual transormation
    % to the identity operation (I). However, any leakage form the driving
    % signal is going to induce Z-rotations on the qubits. Such rotations can
    % be easily cancelled by software, i.e. by not applying any extra
    % electrical pulse. This correction is currently not implemented in this
    % code to make it easier for the students. However, for the simple case in
    % which the qubits ar einitialized in |0>, the fidelity on idel qubits can
    % be simply computed by measuring the qubit along Z. This measurmenet is
    % insensitive with respect to rotazion around Z.
    Question_3_F2_G4(j)=probabilities(3,end,2);
    Question_3_F3_G4(j)=probabilities(3,end,3);
end

% we plot four seperate graphs below
figure;
plot(Gtemp3_range, Question_2_F2_G3);
title('Qubit 2 Fidelity VS G2.2 (Comparison)');
xlabel('G2.2');
ylabel('Q2 Fidelity');
hold on
plot(Gtemp3_range, Question_3_F2_G3);
hold off

figure;
plot(Gtemp3_range, Question_2_F3_G3);
title('Qubit 3 Fidelity VS G2.2');
xlabel('G2.2');
ylabel('Q3 Fidelity');
hold on
plot(Gtemp3_range, Question_3_F3_G3);
hold off

figure;
plot(Gtemp4_range, Question_2_F2_G4);
title('Qubit 2 Fidelity VS G2.3');
xlabel('G2.3');
ylabel('Q2 Fidelity');
hold on
plot(Gtemp4_range, Question_3_F2_G4);
hold off

figure;
plot(Gtemp4_range, Question_2_F3_G4);
title('Qubit 3 Fidelity VS G2.3');
xlabel('G2.3');
ylabel('Q3 Fidelity');
hold on
plot(Gtemp4_range, Question_3_F2_G4);
hold off