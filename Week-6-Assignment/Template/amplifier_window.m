clear;
close all;

% define input signal
% input amplitude 
A=10e-3; %[V]
% input frequency 
fin=1e6; %[Hz]

% number of periods to be simulated
no_periods=20; %[-]
% total simulation time
T=no_periods/fin; %[s]
% points per periond
no_points_per_period=100; %[-]
% simulation time step
dt=1/fin/no_points_per_period; %[s]
% define time vector
t=0:dt:T-dt; %[s]

% input signal
x=A*sin(2*pi*1.02e6*t); %[V]

% define input-referred noise
% equivalent noise temperature
Te=300; %[K]
% reference resistor
R0=50; %[ohm]
%Boltzmann constant
k=1.38e-23; %[m^2kgs^-2K^-1]
% noise spectral density
Sn=4*k*Te*R0; %[V^2/Hz]
noise_rms=sqrt(Sn/dt/2); %[V]
noise=noise_rms*randn(size(t));

% define amplifier
% gain
a1=10; %[-]
a3=-4;%[V^-2]

% compute output signal
x_noise=x+noise;
y=a1*x_noise+a3*x_noise.^3; %[V]

%plot
plot(t,x,t,y)
xlabel('Time [s]')
ylabel('Amplitude [V]')
legend('x','y')

%define window
w=kaiser(length(t),20)';
%w=ones(size(t));

%compute fft
X=abs(fft(x.*w));
XdB=20*log10(X); %[dB]
Y=abs(fft(y.*w));
YdB=20*log10(Y); %[dB]
% FFT resolution
fres=1/T; %[Hz]
% maximum FFT frequency
fmax=1/dt; %[Hz]
% frequency vector
f=0:fres:(fmax-fres); %[Hz]

%compute SNR
% signal bins
maintone = round(fin/fres)+1;
signal_bins = [maintone-7:maintone+7];
%bandwidth for SNR calculation
BW=5e6; %[Hz]
%noise bins
inband_bins=1:round(BW/fres);
noise_bins=setdiff(inband_bins,signal_bins);

%signal power
S=sum(Y(signal_bins).^2);
%noise power
N=sum(Y(noise_bins).^2);
%SNR
SNR=10*log10(S/N);

figure
plot(f,XdB,f,YdB,f(signal_bins),YdB(signal_bins),f(noise_bins),YdB(noise_bins),'Linewidth',2)
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
legend('x','y','y signal','ynoise')
grid
