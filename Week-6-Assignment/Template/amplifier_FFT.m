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
x=A*sin(2*pi*fin*t); %[V]

% define amplifier
% gain
a1=10; %[-]
a3=-40;%[V^-2]

% compute output signal
y=a1*x+a3*x.^3; %[V]

%plot
plot(t,x,t,y)
xlabel('Time [s]')
ylabel('Amplitude [V]')
legend('x','y')


%compute fft
X=abs(fft(x));
XdB=20*log10(X); %[dB]
Y=abs(fft(y));
YdB=20*log10(Y); %[dB]
% FFT resolution
fres=1/T; %[Hz]
% maximum FFT frequency
fmax=1/dt; %[Hz]
% frequency vector
f=0:fres:(fmax-fres); %[Hz]

%plot FFT
figure
hold on
plot(f,XdB,f,YdB)
hold off
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
legend('x','y')
grid
