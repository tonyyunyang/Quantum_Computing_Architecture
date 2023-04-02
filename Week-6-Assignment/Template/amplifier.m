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

% compute output signal
y=a1*x; %[V]

%plot
plot(t,x,t,y)
xlabel('Time [s]')
ylabel('Amplitude [V]')
legend('x','y')