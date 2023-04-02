clear
close all

% fix time axis
dt=1e-9;
T=1e-3;
t=0:dt:T;

% pole (cutoff) frequency
fp=1e6;

%generate filter parameters
a=[1 2*pi*fp*dt-1];
b=[0 2*pi*fp*dt];

% generate white noise
x=randn(size(t));
% apply filter to the noise
y=filter(b,a,x);

%compute fft
XdB=20*log10(abs(fft(x)));
YdB=20*log10(abs(fft(y)));

figure
f=0:1/T:1/dt;
semilogx(f,XdB,f,YdB)
xlabel('f [Hz]')
ylabel('H [dB]')
legend('White noise', 'Filtered white noise')

% test filter response using sinusoid at different frequencies
% genertae frequency array
fin=logspace(3,9,100);
%make sure that sinusoid frequency falls on an FFT bin
fin=round(fin*T)/T;

for i=1:length(fin)
    % generate test inout
    xtest=sin(2*pi*fin(i)*t);
    %apply filter
    ytest=filter(b,a,xtest);
    %compute ratio in amplitude between input and output
    H(i)=max(abs(fft(ytest)))/max(abs(fft(xtest)));
end;

figure
semilogx(fin,20*log10(H))
xlabel('f [Hz]')
ylabel('H [dB]')
title('Frequency response of the filter')