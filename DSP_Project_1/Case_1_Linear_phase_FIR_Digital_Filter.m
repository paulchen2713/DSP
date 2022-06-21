%
% Design of case 1 (odd, symmetric) linear-phase FIR digital filter
%
clear;       % clear workspace
clc;         % clear command window
%
N  = 31;
wp = 0.4*pi;
ws = 0.5*pi;
% Wp = Ws = 1;
%
NH = (N-1)/2;
%
P  = zeros(NH+1, 1);
Qp = zeros(NH+1, NH+1);
Qs = zeros(NH+1, NH+1);
%
for i = 0:NH % 0:1:NH
    if i == 0
        P(i+1) = -2*wp;
    else
        P(i+1) = -2*sin(i*wp)/i;
    end
    %
    for j = 0:NH
       if i == 0 && j == 0
           Qp(i+1, j+1) = wp;
           Qs(i+1, j+1) = pi - ws;
       elseif i==j
           Qp(i+1, j+1) = 0.5*wp + 0.5*sin((i+j)*wp)/(i+j);
           Qs(i+1, j+1) = 0.5*(pi - ws) - 0.5*sin((i+j)*ws)/(i+j);
       else
           Qp(i+1, j+1) = 0.5*sin((i-j)*wp)/(i-j) + 0.5*sin((i+j)*wp)/(i+j);
           Qs(i+1, j+1) = -0.5*sin((i-j)*ws)/(i-j) - 0.5*sin((i+j)*ws)/(i+j);
       end
    end
end
Q = Qp + Qs;
A = -0.5*inv(Q)*P; % inv(Q)
% A = -0.5 * A;
% A = A/Q;
% A = A * P;
%
h = zeros(N,1);
h(NH+1) = A(1);
h(1:NH) = 0.5*A(NH+1:-1:2);
h(NH+2:N) = 0.5*A(2:NH+1);
%
% need to show several plots simulteneously
%
subplot(2,2,1);
stem(0:N-1,h);
xlabel('n');
ylabel('Impulse Response');
%
% Amplitude Response
%
subplot(2,2,2);
AR = abs(freqz(h,1,0:pi/200:pi));
plot(0:1/200:1,AR);
axis([0,1,0,1.1]);
xlabel('Normalized Frequency (\omega/\pi)');
ylabel('Amplitude Response');
%
% Amplitude Response in dB
%
subplot(2,2,4);
plot(0:1/200:1,20*log10(AR));
axis([0,1,-60,10]);
xlabel('Normalized Frequency (\omega/\pi)');
ylabel('Amplitude Response in dB');
pause;
close all; % close formal plots
%
% signal simulation
%
% low frequency input 1
NT = (0:200)';
xL = sin(0.2*pi*NT);
subplot(2,3,1);
plot(NT,xL);
axis([0,200,-3,3]);
%
% low frequency output 2
%
yL = filter(h,1,xL);
subplot(2,3,4);
plot(NT,yL); % plot shows input delay
axis([0,200,-3,3]);
%
% high frequency input 3
%
xH = 1.5*cos(0.7*pi*NT);
subplot(2,3,2);
plot(NT,xH);
axis([0,200,-3,3]);
%
% high frequency output 4
%
yH = filter(h,1,xH);
subplot(2,3,5);
plot(NT,yH); % plot shows transition time
axis([0,200,-3,3]);
%
% high and low combined input 5
%
subplot(2,3,3);
plot(NT,xL + xH);
axis([0,200,-3,3]);
%
% high and low combined output 6
%
y = filter(h,1,xL + xH);
subplot(2,3,6);
plot(NT,y); % plot 6 is pretty the same with 2
axis([0,200,-3,3]);
%
