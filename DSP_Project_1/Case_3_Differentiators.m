%
% Design of case 3 Differentiator
%
clear all;   % clear workspace
clc;         % clear command window
N = 41;
wp = 0.9*pi; % wp = wp1
% no need wp2
% 
%
NH = (N-1)/2;
%
%
P = zeros(NH,1);
Q = zeros(NH,NH);
for i=1:NH % 0:1:NH
    P(i) = 2*wp*cos(i*wp)/i - 2*sin(i*wp)/i^2;
    for j=1:NH
        if i==j
            Q(i,j) = 0.5*wp - 0.5*sin((i+j)*wp)/(i+j);
        else
            Q(i,j) = 0.5*sin((i-j)*wp)/(i-j) - 0.5*sin((i+j)*wp)/(i+j);
       end
    end
end
A = -0.5*inv(Q)*P;
%
%
h = zeros(N,1);
h(1:NH) = 0.5*A(NH:-1:1);
h(NH+2:N) = -0.5*A;
%
% subplot - need to show several plots simulteneously
% -high-pass filter
subplot(1,2,1);
stem(0:N-1,h);
xlabel('n');
ylabel('Impulse Response');
%
% Amplitude Response
subplot(1,2,2);
AR = abs(freqz(h,1,0:pi/200:pi));
plot(0:1/200:1,AR);
axis([0,1,0,3]);
xlabel('Normalized Frequency (\omega/\pi)');
ylabel('Amplitude Response');
%
% case 3 don't need dB plots
%

