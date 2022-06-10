%
% Design of Case 4 Differentiator
%
clear all;   % clear workspace
clc;         % clear command window
N = 30;      % filter length
% wp1 = 0;  因為下限從0開始，故不再需要wp1
% wp2 = pi; 因為上限可到pi，故不再需要wp2
% 
%
NH = N/2;
%
%
P = zeros(NH,1); % declaration 宣告
for i=1:NH % 0:1:NH
    P(i) = -2*sin((i-0.5)*pi)/(i-0.5)^2;
end
Q = 0.5*pi*eye(NH);
A = -0.5*inv(Q)*P;
%
%
h = zeros(N,1); % 怎麼求得h ?
h(1:NH) = 0.5*A(NH:-1:1);
h(NH+1:N) = -0.5*A;
%
% subplot - need to show several plots simulteneously
subplot(1,2,1);
stem(0:N-1,h);
xlabel('n');
ylabel('Impulse Response');
%
% Amplitude Response
subplot(1,2,2);
AR = abs(freqz(h,1,0:pi/200:pi));
plot(0:1/200:1,AR);
axis([0,1,0,3.5]);
xlabel('Normalized Frequency (\omega/\pi)');
ylabel('Amplitude Response');
pause;
close all;
%
%
% signal simulation
%
NT = (0:200)';
x = zeros(1,201);
for in=0:200
   x(in+1) = cos(in*(in/400)*pi); 
end
subplot(2,1,1);
plot(NT,x);
axis([0,200,-3,3]);
%
%
y = filter(h,1,x);
subplot(2,1,2);
plot(NT,y); % plot shows input delay
axis([0,200,-3,3]);
%
