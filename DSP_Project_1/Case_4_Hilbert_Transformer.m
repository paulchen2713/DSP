%
% Design of case 4 Hilbert Transformers
%
clear all;   % clear workspace
clc;         % clear command window
N = 30;
wp1 = 0.05*pi;
% wp2 = pi; 因為上限可到pi，故不再需要wp2
% 
%
NH = N/2;
%
%
P = zeros(NH,1);
Q = zeros(NH,NH);
for i=1:NH % 0:1:NH
    P(i) = 2*cos((i-0.5)*wp1/(i-0.5));
    for j=1:NH
        if i==j
            Q(i,j) = 0.5*(pi-wp1) + 0.5*sin((i+j-1)*wp1)/(i+j-1);
        else
            Q(i,j) = -0.5*sin((i-j)*wp1)/(i-j) + 0.5*sin((i+j-1)*wp1)/(i+j-1);
       end
    end
end
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
axis([0,1,0,1.1]);
xlabel('Normalized Frequency (\omega/\pi)');
ylabel('Amplitude Response');
pause;
%
% case 3 don't need dB plots
