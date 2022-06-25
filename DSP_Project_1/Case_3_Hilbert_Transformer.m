%
% Design of case 3 Hilbert Transformers
%
clear;       % clear workspace
clc;         % clear command window
N = 31;
wp1 = 0.05*pi;
wp2 = 0.95*pi;
% 
%
NH = (N-1)/2;
%
%
P = zeros(NH,1);
Q = zeros(NH,NH);
for i=1:NH % 0:1:NH
    P(i) = -2*(cos(i*wp2) - cos(i*wp1))/i;
    for j=1:NH
        if i==j
            Q(i,j) = 0.5*(wp2-wp1) - 0.5*sin((i+j)*wp2)/(i+j) + 0.5*sin((i+j)*wp1)/(i+j);
        else
            Q(i,j) = 0.5*sin((i-j)*wp2)/(i-j) - 0.5*sin((i-j)*wp1)/(i-j) ....
                     - 0.5*sin((i+j)*wp2)/(i+j) + 0.5*sin((i+j)*wp1)/(i+j);
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
%
% case 3 don't need dB plots
%

