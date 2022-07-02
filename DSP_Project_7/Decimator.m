%
% Experiment_For_Sampling_Rate_Reduction
%
% Decimator
%
clear;   % clear workspace
clc;     % clear command window
%
N = 101;
wp = 0.25*pi;
ws = 0.30*pi;
% Wp = Ws = 1;
%
NH = (N-1)/2;
%
%
P = zeros(NH+1,1);
Qp = zeros(NH+1,NH+1);
Qs = zeros(NH+1,NH+1);
for i=0:NH % 0:1:NH
    if i==0
        P(i+1) = -2*wp;
    else
        P(i+1) = -2*sin(i*wp)/i;
    end
    for j=0:NH
       if i==0 && j==0
           Qp(i+1,j+1) = wp;
           Qs(i+1,j+1) = pi - ws;
       elseif i==j
           Qp(i+1,j+1) = 0.5*wp + 0.5*sin((i+j)*wp)/(i+j);
           Qs(i+1,j+1) = 0.5*(pi - ws) - 0.5*sin((i+j)*ws)/(i+j);
       else
           Qp(i+1,j+1) = 0.5*sin((i-j)*wp)/(i-j) + 0.5*sin((i+j)*wp)/(i+j);
           Qs(i+1,j+1) = -0.5*sin((i-j)*ws)/(i-j) - 0.5*sin((i+j)*ws)/(i+j);
       end
    end
end
Q = Qp + Qs;
A = -0.5*inv(Q)*P;
%
%
h = zeros(N,1);
h(NH+1) = A(1);
h(1:NH) = 0.5*A(NH+1:-1:2);
h(NH+2:N) = 0.5*A(2:NH+1);
%
% need to show several plots simulteneously
%
subplot(2,3,1);
stem(0:N-1,h);
xlabel('n');
ylabel('Impulse Response');
title('Case 1 Lowpass Filter');
%
% Amplitude Response
%
subplot(2,3,4);
AR = abs(freqz(h,1,0:pi/200:pi));
plot(0:1/200:1,AR);
axis([0,1,0,1.1]);
xlabel('Normalized Frequency (\omega/\pi)');
ylabel('Amplitude Response');
title('Case 1 Lowpass Filter');
%
% Amplitude Response in dB
% subplot(2,2,4);
% plot(0:1/200:1,20*log10(AR));
% axis([0,1,-60,10]);
% xlabel('Normalized Frequency (\omega/\pi)');
% ylabel('Amplitude Response in dB');
%
% M = 2
%
h2 = zeros(51, 1);
for i = 0:50
    h2(i + 1) = h(2*i + 1);
end
subplot(2,3,2);
stem(0:(N-1)/2,h2);
xlabel('n');
title('Down Sampling M = 2');
%
subplot(2,3,5);
FR2 = abs(freqz(h2,1,0:pi/200:pi));
plot(0:1/200:1,FR2);
axis([0,1,0,1.1]);
xlabel('Normalized Frequency (\omega/\pi)');
ylabel('Amplitude Response');
title('Down Sampling M = 2');
%
% M = 3
%
h3 = zeros(floor(N/3) + 1, 1);
for i = 0:floor(N/3)
    h3(i + 1) = h(3*i + 1);
end
subplot(2, 3, 3);
stem(0:floor(N/3),h3);
xlabel('n');
title('Down Sampling M = 3');
%
subplot(2, 3, 6);
FR3 = abs(freqz(h3,1,0:pi/200:pi));
plot(0:1/200:1,FR3);
axis([0,1,0,1.1]);
xlabel('Normalized Frequency (\omega/\pi)');
ylabel('Amplitude Response');
title('Down Sampling M = 3');

