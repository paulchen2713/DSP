%
% Design of Two_Channel_Filter_Banks
%
clear all;   % clear workspace
clc;         % clear command window
N = 49;
ws = 0.55*pi;
pointw = 400;
alpha = 1; 
beta = 0.5;
epsilon = 0.00001;
%
%
wp = pi - ws;
NH = (N+1)/2;
deltaw = pi / pointw;
%
%
P0 = zeros(NH, 1);
Qp = zeros(NH, NH);
Qs = zeros(NH, NH);
pointp = 0;
points = 0;
for iw = 0:pointw
    w = iw * deltaw;
    if w <= wp
        pointp = pointp + 1;
        C = zeros(NH, 1);
        for in = 1:NH
            C(in) = cos(w*(in - 0.5));
        end
        P0 = P0 - 2*C;
        Qp = Qp + C*C';
    elseif w >= ws
        points = points + 1;
        C = zeros(NH, 1);
        for in = 1:NH
            C(in) = cos(w*(in - 0.5));
        end
        Qs = Qs + C*C';
    end
end
P0 = wp * P0 / pointp;
Qp = wp * Qp / pointp;
Qs = (pi - ws) * Qs / points;
A = -0.5 * inv(Qp + Qs) * P0;
h = zeros(N+1, 1);
h(1:NH) = 0.5 * A(NH:-1:1);
h(NH+1:N+1) = 0.5 * A;
FR = abs(freqz(h, 1, 0:deltaw:pi));
plot(0:deltaw/pi:pi/pi, FR);
xlabel('Normalized Frequency (\omega/\pi');
ylabel('Magnitude Response');
title('Initial Lowpass Filter');
%
%
deltak = 10000;
iteration = 0;
%
while deltak >= epsilon
    iteration = iteration + 1
    Ak1 = A; % Ak1取前1次的結果
    P = zeros(NH, 1);
    Q1 = zeros(NH, NH);
    Q2 = zeros(NH, NH);
    for iw = 0:pointw
       w = iw * deltaw; 
       C = zeros(NH, 1);
       for in = 1:NH
           C(in) = cos(w*(in - 0.5));
       end
       Cpi = zeros(NH, 1);
       for in = 1:NH
           Cpi(in) = cos((w - pi)*(in - 0.5));
       end
       P = P - 2*((Ak1'*C)*C + (Ak1'*Cpi)*Cpi);
       Q1 = Q1 + ((Ak1'*C)*C + (Ak1'*Cpi)*Cpi) * ((Ak1'*C)*C + (Ak1'*Cpi)*Cpi)';
       if w >= ws
          Q2 = Q2 + C*C'; 
       end
    end
    P = pi * P / pointw;
    Q1 = pi * Q1 / pointw;
    Q2 = alpha * (pi - ws) * Q2 / points;
    A = -0.5 * inv(Q1 + Q2)*P;
    %
    A = beta * A + (1 - beta) *Ak1;
    %
    h = zeros(N+1, 1);
    h(1:NH) = 0.5 * A(NH:-1:1);
    h(NH+1:N+1) = 0.5 * A;
    FR = abs(freqz(h, 1, 0:deltaw:pi));
    subplot(1, 2, 1);
    plot(0:deltaw/pi:pi/pi, FR);
    axis([0, 1, 0, 1.1]);
    xlabel('Normalized Frequency (\omega/\pi');
    ylabel('Magnitude Response');
    title('Lowpass Filter');
    % frequency response of the whole system
    h1 = zeros(N+1, 1);
    for in = 0:N
       h1(in + 1) = (-1)^in * h(in + 1); 
    end
    t = conv(h, h) - conv(h1, h1);
    FRt = abs(freqz(t, 1, 0:deltaw:pi));
    subplot(1, 2, 2);
    plot(0:deltaw/pi:pi/pi, FRt);
    axis([0, 1, 0.8, 1.2]);
    xlabel('Normalized Frequency (\omega/\pi');
    ylabel('Magnitude Response');
    title('Overall System');
    % stop condition
    deltak = norm(A - Ak1) / norm(A)
    pause;
end
