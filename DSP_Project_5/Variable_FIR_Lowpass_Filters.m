%
% Design of Variable_FIR_Lowpass_Filters
%
clear;   % clear workspace
clc;     % clear command window
%
N = 50;
M = 7; 
wp1 = 0.2*pi;
wp2 = 0.6*pi;
wt = 0.2*pi;
pointw = 200; % 取樣點數
pointp = 60;  % 取樣點數
%
%
NH = N/2;
nm = (NH+1) * (M+1);
deltaw = pi/pointw; 
deltap = 1/pointp;
%
%
r = zeros(nm, 1);
Qp = zeros(nm, nm);
Qs = zeros(nm, nm);
PASS = 0; 
STOP = 0;
for ip = 0:pointp
        p = -0.5 + ip*deltap;
        wp = (p+0.5) * (wp2-wp1) + wp1;
        ws = wp + wt;
    for iw = 0:pointw
        w = iw * deltaw;
        c = zeros(nm, 1);
        for inm = 0:nm-1
            n = mod(inm, NH+1);
            m = floor(inm/(NH+1));
            c(inm+1) = p^m * cos(n*w);
        end
        if w <= wp
            r = r - 2*c;
            Qp = Qp + c*c';
            PASS = PASS + 1;
        elseif w >= ws
            Qs = Qs + c*c';
            STOP = STOP + 1;
        end
    end
end
r = 0.5 * (wp1+wp2) * (r/PASS);
Qp = 0.5 * (wp1+wp2) * (Qp/PASS);
Qs = 0.5 * (pi - wp1 - wt) + (pi - wp2 - wt) * (Qs/STOP);
Q = Qp + Qs;
a = -0.5 * inv(Q)*r;
%
%
a2 = reshape(a, NH+1, M+1);
h = zeros(N+1, M+1);
h(NH+1, :) = a2(1, :);
h(1:NH, :) = 0.5 * a2(NH+1:-1:2,:);
h(NH+2:N+1,:) = 0.5 * a2(2:NH+1, :);
%
%
MR = zeros(pointw+1, pointp+1); % Magnitude Response
for ip = 0:pointp
    p = -0.5 + ip*deltap;
    hnp = h(:, 1);
    for im = 1:M
        hnp = hnp + h(:, im+1) * p^(im);
    end
    MR(:, ip+1) = abs(freqz(hnp, 1, 0:deltaw:pi)); % lowpass
end
%
%
XX = zeros(pointw+1, pointp+1);
YY = zeros(pointw+1, pointp+1);
%
for ip = 0:pointp
    XX(:, ip+1) = (0:deltaw:pi)'/pi;
end
for iw = 0:pointw
    YY(iw+1,:) = (-0.5:deltap:0.5);
end
%
%
plot3(XX, YY, MR);
axis([0, 1, -0.5,0.5, 0, 1.1]);
xlabel('Normalized Frequency');
ylabel('Variable p');
zlabel('Magnitude Response');
pause;
%
% print各個子濾波器
%
for im = 0:M
    MRs = abs(freqz(h(:,im+1), 1, 0:pi/200:pi));
    subplot(3, 3, im+1);
    plot(0:1/200:1, MRs);
    axis([0, 1, 0, 100]);
end
%
