%
% Design of 2D FIR filters by McClellan Transformation Techniques
% 2. Circular Type
%
clear;   % clear workspace
clc;     % clear command window
%
N = 20; % -20 ~ 20 共41點
wr = 0.5*pi;
wt = 0.2*pi;
pointw = 200; % 取樣200點
%
%
deltaw = 0.5 * pi / pointw; % 1D
deltaw1 = pi / pointw; % 2D
%
%
r = zeros(2, 1);
Q = zeros(2, 2);
for iw = 0:pointw
   w1 = wr * cos(iw * deltaw);
   w2 = wr * sin(iw * deltaw);
   c = [1, 1-cos(w1)*cos(w2)]';
   r = r -2*(0.5*cos(w1) + 0.5*cos(w2))*c;
   Q = Q + c*c';
end
%
% r = L * r / (pointw + 1);
% Q = L * Q / (pointw + 1);
% L 為線積分長度 但在 inv(Q)*r 會消掉 不影響 T 結果
%
r = r / (pointw + 1);
Q = Q / (pointw + 1);
T = -0.5 * inv(Q) * r;
%
%
wc = acos(T(1));
t11 = T(2);
t00 = -t11;
t10 = 0.5;
t01 = 0.5;
%
% Design of 1D prototype lowpass filters
%
wp = wc;
ws = wp + wt;
P = zeros(N+1, 1);
Qs = zeros(N+1, N+1);
Qp = zeros(N+1, N+1);
pointp = 0;
points = 0;
for iw = 0:pointw
   w = iw * deltaw1;
   if w <= wp
       pointp = pointp + 1;
       c = zeros(N+1, 1);
       for in = 0:N
           c(in+1) = cos(in*w);
       end
       P = P - 2*c;
       Qp = Qp + c*c';
   elseif w >= ws
       points = points + 1;
       c = zeros(N+1, 1);
       for in = 0:N
           c(in+1) = cos(in*w);
       end
       Qs = Qs + c*c';
   end
end
P = wp * P / pointp;
Qp = wp * Qp / pointp;
Qs = (pi - ws) * Qs / points;
A = -0.5 * inv(Qp+Qs) * P;
%
h = zeros(2*N+1, 1);
h(N+1) = A(1);
h(1:N) = 0.5 * flipud(A(2:N+1)); % flip up doum 反轉整個column vector
h(N+2:2*N+1) = 0.5 * A(2:N+1);
FR = abs(freqz(h, 1, 0:deltaw1:pi));
subplot(1, 2, 1);
plot(0:deltaw1/pi:1,FR);
xlabel('Normalized Frequency');
ylabel('Magnitude Response');
title('1D Prototype Filter');
%
% 繪製 a, b係數轉換表格
%
Tab = zeros(N+1, N+1); % Tab a係數轉b係數
Tab(1, 1) = 1;
Tab(2, 2) = 1;
for in = 2:N
    Tab(in+1, 2:N+1) = 2 * Tab(in, 1:N);
    Tab(in+1,:) = Tab(in+1,:) - Tab(in-1,:);
end
B = zeros(N+1, 1);
for in = 0:N
    B = B + A(in+1) * Tab(in+1,:)';
end
F = zeros(3, 3); % 3 x 3 vector(矩陣)
F(1, 1) = 0.25 * t11;
F(1, 2) = 0.5 * t10;
F(1, 3) = 0.25 * t11;
F(2, 1) = 0.5 * t01;
F(2, 2) = t00;
F(2, 3) = 0.5 * t01;
F(3, 1) = 0.25 * t11;
F(3, 2) = 0.5 * t10;
F(3, 3) = 0.25 * t11;
%
%
h2 = zeros(2*N+1, 2*N+1);
h2(N+1, N+1) = B(1); % 中心點
h2(N:N+2, N:N+2) = h2(N:N+2, N:N+2) + B(2)*F;
FN = F;
for in = 2:N
    FN = conv2(FN, F); % conv2() 二維convolution
    h2(N+1-in:N+1+in, N+1-in:N+1+in) = h2(N+1-in:N+1+in, N+1-in:N+1+in) + B(in+1) * FN; 
end
%
%
FR2 = abs(freqz2(h2, -1:2/100:1, -1:2/100:1));
XX = zeros(101, 101);
for i = 0:100
   XX(:, i+1) = (-1:2/100:1)'; 
end
YY = XX';
subplot(1, 2, 2);
plot3(XX, YY, FR2);
axis([-1, 1, -1, 1, 0, 1.1]);
xlabel('\omega_1');
ylabel('\omega_2');
zlabel('Magnitude Response');
title('2D Circularly Symmetric Filter');
%
