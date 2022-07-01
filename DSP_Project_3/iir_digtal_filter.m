%
% Design of IIR Digital Filter
%
clear all;
clc;
%
n=14; % 分母個數
m=14; % 分子個數
ws=0.475*pi; % 高通濾波
wp=0.525*pi; %
delta=0.00001;
epsilon=0.00001;
tau=12;
pointw=10000; % 取樣數
M=100;
%
%
deltaw=pi/pointw;
deltaM=pi/M;
nm=n+m+1; % 係數個數
x=zeros(nm,1);
j=0+i;
%
%
pointpass=0;
pointstop=0;
PASS=zeros(0,1);
STOP=zeros(0,1);
for iw=0:pointw
    w=iw*deltaw;
    if w>=wp
        pointpass=pointpass+1;
        PASS(pointpass)=w;
    elseif w<=ws
        pointstop=pointstop+1;
        STOP(pointstop)=w;
    end
end
Whp=ones(pointpass,1);
Whs=ones(pointstop,1);
%
%
delta_x=100000;
iter=0;
while delta_x>epsilon
    iter=iter+1
    xk_1=x;
    r=zeros(nm,1);
    Qs=zeros(nm,nm);
    Qp=zeros(nm,nm);
    for ip=1:pointpass
        w=PASS(ip);
        ew=zeros(nm,1);
        for in=1:n
            ew(in)=exp(-j*tau*w)*exp(-j*in*w);
        end
        for im=0:m
            ew(in+im+1)=-exp(-j*im*w);
        end
        r=r+Whp(ip)*real(exp(j*tau*w)*ew);
        Qp=Qp+Whp(ip)*ew*ew';
    end
    for is=1:pointstop
        w=STOP(is);
        ew=zeros(nm,1);
        for im=0:m
            ew(in+im+1)=-exp(-j*im*w);
        end
        Qs=Qs+Whs(is)*ew*ew';
    end
    r=(pi-wp)*r/pointpass;
    Qp=real((pi-wp)*Qp/pointpass);
    Qs=real(ws*Qs/pointstop);
    Q=Qp+Qs;
    %
    %
    B=zeros(M+1,nm);
    for im=0:M
         w=im*deltaM;
         for in=1:n
             B(im+1,in)=-cos(in*w);
         end
    end
    d=(1-delta)*ones(M+1,1);
    %
    %
    x=quadprog(Q,r,B,d);
    a=x(1:n); %分母
    b=x(n+1:nm); %分子
    FR=abs(freqz(b,[1;a],0:pi/200:pi));
    plot(0:1/200:1,FR);
    xlabel('Normalized frequency (\omega/\pi)');
    ylabel('Magnitude response');
    pause;
    %
    %
    for ip=1:pointpass
        w=PASS(ip);
        eaw=zeros(n,1);
        for in=1:n
            eaw(in)=exp(-j*in*w);
        end
        Whp(ip)=1/(abs(1+a'*eaw))^2;
    end
    for is=1:pointstop
        w=STOP(is);
        eaw=zeros(n,1);
        for im=1:n
            eaw(im)=exp(-j*im*w);
        end
        Whs(is)=1/(abs(1+a'*eaw))^2;
    end
    delta_x=norm(x-xk_1)/norm(x)
end
%
FR=abs(freqz(b,[1;a],0:pi/200:pi));
subplot(1,2,1);
plot(0:1/200:1,FR);
xlabel('Normalized frequency (\omega/\pi)');
ylabel('Magnitude response');
%
subplot(1,2,2);
zplane(b',[1;a]');
axis([-2,2,-2,2]);
