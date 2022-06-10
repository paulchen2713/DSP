%
% Design of Variable_Fractional_Order_FIR_Differintegrators_Example_1
%
% example 1. (N,  M,  w1,     w2,    p1,   p2)
%          = (40, 5, 0.05pi, 0.95pi, -0.5, +0.5)
clear all;   % clear workspace
clc;         % clear command window
N=40;
M=5;
w1=0.05*pi;
w2=0.95*pi;
p1=-0.5;
p2=0.5;
pointw=200;
pointp=60;
%
%
NH=N/2;
nma=(NH+1)*(M+1);
nmb=NH*(M+1);
deltaw=(w2-w1)/pointw;
deltap=(p2-p1)/pointp;
point=(pointw+1)*(pointp+1);
%
%
ra=zeros(nma,1);
Qa=zeros(nma,nma);
for ip=0:pointp
     p=p1+ip*deltap;
     for iw=0:pointw
          w=w1+iw*deltaw;
          c=zeros(nma,1);
          for i=0:nma-1
               n=mod(i,NH+1);
               m=floor(i/(NH+1));
               c(i+1)=p^m*cos(n*w);               
          end
     ra=ra-2*w^p*cos(p*pi/2)*c;
     Qa=Qa+c*c';
     end
end
ra=ra*(w2-w1)*(p2-p1)/point;
Qa=Qa*(w2-w1)*(p2-p1)/point;
a=-0.5*inv(Qa)*ra;

rb=zeros(nmb,1);
Qb=zeros(nmb,nmb);
for ip=0:pointp
     p=p1+ip*deltap;
     for iw=0:pointw
          w=w1+iw*deltaw;
          s=zeros(nmb,1);
          for i=0:nmb-1
               n=mod(i,NH)+1;
               m=floor(i/NH);
               s(i+1)=p^m*sin(n*w);               
          end
     rb=rb-2*w^p*sin(p*pi/2)*s;
     Qb=Qb+s*s';
     end
end
rb=rb*(w2-w1)*(p2-p1)/point;
Qb=Qb*(w2-w1)*(p2-p1)/point;
b=-0.5*inv(Qb)*rb;
%
%
a2=reshape(a,NH+1,M+1); % even part
he=zeros(N+1,M+1);
he(NH+1,:)=a2(1,:);
he(1:NH,:)=0.5*a2(NH+1:-1:2,:);
he(NH+2:N+1,:)=0.5*a2(2:NH+1,:);
%
b2=reshape(b,NH,M+1); % odd part
ho=zeros(N+1,M+1);
ho(1:NH,:)=0.5*b2(NH:-1:1,:);
ho(NH+2:N+1,:)=-0.5*b2;
%
%
h=he+ho;
%
MR=zeros(pointw+1,pointp+1);
for ip=0:pointp;
     p=p1+ip*deltap;
     hnp=h(:,1);
     for im=1:M
          hnp=hnp+h(:,im+1)*p^im;
     end
     MR(:,ip+1)=abs(freqz(hnp,1,w1:deltaw:w2));
end
%
XX=zeros(pointw+1,pointp+1);
YY=zeros(pointw+1,pointp+1);
for ip=0:pointp
     XX(:,ip+1)=(w1:deltaw:w2)/pi';
end
for iw=0:pointw
     YY(iw+1,:)=p1:deltap:p2;
end
%
plot3(XX,YY,MR);
axis([w1/pi,w2/pi,p1,p2,0,max(max(MR))]);
xlabel('Normalized Frequency');
ylabel('variable P');
zlabel('Magnitude Response');
pause;
%
for im=0:M
     MRs=abs(freqz(h(:,im+1),1,w1:deltaw:w2));
     subplot(3,3,im+1);
     plot(0:1/200:1,MRs);
     axis([0,1,0,3]);
end
