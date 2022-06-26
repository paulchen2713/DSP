%
% Design_of_3D_FIR_Cone-Shaped_Filter
%
clear;
clc;
%
%
A3=zeros(N+1,N+1,N+1);
A3(1,1,1)=h3(N+1,N+1,N+1);
A3(2:N+1,1,1)=2*h3(N+2:2*N+1,N+1,N+1);
A3(1,2:N+1,1)=2*h3(N+1,N+2:2*N+1,N+1);
A3(1,1,2:N+1)=2*h3(N+1,N+1,N+2:2*N+1);
A3(2:N+1,2:N+1,1)=4*h3(N+2:2*N+1,N+2:2*N+1,N+1);
A3(2:N+1,1,2:N+1)=4*h3(N+2:2*N+1,N+1,N+2:2*N+1);
A3(1,2:N+1,2:N+1)=4*h3(N+1,N+2:2*N+1,N+2:2*N+1);
A3(2:N+1,2:N+1,2:N+1)=8*h3(N+2:2*N+1,N+2:2*N+1,N+2:2*N+1);
%
% plot iso_surface for cone filter design
%
point=32;
deltaw=2*pi/32;
HF3=zeros(point+1,point+1,point+1);
XX=zeros(point+1,point+1,point+1);
YY=zeros(point+1,point+1,point+1);
ZZ=zeros(point+1,point+1,point+1);
for i1=0:point
    w1=-pi+i1*deltaw;
    for i2=0:point
        w2=-pi+i2*deltaw;
        for i3=0:point
            w3=-pi+i3*deltaw;
            XX(i1+1,i2+1,i3+1)=w1/pi;
            YY(i1+1,i2+1,i3+1)=w2/pi;
            ZZ(i1+1,i2+1,i3+1)=w3/pi;
            for n1=0:N/2
                for n2=0:N/2
                    for n3=0:N/2
                        HF3(i1+1,i2+1,i3+1)=HF3(i1+1,i2+1,i3+1)+A3(n1+1,n2+1,n3+1)*cos(n1*w1)*cos(n2*w2)*cos(n3*w3);
                    end
                end
            end
        end
    end
end
HF3=abs(HF3);
[F,V]=isosurface(XX,YY,ZZ,HF3,0.99);
[rowV, colV]=size(V);
%
[Z1 ind]=sort(V(:,3));
X1=V(ind,1);
Y1=V(ind,2);
%
Zu=Z1(rowV/2+1:rowV,1);
Xu=X1(rowV/2+1:rowV,1);
Yu=Y1(rowV/2+1:rowV,1);
Zd=Z1(1:rowV/2,1);
Xd=X1(1:rowV/2,1);
Yd=Y1(1:rowV/2,1);
%
ti = -1:.0025:1; 
[XI,YI] = meshgrid(ti,ti);
ZU = griddata(Xu,Yu,Zu,XI,YI);
ZD = griddata(Xd,Yd,Zd,XI,YI);
%
contour3(XI,YI,ZU,35,'k'); hold % with black plot
contour3(XI,YI,ZD,35,'k');
%
view(3);
title('\theta=','FontSize',16);
set( gca,'FontSize', 16);
xlabel('\omega_1 / \pi','FontSize',16);
set(gca,'xtick',linspace(-1,1,5));
ylabel('\omega_2 / \pi','FontSize',16);
set(gca,'ytick',linspace(-1,1,5));
zlabel('\omega_3 / \pi','FontSize',16);
set(gca,'ztick',linspace(-1,1,5));
%
grid on;
hidden off;
%
