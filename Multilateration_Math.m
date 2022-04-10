function [X,Y]=Multilateration_Math(Mode,rx,ry,SNR)

%Initialize Receiver Locations
x1=0;y1=0;
x2=1000;y2=0;
x3=1000;y3=1000;
x4=0;y4=1000;
%Calculate Real distance
r1=sqrt((rx-x1)^2+(ry-y1)^2);
r2=sqrt((rx-x2)^2+(ry-y2)^2);
r3=sqrt((rx-x3)^2+(ry-y3)^2);
r4=sqrt((rx-x4)^2+(ry-y4)^2);
%Initialize Measurement Noise
Var1=r1^2/SNR;
Var2=r2^2/SNR;
Var3=r3^2/SNR;
Var4=r4^2/SNR;
Noise1=sqrt(Var1)*randn();
Noise2=sqrt(Var2)*randn();
Noise3=sqrt(Var3)*randn();
Noise4=sqrt(Var4)*randn();
%Measured distance
d1=r1+Noise1;
d2=r2+Noise2;
d3=r3+Noise3;
d4=r4+Noise4;

% Substract Form LLS
if Mode==1
    H=[x1-x2,y1-y2;
    x1-x3,y1-y3;
    x1-x4,y1-y4];
    x=[d2^2-d1^2-(x2^2+y2^2)+(x1^2+y1^2);
    d3^2-d1^2-(x3^2+y3^2)+(x1^2+y1^2);
    d4^2-d1^2-(x4^2+y4^2)+(x1^2+y1^2)];
    Theta=0.5*inv(H.'*H)*H.'*x;
    X=Theta(1);
    Y=Theta(2);
end
% Direct Form LLS
if Mode==2
    A=[-2*x1, -2*y1, 1;
       -2*x2, -2*y2, 1;
       -2*x3, -2*y3, 1;
       -2*x4, -2*y4, 1];
    b=[d1^2-x1^2-y1^2;
       d2^2-x2^2-y2^2;
       d3^2-x3^2-y3^2;
       d4^2-x4^2-y4^2];
    Theta=inv(A.'*A)*A.'*b;
    X=Theta(1);
    Y=Theta(2);
end
% WLLS
if Mode==3
    A=[-2*x1, -2*y1, 1;
       -2*x2, -2*y2, 1;
       -2*x3, -2*y3, 1;
       -2*x4, -2*y4, 1];
    b=[d1^2-x1^2-y1^2;
       d2^2-x2^2-y2^2;
       d3^2-x3^2-y3^2;
       d4^2-x4^2-y4^2];
    W=(1/4)*[1/(d1^2*Var1),0,0,0;
         0,1/(d2^2*Var2),0,0;
         0,0,1/(d3^2*Var3),0;
         0,0,0,1/(d4^2*Var4)];
    Theta=inv(A.'*W*A)*A.'*W*b;
    X=Theta(1);
    Y=Theta(2);
end
% 2-Step WLLS
if Mode==4
    b=[d1^2-x1^2-y1^2;
       d2^2-x2^2-y2^2;
       d3^2-x3^2-y3^2;
       d4^2-x4^2-y4^2];
    Z=[1,0;0,1;1,1];
    A=[-2*x1, -2*y1, 1;
       -2*x2, -2*y2, 1;
       -2*x3, -2*y3, 1;
       -2*x4, -2*y4, 1];
    W=(1/4)*[1/(d1^2*Var1),0,0,0;
         0,1/(d2^2*Var2),0,0;
         0,0,1/(d3^2*Var3),0;
         0,0,0,1/(d4^2*Var4)];
    Theta1=inv(A.'*W*A)*A.'*W*b;
    X1=Theta1(1);
    Y1=Theta1(2);
    M=[2*X1,0,0;
       0,2*Y1,0;
       0,0,1];
    h=[X1^2,Y1^2,X1^2+Y1^2];
    T=inv(M*(inv(A.'*W*A)*M));
    Theta=inv(Z.'*T*Z)*Z.'*T*h.';
    X=sqrt(Theta(1));
    Y=sqrt(Theta(2));
end
% Cramer-Rao Lower Bound
if Mode==5
    Fisher=zeros(2,2);
    Fisher(1,1)=(rx-x1)^2/(Var1*d1)+(rx-x2)^2/(Var2*d2)+(rx-x3)^2/(Var3*d3)+(rx-x4)^2/(Var4*d4);
    Fisher(1,2)=(rx-x1)*(ry-y1)/(Var1*d1)+(rx-x2)*(ry-y2)/(Var2*d2)+(rx-x3)*(ry-y3)/(Var3*d3)+(rx-x4)*(ry-y4)/(Var4*d4);
    Fisher(2,1)=Fisher(1,2);
    Fisher(2,2)=(ry-y1)^2/(Var1*d1)+(ry-y2)^2/(Var2*d2)+(ry-y3)^2/(Var3*d3)+(ry-y4)^2/(Var4*d4);
    Fisher_inv=inv(Fisher);
    X=Fisher_inv(1,1);
    Y=Fisher_inv(2,2);
end
end