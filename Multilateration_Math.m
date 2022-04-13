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

%Common Matrix
A=[-2*x1, -2*y1, 1;
   -2*x2, -2*y2, 1;
   -2*x3, -2*y3, 1;
   -2*x4, -2*y4, 1];

b=[d1^2-x1^2-y1^2;
   d2^2-x2^2-y2^2;
   d3^2-x3^2-y3^2;
   d4^2-x4^2-y4^2];
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
    Theta=inv(A.'*A)*A.'*b;
    X=Theta(1);
    Y=Theta(2);
end
% WLLS
if Mode==3
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
    Z=[1,0;0,1;1,1];
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
% Constrained WLLS
if Mode==5
    W=(1/4)*[1/(d1^2*Var1),0,0,0;
         0,1/(d2^2*Var2),0,0;
         0,0,1/(d3^2*Var3),0;
         0,0,0,1/(d4^2*Var4)];
    P=[1,0,0;0,1,0;0,0,0];
    q=[0;0;-1];
    %Eigenvalue decomposition, V is diagonal matrix with evalues, U is
    %evector matrix 
    [U,V]=eig(inv(A.'*W*A)*P);
    y1=V(1,1);y2=V(2,2);y3=V(3,3);
    c=q.'*U;c1=c(1);c2=c(2);c3=c(3);
    g=inv(U)*inv(A.'*W*A)*q;g1=g(1);g2=g(2);g3=g(3);
    e=b.'*W*A*U;e1=e(1);e2=e(2);e3=e(3);
    f=inv(U)*inv(A.'*W*A)*A.'*W*b;f1=f(1);f2=f(2);f3=f(3);
    syms m
    if y1==0
        k=vpasolve(c1*f1-0.5*m*c1*g1+(c2*f2)/(1+m*y2)+(c3*f3)/(1+m*y3)-0.5*m*(c2*g2)/(1+m*y2) ...
            -0.5*m*(c3*g3)/(1+m*y3)+e2*f2*y2/(1+m*y2)^2+e3*f3*y3/(1+m*y3)^2- ...
            0.5*m*e2*g2*y2/(1+m*y2)^2-0.5*m*e3*g3*y3/(1+m*y3)^2-...
            0.5*m*c2*f2*y2/(1+m*y2)^2-0.5*m*c3*f3*y3/(1+m*y3)^2+...
            0.25*m^2*c2*g2*y2/(1+m*y2)^2+0.25*m^2*c3*g3*y3/(1+m*y3)^2==0,m);
    elseif y2==0
        k=vpasolve(c2*f2-0.5*m*c2*g2+(c1*f1)/(1+m*y1)+(c3*f3)/(1+m*y3)-0.5*m*(c1*g1)/(1+m*y1) ...
            -0.5*m*(c3*g3)/(1+m*y3)+e1*f1*y1/(1+m*y1)^2+e3*f3*y3/(1+m*y3)^2- ...
            0.5*m*e1*g1*y1/(1+m*y1)^2-0.5*m*e3*g3*y3/(1+m*y3)^2-...
            0.5*m*c1*f1*y1/(1+m*y1)^2-0.5*m*c3*f3*y3/(1+m*y3)^2+...
            0.25*m^2*c1*g1*y1/(1+m*y1)^2+0.25*m^2*c3*g3*y3/(1+m*y3)^2==0,m);
    else
        k=vpasolve(c3*f3-0.5*m*c3*g3+(c2*f2)/(1+m*y2)+(c1*f1)/(1+m*y1)-0.5*m*(c2*g2)/(1+m*y2) ...
            -0.5*m*(c1*g1)/(1+m*y1)+e2*f2*y2/(1+m*y2)^2+e1*f1*y1/(1+m*y1)^2- ...
            0.5*m*e2*g2*y2/(1+m*y2)^2-0.5*m*e1*g1*y1/(1+m*y1)^2-...
            0.5*m*c2*f2*y2/(1+m*y2)^2-0.5*m*c1*f1*y1/(1+m*y1)^2+...
            0.25*m^2*c2*g2*y2/(1+m*y2)^2+0.25*m^2*c1*g1*y1/(1+m*y1)^2==0,m);
    end
    Estimation=zeros(3,1,5);
    for i=1:5
        if imag(k(i))==0
            Estimation(:,:,i)=inv(A.'*W*A+k(i)*P)*(A.'*W*b-0.5*k(i)*q);
        end
    end
    Minimum=1000;
    Min_index=0;
    for j=1:5
        if Estimation(1,1,j)~=0
            J=(A*Estimation(:,:,j)-b).'*W*(A*Estimation(:,:,j)-b);
            if J<=Minimum
                Minimum=J;
                Min_index=j;
            end
        end
    end
    Theta=Estimation(:,:,Min_index);
    X=Theta(1);
    Y=Theta(2);
end

% if Mode==5
%     b=[d1^2-x1^2-y1^2;
%        d2^2-x2^2-y2^2;
%        d3^2-x3^2-y3^2;
%        d4^2-x4^2-y4^2];
%     Z=[1,0;0,1;1,1];
%     A=[-2*x1, -2*y1, 1;
%        -2*x2, -2*y2, 1;
%        -2*x3, -2*y3, 1;
%        -2*x4, -2*y4, 1];
%     W=(1/4)*[1/(d1^2*Var1),0,0,0;
%          0,1/(d2^2*Var2),0,0;
%          0,0,1/(d3^2*Var3),0;
%          0,0,0,1/(d4^2*Var4)];
%     Theta1=inv(A.'*W*A)*A.'*W*b;
%     X1=Theta1(1);
%     Y1=Theta1(2);
%     M=[2*X1,0,0;
%        0,2*Y1,0;
%        0,0,1];
%     h=[X1^2,Y1^2,X1^2+Y1^2];
%     T=inv(M*(inv(A.'*W*A)*M));
%     Theta2=inv(Z.'*T*Z)*Z.'*T*h.';
%     X2=sqrt(Theta2(1));
%     Y2=sqrt(Theta2(2));
%     M2=[2*X2,0,0;
%        0,2*Y2,0;
%        0,0,1];
%     h2=[X2^2,Y2^2,X2^2+Y2^2];
%     T2=inv(M2*(inv(A.'*W*A)*M2));
%     Theta=inv(Z.'*T2*Z)*Z.'*T2*h2.';
%     X=sqrt(Theta(1));
%     Y=sqrt(Theta(2));

% Cramer-Rao Lower Bound
if Mode==6
    Fisher=zeros(2,2);
    Fisher(1,1)=(rx-x1)^2/(Var1*r1)+(rx-x2)^2/(Var2*r2)+(rx-x3)^2/(Var3*r3)+(rx-x4)^2/(Var4*r4);
    Fisher(1,2)=(rx-x1)*(ry-y1)/(Var1*r1)+(rx-x2)*(ry-y2)/(Var2*r2)+(rx-x3)*(ry-y3)/(Var3*r3)+(rx-x4)*(ry-y4)/(Var4*r4);
    Fisher(2,1)=Fisher(1,2);
    Fisher(2,2)=(ry-y1)^2/(Var1*r1)+(ry-y2)^2/(Var2*r2)+(ry-y3)^2/(Var3*r3)+(ry-y4)^2/(Var4*r4);
    Fisher_inv=inv(Fisher);
    X=Fisher_inv(1,1);
    Y=Fisher_inv(2,2);
    NoiseX=sqrt(X)*randn();
    NoiseY=sqrt(Y)*randn();
    X=rx+NoiseX;
    Y=ry+NoiseY;
end
% Approximate Maximum Likelihood 
if Mode==8
    W=(1/4)*[1/(d1^2*Var1),0,0,0;
         0,1/(d2^2*Var2),0,0;
         0,0,1/(d3^2*Var3),0;
         0,0,0,1/(d4^2*Var4)];
    Theta=inv(A.'*W*A)*A.'*W*b;
    X1=Theta(1);
    Y1=Theta(2);
    R1=Theta(3);

end
% Newton Raphson Method
if Mode==7
    Loop=50;
    W=(1/4)*[1/(d1^2*Var1),0,0,0;
         0,1/(d2^2*Var2),0,0;
         0,0,1/(d3^2*Var3),0;
         0,0,0,1/(d4^2*Var4)];
    Theta1=inv(A.'*W*A)*A.'*W*b;
    X1=Theta1(1);
    Y1=Theta1(2);
    Measurements=zeros(2,Loop);
    Theta1=[X1,Y1];
    for i=1:Loop
        r11=sqrt((X1-x1)^2+(Y1-y1)^2);
        r22=sqrt((X1-x2)^2+(Y1-y2)^2);
        r33=sqrt((X1-x3)^2+(Y1-y3)^2);
        r44=sqrt((X1-x4)^2+(Y1-y4)^2);
%         J=[(X1-x1)^2/sqrt((X1-x1)^2+(Y1-y1)^2),(Y1-y1)^2/sqrt((X1-x1)^2+(Y1-y1)^2);
%            (X1-x2)^2/sqrt((X1-x2)^2+(Y1-y2)^2),(Y1-y2)^2/sqrt((X1-x2)^2+(Y1-y2)^2);
%            (X1-x3)^2/sqrt((X1-x3)^2+(Y1-y3)^2),(Y1-y3)^2/sqrt((X1-x3)^2+(Y1-y3)^2);
%            (X1-x4)^2/sqrt((X1-x4)^2+(Y1-y4)^2),(Y1-y4)^2/sqrt((X1-x4)^2+(Y1-y4)^2)];
        G=[(d1-r11)*(X1-x1)/r11+(d2-r22)*(X1-x2)/r22+(d3-r33)*(X1-x3)/r33+(d3-r44)*(X1-x4)/r44;
           (d1-r11)*(Y1-y1)/r11+(d2-r22)*(Y1-y2)/r22+(d3-r33)*(Y1-y3)/r33+(d3-r44)*(Y1-y4)/r44];
        H=[(d1*(Y1-y1)^2-r11^3)/r11^3+(d2*(Y1-y2)^2-r22^3)/r22^3+(d3*(Y1-y3)^2-r33^3)/r33^3+(d4*(Y1-y4)^2-r44^3)/r44^3,...
           (d1-r11)*(X1-x1)*(Y1-y1)/r11^3+(d2-r22)*(X1-x2)*(Y1-y2)/r22^3+(d3-r33)*(X1-x3)*(Y1-y3)/r33^3+(d4-r44)*(X1-x4)*(Y1-y4)/r44^3;
           (d1-r11)*(X1-x1)*(Y1-y1)/r11^3+(d2-r22)*(X1-x2)*(Y1-y2)/r22^3+(d3-r33)*(X1-x3)*(Y1-y3)/r33^3+(d4-r44)*(X1-x4)*(Y1-y4)/r44^3,...
           (d1*(X1-x1)^2-r11^3)/r11^3+(d2*(X1-x2)^2-r22^3)/r22^3+(d3*(X1-x3)^2-r33^3)/r33^3+(d4*(X1-x4)^2-r44^3)/r44^3];
        Theta1=Theta1-inv(H)*G;
        X1=Theta1(1);Y1=Theta1(2);
        Measurements(1,i)=X1;
        Measurements(2,i)=Y1;
    end
    X=X1;Y=Y1;
end

end

