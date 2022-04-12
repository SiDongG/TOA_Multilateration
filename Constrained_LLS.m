
rx=200;ry=600;SNR=1e5;

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
k=vpasolve(c1*f1-0.5*m*c1*g1+(c2*f2)/(1+m*y2)+(c3*f3)/(1+m*y3)-0.5*m*(c2*g2)/(1+m*y2) ...
    -0.5*m*(c3*g3)/(1+m*y3)+e2*f2*y2/(1+m*y2)^2+e3*f3*y3/(1+m*y3)^2- ...
    0.5*m*e2*g2*y2/(1+m*y2)^2-0.5*m*e3*g3*y3/(1+m*y3)^2-...
    0.5*m*c2*f2*y2/(1+m*y2)^2-0.5*m*c3*f3*y3/(1+m*y3)^2+...
    0.25*m^2*c2*g2*y2/(1+m*y2)^2+0.25*m^2*c3*g3*y3/(1+m*y3)^2==0,m);

Estimation=zeros(3,1,5);
for i=1:5
    if imag(k(i))==0
        Estimation(:,:,i)=inv(A.'*W*A+k(i)*P)*(A.'*W*b-0.5*k(i)*q);
    end
end
Minimum=1000;
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
X=Theta(1)
Y=Theta(2)

