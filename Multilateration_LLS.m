
function [X,Y]=Multilateration_LLS(Mode,rx,ry,var)
% Mode: 1=LLS, 2=LLS_Direct, 3=WLLS, 4=2SWLLS
x1=0.1;y1=0.1;z1=0;
x2=-0.1;y2=0.1;z2=0;
x3=-0.1;y3=-0.1;z3=0;
x4=0.1;y4=-0.1;z4=0;

rz=0;
rxULA = phased.OmnidirectionalMicrophoneElement;  %Same location as emitter
rxpos1 = [rx;ry;rz+0.0000001];
rxvel1 = [0;0;0];
rxax1 = azelaxes(0,0);
rxpos2 = [x1;y1;z1];  %Reference Receiver 
rxvel2 = [0;0;0];
rxax2 = azelaxes(0,0);
rxpos3 = [x2;y2;z2];
rxvel3 = [0;0;0];
rxax3 = azelaxes(0,0);
rxpos4 = [x3;y3;z3];
rxvel4 = [0;0;0];
rxax4 = azelaxes(0,0);
rxpos5 = [x4;y4;z4];
rxvel5 = [0;0;0];
rxax5 = azelaxes(0,0);

srcpos = [rx;ry;rz];
srcvel = [0;0;0];
srcax = azelaxes(0,0);
srcULA = phased.OmnidirectionalMicrophoneElement;

Carrier_Fre = 96e3;             % Hz
Propagation_Speed = 340;            % m/s
Operating_range = 68;            % m
pri = (2*Operating_range)/Propagation_Speed;
prf = 1/pri;
bw = 5e3;              % Hz
fs = 9.6e4;            % Hz
waveform = phased.LinearFMWaveform('SampleRate',fs,'SweepBandwidth',bw,'PRF',prf,'PulseWidth',pri/10);
signal = waveform();

radiator = phased.WidebandRadiator('Sensor',srcULA,'PropagationSpeed',Propagation_Speed,'SampleRate',fs,...
    'CarrierFrequency',Carrier_Fre);
collector1 = phased.WidebandCollector('Sensor',rxULA,'PropagationSpeed',Propagation_Speed,'SampleRate',fs,...
    'CarrierFrequency',Carrier_Fre);
collector2 = phased.WidebandCollector('Sensor',rxULA,'PropagationSpeed',Propagation_Speed,'SampleRate',fs,...
    'CarrierFrequency',Carrier_Fre);
collector3 = phased.WidebandCollector('Sensor',rxULA,'PropagationSpeed',Propagation_Speed,'SampleRate',fs,...
    'CarrierFrequency',Carrier_Fre);
collector4 = phased.WidebandCollector('Sensor',rxULA,'PropagationSpeed',Propagation_Speed,'SampleRate',fs,...
    'CarrierFrequency',Carrier_Fre);
collector5 = phased.WidebandCollector('Sensor',rxULA,'PropagationSpeed',Propagation_Speed,'SampleRate',fs,...
    'CarrierFrequency',Carrier_Fre);

channel1 = phased.WidebandFreeSpace('PropagationSpeed',Propagation_Speed,...
    'SampleRate',fs,'OperatingFrequency',Carrier_Fre);
channel2 = phased.WidebandFreeSpace('PropagationSpeed',Propagation_Speed,...
    'SampleRate',fs,'OperatingFrequency',Carrier_Fre);
channel3= phased.WidebandFreeSpace('PropagationSpeed',Propagation_Speed,...
    'SampleRate',fs,'OperatingFrequency',Carrier_Fre);
channel4 = phased.WidebandFreeSpace('PropagationSpeed',Propagation_Speed,...
    'SampleRate',fs,'OperatingFrequency',Carrier_Fre);
channel5 = phased.WidebandFreeSpace('PropagationSpeed',Propagation_Speed,...
    'SampleRate',fs,'OperatingFrequency',Carrier_Fre);

[~,ang1_transmit] = rangeangle(rxpos1,srcpos,srcax);
[~,ang2_transmit] = rangeangle(rxpos2,srcpos,srcax);
[~,ang3_transmit] = rangeangle(rxpos3,srcpos,srcax);
[~,ang4_transmit] = rangeangle(rxpos4,srcpos,srcax);
[~,ang5_transmit] = rangeangle(rxpos5,srcpos,srcax);

sig_transmit = radiator(signal,[ang1_transmit ang2_transmit ang3_transmit ang4_transmit ang5_transmit]);

sigp1 = channel1(sig_transmit(:,1),srcpos,rxpos1,srcvel,rxvel1);
sigp2 = channel2(sig_transmit(:,2),srcpos,rxpos2,srcvel,rxvel2);
sigp3 = channel3(sig_transmit(:,3),srcpos,rxpos3,srcvel,rxvel3);
sigp4 = channel4(sig_transmit(:,4),srcpos,rxpos4,srcvel,rxvel4);
sigp5 = channel5(sig_transmit(:,5),srcpos,rxpos5,srcvel,rxvel5);

[~,ang1_receive] = rangeangle(srcpos,rxpos1,rxax1);
[~,ang2_receive] = rangeangle(srcpos,rxpos2,rxax2);
[~,ang3_receive] = rangeangle(srcpos,rxpos3,rxax3);
[~,ang4_receive] = rangeangle(srcpos,rxpos4,rxax4);
[~,ang5_receive] = rangeangle(srcpos,rxpos5,rxax5);

sigr1 = collector1(sigp1,ang1_receive);
sigr2 = collector2(sigp2,ang2_receive);
sigr3 = collector3(sigp3,ang3_receive);
sigr4 = collector4(sigp4,ang4_receive);
sigr5 = collector5(sigp5,ang5_receive);

%Insert a Pesudo-random noise 
nr=randn(1,length(signal));
ni=randn(1,length(signal));
Noise=(1/40000)*(sqrt(2)/2)*(nr+1i*ni);

sigr2 = sigr2+Noise.';
sigr3 = sigr3+Noise.';
sigr4 = sigr4+Noise.';
sigr5 = sigr5+Noise.';

xcor01_PHAT=abs(fftshift(ifft((fft(sigr1).*conj(fft(sigr2)))./(abs(fft(sigr1)).*abs(fft(sigr2))))));
xcor02_PHAT=abs(fftshift(ifft((fft(sigr1).*conj(fft(sigr3)))./(abs(fft(sigr1)).*abs(fft(sigr3))))));
xcor03_PHAT=abs(fftshift(ifft((fft(sigr1).*conj(fft(sigr4)))./(abs(fft(sigr1)).*abs(fft(sigr4))))));
xcor04_PHAT=abs(fftshift(ifft((fft(sigr1).*conj(fft(sigr5)))./(abs(fft(sigr1)).*abs(fft(sigr5))))));

Peak_Value=zeros(1,4);
Estimated_difference=zeros(1,4);

for i=1:length(xcor01_PHAT)
    if xcor01_PHAT(i)>Peak_Value(1)
        Estimated_difference(1)=(i-0.2*fs-1);
        Peak_Value(1)=xcor01_PHAT(i);
    end
    if xcor02_PHAT(i)>Peak_Value(2)
        Estimated_difference(2)=(i-0.2*fs-1);
        Peak_Value(2)=xcor02_PHAT(i);
    end
    if xcor03_PHAT(i)>Peak_Value(3)
        Estimated_difference(3)=(i-0.2*fs-1);
        Peak_Value(3)=xcor03_PHAT(i);
    end
    if xcor04_PHAT(i)>Peak_Value(4)
        Estimated_difference(4)=(i-0.2*fs-1);
        Peak_Value(4)=xcor04_PHAT(i);
    end
end
%subplot(2,2,1)
%plot(1:length(xcor01_PHAT),xcor01_PHAT)
%subplot(2,2,2)
%plot(1:length(xcor02_PHAT),xcor02_PHAT)
%subplot(2,2,3)
%plot(1:length(xcor03_PHAT),xcor03_PHAT)
%subplot(2,2,4)
%plot(1:length(xcor04_PHAT),xcor04_PHAT)

Estimated_difference=Estimated_difference/fs;
Distance=Estimated_difference*Propagation_Speed;
d1=Distance(1);
d2=Distance(2);
d3=Distance(3);
d4=Distance(4);
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
    W=(1/4)*[1/(d1^2*var),0,0,0;
         0,1/(d2^2*var),0,0;
         0,0,1/(d3^2*var),0;
         0,0,0,1/(d4^2*var)];
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
    W=(1/4)*[1/(d1^2*var),0,0,0;
         0,1/(d2^2*var),0,0;
         0,0,1/(d3^2*var),0;
         0,0,0,1/(d4^2*var)];
    Theta1=inv(A.'*A)*A.'*b;
    X1=Theta1(1);
    Y1=Theta1(2);
    M=[2*X1,0,0;
       0,2*Y1,0;
       0,0,1];
    h=[X1,Y1,X1^2+Y1^2];
    T=inv(M*(inv(A.'*W*A)*M));
    Theta=inv(Z.'*T*Z)*Z.'*T*h.';
    X=sqrt(Theta(1));
    Y=sqrt(Theta(2));
end
% AML 
if Mode==5
    A=[-2*x1, -2*y1, 1;
       -2*x2, -2*y2, 1;
       -2*x3, -2*y3, 1;
       -2*x4, -2*y4, 1];
    b=[d1^2-x1^2-y1^2;
       d2^2-x2^2-y2^2;
       d3^2-x3^2-y3^2;
       d4^2-x4^2-y4^2]; 
    Theta1=inv(A.'*A)*A.'*b;
    X1=Theta1(1);
    Y1=Theta1(2);
    r1=sqrt((X1-x1)^2+(Y1-y1)^2);g1=(X1-x1)/(r1*(r1+d1));h1=(Y1-y1)/(r1*(r1+d1));
    r2=sqrt((X1-x2)^2+(Y1-y2)^2);g2=(X1-x2)/(r2*(r2+d2));h2=(Y1-y2)/(r2*(r2+d2));
    r3=sqrt((X1-x3)^2+(Y1-y3)^2);g3=(X1-x3)/(r3*(r3+d3));h3=(Y1-y3)/(r3*(r3+d3));
    r4=sqrt((X1-x4)^2+(Y1-y4)^2);g4=(X1-x4)/(r4*(r4+d4));h4=(Y1-y4)/(r4*(r4+d4));
    H=[g1*x1+g2*x2+g3*x3+g4*x4,g1*y1+g2*y2+g3*y3+g4*y4;
       h1*x1+h2*x2+h3*x3+h4*x4,h1*y1+h2*y2+h3*y3+h4*y4];
    k=X1^2+Y1^2;k1=x1^2+y1^2;k2=x2^2+y2^2;k3=x3^2+y3^2;k4=x4^2+y4^2;
    x=[g1*(k+k1-d1^2)+g2*(k+k2-d2^2)+g3*(k+k3-d3^2)+g4*(k+k4-d4^2);
       h1*(k+k1-d1^2)+h2*(k+k2-d2^2)+h3*(k+k3-d3^2)+h4*(k+k4-d4^2)];
    Theta=0.5*inv(H.'*H)*H.'*x;
    X=Theta(1);
    Y=Theta(2);
end
% AML Direct Form 
if Mode==6
    A=[-2*x1, -2*y1, 1;
       -2*x2, -2*y2, 1;
       -2*x3, -2*y3, 1;
       -2*x4, -2*y4, 1];
    b=[d1^2-x1^2-y1^2;
       d2^2-x2^2-y2^2;
       d3^2-x3^2-y3^2;
       d4^2-x4^2-y4^2]; 
    Theta1=inv(A.'*A)*A.'*b;
    X1=Theta1(1);
    Y1=Theta1(2);
    r1=sqrt((X1-x1)^2+(Y1-y1)^2);
    r2=sqrt((X1-x2)^2+(Y1-y2)^2);
    r3=sqrt((X1-x3)^2+(Y1-y3)^2);
    r4=sqrt((X1-x4)^2+(Y1-y4)^2);
    Sum1=(-r1*x1+d1*x1)/r1+(-r2*x2+d2*x2)/r2+(-r3*x3+d3*x3)/r3+(-r4*x4+d4*x4)/r4;
    Sum2=(d1-r1)/r1+(d2-r2)/r2+(d3-r3)/r3+(d4-r4)/r4;
    X=Sum1/Sum2;
    Y=0;
end
end
