function [X,Y]=Multilateration_LLS(Mode,rx,ry)
% Mode: 1=LLS, 2=LLS_Direct, 3=WLLS, 4=2SWLLS
x1=-10;y1=-10;z1=0;
x2=10;y2=-10;z2=0;
x3=10;y3=10;z3=0;
x4=-10;y4=10;z4=0;

rz=0;
rxULA = phased.OmnidirectionalMicrophoneElement;  %Same location as emitter
rxpos1 = [rx;ry;rz+0.0001];
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

Carrier_Fre = 96e3;             % 20 kHz
Propagation_Speed = 340;            % 340 m/s
Operating_range = 34;            % 34 m
pri = (2*Operating_range)/Propagation_Speed;
prf = 1/pri;
bw = 5e3;              % 5 kHz
fs = 9.6e4;            % 96 kHz
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
[~,ang5_transmit] = rangeangle(rxpos4,srcpos,srcax);

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


xcor12_PHAT=abs(fftshift(ifft((fft(sigr1).*conj(fft(sigr2)))./(abs(fft(sigr1)).*abs(fft(sigr2))))));
xcor13_PHAT=abs(fftshift(ifft((fft(sigr1).*conj(fft(sigr3)))./(abs(fft(sigr1)).*abs(fft(sigr3))))));
xcor14_PHAT=abs(fftshift(ifft((fft(sigr1).*conj(fft(sigr4)))./(abs(fft(sigr1)).*abs(fft(sigr4))))));
xcor15_PHAT=abs(fftshift(ifft((fft(sigr1).*conj(fft(sigr5)))./(abs(fft(sigr1)).*abs(fft(sigr5))))));

Peak_Value=zeros(1,4);
Estimated_difference=zeros(1,4);

for i=1:length(xcor12_PHAT)
    if xcor12_PHAT(i)>Peak_Value(1)
        Estimated_difference(1)=(i-0.1*fs-1);
        Peak_Value(1)=xcor12_PHAT(i);
    end
    if xcor13_PHAT(i)>Peak_Value(2)
        Estimated_difference(2)=(i-0.1*fs-1);
        Peak_Value(2)=xcor13_PHAT(i);
    end
    if xcor14_PHAT(i)>Peak_Value(3)
        Estimated_difference(3)=(i-0.1*fs-1);
        Peak_Value(3)=xcor14_PHAT(i);
    end
    if xcor15_PHAT(i)>Peak_Value(4)
        Estimated_difference(4)=(i-0.1*fs-1);
        Peak_Value(4)=xcor15_PHAT(i);
    end
end

Estimated_difference=Estimated_difference/fs;
Distance=Estimated_difference*Propagation_Speed;
d1=Distance(1);
d2=Distance(2);
d3=Distance(3);
d4=Distance(4);

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

%if Mode==2
    









end