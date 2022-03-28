x1=0.1;y1=0.1;z1=0;
x2=-0.1;y2=0.1;z2=0;
x3=-0.1;y3=-0.1;z3=0;
x4=0.1;y4=-0.1;z4=0;

rx=5; ry=5; rz=0;
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
subplot(2,2,1)
plot(1:length(xcor01_PHAT),xcor01_PHAT)
subplot(2,2,2)
plot(1:length(xcor02_PHAT),xcor02_PHAT)
subplot(2,2,3)
plot(1:length(xcor03_PHAT),xcor03_PHAT)
subplot(2,2,4)
plot(1:length(xcor04_PHAT),xcor04_PHAT)

Estimated_difference=Estimated_difference/fs;
Distance=Estimated_difference*Propagation_Speed;
d1=Distance(1);
d2=Distance(2);
d3=Distance(3);
d4=Distance(4);
r1=sqrt(4.9^2+4.9^2);
r2=sqrt(5.1^2+4.9^2);
r3=sqrt(5.1^2+5.1^2);
r4=sqrt(4.9^2+5.1^2);

sum=(r1*5.3+d1*5.3-r1*x1-d1*x1)/r1;
sum2=(r2*5.3+d2*5.3-r2*x2-d2*x2)/r2;
sum3=(r3*5.3+d3*5.3-r3*x3-d3*x3)/r3;
sum4=(r4*5.3+d4*5.3-r4*x4-d4*x4)/r4;
total=sum+sum2+sum3+sum4

