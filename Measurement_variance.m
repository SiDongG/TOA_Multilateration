x1=0;y1=0;z1=0;
Carrier_Fre = 96e3;             % Hz
Propagation_Speed = 340;            % m/s
Operating_range = 68;            % m
pri = (2*Operating_range)/Propagation_Speed;
prf = 1/pri;
bw = 5e3;              % Hz
fs = 9.6e4;
var_d=zeros(1,15);
SNR_d=zeros(1,15);
total_d=zeros(1,15);
Num_loop=10;
for loop=1:Num_loop
    for i=2:2:30
        disp(i)
        Measured_distance=zeros(1,100);
        for k=1:200
            rx=i;ry=0;rz=0;
            rxULA = phased.OmnidirectionalMicrophoneElement;
            rxpos1 = [x1;y1;z1];
            rxvel1 = [0;0;0];
            rxax1 = azelaxes(0,0);
            srcpos = [rx;ry;rz];
            srcvel = [0;0;0];
            srcax = azelaxes(0,0);
            srcULA = phased.OmnidirectionalMicrophoneElement;
            waveform = phased.LinearFMWaveform('SampleRate',fs,'SweepBandwidth',bw,'PRF',prf,'PulseWidth',pri/10);
            signal = waveform();

            radiator = phased.WidebandRadiator('Sensor',srcULA,'PropagationSpeed',Propagation_Speed,'SampleRate',fs,...
            'CarrierFrequency',Carrier_Fre);
            collector1 = phased.WidebandCollector('Sensor',rxULA,'PropagationSpeed',Propagation_Speed,'SampleRate',fs,...
                'CarrierFrequency',Carrier_Fre);
            channel1 = phased.WidebandFreeSpace('PropagationSpeed',Propagation_Speed,...
            'SampleRate',fs,'OperatingFrequency',Carrier_Fre);
            [~,ang1_transmit] = rangeangle(rxpos1,srcpos,srcax);
            sig_transmit = radiator(signal,ang1_transmit);
            sigp1 = channel1(sig_transmit(:,1),srcpos,rxpos1,srcvel,rxvel1);
            [~,ang1_receive] = rangeangle(srcpos,rxpos1,rxax1);
            sigr1 = collector1(sigp1,ang1_receive);

            nr=randn(1,length(signal));
            ni=randn(1,length(signal));
            Noise=1/40000*(sqrt(2)/2)*(nr+1i*ni);
            N0=rms(Noise)^2;

            Snr=rms(sigr1)^2/N0;
            SNR_d(i/2)=Snr;

            sigr1 = sigr1+Noise.';
            xcor01_PHAT=abs(fftshift(ifft((fft(sigr1).*conj(fft(signal)))./(abs(fft(sigr1)).*abs(fft(signal))))));

            Position=find(xcor01_PHAT==max(xcor01_PHAT));
            Estimated_difference=Position-0.2*fs-1;
            Distance=Estimated_difference*Propagation_Speed/fs;

            Measured_distance(k)=Distance;
            if k==200
                var_d(i/2)=var(Measured_distance);
            end
        end
        total_d=total_d+var_d;
    end
end
total_d=total_d/Num_loop;

%Result1=1e-4*[0.0043    0.0322    0.0282    0.0754    0.1371    0.1835    0.2139    0.2555    0.2824    0.3665    0.4889
%0.3492    0.4336    0.5999    0.6944]
%Result2=1e-4*[0.0012    0.0331    0.0357    0.0697    0.1533    0.1873    0.1872    0.2854    0.3531    0.4261    0.4658
%0.4938    0.6696    0.7042    0.7251]

