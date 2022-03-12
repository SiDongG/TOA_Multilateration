%% Multipath Cancelation Heuristic Search

%% Definition
clear;clc;close all;
fb=1e3;
bit_period=1/fb; %bit period is 0.001 second 
fc=1e4; %Carrier Frequency
fs=2e4; %Sampling Rate
bits=[1,0,1,0,0,1];
t=0:1/(fs*fb):1/fb;
for i=1:length(bits)
    analog_sig=(bits(i)*2-1)*cos(2*pi*fc*t);
    if i==1
        analog=analog_sig;
    else
        analog=[analog,analog_sig];
    end
end
signal=analog;
path0=signal;
path1=3*[zeros(1,5600),signal];
path2=0.5*[zeros(1,6600),signal];
path3=2*[zeros(1,7200),signal];
path4=0.3*[zeros(1,9400),signal];
signal_r=path0+path1(1:120006)+path2(1:120006)+path3(1:120006)+path4(1:120006);

signal_store=signal_r;
%%
%Amplitude Search range [0:0.1:5]
%Delay Search range [0:10us:500us]
%Delay Search sample [0:200:10000]
% multipath_delay(misec) [0,280us,330us,360us,470us]
% delay_samples [0,5600,6600,7200,9400]

%% Testing
Amp_range=0:0.1:5;
Delay_range=0:200:10000;
Delay_list=zeros(1,4);
Amp_list=zeros(1,4);
for loop=1:4
    Min_rms=10;Min_delay=0;Min_amp=0;
    for i=2:length(Amp_range)
        disp(i)
        amp=Amp_range(i);
        for j=2:length(Delay_range)
            delay=Delay_range(j);
            S=amp*[zeros(1,delay),signal];
            signal_test=signal_r-S(1:120006);
            if rms(signal_test)<Min_rms
                Min_rms=rms(signal_test);
                Min_delay=delay;
                Min_amp=amp;
            end
        end
    end
    Delay_list(loop)=Min_delay;
    Amp_list(loop)=Min_amp;
    S=Min_amp*[zeros(1,Min_delay),signal];
    signal_r=signal_r-S(1:120006);
end
K=3*[zeros(1,5600),signal];
signal_k=signal_store-K(1:120006);
cross=fftshift(ifft(fft(signal).*conj(fft(signal_store))));
cross2=fftshift(ifft(fft(signal).*conj(fft(signal_r))));
subplot(2,3,1)
plot(1:length(cross),cross);
subplot(2,3,2)
plot(1:length(cross2),cross2);
subplot(2,3,3)
plot(1:length(signal_store),signal_store);
subplot(2,3,4)
plot(1:length(signal_r),signal_r);
subplot(2,3,5)
plot(1:length(signal_k),signal_k);

cross3=fftshift(ifft(fft(signal).*conj(fft(signal_r))));

