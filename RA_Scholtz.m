%% Robert A Scholtz MULTIpath characteristic paper method
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

path0=[signal,zeros(1,9400)];
path1=3*[zeros(1,5600),signal,zeros(1,3800)];
path2=0.8*[zeros(1,6600),signal,zeros(1,2800)];
path3=2*[zeros(1,7200),signal,zeros(1,2200)];
path4=0.5*[zeros(1,9400),signal];
signal_r=path0+path1+path2+path3+path4;
signal_store=signal_r;

%% Amplitude Estimation Test
%w=[zeros(1,5600),signal,zeros(1,3800)].';

%c=inv(w.'*w)*w.'*signal_r.';

Delay_list=zeros(1,4);
Amp_list=zeros(1,4);
subplot(1,2,1)
Cross_ref=abs(fftshift(ifft((fft(path0).*conj(fft(path0)))./(abs(fft(path0)).*abs(fft(path0))))));
plot(1:length(Cross_ref),Cross_ref)
Start=find(Cross_ref==max(Cross_ref));
for loop=1:4
    Cross=abs(fftshift(ifft((fft(signal_r).*conj(fft(path0)))./(abs(fft(signal_r)).*abs(fft(path0))))));
    x=find(Cross==max(Cross));
    delay=x-Start;
    if delay==0
        break
    end
    Delay_list(loop)=delay;
    w=[zeros(1,delay),signal,zeros(1,9400-delay)].';
    A=inv(w.'*w)*w.'*signal_r.';
    Amp_list(loop)=A;
    S=A*w.';
    signal_r=signal_r-S;
end

Cross_final=abs(fftshift(ifft((fft(signal_r).*conj(fft(path0)))./(abs(fft(signal_r)).*abs(fft(path0))))));
subplot(1,2,2)
plot(1:length(Cross_final),Cross_final)


