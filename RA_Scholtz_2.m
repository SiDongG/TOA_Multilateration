%% Robert A Scholtz MULTIpath characteristic paper method
%% Using direct correlation matrix across all delays instead of single serach 

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

path0=[signal,zeros(1,9400)];
path1=3*[zeros(1,5600),signal,zeros(1,3800)];
path2=0.8*[zeros(1,6600),signal,zeros(1,2800)];
path3=2*[zeros(1,7200),signal,zeros(1,2200)];
path4=0.5*[zeros(1,9400),signal];
signal_r=path0+path1+path2+path3+path4;
signal_store=signal_r;

%% Cancellation

W=[path1.'/3,path2.'/0.8,path3.'/2,path4.'/0.5];
C=inv(W.'*W)*W.'*signal_r.';


