%clear;clc;close all;
SNRdb_range=0:5:45;
SNR_List=10.^(SNRdb_range/10);
rx=300;
ry=600;
Error=zeros(8,length(SNR_List));
loop=100;
%Measurement=zeros(2,loop);
for Mode=1:6
    
    disp(Mode)
    for v=1:length(SNR_List)
        SNR=SNR_List(v);
        sum=0;
%         Measurement=zeros(1,loop);
        for i=1:loop
            [X,Y]=Multilateration_Math(Mode,rx,ry,SNR);
            sum=sum+(X-rx)^2+(Y-ry)^2;
%            Measurement(1,i)=X-rx;
%            Measurement(2,i)=Y-ry;
        end
%        Error(Mode,v)=var(Measurement(1,:))+var(Measurement(2,:));
        Error(Mode,v)=sum/loop;
    end

end  
for v=1:length(SNR_List)
    SNR=SNR_List(v);
    [X,Y]=Multilateration_Math(8,rx,ry,SNR);
    Error(8,v)=X+Y;
end

figure;
hold on; box on;
semilogy(SNRdb_range,log10(Error(1,:)),'k-','LineWidth',1);
semilogy(SNRdb_range,log10(Error(2,:)),'g-','LineWidth',1);
semilogy(SNRdb_range,log10(Error(3,:)),'r-','LineWidth',1);
semilogy(SNRdb_range,log10(Error(4,:)),'b-','LineWidth',1);
semilogy(SNRdb_range,log10(Error(5,:)),'LineWidth',1);
semilogy(SNRdb_range,log10(Error(6,:)),'LineWidth',1);
%semilogy(SNRdb_range,log10(Error(7,:)),'LineWidth',1);
semilogy(SNRdb_range,log10(Error(8,:)),'b--','LineWidth',1);



legend('LLS1','LLS2','WLLS','2SWLLS','Grid','CWLLLS','CRLB')