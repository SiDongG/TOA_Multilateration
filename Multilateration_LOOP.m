
Error=zeros(4,10);
y=0;
for Mode=1:4
    index=1;
    disp(Mode)
    for x=2:2:20
        sum=0;
        for i=1:10
            [X,Y]=Multilateration_LLS(Mode,x,y);
            sum=(X-x)^2+(Y-y)^2;
        end
        Error(Mode,index)=sum;
        index=index+1;
    end
end        

figure;
hold on; box on;
plot(2:2:20,Error(1,:),'k-','LineWidth',1);
plot(2:2:20,Error(2,:),'g-','LineWidth',1);
plot(2:2:20,Error(3,:),'r-','LineWidth',1);
plot(2:2:20,Error(4,:),'b-','LineWidth',1);

legend('LLS1','LLS2','WLLS','2SWLLS')

        