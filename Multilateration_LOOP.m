
Error=zeros(3,9);
for Mode=1:3
    index=1;
    for r=5:5:45
        sum=0;
        for x=-r:r/10:r
            y=sqrt(r^2-x^2);
            [X,Y]=Multilateration_LLS(Mode,x,y);
            sum=(X-x)^2+(Y-y)^2;
        end
        Error(Mode,index)=sum;
        index=index+1;
    end
end        
        