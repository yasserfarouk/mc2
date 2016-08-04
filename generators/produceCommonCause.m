function [X,Y,P,t,d]=produceCommonCause(nPatterns,patternLength,M,a,b)
% prodcues example signals with patterns that are activated based on a
% a single input can cause multiple outputs 
% causal model
% nPatterns     number of patterns to use in X and Y
% patternLength the length of a single pattern
% M             the model used as a matrix of d*n*2 elements
%               if M(:,i,1)=j this means that i in X
%               activates pattern j in Y. if j<0 then this one is ignored
%               if M(:,i,2)=M(:,i,3)=x then pattern j will be activated in Y after a
%               delay of x steps
%               if M(:,i,2)!=M(:,i,3) then pattern j will be activated in Y after a
%               delay of between M(:,i,2) and M(:,i,3) steps
%
%
    X=zeros(1,(nPatterns+1)*patternLength);
    Y=zeros(1,(nPatterns+10)*patternLength);
    P=zeros(1,nPatterns);
    for i=1:nPatterns
        P(i)=(round(5*rand(1,1)));
        X(1+(i-1)*patternLength:i*patternLength)=generatePattern(P(i),patternLength);        
    end
    [t,d]=model(P,M);
    for i=1:nPatterns
        s=+sum(d(1:i));
        Y(1+(i-1)*patternLength+s:i*patternLength+s)=generatePattern(t(i),patternLength);
    end
end
function [YP,delay]=model(XP,M)
    b=zeros(size(XP));
    a=zeros(size(XP));
    YP=zeros(size(XP));
    for d=1:size(M,1)
        for i=1:numel(XP)
            if(XP(i) ~=0) 
                YP(i)=M(XP(i),1);
                a(i)=M(XP(i),2);
                b(i)=M(XP(i),3);
            else
                YP(i)=0;
                a(i)=0;
                b(i)=0;
            end
    %         switch(XP(i))
    %             case 0
    %                 YP(i)=0;
    %                 a(i)=4;
    %                 b(i)=4;
    %             case 1 
    %                 YP(i)=3;
    %                 a(i)=4;
    %                 b(i)=4;
    %             case 2
    %                 YP(i)=5;
    %                 a(i)=1;
    %                 b(i)=1;
    %             case 3
    %                 YP(i)=4;
    %                 a(i)=7;
    %                 b(i)=7;
    %             case 4
    %                 YP(i)=5;
    %                 a(i)=16;
    %                 b(i)=16;
    %             case 5
    %                 YP(i)=2;
    %                 a(i)=1;
    %                 b(i)=1;
    %         end
        end
    end
     delay=a + (b-a) .* rand(size(XP));
end
function X=generatePattern(type,length)
    X=zeros(length,1);
    length=length-1;
    T=length/3;
    switch(type)
        case 0
            X=zeros(1,length+1);
        case 1
            n=length/T;
            X=sin(0:2*pi/T:2*pi*n);
        case 4
            n=length/T;
            X=atan(0:2*pi/T:2*pi*n);
        case 2
            X=-1:2/(length):1;
        case 5
            X=1:-2/(length):-1;
        case 3
            X=0.8*ones(1,length+1);
    end
end