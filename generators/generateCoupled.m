function x = generateCoupled(system,params,L)
%GENERATECOUPLED Summary of this function goes here
%   Detailed explanation goes here
switch(system)
    case 'two'
        x=zeros(L,2);
        x(1,1)=params(1);
        x(1,2)=params(2);
        rx=params(3); ry=params(4);
        betaxy=params(5); betayx=params(6);
        for t=1:L-1
            x(t+1,1)=x(t,1)*(rx-rx*x(t,1)-betaxy*x(t,2));
            x(t+1,2)=x(t,2)*(ry-ry*x(t,2)-betayx*x(t,1));
        end
        
    case 'five'
        % diverges to infinity. TODO: find appropriate parameters
        x=zeros(5,L);
        x(:,1)=0.01*randn(5,1);
        for t=1:L-1
            x(1,t+1) = x(1,t) *(4 - 4*x(1,t) - 2* x(2,t) - 0.4* x(3,t));
            x(2,t+1) = x(2,t) *(3.1 - 0.31* x(1,t) - 3.1* x(2,t) - 0.93* x(3,t));
            x(3,t+1) = x(3,t) *(2.12 + 0.636* x(1,t) + 0.636* x(2,t) - 2.12* x(3,t));
            x(4,t+1) = x(4,t) *(3.8 - 0.111* x(1,t) - 0.011* x(2,t) + 0.131* x(3,t) - 3.8* x(4,t));
            x(5,t+1) = x(5,t) *(4.1 - 0.082* x(1,t) - 0.111* x(2,t) - 0.125* x(3,t) - 4.1* x(5,t));
        end
        x=x';
end
end

