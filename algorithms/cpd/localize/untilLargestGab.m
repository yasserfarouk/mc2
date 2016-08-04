function i=untilLargestGab(x,argvarin)
% returns the index after which the largest gab in x appears. x should be
% sorted in descending order. argvarin is not used and is here just as a place holder for the
% function to be the same as all other until* functions
n=int32(numel(x));
if(n==0)
    i=0;
    return;
end
if(n==1)
    i=1;
    return;
end
i=1; m=x(1)-x(2);
for j=2:n-1
    d=x(j)-x(j+1);
    if(d)>m
        m=d;
        i=j;
    end
end