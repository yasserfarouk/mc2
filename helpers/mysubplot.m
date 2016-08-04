function [ f ] = mysubplot( i,j,k,eachAlone )
%MYSUBPLOT Summary of this function goes here
%   Detailed explanation goes here
if eachAlone
    f=figure;
else
    f=subplot(i,j,k);
end

end

