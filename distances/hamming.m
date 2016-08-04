function [ d ] = hamming( s1,s2 )
%HAMMING Summary of this function goes here
%   Detailed explanation goes here
s=(s1-s2);
d=numel(find(s~=0));
end

