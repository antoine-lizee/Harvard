function [ P ] = simuprob( N, f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
P=[];
n=0;
while n<N
    P0=rand(N-n,1)*12-6;
    p=rand(size(P0));
    P0(p>(-cos(P0*f)/2+1/2))=[];
    P=[P; P0];
    n=n+length(P0);
end
%hist(P,100);



end


