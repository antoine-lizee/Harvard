function [ x ] = verifchamp( results, field )
% UNTITLED2 Summary of this function goes here
% Vérification de la redondance de la structure results
%   Detailed explanation goes here

y=results(1).field; 

try
    f=fieldnames('y')

for i=1:length(coleochaeteresults)
   x=(coleochaeteresults(i).planes(1).cell==i)*x
   end;
x

end

