function [ area ] = getareageom( rho, theta, Ri, beta)
%GETAREAGEOM Summary of this function goes here
%   This function computes the area of a solution defined by its points in
%   the 2D plan.
%   One must be cautious of the order with which rho, theta, Ri and beta is
%   provided. There are 2*n subquarter, and each argument is a 2*n vector.
%   
%   See 'findfixedpoint', the function that uses 'getareageom'.
%
%Copyright Antoine Lizée @ Dumais Laboratory, Harvard University
%July 2011

polyarea=rho(1:end-1).*rho(2:end).*sin(theta)/2;
arcarea=Ri.^2.*(beta-sin(beta).*cos(beta)); 

area=sum(polyarea+arcarea);

        
        

end

