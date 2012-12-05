function [beta  rho theta_b]=getparam(delta, DELTA, theta)
%GETANGLES Summary of this function goes here
%This simple function provide us with new, easy-to-use parameters for the
%figure from the set of parameters (hi,Ri,thetai) we use to get the solution.
%In particular, these new parameters are needed to compute area and
%perimeter of the cell solution. 
%N.B. : rho is calculated for Ri=1.
%
%As many functions used in the fixed-point calculation they are quite
%simple in term of what is computed, by cannot be understood at all without
%the suitable drawing provided in the documentation.
%
%For concern of speed in calculation, the verifications of the values can
%be each time commented. Nevertheless, it is strongly recommended to let
%the verification be computed for testing : it is at this point of the
%determination of the solution that some "end of solution set" is exprimed.
%
%Copyright Antoine Lizée @ Dumais Laboratory, Harvard University
%July 2011

    %Initialisation
    theta=theta*pi/180;
    beta=zeros(1,2);
      
    %Calculation of the angles of the arc circles
    R=sqrt(1+DELTA^2);
    Xi=(1-delta);
    Xj=Xi*DELTA;
    ratio=(Xi^2+R^2-Xj^2)/(2*Xi*R);
    beta(1)=acos(1/R)-acos(ratio);
    beta(2)=theta-pi/2-beta(1);
    %Verification
    beta2verif=acos(DELTA/R)-acos((Xj^2+R^2-Xi^2)/(2*Xj*R));
    if abs(beta2verif-beta(2))>0.0000001
        error('calculation of the beta angles is wrong....')
    end
    beta=beta/2; %From the angles at the center of the circle to the angle between the arc and the edge
    
    %Calculation of the legnth rho between the center and the intersection
    %of the wals:
    hi=delta;
    ei=2*sin(beta(1)); %Cordlength of the edge in arc of circle
    rho=sqrt(hi^2+ei^2-2*hi*ei*sin(beta(1)));
    %Verification
    hj=delta*DELTA;
    ej=2*sin(beta(2))*DELTA; 
    rhoverif=sqrt(hj^2+ej^2-2*hj*ej*sin(beta(2)));    
    if abs(rho-rhoverif)>0.000000001
        error('calculation of the rho length is wrong....')
    end
    
    %Calculation of the angle between the segment hi and the connection to
    %the intersection.
    theta_b(1)=acos((rho^2+hi^2-ei^2)/(2*abs(hi)*rho));
    theta_b(2)=theta-theta_b(1);
    %Verification
    theta_b2verif=acos((rho^2+hj^2-ej^2)/(2*abs(hj)*rho));
    if abs(abs(theta_b(2))-abs(theta_b2verif))>0.0000001
        error('calculation of the theta_b angle is wrong....')
    end
    
end
