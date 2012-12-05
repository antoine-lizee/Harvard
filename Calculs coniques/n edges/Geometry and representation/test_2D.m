% To test and visualize the solutions give by 'find fixedpoint'
% Very useful to draw figures
%
%   Copyright Antoine Lizée July 2011 
%   @Dumais lab, OEB dpt, Harvard University.

div=1; % With or without division
n=7;
type=2;

% figure(20);
% clf
hold all

for gamma=-70
    
    disp(gamma)
    switch type
        case 1
            [ V beta_final perim ]=findfixedpoint(n,gamma,1);
        case 2
            [ V beta_final perim ]=findfixedpoint(n+0.2,gamma,1);
    end
    [~, cell]=draw_2D(V,beta_final);
    
    if div
        M=findcarac(cell,1);
        switch type
            case 1
            [ ~,~, Vi, Ai, ~, lengthi] = findplane( M,n,n+2); % For type 1
            case 2
            [ ~,~, Vi, Ai, ~, lengthi] = findplane( M,n-1,n+1); % For type 2
        end
        P=findP(Vi,Ai(3),0:0.01:1);
        plot(P(:,1),P(:,2),'-b','Linewidth',1);
    end
    
end  

set(gca, 'YTick', [], 'XTick', [], 'XTickMode', 'manual', 'YTickMode', 'manual','Box', 'on');

