function [ M ] = findcarac( N,n )
%FINDCARAC
%   This function computes many useful values for the edges of one cell "n"
%   of a cellNetwork "N".
%
%See documentation to understand further the role of the function in the
%entire program divisionplanes. This documentation is available from the
%author or the corresponding professor Jacques Dumais.
%
%Copyright Antoine Lizée @ Dumais Laboratory, OEB Dpt of Harvard University
%June 2011 email :antoine.lizee@polytechnique.edu

M=struct('tissue',inputname(1),'cell',n,'area',cellArea(N,n));

% Flipping the cells that are not stored counter clockwise :
for i=n
    if cellArea(N,n)<0;
        N.c{i}=-N.c{i}(end:-1:1);
        disp(['Warning : the cell number ' num2str(i) ' has a clockwise orientation in the indexation of edges']);
    end
end

% Extraction of the edge information :
M.I_e=N.c{n}'; % Index of the edges
M.I_v=N.e(abs(M.I_e),1:2); % Index of the two vertices of the edge
M.a=N.e(abs(M.I_e),3); % Angle information of the edge
% Sorting of the edges
boole=M.I_e<=0;
M.a(boole)=-M.a(boole);
M.I_v(boole,:)=circshift(M.I_v(boole,:),[0,1]);

M.v1=N.v(M.I_v(:,1),:); % Coordinates of the first vertex of the edge
M.v2=N.v(M.I_v(:,2),:); % Coordinates of the second vertex of the edge


%Calculation of edge properties
M.xy=M.v2-M.v1; % Coordinates of the edge vector
M.a_e=atan2(M.xy(:,2),M.xy(:,1)); % Angle of the edge vector
M.l_e=sqrt(sum(M.xy.^2,2)); % Length of the edge
M.T=M.xy./(M.l_e*[1,1]); % Tangent vector of the edge
M.N=[-M.T(:,2) M.T(:,1)]; % normal vector, "left"-rotated (positively)-rotated


% Computation of the arc properties
M.R_e=M.l_e./(2*sin(M.a)); % curvature radius. Can take the "Inf" value
M.r_e=M.l_e./(2*tan(M.a)); % Distance from the edge to the center of the circle
M.VC=M.T.*(M.l_e/2*[1,1]) + M.N.*(M.r_e*[1,1]); % coordinates of the distance from the edge to the center
M.C=N.v(M.I_v(:,1),:) + M.VC; % Coordinates of the center
M.alpha(:,1)=pi + atan2(M.VC(:,2),M.VC(:,1)); % first angle of the arc, in [0 2pi]
VC2=M.VC-M.xy;
M.alpha(:,2)=pi + atan2(VC2(:,2),VC2(:,1)); % second angle of the arc
M.R_e=abs(M.R_e); % Forget the signum of the Radius.
M.r_e=abs(M.r_e);
M.s=sign(M.a);

%---------------------------------------------------------------------------

% This is only useful for consideration
% of arc plot, and for debugging. For speed issues, the following lines are
% commented. Indeed, it prevents from doing some useless calculation in finda 
% for having the [0, 2pi] desired interval)

% % Inversion of the angles which are in the wrong way. It occurs when alpha
% % 1 and 2 are on each side of the pi river. 
% indexalpha1=find((M.alpha(:,2)-M.alpha(:,1))>3/2*pi); % In this case, the difference is greater than an arbitrary angle, let say 3/2*pi
% indexalpha2=find((M.alpha(:,2)-M.alpha(:,1))<-3/2*pi);
% M.alpha(indexalpha1,1)=M.alpha(indexalpha1,1)+2*pi;
% M.alpha(indexalpha2,2)=M.alpha(indexalpha2,2)+2*pi;
% % Therefore, M.alpha is in the [0,4pi]

%---------------------------------------------------------------------------
    
end

% % Debug : copy in a script file and uncomment the following lines 
% %after having uncommented the previous sction of the function.
% figure(2)
% clf
% load('Microsorum1_1-treated')
% plot(originaltissue);
% axis equal
% T=[];
% for i=1:length(originaltissue.c)
%     M=findcarac(originaltissue,i);
% 
%     %disp([M.I_e M.I_v M.v1 M.v2 M.a_e M.xy M.l_e M.R_e M.r_e M.T M.N M.VC M.C])
% 
%     testp=[M.C, M.R_e, M.r_e, M.alpha];
%     testp(isinf(testp(:,1)),:)=[];
%     testp(isnan(testp(:,1)),:)=[];
%     hold on
%     plot(testp(:,1),testp(:,2),'or');
%     T=[T;testp];
%     for j=1:size(testp,1)
%         fnplt(rsmak('arc',abs(testp(j,3)),testp(j,1:2)',testp(j,5:6)));
%     end
% end



