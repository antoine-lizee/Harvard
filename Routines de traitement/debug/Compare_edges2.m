
% Compare different normalisations
s=zeros(5,3);
for i=1:5
    display(['Microsorum' int2str(i) '_1-treated'])
    load(['Microsorum' int2str(i) '_1-treated']);
    figure(i);
    p(i)=plotedges(originaltissue,'no','normmean', '3d');
    title(['Histogramme 3d de Microsorum' int2str(i)])
    s(i,1:2)=[p(i).z.N_h3d sum(sum(p(i).z.totz))*range(p(i).z.x)*range(p(i).z.y)];
    s(i,3)=s(i,2)/s(i,1);
end
display(s)
figure(99); plot(s(:,1),s(:,2),'.')
deb=[p(1).z.deb p(2).z.deb p(3).z.deb p(4).z.deb p(5).z.deb];


% Tips I want to store : (Previous versions of this file)

% Plot in different figures
%for i=1:5
%     load(['Microsorum' int2str(i) '_1-treated']);
%     p=plotedges(originaltissue,'histsize','normmean', 'windows');
%     current_figures=[1,4];
%     copyobj(allchild(1),figure(9));
%     for j=current_figures
%         h=figure(j*10+i)
%         copyobj(allchild(j),h);
%         clf(j);
%     end
%     
% end


% Plot in different colors
% set(0,'DefaultLineMarkerSize',2);
% colr=['y', 'r', 'c', 'm', 'g']
% for i=1:5
%     load(['Microsorum' int2str(i) '_1-treated']);   
%     figure(1);
%     p=plotedges(originaltissue,'normmean', 'windows');
%     figure(2);
%     polar(p.t,p.r,['.' colr(i)]);
%     hold on
% end
% set(0,'DefaultLineMarkerSize',6);

