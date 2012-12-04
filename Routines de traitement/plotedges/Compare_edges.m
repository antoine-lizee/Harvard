% Plot several contour plots :
clear p
for i=1:5
    filename=['Microsorum' int2str(i) '_1-treated'];
%     filename='zinnia1_5_22_A_ssd_010-treated';
    display(filename)
    load(filename);
    figure(i);
    set(gcf,'Position', [1 38 1120 587]);
    p(i)=plotedges(originaltissue,'box','hist','histsize','sub','no','normmean', '3d','stat');
    set(i,'Name',['Données d''analyse de Microsorum' int2str(i)]);
%     set(i,'Name','Données d''analyse de Zinnia')
    print(i,'-r300','-dmeta',['M' int2str(i) 'datanalysis']);
    figure(10+i);
    set(gcf,'Position', [1 38 1120 587]);
    
    contour3d(p(i));
    
    title(['Lignes de niveau de Microsorum' int2str(i)]);
    set(gcf,'Name',['Lignes de niveau de Microsorum' int2str(i)]);
%     title('Lignes de niveau de Zinnia');
%     set(gcf,'Name','Lignes de niveau de Zinnia');
%    print('-r300','-dmeta',['M' int2str(i) 'contourplot']);
end
