load('tissue_test')
n=40;
k=2;
l=5;
new=1;
for h=1:n
    disp(h);
    f=ceil(h/4);
    figure(f);
    subplot(2,2,h-(f-1)*4);
    %clf;
    if new
    tissue=originaltissue3;
    tissue.e([1,2,5],3)=rand(3,1)*2-1;
    tissue.v([1,2],1)=rand(2,1)+[-0.5; 0.5];
    eval(['tissue' num2str(h) '=tissue ;']);
    end
    name=['tissue' num2str(h) ];
    eval(['tissue=' name ';']);
    M=findcarac(tissue,1);
    cla;
    plot(tissue,'-r*','Linewidth',2);
    axis equal
    [ ~, s1(1) ]=finds2(1,M,l,k);
    [ ~, s1(2) ]=finds2(0,M,l,k);
    for i=[0,1,s1]
        [ s0 ~, a V A]=finda(i,M,k,l);
        P=findP(V,A(3),0:0.01:1);
        plot(P(:,1),P(:,2),'-g','Linewidth',1);
%         text(V(1:2,1),V(1:2,2),{num2str(s1,'\\color{black}  % .1f'); num2str(s0,'\\color{black}  % .1f')});
         text(V(2,1),V(2,2),num2str(a,'\\color{green} \\bf % .2f  \\rm'),'VerticalAlignment','Top');
    end
    [ sa sb V A af]=findplane(M,k,l);
    if isempty(sa), disp(['no solution for the area constraint for tissue: ' num2str(h)]); title('no sol'); continue; end
    if isempty(V), disp(['no intersection for tissue: ' num2str(h)]); title('no int'); continue; end
    P=findP(V,A(3),0:0.01:1);
    plot(P(:,1),P(:,2),'-b','Linewidth',2);
    text(V(1:2,1),V(1:2,2),{num2str(sb,'\\color{blue}  % .1f'); num2str(sa,'\\color{blue}  % .1f')});
    text(V(2,1),V(2,2),num2str((af-0.5)*10^3,'\\color{blue} \\bf % .2f  \\rm'),'VerticalAlignment','Top');
end
maximize(1:ceil(n/4));
% delete(findobj('Type','Figure'));
%clear('-regexp','^tissue\w*');