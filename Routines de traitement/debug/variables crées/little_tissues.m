% from microsorum 4
originaltissue2=originaltissue;
originaltissue2.c(1:7)=[];
originaltissue2.c(4:end)=[];
originaltissue2=clean(originaltissue2);
clf;
plot(originaltissue2);
axis equal;
figure(gcf);