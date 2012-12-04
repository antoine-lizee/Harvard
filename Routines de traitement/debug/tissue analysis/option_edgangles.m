

figure(1);
clf
plot(originaltissue2_1,'-x','Linewidth',2,'num', 'MarkerSize',4);
axis equal;
figure(2); clf
p11=plotedges(originaltissue2_1,'low');
figure(3); clf
p12=plotedges(originaltissue2_1,'low','edgeswithangles');
