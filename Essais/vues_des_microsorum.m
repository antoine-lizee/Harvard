for i=1:5, 
    figure(i)
    load(['D:\Labwork\Programme S�bastien Besson\Microsorum' mat2str(i) '_1-treated.mat'])
    plot(idealtissue,'b'), plot(originaltissue)
end
    