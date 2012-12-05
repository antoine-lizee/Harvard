%This script creates example figures that present the final results.
G=[330 260 190 100 20 ]; % Gamma angles, approximately chosen when the dividing length is the most optimal compared to other n values.

for i=2:4
    [~, h]=test_3D_function(G(i),i,100+i,i);
end

[~, h]=test_3D_function(G(5),5.2,105,5);

maximize(2:5);
