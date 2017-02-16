function [colormap1, colormap2] = assignColors(n, Gdist1, Gdist2, Xi, Xj)

colors = sort(randperm(256, length(Xi)));

colormap1 = zeros(n,1);
colormap2 = zeros(n,1);
colormap1(Xi) = colors;
colormap2(Xj) = colors;

for i=1:n
    if(colormap1(i) == 0)
        [m, index] = min(Gdist1(i,Xi));
        colormap1(i) = colormap1(Xi(index));
    end
    if(colormap2(i) == 0)
        [m, index] = min(Gdist2(i,Xj));
        colormap2(i) = colormap2(Xj(index));
    end
end
