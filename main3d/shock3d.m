clc
clear
M=dlmread('anidata.txt');
x = M(1, :);
grids = numel(x);
[X, Y] = meshgrid(x, x);
i = 2;
j = 2;
hfig = figure
F(1)=getframe(hfig);
while(1)
    clf
    tempz = M(i:(i+grids-1), :);
    surf(X, Y, tempz);
    i = i+grids;
    view([60, 30]);
    xlabel('x')
    ylabel('y')
    zlabel('pressure')
    drawnow
    F(j) = getframe(hfig);
    j = j+1;
%     pause(0.2)
    if i>length(M)-grids
        break
    end
end
movie2avi(F,'pressure.avi','FPS',10)
    