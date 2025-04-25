function [hf,hp]=plotsem_alpha(x,y,semL,semU,color1,color2, a, lineStyle, lineWidth)
% has separate upper and lower sems
x=x(:)';
y=y(:)';
semL=semL(:)';
semU=semU(:)';

hf=fill([x,x(end:-1:1)],[y+semU,y(end:-1:1)-semL(end:-1:1)]',color2, 'FaceAlpha', a);
set(hf,'edgec','none'); % color of the edge of patches
hold on
hp=plot(x,y,'color',color1,'linew',lineWidth,'LineStyle',lineStyle);


