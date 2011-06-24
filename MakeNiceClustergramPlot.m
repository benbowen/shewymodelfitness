function MakeNiceClustergramPlot(z,xt,yt)
imagesc(z)
ax=gca;
view(2)
colormap hot
axis equal
colorbar
xt=regexprep(xt,'_',' ')
set(gca,'tickdir','out','yticklabel',yt,'ytick',-0.5:length(yt)+0.5,'xtick',0.5:length(xt)-0.5,'xticklabel',xt,'ydir','normal');

% fix the x position
usethesevalues=get(ax,'XTickLabel');
set(ax,'XTickLabel',[]);
putthemhere=get(ax,'XTick');
ylimits=get(ax,'Ylim');
ypos=ylimits(1)-.002*(ylimits(2)-ylimits(1));
th=text(putthemhere+0.5,ypos*ones(1,length(putthemhere)),usethesevalues,... 
'HorizontalAlignment','right','rotation',90,'parent',ax);

% fix the y position
usethesevalues=get(ax,'YTickLabel');
set(ax,'YTickLabel',[]);
putthemhere=get(ax,'YTick');
putthemhere=putthemhere(2:end-1);
xlimits=get(ax,'Xlim');
xpos=0;%xlimits(1);%*(xlimits(2)-xlimits(1));
th=text(xpos*ones(1,length(putthemhere)),putthemhere+0.5,usethesevalues,... 
'HorizontalAlignment','right','rotation',0,'parent',ax);


axis tight
grid on