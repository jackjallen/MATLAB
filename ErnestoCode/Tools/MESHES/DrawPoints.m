function p= DrawPoints
% 
% left click to add a new point
% middle click to remove the closest point
% right click to end
% 
% output= [ x1 y1
%           x2 y2
%           ...
%           xN yN]
%         

p=[];

figure(gcf);
a= axis;

while 1
  cla
  for i=1:size(p,1)
    hold on;
    plot( p(i,1) , p(i,2) , '*');
    text( p(i,1) , p(i,2) ,...
            sprintf('%d',i),'Color',[0 0 0],...
            'HorizontalAlignment', 'left','VerticalAlignment', 'bottom');
  end
  hold off;
  axis( a );

  [xc,yc,button]=ginput(1);
  if     button==3 %finish
    return
  elseif button==2 %delete the closest
    [d,i]= min( (xc-p(:,1) ).^2  + (yc-p(:,2)).^2);
    p=p([1:i-1 i+1:end],:);
  elseif button==1 %landmark
    p= [p; xc yc];
  end
end
