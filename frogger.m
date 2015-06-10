%% frogger

clear all
close all
clc

road = ones(100,50);
[road_xdim road_ydim ] =size(road)
car1_length = 10
car1 = 1
while car1 + car1_length + 1 > 0;
        car1 = car1+1;
        road(50:60,car1:car1 + car1_length) = 0; % Car 1
        imagesc(road)
end
% while car1 + car1_length + 1 == 0; 
%         car1 = i-1
%         road(50:60,car1:car1+car1_length) = 0;
%     
%         imagesc(road);
% end
%     

%%
   
    
width=15;
%map=zeros(width);
window=figure('position',[100 100 500 500],'menubar','none');
set(gcf,'Renderer','OpenGL');
supercolourmap=colormap(jet(width.^2));
cmap=zeros(width^2+2,3);

while 1==1
figure(window)
map=zeros(width);
frontx=10;
fronty=10;
cont=1;
slength=4;


dx=1;
dy=0;


fronty=fronty+dy;
        frontx=frontx+dx;
        
        if fronty>width
            fronty=1;
        elseif fronty<1
            fronty=width;
        end
        
        if frontx>width
            frontx=1;
        elseif frontx<1;
            frontx=width;
        end
end

