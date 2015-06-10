% James Smith & Martin Hailstone
% 02/03/2015
% Snake Game

clear all
close all
clc

%global dx
%global dy

load('ONBI_h_score.mat')

width=15;

%map=zeros(width);
window=figure('position',[100 100 500 500],'menubar','none');
set(gcf,'Renderer','OpenGL');
supercolourmap=colormap(jet(width.^2));
cmap=zeros(width^2+2,3);
while 1==1
figure(window)
map=zeros(width);
car1_x=1;
car1_y=10;
car2_x=10;
car2_y=3;

car1_length=4;
car2_length=8;

cont=1;


foodno=width^2+2;


car1_dx = 1;
car1_dy = 0;
car2_dx = 1;
car2_dy = 0;

    while cont==1;
        t_s=tic;
        
        set(window, 'KeyPressFcn', @(x,y)disp(get(window,'CurrentCharacter')));
        a = get(window,'CurrentCharacter');
        
       
        car1_y=car1_y+car2_dy;
        car1_x=car1_x+car2_dx;
        
        if car1_x>width
            car1_x=1;
        elseif car1_x<1
            car1_x=width;
        end
        
        car2_y=car2_y+car2_dy;
        car2_x=car2_x+car2_dx;
        
        if car2_x>width
            car2_x=1;
        elseif car2_x<1
            car2_x=width;
        end
       
        
        if map(car1_y,car1_x)>=1 && map(car1_y,car1_x)<foodno
            uicontrol('Position',[165 225 200 50],'style','text','string','YOU DIED','FontSize',25)
            cont=0;
        
        else
            map(car1_y,car1_x)=car1_length;
        end
        
        if map(car2_y,car2_x)>=1 && map(car2_y,car2_x)<foodno
            uicontrol('Position',[165 225 200 50],'style','text','string','YOU DIED','FontSize',25)
            cont=0;
        
        else
            map(car2_y,car2_x)=car2_length;
        end
        map(map~=foodno)=map(map~=foodno)-1;
        map(map<0)=0;
        
        %%
        cmap=zeros(width^2,3);
        % cmap(2:1+slength,:)=colormap(jet(slength));
        % cmap(2:1+slength,:)=colormap(jet(slength));
        xq=mat2gray(0:(car1_length-1))*width^2+1;
        minicolourmap=interp1(1:width^2,supercolourmap,xq,'nearest');
        cmap(2:1+car1_length,:)=minicolourmap;
        
        xq=mat2gray(0:(car2_length-1))*width^2+1;
        minicolourmap=interp1(1:width^2,supercolourmap,xq,'nearest');
        cmap(2:1+car2_length,:)=minicolourmap;
        
        cmap(end,:)=[1;1;1];
        
        %%
        
        
        image(map);colormap(cmap)
        score=(car1_length-4)*10;
         score=(car2_length-4)*10;
        text(1,1, ['Score ',num2str(score)],'color',[1 1 1])
        axis off
        pt=0.07-toc(t_s);
        pt(pt<0)=0;
        pause(pt)
        
       
    end
    
    
   
end