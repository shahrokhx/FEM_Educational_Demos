%%               E D U C A T I O N A L      S N I P P E T S                     
%__________________________________________________________________________
% 
%                       Finite Element Methods
%                     Developed by SHAHROKH SHAHI 
%                           (www.sshahi.com)
%
%                   Georgia Institute of Technology
%__________________________________________________________________________
%
% Iso-parametric formulation    
% Quad elements shape functions; higher-Order Shape Functions (Q8)
%
%
%% Initial Setup 
clc
clear
% close all

format short g
format compact

% adding library folder(s) to the path
path = mfilename('fullpath');
path(end-length(mfilename):end)=[];
addpath(fullfile(path,'lib'));

%% General Information


                   %   [1]   [2]   [3]   [4]   [5]   [6]   [7]   [8] 
naturalCoordinates = [ -1    -1     1     1    -1     0     1     0
                       -1     1     1    -1     0     1     0    -1];


[nDim,nElemNode] = size(naturalCoordinates);                   
                   
xi = naturalCoordinates(1,:);
eta= naturalCoordinates(2,:);

convIndex = convhull(xi,eta);


fig = figure(1); 
clf(fig);


axes(1) = subplot(3,3,1);
hold(axes(1),'on');
grid(axes(1),'on')

patch(axes(1),xi(convIndex),eta(convIndex),'g','FaceAlpha',0.5)
for i = 1 : size(naturalCoordinates,2)
    txt = text(axes(1),xi(i),eta(i),num2str(i));
    txt.BackgroundColor = 'y';
    txt.FontSize = 10;
    txt.FontWeight = 'normal';
end

axis(axes(1),[-2 2 -2 2])
axis(axes(1),'equal')
ax(1) = line(axes(1),axes(1).XLim,[0 0]);
ax(2) = line(axes(1),[0 0],axes(1).YLim);
set(ax,'LineWidth',2,'Color','k')
axLabel(1) = text(axes(1),axes(1).XLim(2),0,'\xi');
axLabel(2) = text(axes(1),0,axes(1).YLim(2),'\eta');
set(axLabel,'FontSize',20)
set(axLabel,'BackgroundColor','y')

%% Working with Shape Functions


for iNode = 1 : nElemNode
    N = shapefunctions(nDim,nElemNode,[xi(iNode),eta(iNode)]);
end

xiRange  = -1 : .05 : 1;
etaRange = -1 : .05 : 1; 

[XI,ETA] = meshgrid(xiRange,etaRange);
Z(:,:,nElemNode) = zeros(size(XI));
for i = 1 : length(xiRange)
    for j = 1 : length(etaRange)
        N = shapefunctions(nDim,nElemNode,[xiRange(i),etaRange(j)]);
        for k = 1 : nElemNode
            Z(i,j,k) = N(k);
        end
    end
end


for k = 1 : nElemNode
    idx = k + 1;
    axes(idx) = subplot(3,3,idx);
    hold(axes(idx),'on');
    grid(axes(idx),'on');
    surfc(axes(idx),XI,ETA,Z(:,:,k),'FaceAlpha',0.8);
    
    view(axes(idx),[30,27])
    axis(axes(idx),'equal')
    
    ax(1) = line(axes(idx),axes(idx).XLim,[0 0],[0,0]);
    ax(2) = line(axes(idx),[0 0],axes(idx).YLim,[0,0]);
    ax(3) = line(axes(idx),[0 0],[0,0],axes(idx).ZLim);
    set(ax,'LineWidth',2,'Color','k')
    
    title(axes(idx),['N_',num2str(k)])
    
    patch(axes(idx),xi(convIndex),eta(convIndex),'g','FaceAlpha',0.3)
    for i = 1 :nElemNode
        txt = text(axes(idx),xi(i),eta(i),num2str(i));
        txt.BackgroundColor = 'y';
        txt.FontSize = 10;
        txt.FontWeight = 'Bold';
    end
end

%% export figure
fig.Position = [50, 50, 1400,800];
% exportgraphics(fig, [mfilename,'.pdf'], 'ContentType','vector')
% exportgraphics(fig, [mfilename,'.png'], 'ContentType','image')
