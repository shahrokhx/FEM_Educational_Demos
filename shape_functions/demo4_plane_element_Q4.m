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
% Quad elements shape functions (4-nodes)
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
addpath(fullfile(path,'lib'))

%% General Information


nDim = 2;      % number of dimensions (2D)
nDof = 2;      % number of degrees of freedom
nElemNode = 4; % number of nodes per elements

%              [1]   [2]   [3]   [4]
coordinates = [ 2     6     5     3.5     %  x
                1     2     4     4];     %  y
            
x = coordinates(1,:);
y = coordinates(2,:);


                   %   [1]    [2]    [3]    [4]  
naturalCoordinates = [ -1     -1      1      1
                       -1      1      1     -1 ];

xi = naturalCoordinates(1,:);
eta= naturalCoordinates(2,:);


fig = figure(1); 
clf(fig);

axes(1) = subplot(2,4,[1,2]);
set(axes(1),'Tag','global')
hold(axes(1),'on');
grid(axes(1),'on')
% plot(axes(1),x,y,'O')
patch(axes(1),x,y,'b','FaceAlpha',0.5)
nodeNumbers = [1, 4, 3, 2];
for i = 1 : size(coordinates,2)
    txt = text(axes(1),x(i),y(i),num2str(nodeNumbers(i)));
    txt.BackgroundColor = 'c';
    txt.FontSize = 15;
    txt.FontWeight = 'Bold';
end
axis(axes(1),[2 7 0 6])
axis(axes(1),'equal')

axes(2) = subplot(2,4,[3,4]);
set(axes(2),'Tag','local')
hold(axes(2),'on');
grid(axes(2),'on')
% plot(axes(1),x,y,'O')
patch(axes(2),xi,eta,'g','FaceAlpha',0.5)
for i = 1 : size(coordinates,2)
    txt = text(axes(2),xi(i),eta(i),num2str(i));
    txt.BackgroundColor = 'y';
    txt.FontSize = 15;
    txt.FontWeight = 'Bold';
end

axis(axes(2),[-2 2 -2 2])
axis(axes(2),'equal')
ax(1) = line(axes(2),axes(2).XLim,[0 0]);
ax(2) = line(axes(2),[0 0],axes(2).YLim);
set(ax,'LineWidth',2,'Color','k')
axLabel(1) = text(axes(2),axes(2).XLim(2),0,'\xi');
axLabel(2) = text(axes(2),0,axes(2).YLim(2),'\eta');
set(axLabel,'FontSize',20)
set(axLabel,'BackgroundColor','y')

%% Working with Shape Functions


for iNode = 1 : nElemNode
    N = shapefunctions(nDim,nElemNode,[xi(iNode),eta(iNode)]);
end

xiRange  = -1 : .1 : 1;
etaRange = -1 : .1 : 1; 

[XI,ETA] = meshgrid(xiRange,etaRange);
Z(:,:,nElemNode) = zeros(size(XI));
for i = 1 : length(xiRange)
    for j = 1 : length(etaRange)
        for k = 1 : nElemNode
            N = shapefunctions(nDim,nElemNode,[xiRange(i),etaRange(j)]);
            Z(i,j,k) = N(k);
        end
    end
end


for k = 1 : nElemNode
    idx = k + 2;
    axes(idx) = subplot(2,4,4+k);
    set(axes(idx),'Tag','N');
    set(axes(idx),'UserData',k);
    hold(axes(idx),'on');
    grid(axes(idx),'on');
    surfc(axes(idx),XI,ETA,Z(:,:,k),'FaceAlpha',0.5);
    view(axes(idx),[30,27])
    axis(axes(idx),'equal')
    
    ax(1) = line(axes(idx),axes(idx).XLim,[0 0],[0,0]);
    ax(2) = line(axes(idx),[0 0],axes(idx).YLim,[0,0]);
    ax(3) = line(axes(idx),[0 0],[0,0],axes(idx).ZLim);
    set(ax,'LineWidth',2,'Color','k')
    axLabel(1) = text(axes(idx),axes(idx).XLim(2),0,'\xi');
    axLabel(2) = text(axes(idx),0,axes(idx).YLim(2),'\eta');
    set(axLabel,'FontSize',10)
    set(axLabel,'BackgroundColor','none')

    patch(axes(idx),xi,eta,'o','FaceAlpha',0.2)
    
    for i = 1 : size(coordinates,2)
        txt = text(axes(idx),xi(i),eta(i),num2str(i));
        txt.BackgroundColor = 'y';
        txt.FontSize = 10;
        txt.FontWeight = 'Bold';
    end
end    

%% Pick a point in natural coordinates:
xi_pick  = 0.5;
eta_pick = 0.5;


plot(axes(2),xi_pick,eta_pick,'r*','Tag','mapPoint')

N = shapefunctions(nDim,nElemNode,[xi_pick,eta_pick]);
xyMap = coordinates * N;
h = plot(axes(1),xyMap(1),xyMap(2),'r*','Tag','mapPoint');



fig.UserData = coordinates;
set (gcf, 'WindowButtonMotionFcn', @mouseMove);


%%
function mouseMove (object, eventdata)

axes = object.Children;
coordinates = object.UserData;

axesGlobal = findobj(axes,'Tag','global');
axesLocal  = findobj(axes,'Tag','local');
axesN      = findobj(axes,'Tag','N');

C2 = get (axesLocal, 'CurrentPoint');

xi = C2(1,1);
eta= C2(1,2);

C1 = get (axes(3), 'CurrentPoint');
x = C1(1,1);
y = C1(1,2);


if abs(xi)<=1 && abs(eta)<=1
    titleTxt = sprintf('(%f , %f) \n',xi,eta);
    title(axesLocal, ['(\xi,\eta)  =  ',titleTxt]);

    N = shapefunctions(2,4,[xi,eta]);
    xyMap = coordinates * N;
    
    titleTxt = sprintf('(X,Y) = (%f , %f) \n',xyMap(1),xyMap(2));
    title(axesGlobal, titleTxt);
    
    points = findobj(object,'Tag','mapPoint');
    set(points,'Visible','off')
    h(1) = plot(axesLocal, xi,eta,'r*','Tag','mapPoint');
    h(2) = plot(axesGlobal,xyMap(1),xyMap(2),'r*','Tag','mapPoint');
    
    
    tmp = N(2); N(2) = N(4); N(4)=tmp;
    for k = 1 : length(axesN)
        
        axx = findobj(axesN,'UserData',k);
        h(k) = plot(axx, xi,eta,'rO','Tag','mapPoint');
        title (axx,['N_',num2str(k),'=',num2str(N(k))])
        
        plot3(axx,xi,eta,N(k),'r.','Tag','mapPoint','MarkerSize',10)
    end
    set(h,'MarkerSize',6)
    set(h,'MarkerFaceColor','r')
end


if (x >= axes(3).XLim(1))&&(x <= axes(3).XLim(2))&&...
   (y >= axes(3).YLim(1))&&(y <= axes(3).YLim(2))


end

end



%%
function mouseMove1 (object, eventdata)

axes = object.Children;
C = get (axes(2), 'CurrentPoint');
title(axes(2), ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);


coordinates = object.UserData;

xi = C(1,1);
eta= C(1,2);

if abs(xi)<=1 && abs(eta)<=1
    N = shapefunctions(2,4,[xi,eta]);
    xyMap = coordinates * N;
    title(axes(3), ['(X,Y) = (',num2str(xyMap(1)),',',num2str(xyMap(2)),')']);
    
    points = findobj(axes(3),'Tag','mapPoint');
    set(points,'Visible','off')
    h = plot(axes(3),xyMap(1),xyMap(2),'r*');
    set(h,'Tag','mapPoint')
end

end
