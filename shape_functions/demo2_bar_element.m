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
% Bar element shape function (2-node bar element)
% Interactive demonstration + Displacement demonstration

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
%              [1]   [2]   
coordinates = [ 2     6       %  x
                1     4  ];   %  y

[nDim,nElemNode] = size(coordinates);   


x = coordinates(1,:);
y = coordinates(2,:);


                   %   [1]    [2]  
naturalCoordinates = [ -1      1     
                        0      0   ];

xi = naturalCoordinates(1,:);
eta= naturalCoordinates(2,:);


fig = figure(1); 
clf(fig);

axes(1) = subplot(3,4,[1,2]);
set(axes(1),'Tag','global')
hold(axes(1),'on');
grid(axes(1),'on')

patch(axes(1),x,y,'b','FaceAlpha',0.5)
for i = 1 : size(coordinates,2)
    txt = text(axes(1),x(i),y(i),num2str(i));
    txt.BackgroundColor = 'c';
    txt.FontSize = 15;
    txt.FontWeight = 'Bold';
end
axis(axes(1),[2 7 0 6])
axis(axes(1),'equal')




axes(2) = subplot(3,4,[3,4]);
set(axes(2),'Tag','local')
hold(axes(2),'on');
grid(axes(2),'on')
patch(axes(2),xi,eta,'g','FaceAlpha',0.5)
plot(axes(2),xi,eta,'S','MarkerSize',10,'MarkerFaceColor','k')
for i = 1 : size(coordinates,2)
    txt = text(axes(2),xi(i),eta(i)+.5,num2str(i));
    txt.BackgroundColor = 'y';
    txt.FontSize = 15;
    txt.FontWeight = 'Bold';
end

axis(axes(2),[-2 2 -2 2])
axis(axes(2),'equal')
ax(1) = line(axes(2),axes(2).XLim,[0 0]);
% ax(2) = line(axes(2),[0 0],axes(2).YLim);
set(ax,'LineWidth',.5,'Color','k','LineStyle','--')
axLabel(1) = text(axes(2),axes(2).XLim(2),0,'\xi');
% axLabel(2) = text(axes(2),0,axes(2).YLim(2),'\eta');
set(axLabel,'FontSize',15)
% set(axLabel,'BackgroundColor','y')

%% Working with Shape Functions


for iNode = 1 : nElemNode
    N = shapefunctions(nDim-1,nElemNode,xi(iNode));
end

xiRange  = -1 : .1 : 1;
Z = zeros(length(xiRange),nElemNode);
for i = 1 : length(xiRange)
    N = shapefunctions(nDim-1,nElemNode,xiRange(i))';
    Z(i,:) = [N(2) N(1)];
end

for k = 1 : nElemNode
    idx = k + 2;
    if nElemNode == 2
        axes(idx) = subplot(3,4,[4+2*k-1,4+2*k]);
    elseif nElemNode == 3 
        axes(idx) = subplot(3,4,4+k);
    end
    set(axes(idx),'Tag','N');
    set(axes(idx),'UserData',k);
    hold(axes(idx),'on');
    grid(axes(idx),'on');
    plot(xiRange,Z(:,k))
    axis(axes(idx),'equal')
    axis(axes(idx),[-2 2 -.5 1])
    
    ax(1) = line(axes(idx),axes(idx).XLim,[0 0]);
    set(ax,'LineWidth',.5,'Color','k','LineStyle','--')
    axLabel(1) = text(axes(idx),axes(idx).XLim(2),0,'\xi');
    set(axLabel,'FontSize',10)
    set(axLabel,'BackgroundColor','none')
    
     patch(axes(idx),xi,eta,'o')
     
     for i = 1 : size(coordinates,2)
        txt = text(axes(idx),xi(i),eta(i),num2str(i));
        txt.BackgroundColor = 'y';
        txt.FontSize = 10;
        txt.FontWeight = 'Bold';
    end
end

%% Deformation Demonstration
u1 = 1;
u2 = 2;

idx = idx + 1;
axes(idx) = subplot(3,4,[9,12]);
set(axes(idx),'Tag','U');
set(axes(idx),'UserData',[u1,u2]);
hold(axes(idx),'on');
grid(axes(idx),'on');
axis(axes(idx),'equal')


uMat = zeros(length(xiRange),2);
for i = 1 : length(xiRange)
    N = shapefunctions(nDim-1,nElemNode,xiRange(i))';
    uMat(i,:) = [N(2)*u1 N(1)*u2];
end

p(1) = plot(axes(idx),xiRange,uMat(:,1));
p(2) = plot(axes(idx),xiRange,uMat(:,2));
p(3) = plot(axes(idx),xiRange,sum(uMat,2));
set(p,'LineWidth',1.5)
legend(axes(idx),{'N_1\timesu_1','N_2\timesu_2','U = \Sigma(N_i\timesu_i)'},...
                 'Location','best','AutoUpdate','off')
lin(1) = line([-1,-1],[0,uMat(1,1)]);
lin(2) = line([ 1, 1],[0,uMat(end,end)]);
set(lin,'LineStyle','-.','LineWidth',1,'Color','k')

u1s = sprintf('%g',u1);
u2s = sprintf('%g',u2);
text(axes(idx),-1.5,0.5,['u_1 = ',u1s])
text(axes(idx),1.1,1,['u_2 = ',u2s])


%% export figure
fig.Position = [50, 50, 1400,800];
% exportgraphics(fig, [mfilename,'.pdf'], 'ContentType','vector')
% exportgraphics(fig, [mfilename,'.png'], 'ContentType','image')

%% Pick a point in natural coordinates:
% xi_pick  = 0.5;
% eta_pick = 0.5;
% 
% 
% plot(axes(2),xi_pick,eta_pick,'r*','Tag','mapPoint')
% 
% N = shapefunctions(nDim,nElemNode,[xi_pick,eta_pick]);
% xyMap = coordinates * N;
% h = plot(axes(1),xyMap(1),xyMap(2),'r*','Tag','mapPoint');
% 
% 
% 
fig.UserData = coordinates;
set (gcf, 'WindowButtonMotionFcn', @mouseMove);


%%
function mouseMove (object, eventdata)

axes = object.Children;
coordinates = object.UserData;
node1 = coordinates(:,1);

dCoord = diff(coordinates');
L      = norm(dCoord);
coordinates = [0 L];
c = dCoord(2)/L;
s = dCoord(1)/L;
T = [c s; -s c];

axesGlobal = findobj(axes,'Tag','global');
axesLocal  = findobj(axes,'Tag','local');
axesN      = findobj(axes,'Tag','N');

C1 = get (axesLocal, 'CurrentPoint');
xi = C1(1,1);
eta= C1(1,2);


if abs(xi)<=1 && abs(eta)<=1

    title(axesLocal, ['\xi = ',num2str(C1(1,1))]);

    N = shapefunctions(1,2,xi);
    tmp = N(1); N(1) = N(2); N(2) = tmp; 
    xyMap = coordinates * N;
    xyMap = (T * [0 xyMap]') + node1;
    titleTxt = sprintf('(X,Y) = (%f , %f) \n',xyMap(1),xyMap(2));
    title(axesGlobal, titleTxt);
    
    points = findobj(object,'Tag','mapPoint');
    set(points,'Visible','off')
    h(1) = plot(axesLocal, xi,0,'r*','Tag','mapPoint');
    h(2) = plot(axesGlobal,xyMap(1),xyMap(2),'r*','Tag','mapPoint');
    
    
    for k = 1 : length(axesN)
        axx = findobj(axesN,'UserData',k);
        h(k) = plot(axx, xi,0,'rO','Tag','mapPoint');
        title (axx,['N_',num2str(k),' = ',num2str(N(k))])
        
        plot(axx,xi,N(k),'r.','Tag','mapPoint','MarkerSize',10)
    end
    set(h,'MarkerFaceColor','r')
    set(h,'MarkerSize',6)
    
    axesU = findobj(axes,'Tag','U');
    u = axesU.UserData;
    h(3) = plot(axesU, xi,N(1)*u(1),'r^','Tag','mapPoint');
    h(4) = plot(axesU, xi,N(2)*u(2),'rv','Tag','mapPoint');
    h(5) = plot(axesU, xi,N'*u'    ,'rH','Tag','mapPoint');
    titleTxt = sprintf ('(%6.4f * %g) + (%6.4f * %g)   =   %f',...
                        N(1),u(1),N(2),u(2),N'*u');
    title(axesU,['U   =   (N_1\timesu_1) + (N_2\timesu_2)   =   ',titleTxt])
    set(h(5),'MarkerSize',8)
    
    
    
end

end

