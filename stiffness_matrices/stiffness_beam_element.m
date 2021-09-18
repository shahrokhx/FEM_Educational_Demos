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
% Calculating a beam element stiffness matrix by numerical integration
% Ref: Example 5.1 Introduction to Finite Elements in Engineering 
%

%% Initial Setup 
clc
clear
close all

format short g
format compact

% adding library folder(s) to the path
path = mfilename('fullpath');
path(end-length(mfilename):end)=[];
addpath(fullfile(path,'lib_integration'))
%% General Information

nDim = 2;      % number of dimensions (2D)
nDof = 2;      % number of degrees of freedom
nElemNode = 2; % number of nodes per elements


coordinates = [0 1 ];

% material properties
E = 200e9;
I = 4e6 * 1e-6;

%% Numerical Integration Setup
xiList = [-0.5773502692,  0.5773502692];
wList  = [1.0          ,  1.0         ];

%% Calculation of Element Stiffness Matrix

L = coordinates(2) - coordinates(1);
k = zeros(nElemNode * nDof); % initializing element local stiffness matrix

for intpt = 1 : length(xiList)
    xi = xiList(intpt);
    w  = wList(intpt);
    
    d2Hdxi2 = shapefunctionderivs(xi,L);
    
    k = k + (8*E*I/L^3) * d2Hdxi2;
end

k


%% Helper functions
function d2Hdxi2 = shapefunctionderivs(xi,le)
    
    d2Hdxi2 = [(3/2)*xi, (-1 + 3*xi)/2 * (le/2), ...
              -(3/2)*xi, ( 1 + 3*xi)/2 * (le/2)];

    d2Hdxi2 = d2Hdxi2' * d2Hdxi2;
end