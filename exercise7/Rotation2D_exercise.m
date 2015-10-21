%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostic Medical Image Processing                      
% WS 2014/15                                           
% Exercise: 2D rotation using complex number
% NOTE: Complete the '???' lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

% Point correspondences
qk = [1 2;
      3 6;
      4 6;
      3 4];
  
pk = [-0.7 2.1;
      -2.1 6.4;
      -1.4 7.1;
      -0.7 4.9];

% To avoid solving the nonlinear problem, we can make use of complex
% number to construct a linear system equation, i.e., A*transf = b and 
% solve for transf (Check the slides how to build the system matrix A)
n = size(qk,1);

%A = [[qk(:,1);qk(:,2)],[-qk(:,2);qk(:,1)],[ones(n,1);zeros(n,1)],[zeros(n,1);ones(n,1)]];
A = [qk(:,1),-qk(:,2),ones(n,1),zeros(n,1);qk(:,2),qk(:,1),zeros(n,1),ones(n,1)];
b = [pk(:,1);pk(:,2)];

% Solve the linear equation to obtain transf, i.e., rotation and translation in
% complex number form
transf = pinv(A) * b;
r = transf(1:2)/sqrt(transf(1)^2+transf(2)^2);
t = transf(3:4);

% Rotation angle
angle = atan2(r(2),r(1));
disp('Rotation Angle (degree)')
angle*180/pi

% Build an Euler rotation matrix such that we can rotate qk, followed by 
% translation of qk to get evaluated pk
R = [cos(angle),-sin(angle);sin(angle),cos(angle)];
evalpk = qk*R'+repmat(t',size(qk,1),1);

% Compute the error between real pk and evaluated pk
error1 = norm(mean(pk-evalpk,1))


% Plot properties
markerSize = 15;
lineWidth = 3;
fontSize = 20; 

% Plot qk, pk, and evalpk
figure(1);
hold on;
grid on;
plot(qk(:,1)',qk(:,2)','b+','LineWidth',lineWidth,'MarkerSize',markerSize);
plot(pk(:,1)',pk(:,2)','g+','LineWidth',lineWidth,'MarkerSize',markerSize);
plot(evalpk(:,1)',evalpk(:,2)','r+','LineWidth',lineWidth,'MarkerSize',markerSize);
plot([min([qk(:,1)' pk(:,1)'])-5 max([qk(:,1)' pk(:,1)'])+5],[0 0],'k-','LineWidth',lineWidth);
plot([0 0],[min([qk(:,2)' pk(:,2)'])-5 max([qk(:,2)' pk(:,2)'])+5],'k-','LineWidth',lineWidth);
set(gca,'FontSize',fontSize);
axis image;
axis([min([qk(:,1)' pk(:,1)'])-5 max([qk(:,1)' pk(:,1)'])+5 min([qk(:,2)' pk(:,2)'])-5 max([qk(:,2)' pk(:,2)'])+5]);

legend('Original points', 'Rotated points', 'Calculated points');
  
