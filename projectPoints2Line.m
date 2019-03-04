function [t, PNT] = projectPoints2Line(P,L)
% t = projectPoints2Line(P,L)
% Calculates the parametric variable of the points P after they have
% been projected onto line L
% P : is a list of points [X Y]
% L: is the line defined as [X1 Y1 X2 Y2]
% t: is the parametric variable. When t = 0 the point is projected onto
%    (X1, Y1) point of line. When t = 1 the projection of the point is on
%    (X2, Y2). Therefore the projection of the points with t>=0 and t<=1
%    lay inside the line.
% PNT: The coordinates of the projected points


% https://gamedev.stackexchange.com/questions/72528/how-can-i-project-a-3d-point-onto-a-3d-line
% Project a point P onto line defined by A and B
% AB = B - A
% AP = P - A
% t = dot(AP,AB)/dot(AB,AB)
% projected point A + t*AB

% The code is actually a vectorized way to do the following
% AB = L([3 4]) - L([1 2]); 
% AP = P -L([1 2]);
% t = dot(AP,AB)/dot(AB,AB)

Ldst=sqrt((L(:,3)-L(:,1)).^2+(L(:,4)-L(:,2)).^2);
t = ((L(:,2)-P(:,2)).*(L(:,2)-L(:,4))-(L(:,1)-P(:,1)).*(L(:,3)-L(:,1)))./Ldst.^2;

pntsX = L(:,1) + t.*(L(:,3) - L(:,1));
pntsY = L(:,2) + t.*(L(:,4) - L(:,2));

PNT = [pntsX pntsY];
