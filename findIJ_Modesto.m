function [I, J] = findIJ_Modesto(X,Y)

% define the coordinate system
xcoord = @(ii)200 + (ii-1).*400;
ycoord = @(ii)61200 - (200 + (ii-1).*400);

for ii = 1:137
   x = xcoord(ii);
   xmin = x-200;
   xmax = x+200;
   if (X >=xmin && X <=xmax) || (ii == 1 && X < xmin) || (ii == 137 && X > xmax)
       J = ii;
       break;
   end
end

for ii = 1:153
   y = ycoord(ii);
   ymin = y-200;
   ymax = y+200;
   if (Y >=ymin && Y<= ymax) || (ii == 1 && Y > ymax) || (ii == 153 && Y < ymin)
       I = ii;
       break;
   end
end