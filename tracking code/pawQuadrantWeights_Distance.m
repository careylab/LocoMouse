function weights = pawQuadrantWeights_Distance(x, y, Mx, My, Px, Py, max_distance)
% max_distance is optional and always defined on the normalized coordinate
% system, meaning that the values range from 0 to 1;
if ~isvector(x) || ~isvector(y)
    error('x and y must be vectors');
end

if ~exist('max_distance')
    max_distance = 1;
else
    max_distance = min(max(max_distance,0),1);
end
   
Nx = length(x);
Ny = length(y);

if Nx ~= Ny
    error('x and y must have the same number of elements!');
end

x = (x-Mx(1))./(Mx(2)-Mx(1)+1);
y = (y-My(1))./(My(2)-My(1)+1);

weights = 1 - sqrt((x - Px).^2 + (y - Py).^2)/sqrt(2);

if max_distance ~= 1
    weights(weights < (1-max_distance)) = 0;
end



