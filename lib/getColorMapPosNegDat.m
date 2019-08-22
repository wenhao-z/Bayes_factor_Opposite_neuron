function [cmap, cIdxZero, colortick]= getColorMapPosNegDat(cAxis, cLength)

if nargin < 2
   cLength = 64; 
end

colorLim = cAxis;
colortick = linspace(colorLim(1), colorLim(2), cLength);

% cSpec = lines(2);
cSpec = cool(2);
% cSpec = jet(2);

cSpec = flip(cSpec, 1);
% cSpec = [1,0,0; ...
%     0,0,1];

% Adjust the intensity of the color limit according to the max of value.
% Ensure that the gradients of color bar towards two directions are the
% same.

colorRatio = 3;
if abs(colorLim(1)) < abs(colorLim(2))
    cSpec(1,:) = (cSpec(1,:) * abs(colorLim(1))*colorRatio + ones(1,3) * abs(colorLim(2)))./ ...
        (abs(colorLim(1))*colorRatio + abs(colorLim(2)));
else
    cSpec(2,:) = (cSpec(2,:) * abs(colorLim(2))*colorRatio + ones(1,3) * abs(colorLim(1)))./ ...
        (abs(colorLim(1)) + abs(colorLim(2))*colorRatio );
end

[~, cIdxZero] = min(abs(colortick));

cmap = zeros(cLength, 3);

cmap(cIdxZero,:) = ones(1,3);
for iter = 1:3
    cmap(1:cIdxZero, iter) = linspace(cSpec(1,iter), 1, cIdxZero);
    cmap(cIdxZero:end, iter) = linspace(1, cSpec(2,iter), cLength+1-cIdxZero);
end
