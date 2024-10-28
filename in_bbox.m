function [YN,rx,ry] = in_bbox(xp,yp,bbox)

% xp yp define the centroid of a particle
% bbox is a bounding box based on regionprops
% YN --> yes or no
% rx (0-1) if yes, 0 means at center, 1 means at edges on x axis
% ry (0-1) if yes, 0 means at center, 1 means at edges on y axis
% rx = ry = nan if no

bbox_x1 = bbox(1,1) + 0.5;
bbox_y1 = bbox(1,2) + 0.5;
bbox_x2 = bbox(1,1) + bbox(1,3) - 0.5;
bbox_y2 = bbox(1,2) + bbox(1,4) - 0.5;

if xp >= bbox_x1 && xp <= bbox_x2 && ...
        yp >= bbox_y1 && yp <= bbox_y2
    YN = true;
    rx = 2 * abs(xp - 0.5*(bbox_x1 + bbox_x2))/bbox(1,3);
    ry = 2 * abs(yp - 0.5*(bbox_y1 + bbox_y2))/bbox(1,4);
else
    YN = false;
    rx = nan;
    ry = nan;
end

end