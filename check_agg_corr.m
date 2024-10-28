function agg = check_agg_corr(bbox1,bbox2,bbox3,im1,im2)

% agg = 1 means only the first particle should be kept
% agg = 2 means only the 2nd particle should be kept
% agg = 3 means it is a successful aggregate
% agg = 0 means there is issue
bbox1_x1 = bbox1(1,1) + 0.5;
bbox1_y1 = bbox1(1,2) + 0.5;
bbox1_x2 = bbox1(1,1) + bbox1(1,3) - 0.5;
bbox1_y2 = bbox1(1,2) + bbox1(1,4) - 0.5;
x1d = double(im1(bbox1_y1:bbox1_y2,bbox1_x1:bbox1_x2));

bbox2_x1 = bbox2(1,1) + 0.5;
bbox2_y1 = bbox2(1,2) + 0.5;
bbox2_x2 = bbox2(1,1) + bbox2(1,3) - 0.5;
bbox2_y2 = bbox2(1,2) + bbox2(1,4) - 0.5;
x2d = double(im1(bbox2_y1:bbox2_y2,bbox2_x1:bbox2_x2));

bbox3_x1 = bbox3(1,1) + 0.5;
bbox3_y1 = bbox3(1,2) + 0.5;
bbox3_x2 = bbox3(1,1) + bbox3(1,3) - 0.5;
bbox3_y2 = bbox3(1,2) + bbox3(1,4) - 0.5;
x3d = double(im2(bbox3_y1:bbox3_y2,bbox3_x1:bbox3_x2));
x3d_pad = zeros(size(x3d)*2);
x3d_pad(1:size(x3d,1),1:size(x3d,2)) = x3d;

bbox12_x1 = min(bbox1_x1,bbox2_x1);
bbox12_x2 = max(bbox1_x2,bbox2_x2);
bbox12_y1 = min(bbox1_y1,bbox2_y1);
bbox12_y2 = max(bbox1_y2,bbox2_y2);
x12d = double(im1(bbox12_y1:bbox12_y2,bbox12_x1:bbox12_x2));


% if the sub-aggregate is way too large, deny
% zero pad x1d and x2d
if size(x1d,1) > size(x3d_pad,1) || size(x1d,2) > size(x3d_pad,2)
    agg = 2;
    if size(x2d,1) > size(x3d_pad,1) || size(x2d,2) > size(x3d_pad,2)
        agg = 0;
    end
end
if size(x2d,1) > size(x3d_pad,1) || size(x2d,2) > size(x3d_pad,2)
    agg = 1;
    if size(x1d,1) > size(x3d_pad,1) || size(x1d,2) > size(x3d_pad,2)
        agg = 0;
    end
end


if ~exist('agg','var')
    % if bbox1 is within bbox2 or bbox2 within bbox1
    if in_bbox(bbox1_x1,bbox1_y1,bbox2) || in_bbox(bbox1_x2,bbox1_y1,bbox2) || ...
            in_bbox(bbox1_x1,bbox1_y2,bbox2) || in_bbox(bbox1_x2,bbox1_y2,bbox2)
        agg = 3;
    elseif in_bbox(bbox2_x1,bbox2_y1,bbox1) || in_bbox(bbox2_x2,bbox2_y1,bbox1) || ...
            in_bbox(bbox2_x1,bbox2_y2,bbox1) || in_bbox(bbox2_x2,bbox2_y2,bbox1)
        agg = 3;
    else
        if size(x12d,1) > size(x3d_pad,1) || size(x12d,2) > size(x3d_pad,2)
            [~,~,corr1,~] = best_corr(x1d,x3d_pad);
            [~,~,corr2,~] = best_corr(x2d,x3d_pad);

            if corr1 > corr2
                agg = 1;
            else
                agg = 2;
            end

        else
            % compare the best correlations for the 3 cases
            [~,~,corr1,~] = best_corr(x1d,x3d_pad);
            [~,~,corr2,~] = best_corr(x2d,x3d_pad);
            [~,~,corr12,~] = best_corr(x12d,x3d_pad);
            if corr12 == max([corr1,corr2,corr12]) || corr12 > 0.95
                agg = 3;
            elseif corr1 == max([corr1,corr2,corr12])
                agg = 1;
            elseif corr2 == max([corr1,corr2,corr12])
                agg = 2;
            else
                agg = 0;
            end
        end
    end
end

end