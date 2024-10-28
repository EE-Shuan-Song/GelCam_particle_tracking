function YN = in_mask(g_x,g_y,mask)
% make the mask logical
mask = mask > 0;
if g_x >= size(mask,2) || g_y >= size(mask,1) || ...
        g_x < 1 || g_y < 1
    YN = false;
elseif mask(ceil(g_y),ceil(g_x)) > 0 && mask(floor(g_y),ceil(g_x)) > 0 && ...
        mask(ceil(g_y),floor(g_x)) > 0 && mask(floor(g_y),floor(g_x)) > 0
    YN = true;
else
    YN = false;
end

end