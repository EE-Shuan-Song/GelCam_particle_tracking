function [x2d,x2d_bb] = get_x2d(im2,gcx,gcy,size_x1d,incSize)
% x2d is the second interrogation window
% x2d_bb is the bounding box location with x_min, x_max, y_min, y_max
% im2 is the raw 2nd image
% gcx, gcy is the centroid from the guess
% size_x1d is the size of x1d
% incSize is 32 by default
size_x2d = size_x1d + incSize;
x2d_bbx = round(gcx - size_x2d(2)/2);
x2d_bby = round(gcy - size_x2d(1)/2);
x2d_bbx2 = x2d_bbx + size_x2d(2);
x2d_bby2 = x2d_bby + size_x2d(1);
[im_h,im_w] = size(im2);
if x2d_bbx < 1
    x2d_bbx = 1;
    x2d_bbx2 = size_x2d(2);
end
if x2d_bby < 1
    x2d_bby = 1;
    x2d_bby2 = size_x2d(1);
end
if x2d_bbx2 > im_w
    x2d_bbx2 = im_w;
    x2d_bbx = im_w - size_x2d(2) + 1;
end
if x2d_bby2 > im_h
    x2d_bby2 = im_h;
    x2d_bby = im_h - size_x2d(1) + 1;
end
x2d = double(im2(x2d_bby:x2d_bby2,x2d_bbx:x2d_bbx2));
x2d_bb = [x2d_bbx,x2d_bby,x2d_bbx2,x2d_bby2];

end