function table = particle_props(framenum,im_bw,im_g,im_c)
% use bw images + gray images + color images
% gray images come after the background removal
% color images is the raw image
[h,w] = size(im_bw);

props = regionprops('table',im_bw,...
    'Area','BoundingBox',...
    'Circularity','Perimeter',...
    'MajorAxisLength','MinorAxisLength',...
    'Centroid','Image');

num_particles = height(props);
box = props.BoundingBox;
% index in each single frame
Index = (1:num_particles)';
% frame number
FrameNum = framenum*ones(num_particles,1);
% define two cells at first
GrayImage = cell(num_particles,1);
ColorImage = cell(num_particles,1);


for i = 1:num_particles
    bx = box(i,1) + 0.5;
    by = box(i,2) + 0.5;
    bx2 = bx + box(i,3) - 1;
    by2 = by + box(i,4) - 1;
    % find interrogation window properties
    if bx - 8 < 1
        ix = 1;
    else
        ix = bx - 8;
    end
    if by - 8 < 1
        iy = 1;
    else
        iy = by - 8;
    end
    if bx2 + 8 > w
        ix2 = w;
    else
        ix2 = bx2 + 8;
    end
    if by2 + 8 > h
        iy2 = h;
    else
        iy2 = by2 + 8;
    end
    % width and height of the interrogation window
    iw = ix2 - ix + 1;
    ih = iy2 - iy + 1;

    GrayImage(i,1) = mat2cell(im_g(iy:iy2,ix:ix2),ih,iw);
    ColorImage(i,1) = mat2cell(im_c(iy:iy2,ix:ix2,1:3),ih,iw,3);

end


table = addvars(props,FrameNum,Index,'Before','Area');
table = addvars(table,GrayImage,ColorImage);

end