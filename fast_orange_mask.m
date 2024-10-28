function mask = fast_orange_mask(im)

hsvImage = rgb2hsv(im);
orangeMask = (hsvImage(:,:,1) >= 0) & (hsvImage(:,:,1) <= 0.25) ...
    & (hsvImage(:,:,3) > 0.3) & (hsvImage(:,:,2) > 0.1);
mask = bwareafilt(orangeMask,[10000 inf]);
mask = imdilate(mask,strel('disk',5));
mask = ~mask;

end