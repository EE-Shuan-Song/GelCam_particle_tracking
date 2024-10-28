function mask = fast_dynamic_mask(im)

r_im = im(:,:,1);
T = otsuthresh(nonzeros(r_im(:)));
im_bw = imbinarize(r_im,T);
im_bw2 = bwareafilt(im_bw,[40000 inf]);
im_bw2 = imdilate(im_bw2,strel('disk',9));
mask = ~im_bw2;


end