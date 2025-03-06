clear
close all

dir = fullfile(pwd,'pics/raw/Fri_Aug_31_16-52-42_2018.jpg');
dir_test = fullfile(pwd,'pics/raw/Wed_Sep__5_00-37-51_2018.jpg');

bar_exist = 1;

im = imread(dir);
im_test = imread(dir_test);
imageSize = size(im);
mask0 = zeros(imageSize(1),imageSize(2));
offset = int32((imageSize(2) - imageSize(1)/3*4)/2);
mask0(:,[1:offset,imageSize(2)-offset+1:imageSize(2)])=1;

all_black = zeros(imageSize(1),imageSize(2),'uint8');
r_im = cat(3,im(:,:,1),all_black,all_black);
g_im = cat(3,all_black,im(:,:,2),all_black);
b_im = cat(3,all_black,all_black,im(:,:,3));

figure;
subplot(2,2,1);
imshow(im)
title('raw')
set(gca,'fontsize',16)
subplot(2,2,2);
imshow(r_im)
title('red')
set(gca,'fontsize',16)
subplot(2,2,3);
imshow(g_im)
title('green')
set(gca,'fontsize',16)
subplot(2,2,4);
imshow(b_im)
title('blue')
set(gca,'fontsize',16)

% get a mask

% im_t = rgb2gray(im);

im_r = im(:,:,1);
T = otsuthresh(im_r(:));
if T < 0.8
    T = T + 0.05;
end
im_bw0 = imbinarize(im_r,T);
im_bw1 = (im_bw0 + mask0)>0;
im_bw2 = bwareafilt(im_bw1,[5000 inf]);
im_bw2(:,[1,imageSize(2)]) = 1;
im_bw2([1,imageSize(1)],:) = 1;
im_bw2(imageSize(1),imageSize(2)/2) = 0;
im_bw2(1,imageSize(2)/2) = 0;
im_bw3 = imfill(im_bw2,'holes');

if bar_exist == 1

    figure;
    imshow(im_bw3);hold on;

    [H, theta, rho] = hough(im_bw0);
    peaks = houghpeaks(H, 2);
    lines = houghlines(im_bw3, theta, rho, peaks);
    extensionLength = 100;
    xy1 = zeros(2,2);
    xy2 = zeros(2,2);
    m = 1; n = 1;
    for k = 1:length(lines)
        % Get the line endpoints
        % xy1 is the starting point
        % xy2 is the ending point
        xy1(m,:) = lines(k).point1;
        m = m + 1;
        xy2(n,:) = lines(k).point2;
        n = n + 1;
    end

    line([xy1(1,1), xy2(1,1)], [xy1(1,2), xy2(1,2)], 'LineWidth', 1, 'Color', 'green');
    line([xy1(2,1), xy2(2,1)], [xy1(2,2), xy2(2,2)], 'LineWidth', 1, 'Color', 'green');

    vector1 = (xy1(1,2) - xy2(1,2))/(xy1(1,1) - xy2(1,1));
    vector2 = (xy1(2,2) - xy2(2,2))/(xy1(2,1) - xy2(2,1));
    im_bw4 = zeros(imageSize(1),imageSize(2));
    for i = 1:imageSize(2)
        linePixel_1y = round(xy1(1,2) - (xy1(1,1) - i)*vector1);
        if linePixel_1y <= imageSize(1) && linePixel_1y >= 1
            im_bw4(linePixel_1y,i) = 1;
        end
        
        linePixel_2y = round(xy1(2,2) - (xy1(2,1) - i)*vector1);
        if linePixel_2y <= imageSize(1) && linePixel_2y >= 1
            im_bw4(linePixel_2y,i) = 1;
        end
    end
    im_bw4 = imdilate(im_bw4,strel('disk',9));
    im_bw4 = (im_bw4 + im_bw3) > 0;
    im_bw4 = imfill(im_bw4,'holes');
    im_bw4 = imdilate(im_bw4,strel('disk',7));
    % figure;
    % imshow(im_bw4);


end

figure;
subplot(2,3,1)
imshow(im_bw0);
subplot(2,3,2)
imshow(im_bw1)
subplot(2,3,3)
imshow(im_bw2)
subplot(2,3,4)
imshow(im_bw3)
if bar_exist == 1
    subplot(2,3,5)
    imshow(im_bw4)
end

im_masked = uint8(double(im_test).*(~im_bw4));

figure('name','mask performance');
imshowpair(im_test,im_masked,'montage')
