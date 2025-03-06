clear
close all

% test particle tracking
% + aggregation
% + fragmentation

load('Epoch3_STT3_flux01.mat','T_g','T_l','bad_frame','ind_start','ind_end');
load time_series.mat
load('PIVlab_results.mat')
static_mask = imread('static_mask_Epoch3_STT3.png');
test_dir = 'Epoch3_103_STT3';

frameNumbers = ind_start:ind_end;
frameNumbers(ismember(frameNumbers, bad_frame)) = [];
frames = length(frameNumbers);
draw_boundary = 0;
%% obtain particle stats
imageCell = cell(1, frames);
imagePath = sprintf("bgrem/bgrem_%05d.png", frameNumbers(1));
image = imread(imagePath);
[im_h,im_w] = size(image);

for i = 1:frames
    frameNumber = frameNumbers(i);
    imagePath = sprintf("bgrem/bgrem_%05d.png", frameNumber);
    image = imread(imagePath);
    imageCell{1, i} = image;
    bwImage = imbinarize(image, T_g);
    imageCell{2, i} = bwImage;
    bw_A = bwareafilt(imageCell{2,i},[50 inf]);
    bwImage2 = imbinarize(image, 0.15);
    bwImage2 = imdilate(bwImage2,strel('disk', 2));
    bwImage2 = imerode(bwImage2,strel('disk', 2));
    bw_hysteresis = imreconstruct(bw_A,bwImage2);
    imageCell{3,i} = bw_hysteresis;
end
Particles = table;
for i = 1:frames
    im_g = imread(fullfile(test_dir,'bgrem',['bgrem_' num2str(frameNumbers(i),'%05d') '.png']));
    im_c = imread(fullfile(test_dir,'pics/raw',list(time_index(frameNumbers(i))).name));

    if draw_boundary == 1
        [B,L] = bwboundaries(imageCell{3,i},'noholes');
        for k = 1:length(B)
            boundary = B{k};
            for j = 1:length(boundary)
                im_c(boundary(j,1), boundary(j,2), 1) = 255;im_c(boundary(j,1), boundary(j,2), 2) = 0;im_c(boundary(j,1), boundary(j,2), 3) = 0;
            end
        end
    end

    props1 = particle_props(frameNumbers(i),imageCell{3,i},im_g,im_c);
    Particles = [Particles;props1];
end

save('particle_id.mat','Particles')
