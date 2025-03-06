clear
close all

% file to save images after
% 1. red channel images
% 2. applying static masks
% 3. background removal after svd
% 4. applying dyanmic masks (if needed)

%% initialization
load time_series.mat
read_dir = list(1).folder;
if ~exist(fullfile(pwd,'bgrem/'),'dir')
    mkdir(fullfile(pwd,'bgrem/'));
end
write_dir = fullfile(pwd,'bgrem/');
static_mask = imread('static_mask_Epoch3_STT3.png');
im0 = imread(fullfile(read_dir,list(1).name));
r_im0 = im0(:,:,1);
[H,W]  = size(r_im0);
save_bg = 0;

% this requires manual determination
ind_start = 43;
ind_int1 = 255; % at 263 the bungee comes off
ind_end = 455; % last frame
bad_frame = [];

N_modes = 1;
N_img1 = 0;
N_img2 = 0;
for i = ind_start:ind_int1
    if ismember(i,bad_frame) ~= 1
        N_img1 = N_img1 + 1;
    end
end
M1 = zeros(N_img1,H*W,'single');


m = 1;
for i = ind_start:ind_int1
    if ismember(i,bad_frame) ~= 1
        j = time_index(i);
        im = imread(fullfile(read_dir,list(j).name));
        r_im = im(:,:,1); % select red channel
        r_im_masked = single(r_im).*static_mask; % apply the static masks
        M1(m,:) = single(r_im_masked(:));
        m = m + 1;
    end
end

% pod based background removal
C1 = M1*transpose(M1);
[U,S,~] = svd(C1);
mod_i = N_modes;
index = 1:N_img1;
PHI = transpose(M1)*U(index,mod_i);
PHI = PHI/sqrt(sum(PHI.^2));
TCoeff = M1*PHI;

m = 1;
for i = ind_start:ind_int1
    if ismember(i,bad_frame) ~= 1
        j = time_index(i);
        im = imread(fullfile(read_dir,list(j).name));
        % get a dynamic mask
        im = uint8(static_mask.*(double(im)));
        dynamic_mask = fast_dynamic_mask(im);
        % remove background
        A = im(:,:,1);
        tmp_im = single(A(:))';
        tmp_im = tmp_im - TCoeff(m,1) * PHI(:,1)';
        m = m + 1;
        c_img = reshape(tmp_im,H,W);
        % apply the dyanmic mask
        c_img_masked = c_img.*dynamic_mask;
        % save file
        imwrite(uint8(c_img_masked),fullfile(write_dir,['bgrem_00' num2str(i,'%03d') '.png']))


        if save_bg
            bf = zeros(size(A),'like',A);
            bf(:) = c_img(:);
            bg = A-bf;
            imwrite(uint8(bg),fullfile(write_dir,['bg_00' num2str(i,'%03d') '.png']))
        end
    end
end

clearvars M1 PHI tem_im

%% second part if the bungee problem occurs
if ind_int1 < ind_end
    for i = ind_int1+1:ind_end
        if ismember(i,bad_frame) ~= 1
            N_img2 = N_img2 + 1;
        end
    end
    M2 = zeros(N_img2,H*W,'single');

    m = 1;
    for i = ind_int1+1:ind_end
        if ismember(i,bad_frame) ~= 1
            j = time_index(i);
            im = imread(fullfile(read_dir,list(j).name));
            r_im = im(:,:,1); % select red channel
            r_im_masked = single(r_im).*static_mask; % apply the static masks
            M2(m,:) = single(r_im_masked(:));
            m = m + 1;
        end
    end

    % pod based background removal
    C2 = M2*transpose(M2);
    [U,S,~] = svd(C2);
    mod_i = N_modes;
    index = 1:N_img2;
    PHI = transpose(M2)*U(index,mod_i);
    PHI = PHI/sqrt(sum(PHI.^2));
    TCoeff = M2*PHI;

    m = 1;
    for i = ind_int1+1:ind_end
        if ismember(i,bad_frame) ~= 1
            j = time_index(i);
            im = imread(fullfile(read_dir,list(j).name));
            % get a dynamic mask
            im = uint8(static_mask.*(double(im)));
            dynamic_mask = fast_dynamic_mask(im);
            % remove background
            A = im(:,:,1);
            tmp_im = single(A(:))';
            tmp_im = tmp_im - TCoeff(m,1) * PHI(:,1)';
            m = m + 1;
            c_img = reshape(tmp_im,H,W);
            % apply the dyanmic mask
            c_img_masked = c_img.*dynamic_mask;
            % save file
            imwrite(uint8(c_img_masked),fullfile(write_dir,['bgrem_00' num2str(i,'%03d') '.png']))


            if save_bg
                bf = zeros(size(A),'like',A);
                bf(:) = c_img(:);
                bg = A-bf;
                imwrite(uint8(bg),fullfile(write_dir,['bg_00' num2str(i,'%03d') '.png']))
            end
        end
    end
end
