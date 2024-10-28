clear
close all

% test particle tracking
% + aggregation
% + fragmentation

load('Epoch3_STT3_flux01.mat','T_g','T_l','bad_frame','ind_start','ind_end');
load time_series.mat
load('/Users/yxs/Documents/Research/Marine Snow Camera/GelCam/Datasets/EXPORTS 2018 (updated)/Epoch3_103_STT3/PIV/PIVlab_results.mat')
static_mask = imread('/Users/yxs/Documents/Research/Marine Snow Camera/GelCam/Datasets/EXPORTS 2018 (updated)/Epoch3_103_STT3/static_mask_Epoch3_STT3.png');
test_dir = '/Users/yxs/Documents/Research/Marine Snow Camera/GelCam/Datasets/EXPORTS 2018 (updated)/Epoch3_103_STT3';

frameNumbers = ind_start:ind_end;
frameNumbers(ismember(frameNumbers, bad_frame)) = [];
frames = length(frameNumbers);
draw_boundary = 0;
%% obtain particle stats
imageCell = cell(1, frames);
imagePath = sprintf("/Users/yxs/Documents/Research/Marine Snow Camera/GelCam/Datasets/EXPORTS 2018 (updated)/Epoch3_103_STT3/bgrem/bgrem_%05d.png", frameNumbers(1));
image = imread(imagePath);
[im_h,im_w] = size(image);

for i = 1:frames
    frameNumber = frameNumbers(i);
    imagePath = sprintf("/Users/yxs/Documents/Research/Marine Snow Camera/GelCam/Datasets/EXPORTS 2018 (updated)/Epoch3_103_STT3/bgrem/bgrem_%05d.png", frameNumber);
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

%% try to obtain the interpolated velocity
PIV_ind = ind_start:ind_end-1;
PIV_ind(ismember(PIV_ind, bad_frame)) = [];

velocityCell = cell(2,frames-1);
guessC = cell(2,frames-1);
for i = 1:frames-1 % or i = 1 -> the 2nd last frame
    % particle index in the entire table "Particles"
    ind_p = find(Particles.FrameNum == frameNumbers(i));
    num_p = length(ind_p);
    est_u = zeros(num_p,1);
    est_v = zeros(num_p,1);
    meshX = x{PIV_ind == frameNumbers(i),1};
    meshY = y{PIV_ind == frameNumbers(i),1};
    meshU = u_filtered{PIV_ind == frameNumbers(i),1};
    meshV = v_filtered{PIV_ind == frameNumbers(i),1};

    for j = ind_p(1):ind_p(end)
        x_p = Particles.Centroid(j,1);
        y_p = Particles.Centroid(j,2);

        % region 1
        if x_p > min(meshX(:)) && x_p < max(meshX(:)) && ...
                y_p > min(meshY(:)) && y_p < max(meshY(:))
            est_u(j-ind_p(1)+1,1) = interp2(meshX,meshY,meshU,x_p,y_p,'spline');
            est_v(j-ind_p(1)+1,1) = interp2(meshX,meshY,meshV,x_p,y_p,'spline');

            % region 2
        elseif x_p <= min(meshX(:)) && y_p > min(meshY(:)) && y_p < max(meshY(:))
            est_u(j-ind_p(1)+1,1) = interp1(meshY(:,1),meshU(:,1),y_p,'linear');
            est_v(j-ind_p(1)+1,1) = interp1(meshY(:,1),meshV(:,1),y_p,'linear');
        elseif x_p >= max(meshX(:)) && y_p > min(meshY(:)) && y_p < max(meshY(:))
            est_u(j-ind_p(1)+1,1) = interp1(meshY(:,end),meshU(:,end),y_p,'linear');
            est_v(j-ind_p(1)+1,1) = interp1(meshY(:,end),meshV(:,end),y_p,'linear');

            % region 3
        elseif y_p <= min(meshY(:)) && x_p > min(meshX(:)) && x_p < max(meshX(:))
            est_u(j-ind_p(1)+1,1) = interp1(meshX(1,:),meshU(1,:),x_p,'linear');
            est_v(j-ind_p(1)+1,1) = interp1(meshX(1,:),meshV(1,:),x_p,'linear');
        elseif y_p >= max(meshY(:)) && x_p > min(meshX(:)) && x_p < max(meshX(:))
            est_u(j-ind_p(1)+1,1) = interp1(meshX(1,:),meshU(1,:),x_p,'linear');
            est_v(j-ind_p(1)+1,1) = interp1(meshX(1,:),meshV(1,:),x_p,'linear');

            % region 4
        elseif x_p <= min(meshX(:)) && y_p <= min(meshY(:))
            est_u(j-ind_p(1)+1,1) = meshU(1,1);
            est_v(j-ind_p(1)+1,1) = meshV(1,1);
        elseif x_p >= max(meshX(:)) && y_p <= min(meshY(:))
            est_u(j-ind_p(1)+1,1) = meshU(end,1);
            est_v(j-ind_p(1)+1,1) = meshV(end,1);
        elseif x_p <= min(meshX(:)) && y_p >= max(meshY(:))
            est_u(j-ind_p(1)+1,1) = meshU(1,end);
            est_v(j-ind_p(1)+1,1) = meshV(1,end);
        elseif x_p >= max(meshX(:)) && y_p >= max(meshY(:))
            est_u(j-ind_p(1)+1,1) = meshU(end,end);
            est_v(j-ind_p(1)+1,1) = meshV(end,end);
        end
    end

    velocityCell{1,i} = est_u;
    velocityCell{2,i} = est_v;
    % guessed centroids in the next frame
    guessC{1,i} = Particles.Centroid(ind_p,1) + velocityCell{1,i};
    guessC{2,i} = Particles.Centroid(ind_p,2) + velocityCell{2,i};
end
%% start tracking
% step 1: initialization
% ==> here it is frame #36
ind_p = find(Particles.FrameNum == frameNumbers(1));

% easy initialization of A
A = num2cell(ind_p);

% % apply the same area threshold >= 100 to initialize
% A = {};
% for i = 1:ind_p(end)
%     if Particles.Area(i) >= 100
%         A{length(A)+1,1} = i;
%     end
% end

% define a dist that is the distance from guessC to exact C
% dist = []; m = 1;
% n6 = 1;

% step 2: run the code
% -1: cannot find
% -2: too small
% -3: zooplankton
% -4: entering non-masked areas
% -5: no strong correlation
for frame_i = 2:frames

    % get the next bgrem gray image before the particle for loop
    im2 = imread(fullfile(test_dir,'bgrem',['bgrem_' num2str(frameNumbers(frame_i),'%05d') '.png']));

    % get the particle id in 2nd intterogation window
    ind_p = find(Particles.FrameNum == frameNumbers(frame_i));

    % only try to track the last particle id
    for p_id = 1:length(A)

        % let's say -1 is the ending of the particle tracking
        if A{p_id}(end) > 0

            % x1d_id is the main index from table Particles
            x1d_id = A{p_id}(end);

            % now we have the 1st interrogation window
            x1d = double(cell2mat(Particles.GrayImage(x1d_id)));

            % guessed centroid in the next frame
            gcx = guessC{1,frame_i-1}(Particles.Index(x1d_id));
            gcy = guessC{2,frame_i-1}(Particles.Index(x1d_id));

            size_x1d = size(x1d);

            % get x2d
            if max(size_x1d) <= 48 % bounding box < 32
                [x2d,x2d_bb] = get_x2d(im2,gcx,gcy,size_x1d,16);
            else
                [x2d,x2d_bb] = get_x2d(im2,gcx,gcy,size_x1d,32);
            end

            % fft & ifft method to get the best correlation
            [maxy,maxx,corr,R_fft] = best_corr(x1d,x2d);

           
            % get the tracked center
            xp = floor(maxx + x2d_bb(1) + size_x1d(2)/2);
            yp = floor(maxy + x2d_bb(2) + size_x1d(1)/2);

            % initialize the number of successful trackings
            YN_tracked = 0;
            p_tracked = [];
            r_tracked = [];

            if ~isnan(corr)
                % recall the particle id from the 2nd bgrem image
                for p2_id = ind_p(1):ind_p(end)

                    % check if tracked center is in any bbox
                    [YN,rx,ry] = in_bbox(xp,yp,Particles.BoundingBox(p2_id,:));

                    if YN == 1
                        YN_tracked = YN_tracked + 1;
                        p_tracked = [p_tracked p2_id];
                        r_tracked = [r_tracked rx+ry];
                        % dist(m) =  sqrt((Particles.Centroid(p2_id,1) - gcx)^2 + ...
                        %     (Particles.Centroid(p2_id,2) - gcy)^2);
                        % pxl_g = double(cell2mat(Particles.GrayImage(A{p_id}(end))));
                        % pxl = sort(pxl_g(:),'descend');
                        %
                        % % only obtain the mean gray pixel values from
                        % % 51 - area
                        % mean_g_1(m) = mean(pxl(51:Particles.Area(A{p_id}(end))));
                        % m = m + 1;


                    end
                end
            end

            % only 1 particle tracked
            if YN_tracked == 1
                A{p_id} = [A{p_id} p_tracked(1)];

                % more than 1 particles tracked
            elseif YN_tracked > 1
                A{p_id} = [A{p_id} p_tracked(r_tracked == min(r_tracked))];

                % if nothing gets tracked
            elseif YN_tracked == 0

                % check if the guessed center is in the mask or not
                YN_mask = in_mask(gcx,gcy,static_mask);
                if YN_mask ~= 1

                    % end of the tracking
                    % -1 means out of mask
                    A{p_id} = [A{p_id} -1];
                else
                    % obtain the coordinates of 4 corners of the bbox
                    bbox = Particles.BoundingBox(A{p_id}(end),:);
                    bbox_x1 = bbox(1,1) + 0.5;
                    bbox_y1 = bbox(1,2) + 0.5;
                    bbox_x2 = bbox(1,1) + bbox(1,3) - 0.5;
                    bbox_y2 = bbox(1,2) + bbox(1,4) - 0.5;
                    YN_bbox_mask = in_mask(bbox_x1,bbox_y1,static_mask) + ...
                        in_mask(bbox_x1,bbox_y2,static_mask) + ...
                        in_mask(bbox_x2,bbox_y1,static_mask) + ...
                        in_mask(bbox_x2,bbox_y2,static_mask);

                    % check if 4 corners are all within the mask
                    if YN_bbox_mask < 4
                        A{p_id} = [A{p_id} -1];

                        % weak correlation
                    elseif corr < 0.5
                        A{p_id} = [A{p_id} -2];

                        % if particle is too small in area
                    elseif Particles.Area(A{p_id}(end)) <= 100
                        A{p_id} = [A{p_id} -3];

                    elseif Particles.Area(A{p_id}(end)) > 2000
                        A{p_id} = [A{p_id} -5];

                    else
                        % obtain pixel grayscale values and sort
                        pxl_g = double(cell2mat(Particles.GrayImage(A{p_id}(end))));
                        pxl = sort(pxl_g(:),'descend');

                        % only obtain the mean gray pixel values from
                        % 51 - area
                        mean_g_int = mean(pxl(51:Particles.Area(A{p_id}(end))));
                        if mean_g_int/255 <= T_g && corr >= 0.7
                            % mean_g_6(n6) = mean_g_int;
                            % n6 = n6 + 1;
                            A{p_id} = [A{p_id} -4];
                        else

                            % all other situations
                            A{p_id} = [A{p_id} -6];
                        end
                    end
                end
            end
        end
    end

    % we have to update cell A every time
    % find remaining p_id in 2nd bgrem image
    for p2_id = ind_p(1):ind_p(end)

        % (the 1st area >= 100)
        % && Particles.Area(p2_id) >= 100
        if ~any(cellfun(@(x) ismember(p2_id, x), A))
            A{end+1,1} = p2_id;
        end
    end
    disp(frame_i)
end


%% overlap sanity check
% last element of the tracking
A_end = zeros(length(A),1);
for i = 1:length(A)
    if A{i,1}(end) < 0
        A_end(i) = A{i,1}(end-1);
    else
        A_end(i) = A{i,1}(end);
    end
end

% find the same row
ovlp_row0 = {};
m = 1;
for i = 1:length(A_end)
    if size(find(A_end == A_end(i)),1) > 1
        if ismember(i,[ovlp_row0{:};]) ~= 1
            ovlp_row0{m,1} = (find(A_end == A_end(i)))';
            m = m + 1;
        end
    end
end
ovlp_row = zeros(6,1); m = 1;
for i = 1:length(ovlp_row0)
    if size(ovlp_row0{i},2) == 2
        ovlp_row(m,1:2) = ovlp_row0{i};
        m = m + 1;
    else
        for j = 2:size(ovlp_row0{i},2)
            ovlp_row(m,1:2) = [ovlp_row0{i}(1),ovlp_row0{i}(j)];
            m = m + 1;
        end
    end
end

% find the element in these rows
for i = 1:length(ovlp_row)
    same_element = intersect(A{ovlp_row(i,1),1},A{ovlp_row(i,2),1});
    same_element = same_element(same_element>0);
    % #1 & #2 sub-aggregates
    % #3 aggregate
    ovlp_row(i,3) = A{ovlp_row(i,1),1}(find(A{ovlp_row(i,1),1} == min(same_element)) - 1);
    ovlp_row(i,4) = A{ovlp_row(i,2),1}(find(A{ovlp_row(i,2),1} == min(same_element)) - 1);
    ovlp_row(i,5) = min(same_element);

    area1 = Particles.Area(ovlp_row(i,3));
    area2 = Particles.Area(ovlp_row(i,4));
    area3 = Particles.Area(ovlp_row(i,5));

    ovlp = check_agg_area(area1,area2,area3);

    if ovlp == 0
        bbox1(1,1:4) = Particles.BoundingBox(ovlp_row(i,3),:);
        bbox2(1,1:4) = Particles.BoundingBox(ovlp_row(i,4),:);
        bbox3(1,1:4) = Particles.BoundingBox(ovlp_row(i,5),:);
        im1 = imread(fullfile(test_dir,'bgrem',['bgrem_' num2str(Particles.FrameNum(ovlp_row(i,3)),'%05d') '.png']));
        im2 = imread(fullfile(test_dir,'bgrem',['bgrem_' num2str(Particles.FrameNum(ovlp_row(i,5)),'%05d') '.png']));
        % this code: agg = 3 means strong evidence of aggregation
        ovlp = check_agg_corr(bbox1,bbox2,bbox3,im1,im2);
    end


    ovlp_row(i,6) = ovlp;
end

remove_row = [];
for i = 1:length(ovlp_row)
    switch ovlp_row(i,6)
        case 1
            A_row = ovlp_row(i,2);
            % clean all identical elements from that row
            % and run the particle matching for that particle
            % excluding the aggregate
            ind_clean = find(A{A_row} == ovlp_row(i,4));
            if ~isempty(ind_clean)
                A{A_row} = A{A_row}(1:ind_clean);
                % run the particle tracking again
                % but exclude the old wrong match
                frame_i = find(frameNumbers == Particles.FrameNum(ovlp_row(i,5)));
                im2 = imread(fullfile(test_dir,'bgrem',['bgrem_' num2str(frameNumbers(frame_i),'%05d') '.png']));
                bbox(1,1:4) = Particles.BoundingBox(ovlp_row(i,5),:);
                bbox_x1 = bbox(1,1) + 0.5;
                bbox_y1 = bbox(1,2) + 0.5;
                bbox_x2 = bbox(1,1) + bbox(1,3) - 0.5;
                bbox_y2 = bbox(1,2) + bbox(1,4) - 0.5;
                mask = cell2mat(Particles.Image(ovlp_row(i,5)));
                mask = ~mask;
                im2(bbox_y1:bbox_y2,bbox_x1:bbox_x2) = uint8(double(im2(bbox_y1:bbox_y2,bbox_x1:bbox_x2)).*mask);
                ind_p = find(Particles.FrameNum == frameNumbers(frame_i));
                p_id = A_row;
                x1d_id = A{p_id}(end);
                x1d = double(cell2mat(Particles.GrayImage(x1d_id)));
                gcx = guessC{1,frame_i-1}(Particles.Index(x1d_id));
                gcy = guessC{2,frame_i-1}(Particles.Index(x1d_id));
                size_x1d = size(x1d);
                if max(size_x1d) <= 48
                    [x2d,x2d_bb] = get_x2d(im2,gcx,gcy,size_x1d,16);
                else
                    [x2d,x2d_bb] = get_x2d(im2,gcx,gcy,size_x1d,32);
                end
                [maxy,maxx,corr,R_fft] = best_corr(x1d,x2d);
                xp = floor(maxx + x2d_bb(1) + size_x1d(2)/2);
                yp = floor(maxy + x2d_bb(2) + size_x1d(1)/2);
                YN_tracked = 0;
                p_tracked = [];
                r_tracked = [];
                if ~isnan(corr)
                    for p2_id = ind_p(1):ind_p(end)
                        if p2_id ~= ovlp_row(i,5)
                            [YN,rx,ry] = in_bbox(xp,yp,Particles.BoundingBox(p2_id,:));
                            if YN == 1
                                YN_tracked = YN_tracked + 1;
                                p_tracked = [p_tracked p2_id];
                                r_tracked = [r_tracked rx+ry];
                            end
                        end
                    end
                end
                if YN_tracked == 1
                    A{p_id} = [A{p_id} p_tracked(1)];
                elseif YN_tracked > 1
                    A{p_id} = [A{p_id} p_tracked(r_tracked == min(r_tracked))];
                elseif YN_tracked == 0
                    YN_mask = in_mask(gcx,gcy,static_mask);
                    if YN_mask ~= 1
                        A{p_id} = [A{p_id} -1];
                    else
                        bbox = Particles.BoundingBox(A{p_id}(end),:);
                        bbox_x1 = bbox(1,1) + 0.5;
                        bbox_y1 = bbox(1,2) + 0.5;
                        bbox_x2 = bbox(1,1) + bbox(1,3) - 0.5;
                        bbox_y2 = bbox(1,2) + bbox(1,4) - 0.5;
                        YN_bbox_mask = in_mask(bbox_x1,bbox_y1,static_mask) + ...
                            in_mask(bbox_x1,bbox_y2,static_mask) + ...
                            in_mask(bbox_x2,bbox_y1,static_mask) + ...
                            in_mask(bbox_x2,bbox_y2,static_mask);
                        if YN_bbox_mask < 4
                            A{p_id} = [A{p_id} -1];
                        elseif corr < 0.5
                            A{p_id} = [A{p_id} -2];
                        elseif Particles.Area(A{p_id}(end)) <= 100
                            A{p_id} = [A{p_id} -3];
                        elseif Particles.Area(A{p_id}(end)) > 2000
                            A{p_id} = [A{p_id} -5];
                        else
                            pxl_g = double(cell2mat(Particles.GrayImage(A{p_id}(end))));
                            pxl = sort(pxl_g(:),'descend');
                            mean_g_int = mean(pxl(51:Particles.Area(A{p_id}(end))));
                            if mean_g_int/255 <= T_g && corr >= 0.7
                                A{p_id} = [A{p_id} -4];
                            else
                                A{p_id} = [A{p_id} -6];
                            end
                        end
                    end
                end
                % re-write that row and remove old trackings
                if A{p_id}(end) > 0

                    % check if it is getting frames from a new line
                    m = 0;
                    for j = 1:length(A)
                        if A{p_id}(end) == A{j,1}(1)
                            A{p_id} = [A{p_id}(1:end-1) A{j}(1:end)];
                            remove_row = [remove_row j];
                            m = m + 1;
                        end
                    end
                    if m == 0
                        for j = 1:length(A)
                            if ismember(A{p_id,1}(end),A{j})
                                ind = find(A{j,1} == A{p_id,1}(end));
                                A{p_id} = [A{p_id}(1:end-1) A{j}(ind:end)];
                            end
                        end
                    end
                end
            end


        case 2
            A_row = ovlp_row(i,1);
            ind_clean = find(A{A_row} == ovlp_row(i,3));
            if ~isempty(ind_clean)
                A{A_row} = A{A_row}(1:ind_clean);
                % run the particle tracking again
                % but exclude the old wrong match
                frame_i = find(frameNumbers == Particles.FrameNum(ovlp_row(i,5)));
                im2 = imread(fullfile(test_dir,'bgrem',['bgrem_' num2str(frameNumbers(frame_i),'%05d') '.png']));
                bbox(1,1:4) = Particles.BoundingBox(ovlp_row(i,5),:);
                bbox_x1 = bbox(1,1) + 0.5;
                bbox_y1 = bbox(1,2) + 0.5;
                bbox_x2 = bbox(1,1) + bbox(1,3) - 0.5;
                bbox_y2 = bbox(1,2) + bbox(1,4) - 0.5;
                mask = cell2mat(Particles.Image(ovlp_row(i,5)));
                mask = ~mask;
                im2(bbox_y1:bbox_y2,bbox_x1:bbox_x2) = uint8(double(im2(bbox_y1:bbox_y2,bbox_x1:bbox_x2)).*mask);
                ind_p = find(Particles.FrameNum == frameNumbers(frame_i));
                ind_p(ind_p == ovlp_row(i,5)) =  [];
                p_id = A_row;
                x1d_id = A{p_id}(end);
                x1d = double(cell2mat(Particles.GrayImage(x1d_id)));
                gcx = guessC{1,frame_i-1}(Particles.Index(x1d_id));
                gcy = guessC{2,frame_i-1}(Particles.Index(x1d_id));
                size_x1d = size(x1d);
                if max(size_x1d) <= 48
                    [x2d,x2d_bb] = get_x2d(im2,gcx,gcy,size_x1d,16);
                else
                    [x2d,x2d_bb] = get_x2d(im2,gcx,gcy,size_x1d,32);
                end
                [maxy,maxx,corr,R_fft] = best_corr(x1d,x2d);
                xp = floor(maxx + x2d_bb(1) + size_x1d(2)/2);
                yp = floor(maxy + x2d_bb(2) + size_x1d(1)/2);
                YN_tracked = 0;
                p_tracked = [];
                r_tracked = [];
                if ~isnan(corr)
                    for p2_id = ind_p(1):ind_p(end)
                        if p2_id ~= ovlp_row(i,5)
                            [YN,rx,ry] = in_bbox(xp,yp,Particles.BoundingBox(p2_id,:));
                            if YN == 1
                                YN_tracked = YN_tracked + 1;
                                p_tracked = [p_tracked p2_id];
                                r_tracked = [r_tracked rx+ry];
                            end
                        end
                    end
                end
                if YN_tracked == 1
                    A{p_id} = [A{p_id} p_tracked(1)];
                elseif YN_tracked > 1
                    A{p_id} = [A{p_id} p_tracked(r_tracked == min(r_tracked))];
                elseif YN_tracked == 0
                    YN_mask = in_mask(gcx,gcy,static_mask);
                    if YN_mask ~= 1
                        A{p_id} = [A{p_id} -1];
                    else
                        bbox = Particles.BoundingBox(A{p_id}(end),:);
                        bbox_x1 = bbox(1,1) + 0.5;
                        bbox_y1 = bbox(1,2) + 0.5;
                        bbox_x2 = bbox(1,1) + bbox(1,3) - 0.5;
                        bbox_y2 = bbox(1,2) + bbox(1,4) - 0.5;
                        YN_bbox_mask = in_mask(bbox_x1,bbox_y1,static_mask) + ...
                            in_mask(bbox_x1,bbox_y2,static_mask) + ...
                            in_mask(bbox_x2,bbox_y1,static_mask) + ...
                            in_mask(bbox_x2,bbox_y2,static_mask);
                        if YN_bbox_mask < 4
                            A{p_id} = [A{p_id} -1];
                        elseif corr < 0.5
                            A{p_id} = [A{p_id} -2];
                        elseif Particles.Area(A{p_id}(end)) <= 100
                            A{p_id} = [A{p_id} -3];
                        elseif Particles.Area(A{p_id}(end)) > 2000
                            A{p_id} = [A{p_id} -5];
                        else
                            pxl_g = double(cell2mat(Particles.GrayImage(A{p_id}(end))));
                            pxl = sort(pxl_g(:),'descend');
                            mean_g_int = mean(pxl(51:Particles.Area(A{p_id}(end))));
                            if mean_g_int/255 <= T_g && corr >= 0.7
                                A{p_id} = [A{p_id} -4];
                            else
                                A{p_id} = [A{p_id} -6];
                            end
                        end
                    end
                end
                % re-write that row and remove old trackings
                if A{p_id}(end) > 0

                    % check if it is getting frames from a new line
                    m = 0;
                    for j = 1:length(A)
                        if A{p_id}(end) == A{j,1}(1)
                            A{p_id} = [A{p_id}(1:end-1) A{j}(1:end)];
                            remove_row = [remove_row j];
                            m = m + 1;
                        end
                    end
                    if m == 0
                        for j = 1:length(A)
                            if ismember(A{p_id,1}(end),A{j})
                                ind = find(A{j,1} == A{p_id,1}(end));
                                A{p_id} = [A{p_id}(1:end-1) A{j}(ind:end)];
                            end
                        end
                    end
                end
            end

        case 3
            % do nothing

        otherwise
            % should not have these cases
            % if this case occurs, we stop the tracking for both particles
            A_row = ovlp_row(i,1);
            ind_clean = find(A{A_row} == ovlp_row(i,3));
            if ~isempty(ind_clean)
                A{A_row} = A{A_row}(1:ind_clean);
                % run the particle tracking again
                % but exclude the old wrong match
                frame_i = find(frameNumbers == Particles.FrameNum(ovlp_row(i,5)));
                im2 = imread(fullfile(test_dir,'bgrem',['bgrem_' num2str(frameNumbers(frame_i),'%05d') '.png']));
                bbox(1,1:4) = Particles.BoundingBox(ovlp_row(i,5),:);
                bbox_x1 = bbox(1,1) + 0.5;
                bbox_y1 = bbox(1,2) + 0.5;
                bbox_x2 = bbox(1,1) + bbox(1,3) - 0.5;
                bbox_y2 = bbox(1,2) + bbox(1,4) - 0.5;
                mask = cell2mat(Particles.Image(ovlp_row(i,5)));
                mask = ~mask;
                im2(bbox_y1:bbox_y2,bbox_x1:bbox_x2) = uint8(double(im2(bbox_y1:bbox_y2,bbox_x1:bbox_x2)).*mask);
                ind_p = find(Particles.FrameNum == frameNumbers(frame_i));
                ind_p(ind_p == ovlp_row(i,5)) =  [];
                p_id = A_row;
                x1d_id = A{p_id}(end);
                x1d = double(cell2mat(Particles.GrayImage(x1d_id)));
                gcx = guessC{1,frame_i-1}(Particles.Index(x1d_id));
                gcy = guessC{2,frame_i-1}(Particles.Index(x1d_id));
                size_x1d = size(x1d);
                if max(size_x1d) <= 48
                    [x2d,x2d_bb] = get_x2d(im2,gcx,gcy,size_x1d,16);
                else
                    [x2d,x2d_bb] = get_x2d(im2,gcx,gcy,size_x1d,32);
                end
                [maxy,maxx,corr,R_fft] = best_corr(x1d,x2d);
                xp = floor(maxx + x2d_bb(1) + size_x1d(2)/2);
                yp = floor(maxy + x2d_bb(2) + size_x1d(1)/2);
                YN_tracked = 0;
                p_tracked = [];
                r_tracked = [];
                if ~isnan(corr)
                    for p2_id = ind_p(1):ind_p(end)
                        [YN,rx,ry] = in_bbox(xp,yp,Particles.BoundingBox(p2_id,:));
                        if YN == 1
                            YN_tracked = YN_tracked + 1;
                            p_tracked = [p_tracked p2_id];
                            r_tracked = [r_tracked rx+ry];
                        end
                    end
                end
                if YN_tracked == 1
                    A{p_id} = [A{p_id} p_tracked(1)];
                elseif YN_tracked > 1
                    A{p_id} = [A{p_id} p_tracked(r_tracked == min(r_tracked))];
                elseif YN_tracked == 0
                    YN_mask = in_mask(gcx,gcy,static_mask);
                    if YN_mask ~= 1
                        A{p_id} = [A{p_id} -1];
                    else
                        bbox = Particles.BoundingBox(A{p_id}(end),:);
                        bbox_x1 = bbox(1,1) + 0.5;
                        bbox_y1 = bbox(1,2) + 0.5;
                        bbox_x2 = bbox(1,1) + bbox(1,3) - 0.5;
                        bbox_y2 = bbox(1,2) + bbox(1,4) - 0.5;
                        YN_bbox_mask = in_mask(bbox_x1,bbox_y1,static_mask) + ...
                            in_mask(bbox_x1,bbox_y2,static_mask) + ...
                            in_mask(bbox_x2,bbox_y1,static_mask) + ...
                            in_mask(bbox_x2,bbox_y2,static_mask);
                        if YN_bbox_mask < 4
                            A{p_id} = [A{p_id} -1];
                        elseif corr < 0.5
                            A{p_id} = [A{p_id} -2];
                        elseif Particles.Area(A{p_id}(end)) <= 100
                            A{p_id} = [A{p_id} -3];
                        elseif Particles.Area(A{p_id}(end)) > 2000
                            A{p_id} = [A{p_id} -5];
                        else
                            pxl_g = double(cell2mat(Particles.GrayImage(A{p_id}(end))));
                            pxl = sort(pxl_g(:),'descend');
                            mean_g_int = mean(pxl(51:Particles.Area(A{p_id}(end))));
                            if mean_g_int/255 <= T_g && corr >= 0.7
                                A{p_id} = [A{p_id} -4];
                            else
                                A{p_id} = [A{p_id} -6];
                            end
                        end
                    end
                end
                % re-write that row and remove old trackings
                if A{p_id}(end) > 0
                    m = 0;
                    for j = 1:length(A)
                        if A{p_id}(end) == A{j,1}(1)
                            A{p_id} = [A{p_id}(1:end-1) A{j}(1:end)];
                            remove_row = [remove_row j];
                            m = m+1;
                        end
                    end
                    if m == 0
                        for j = 1:length(A)
                            if ismember(A{p_id,1}(end),A{j})
                                ind = find(A{j,1} == A{p_id,1}(end));
                                A{p_id} = [A{p_id}(1:end-1) A{j}(ind:end)];
                            end
                        end
                    end
                end
            end

            A_row = ovlp_row(i,2);
            ind_clean = find(A{A_row} == ovlp_row(i,4));
            if ~isempty(ind_clean)
                A{A_row} = A{A_row}(1:ind_clean);
                % run the particle tracking again
                % but exclude the old wrong match
                frame_i = find(frameNumbers == Particles.FrameNum(ovlp_row(i,5)));
                im2 = imread(fullfile(test_dir,'bgrem',['bgrem_' num2str(frameNumbers(frame_i),'%05d') '.png']));
                bbox(1,1:4) = Particles.BoundingBox(ovlp_row(i,5),:);
                bbox_x1 = bbox(1,1) + 0.5;
                bbox_y1 = bbox(1,2) + 0.5;
                bbox_x2 = bbox(1,1) + bbox(1,3) - 0.5;
                bbox_y2 = bbox(1,2) + bbox(1,4) - 0.5;
                mask = cell2mat(Particles.Image(ovlp_row(i,5)));
                mask = ~mask;
                im2(bbox_y1:bbox_y2,bbox_x1:bbox_x2) = uint8(double(im2(bbox_y1:bbox_y2,bbox_x1:bbox_x2)).*mask);
                ind_p = find(Particles.FrameNum == frameNumbers(frame_i));
                ind_p(ind_p == ovlp_row(i,5)) =  [];
                p_id = A_row;
                x1d_id = A{p_id}(end);
                x1d = double(cell2mat(Particles.GrayImage(x1d_id)));
                gcx = guessC{1,frame_i-1}(Particles.Index(x1d_id));
                gcy = guessC{2,frame_i-1}(Particles.Index(x1d_id));
                size_x1d = size(x1d);
                if max(size_x1d) <= 48
                    [x2d,x2d_bb] = get_x2d(im2,gcx,gcy,size_x1d,16);
                else
                    [x2d,x2d_bb] = get_x2d(im2,gcx,gcy,size_x1d,32);
                end
                [maxy,maxx,corr,R_fft] = best_corr(x1d,x2d);
                xp = floor(maxx + x2d_bb(1) + size_x1d(2)/2);
                yp = floor(maxy + x2d_bb(2) + size_x1d(1)/2);
                YN_tracked = 0;
                p_tracked = [];
                r_tracked = [];
                if ~isnan(corr)
                    for p2_id = ind_p(1):ind_p(end)
                        [YN,rx,ry] = in_bbox(xp,yp,Particles.BoundingBox(p2_id,:));
                        if YN == 1
                            YN_tracked = YN_tracked + 1;
                            p_tracked = [p_tracked p2_id];
                            r_tracked = [r_tracked rx+ry];
                        end
                    end
                end
                if YN_tracked == 1
                    A{p_id} = [A{p_id} p_tracked(1)];
                elseif YN_tracked > 1
                    A{p_id} = [A{p_id} p_tracked(r_tracked == min(r_tracked))];
                elseif YN_tracked == 0
                    YN_mask = in_mask(gcx,gcy,static_mask);
                    if YN_mask ~= 1
                        A{p_id} = [A{p_id} -1];
                    else
                        bbox = Particles.BoundingBox(A{p_id}(end),:);
                        bbox_x1 = bbox(1,1) + 0.5;
                        bbox_y1 = bbox(1,2) + 0.5;
                        bbox_x2 = bbox(1,1) + bbox(1,3) - 0.5;
                        bbox_y2 = bbox(1,2) + bbox(1,4) - 0.5;
                        YN_bbox_mask = in_mask(bbox_x1,bbox_y1,static_mask) + ...
                            in_mask(bbox_x1,bbox_y2,static_mask) + ...
                            in_mask(bbox_x2,bbox_y1,static_mask) + ...
                            in_mask(bbox_x2,bbox_y2,static_mask);
                        if YN_bbox_mask < 4
                            A{p_id} = [A{p_id} -1];
                        elseif corr < 0.5
                            A{p_id} = [A{p_id} -2];
                        elseif Particles.Area(A{p_id}(end)) <= 100
                            A{p_id} = [A{p_id} -3];
                        elseif Particles.Area(A{p_id}(end)) > 2000
                            A{p_id} = [A{p_id} -5];
                        else
                            pxl_g = double(cell2mat(Particles.GrayImage(A{p_id}(end))));
                            pxl = sort(pxl_g(:),'descend');
                            mean_g_int = mean(pxl(51:Particles.Area(A{p_id}(end))));
                            if mean_g_int/255 <= T_g && corr >= 0.7
                                A{p_id} = [A{p_id} -4];
                            else
                                A{p_id} = [A{p_id} -6];
                            end
                        end
                    end
                end
                % re-write that row and remove old trackings
                if A{p_id}(end) > 0

                    m = 0;
                    for j = 1:length(A)
                        if A{p_id}(end) == A{j,1}(1)
                            A{p_id} = [A{p_id}(1:end-1) A{j}(1:end)];
                            m = m+1;
                            remove_row = [remove_row j];
                        end
                    end
                    if m == 0
                        for j = 1:length(A)
                            if ismember(A{p_id,1}(end),A{j})
                                ind = find(A{j,1} == A{p_id,1}(end));
                                A{p_id} = [A{p_id}(1:end-1) A{j}(ind:end)];
                            end
                        end
                    end
                end
                
            end

    end

end
A(remove_row) = [];


%% fragmentation sanity check
velocityCell = cell(2,frames-1);
% this guess should be one frame ago
guessC_2 = cell(2,frames-1);
for i = 1:frames-1

    % particle index in the entire table "Particles"
    % find the particles in the "next" frame
    ind_p = find(Particles.FrameNum == (frameNumbers(i+1)));
    num_p = length(ind_p);
    est_u = zeros(num_p,1);
    est_v = zeros(num_p,1);
    meshX = x{PIV_ind == frameNumbers(i),1};
    meshY = y{PIV_ind == frameNumbers(i),1};
    meshU = u_filtered{PIV_ind == frameNumbers(i),1};
    meshV = v_filtered{PIV_ind == frameNumbers(i),1};

    for j = ind_p(1):ind_p(end)
        x_p = Particles.Centroid(j,1);
        y_p = Particles.Centroid(j,2);

        % region 1
        if x_p > min(meshX(:)) && x_p < max(meshX(:)) && ...
                y_p > min(meshY(:)) && y_p < max(meshY(:))
            est_u(j-ind_p(1)+1,1) = interp2(meshX,meshY,meshU,x_p,y_p,'spline');
            est_v(j-ind_p(1)+1,1) = interp2(meshX,meshY,meshV,x_p,y_p,'spline');

            % region 2
        elseif x_p <= min(meshX(:)) && y_p > min(meshY(:)) && y_p < max(meshY(:))
            est_u(j-ind_p(1)+1,1) = interp1(meshY(:,1),meshU(:,1),y_p,'linear');
            est_v(j-ind_p(1)+1,1) = interp1(meshY(:,1),meshV(:,1),y_p,'linear');
        elseif x_p >= max(meshX(:)) && y_p > min(meshY(:)) && y_p < max(meshY(:))
            est_u(j-ind_p(1)+1,1) = interp1(meshY(:,end),meshU(:,end),y_p,'linear');
            est_v(j-ind_p(1)+1,1) = interp1(meshY(:,end),meshV(:,end),y_p,'linear');

            % region 3
        elseif y_p <= min(meshY(:)) && x_p > min(meshX(:)) && x_p < max(meshX(:))
            est_u(j-ind_p(1)+1,1) = interp1(meshX(1,:),meshU(1,:),x_p,'linear');
            est_v(j-ind_p(1)+1,1) = interp1(meshX(1,:),meshV(1,:),x_p,'linear');
        elseif y_p >= max(meshY(:)) && x_p > min(meshX(:)) && x_p < max(meshX(:))
            est_u(j-ind_p(1)+1,1) = interp1(meshX(1,:),meshU(1,:),x_p,'linear');
            est_v(j-ind_p(1)+1,1) = interp1(meshX(1,:),meshV(1,:),x_p,'linear');

            % region 4
        elseif x_p <= min(meshX(:)) && y_p <= min(meshY(:))
            est_u(j-ind_p(1)+1,1) = meshU(1,1);
            est_v(j-ind_p(1)+1,1) = meshV(1,1);
        elseif x_p >= max(meshX(:)) && y_p <= min(meshY(:))
            est_u(j-ind_p(1)+1,1) = meshU(end,1);
            est_v(j-ind_p(1)+1,1) = meshV(end,1);
        elseif x_p <= min(meshX(:)) && y_p >= max(meshY(:))
            est_u(j-ind_p(1)+1,1) = meshU(1,end);
            est_v(j-ind_p(1)+1,1) = meshV(1,end);
        elseif x_p >= max(meshX(:)) && y_p >= max(meshY(:))
            est_u(j-ind_p(1)+1,1) = meshU(end,end);
            est_v(j-ind_p(1)+1,1) = meshV(end,end);
        end
    end

    velocityCell{1,i} = est_u;
    velocityCell{2,i} = est_v;
    % guessed centroids in the previous frame
    guessC_2{1,i} = Particles.Centroid(ind_p,1) - velocityCell{1,i};
    guessC_2{2,i} = Particles.Centroid(ind_p,2) - velocityCell{2,i};
end

% initialize a matrix called B for fragmentation
% flag = 0 means the 1st frame
B = zeros(length(A),2);
A_end = zeros(length(A),1);
% obtain the last element from A
for row_A = 1:length(A)
    if A{row_A,1}(end) < 0
        A_end(row_A,1) = A{row_A,1}(end-1);
    else
        A_end(row_A,1) = A{row_A,1}(end);
    end
end

%
for row_A = 1:length(A)

    % p_id for 2nd particle
    p2_id = A{row_A,1}(1);
    B(row_A,1) = p2_id;
    frame_i2 = find(frameNumbers == Particles.FrameNum(p2_id));

    % exclude the frame #1
    if frame_i2 ~= 1

        % interrogation window for 2nd particle
        x2d = double(cell2mat(Particles.GrayImage(p2_id)));

        % guessed centroid in the previous frame
        gcx = guessC_2{1,frame_i2-1}(Particles.Index(p2_id));
        gcy = guessC_2{2,frame_i2-1}(Particles.Index(p2_id));

        % check if the guessed center is in the mask or not
        YN_mask = in_mask(gcx,gcy,static_mask);

        % im1 is always the previous frame
        im1 = imread(fullfile(test_dir,'bgrem',['bgrem_' num2str(frameNumbers(frame_i2-1),'%05d') '.png']));

        size_x2d = size(x2d);
        % get x1d
        if max(size_x2d) <= 48 % bounding box < 32
            [x1d,x1d_bb] = get_x2d(im1,gcx,gcy,size_x2d,16);
        else
            [x1d,x1d_bb] = get_x2d(im1,gcx,gcy,size_x2d,32);
        end

        % fft & ifft method to get the best correlation
        [maxy,maxx,corr,R_fft] = best_corr(x2d,x1d);

        % get the tracked center
        xp = floor(maxx + x1d_bb(1) + size_x2d(2)/2);
        yp = floor(maxy + x1d_bb(2) + size_x2d(1)/2);

        % obtain indices of all particles in the previous frame
        ind_p = find(Particles.FrameNum == frameNumbers(frame_i2-1));

        % initialize the number of successful trackings
        YN_tracked = 0;
        p_tracked = [];
        r_tracked = [];
        if ~isnan(corr)
            for p1_id = ind_p(1):ind_p(end)
                % check if tracked center is in any bbox
                [YN,rx,ry] = in_bbox(xp,yp,Particles.BoundingBox(p1_id,:));

                if YN == 1
                    YN_tracked = YN_tracked + 1;
                    p_tracked = [p_tracked p1_id];
                    r_tracked = [r_tracked rx+ry];

                end
            end
        end

        % only 1 particle tracked
        if YN_tracked == 1 && corr > 0.7
            B(row_A,2) = p_tracked(1);
        elseif YN_tracked > 1
            B(row_A,2) = p_tracked(r_tracked == min(r_tracked));
        end

    end

end

% this is combining two rows into one due to the mismatch for some reason
mis_match = 0;
% A2 is the tracking matrix after fragmentation check
A2 = A;
B2 = B;
remove_row2 = [];
% have to do it reservely
% otherwise we may lose some unique particles
% for example line 2001 -> line 1001 -> line 1
% using the ascending order will lose particles in line 2001 and we will
% only have particles in line 1 & 1001 because particles in line 2001 is
% going to line 1001. And line 1001 is deleted!
for i = length(B):-1:1
    if B(i,2) ~= 0
        for j = length(A_end):-1:1
            if B(i,2) == A_end(j)
                mis_match = mis_match + 1;
                if A2{j,1}(end) < 0
                    A2{j,1} = [A2{j,1}(1:end-1) A2{i,1}(1:end)];
                    % update A_end every time
                    if A2{i,1}(end) < 0
                        A_end(j) = A2{i,1}(end-1);
                    else
                        A_end(j) = A2{i,1}(end);
                    end
                else
                    A2{j,1} = [A2{j,1}(1:end) A2{i,1}(1:end)];
                    % update A_end every time
                    if A2{i,1}(end) < 0
                        A_end(j) = A2{i,1}(end-1);
                    else
                        A_end(j) = A2{i,1}(end);
                    end
                end
                remove_row2 = [remove_row2 i];
                B2(i,2) = 0;
            end
        end 
    end
end


%% now we have to rearrange everything

% A2 is the new tracking matrix
A2(remove_row2) = [];

% B2 now has the match particle id
B2(remove_row2,:) = [];
% A3 = [A2{:};];
% A3 = A3(A3>0);

% % check whether all particles are detected or not
% A4 = unique(A3); 

% #1 #2 sub-aggregate 1 id & row
% #3 #4 sub-aggregate 2 id & row
% #5    aggregate 3 id
ovlp_in = zeros(sum(ovlp_row(:,6)==3),5);
m = 1;
for i = 1:length(ovlp_row)
    if ovlp_row(i,6) == 3
        ovlp_in(m,1) = ovlp_row(i,3);
        for j = 1:length(A2)
            if ismember(ovlp_row(i,3),A2{j,1})
                ovlp_in(m,2) = j;
            end
        end
        ovlp_in(m,3) = ovlp_row(i,4);
        for j = 1:length(A2)
            if ismember(ovlp_row(i,4),A2{j,1})
                ovlp_in(m,4) = j;
            end
        end
        ovlp_in(m,5) = ovlp_row(i,5);
        m = m + 1;
    end
end
% #1 #2 aggregate 3 id & row
% #3 #4 sub-aggregate 4 id & row
ovlp_out = zeros(sum(B2(:,2)>0),4);
m = 1;
for i = 1:length(B2)
    if B2(i,2) ~= 0
        ovlp_out(m,1) = B2(i,2);
        for j = 1:length(A2)
            if ismember(B2(i,2),A2{j,1})
                ovlp_out(m,2) = j;
            end
        end
        ovlp_out(m,3) = B2(i,1);
        for j = 1:length(A2)
            if ismember(B2(i,1),A2{j,1})
                ovlp_out(m,4) = j;
            end
        end
        m = m + 1;
    end
end

save('particle_id.mat','Particles')
save('tracking.mat','A2','ovlp_in','ovlp_out')


%% Part 2: see if we can get matched between two images over a couple of frames
% get the index for PIV ready
PIV_ind = ind_start:ind_end-1;
PIV_ind(ismember(PIV_ind,bad_frame)) = [];

% initialize the start & end particle id
% possible start means particles that potentially miss tracking due to
% various reasons.
% possible end means particles that cannot find their matches in the next
% frame.
possible_start = [];
possible_end = [];
start_count = 0;
end_count = 0;
for i = 1:length(A2)
    if Particles.FrameNum(A2{i,1}(1)) ~= Particles.FrameNum(1)
        possible_start = [possible_start A2{i,1}(1)];
    end
    if ismember(A2{i,1}(end),[-1 -2 -3 -4 -6])
        possible_end = [possible_end A2{i,1}(end-1)];
    end
end
possible_end = unique(possible_end);
possible_end(2,:) = 0;

for i = 1:length(possible_end)
    p1_id = possible_end(1,i);

    % initialize the number of successful trackings
    YN_tracked = 0;
    p_tracked = [];
    % corr_tracked = [];
    % dist_tracked = [];

    x_p = zeros(frames_max_interval+1,1);
    y_p = zeros(frames_max_interval+1,1);

    x_p(1,1) = Particles.Centroid(p1_id,1);
    y_p(1,1) = Particles.Centroid(p1_id,2);

    for k = 1:min(frames_max_interval,find(frameNumbers == ind_end)-find(frameNumbers == Particles.FrameNum(p1_id)))

        % obtain the frame for one step in PIV
        % avoid bad frame
        p_frame = frameNumbers(find(frameNumbers == Particles.FrameNum(p1_id))-1+k);

        % obtain the mesh grid for X Y U V
        meshX = x{PIV_ind == p_frame,1};
        meshY = y{PIV_ind == p_frame,1};
        meshU = u_filtered{PIV_ind == p_frame,1};
        meshV = v_filtered{PIV_ind == p_frame,1};

        % now do the interpolations based on PIV
        % region 1
        if x_p(k,1) > min(meshX(:)) && x_p(k,1) < max(meshX(:)) && ...
                y_p(k,1) > min(meshY(:)) && y_p(k,1) < max(meshY(:))
            est_u = interp2(meshX,meshY,meshU,x_p(k,1),y_p(k,1),'spline');
            est_v = interp2(meshX,meshY,meshV,x_p(k,1),y_p(k,1),'spline');
            % region 2
        elseif x_p(k,1) <= min(meshX(:)) && y_p(k,1) > min(meshY(:)) && y_p(k,1) < max(meshY(:))
            est_u = interp1(meshY(:,1),meshU(:,1),y_p(k,1),'linear');
            est_v = interp1(meshY(:,1),meshV(:,1),y_p(k,1),'linear');
        elseif x_p(k,1) >= max(meshX(:)) && y_p(k,1) > min(meshY(:)) && y_p(k,1) < max(meshY(:))
            est_u = interp1(meshY(:,end),meshU(:,end),y_p(k,1),'linear');
            est_v = interp1(meshY(:,end),meshV(:,end),y_p(k,1),'linear');
            % region 3
        elseif y_p(k,1) <= min(meshY(:)) && x_p(k,1) > min(meshX(:)) && x_p(k,1) < max(meshX(:))
            est_u = interp1(meshX(1,:),meshU(1,:),x_p(k,1),'linear');
            est_v = interp1(meshX(1,:),meshV(1,:),x_p(k,1),'linear');
        elseif y_p(k,1) >= max(meshY(:)) && x_p(k,1) > min(meshX(:)) && x_p(k,1) < max(meshX(:))
            est_u = interp1(meshX(1,:),meshU(1,:),x_p(k,1),'linear');
            est_v = interp1(meshX(1,:),meshV(1,:),x_p(k,1),'linear');
            % region 4
        elseif x_p(k,1) <= min(meshX(:)) && y_p(k,1) <= min(meshY(:))
            est_u = meshU(1,1);
            est_v = meshV(1,1);
        elseif x_p(k,1) >= max(meshX(:)) && y_p(k,1) <= min(meshY(:))
            est_u = meshU(end,1);
            est_v = meshV(end,1);
        elseif x_p(k,1) <= min(meshX(:)) && y_p(k,1) >= max(meshY(:))
            est_u = meshU(1,end);
            est_v = meshV(1,end);
        elseif x_p(k,1) >= max(meshX(:)) && y_p(k,1) >= max(meshY(:))
            est_u = meshU(end,end);
            est_v = meshV(end,end);
        end

        % update the estimated centroids based on PIV
        x_p(1+k,1) = x_p(k,1) + est_u;
        y_p(1+k,1) = y_p(k,1) + est_v;
    end

    for j = 1:length(possible_start)

        p2_id = possible_start(j);

        % we cannot search for discontinuous tracking for > 16 hours
        % first condition uses "find" to avoid bad frames
        % second condition is based on time (16 hours)
        if find(frameNumbers == Particles.FrameNum(p2_id)) - find(frameNumbers == Particles.FrameNum(p1_id)) > 0 && ...
                Particles.FrameNum(p2_id) - Particles.FrameNum(p1_id) < frames_max_interval

            % calculate the steps            
            % the step number determines the interrogation window & PIV use
            % use find function to avoid bad frames
            steps = find(frameNumbers == Particles.FrameNum(p2_id)) - find(frameNumbers == Particles.FrameNum(p1_id));

            x_p_guess = x_p(1+steps,1);
            y_p_guess = y_p(1+steps,1);

            % only necessary to calculate p2 if conditions are met
            x_p2 = Particles.Centroid(p2_id,1);
            y_p2 = Particles.Centroid(p2_id,2);

            % calulate the distance between two points
            % dist0 = sqrt((x_p2-x_p(1))^2 + (y_p2-y_p(1))^2);
            dist = sqrt((x_p2-x_p_guess)^2 + (y_p2-y_p_guess)^2);

            % determine if the distance is within the range
            % if the area is similar we get it approved
            if dist <= 8 &&...
                    Particles.Area(p1_id) < 1.5*Particles.Area(p2_id) &&...
                    Particles.Area(p2_id) < 1.5*Particles.Area(p1_id)

                YN_tracked = YN_tracked + 1;
                % dist_tracked = [dist_tracked dist];
                p_tracked = [p_tracked p2_id];
                % corr_tracked = [corr_tracked 1];

                % according to the histogram, the distance between estimations
                % and real match is <6.82 (95%). 32 pixels is safe.
            elseif dist <= 16

                % get 1st interrogation window
                x1d = double(cell2mat(Particles.GrayImage(p1_id)));
                size_x1d = size(x1d);

                im2_x2d = double(cell2mat(Particles.GrayImage(p2_id)));
                x2d = zeros(size(im2_x2d)+16);
                x2d(1:size(im2_x2d,1),1:size(im2_x2d,2)) = im2_x2d;

                % % load 2nd image
                % im2_id = frameNumbers(frameNumbers == Particles.FrameNum(p2_id));
                % % get the next bgrem gray image before the particle for loop
                % im2 = imread(fullfile(bgrem_dir,'bgrem',['bgrem_' num2str(im2_id,'%05d') '.png']));
                % % get 2nd interrogation window
                % [x2d,~] = get_x2d(im2,x_p2,y_p2,size_x1d,8);

                % compare the bbox and area
                if size_x1d(1) <= size(x2d,1) && size_x1d(2) <= size(x2d,2) &&...
                        Particles.Area(p1_id) < 2*Particles.Area(p2_id) &&...
                        Particles.Area(p2_id) < 2*Particles.Area(p1_id)
                    % cross correlation
                    [~,~,corr,~] = best_corr(x1d,x2d);

                    if ~isnan(corr)
                        if corr >= 0.5
                            YN_tracked = YN_tracked + 1;
                            % dist_tracked = [dist_tracked dist];
                            p_tracked = [p_tracked p2_id];
                            % corr_tracked = [corr_tracked corr];
                        end
                    end
                end
            end
            clearvars corr
        end

    end

    if YN_tracked == 1
        possible_end(2,i) = p_tracked;
    elseif YN_tracked > 1
        possible_end(2,i) = min(p_tracked);
        % possible_end(2,i) = p_tracked(dist_tracked == min(dist_tracked));
        % possible_end(2,i) = p_tracked(corr_tracked == max(corr_tracked));
    end


    disp(i)
end


num_match = sum(possible_end(2,:)>0);
matches = zeros(num_match,2);
m = 1;
for i = 1:length(possible_end)
    if possible_end(2,i) > 0
        matches(m,1) = possible_end(1,i);
        matches(m,2) = possible_end(2,i);
        m = m + 1;
    end
end

A3 = A2;
for i = 1:length(A2)
    A2_start(i) = A2{i,1}(1); 
    if A2{i,1}(end) < 0
        A2_end(i) = A2{i,1}(end-1);
    else
        A2_end(i) = A2{i,1}(end);
    end
end
m = 0;
for i = length(possible_end):-1:1
    if possible_end(2,i) > 0
        if length(find(possible_end(2,:) == possible_end(2,i))) == 1
            A3_ind1 = find(A2_start == possible_end(2,i));
            A3_ind2 = find(A2_end == possible_end(1,i));
            for j = 1:length(A3_ind2)
                A3{A3_ind2(j),1} = [A3{A3_ind2(j),1}(1:end-1) A3{A3_ind1}(1:end)];
            end
            A3{A3_ind1,1} = [];
            m = m +1;
        else
            A3_ind1 = find(A2_start == possible_end(2,i));
            A3_inds2 = find(possible_end(2,:) == possible_end(2,i));
            A3_ind2 = find(A2_end == max(possible_end(1,A3_inds2)),1);
            
            possible_end(2,A3_inds2(A3_inds2~=A3_ind2)) = 0;
            for j = 1:length(A3_ind2)
                A3{A3_ind2(j),1} = [A3{A3_ind2(j),1}(1:end-1) A3{A3_ind1}(1:end)];
            end
            A3{A3_ind1,1} = [];
        end
    end
end
A4 = A3;
isEmptyCell = cellfun(@isempty, A3);
A4(isEmptyCell) = [];


%% Part 3: add one more feature to see if any particle comes from the unmasked area
A4_start = zeros(length(A4),2);
for i = 1:length(A4)
    A4_start(i,1) = A4{i,1}(1);
end

frameNumbers = ind_start:ind_end;
frameNumbers(ismember(frameNumbers,bad_frame)) = [];
frames = length(frameNumbers);
PIV_ind = ind_start:ind_end-1;
PIV_ind(ismember(PIV_ind,bad_frame)) = [];

for i = 1:length(A4)
    p_id = A4_start(i,1);
    if Particles.FrameNum(p_id) > Particles.FrameNum(1)

        x_p = Particles.Centroid(p_id,1);
        y_p = Particles.Centroid(p_id,2);

        % avoid bad frame
        p_frame = frameNumbers(find(frameNumbers == Particles.FrameNum(p_id))-1);

        % obtain the mesh grid for X Y U V
        meshX = x{PIV_ind == p_frame,1};
        meshY = y{PIV_ind == p_frame,1};
        meshU = u_filtered{PIV_ind == p_frame,1};
        meshV = v_filtered{PIV_ind == p_frame,1};

        if x_p > min(meshX(:)) && x_p < max(meshX(:)) && ...
                y_p > min(meshY(:)) && y_p < max(meshY(:))
            est_u = interp2(meshX,meshY,meshU,x_p,y_p,'spline');
            est_v = interp2(meshX,meshY,meshV,x_p,y_p,'spline');
            % region 2
        elseif x_p <= min(meshX(:)) && y_p > min(meshY(:)) && y_p < max(meshY(:))
            est_u = interp1(meshY(:,1),meshU(:,1),y_p,'linear');
            est_v = interp1(meshY(:,1),meshV(:,1),y_p,'linear');
        elseif x_p >= max(meshX(:)) && y_p > min(meshY(:)) && y_p < max(meshY(:))
            est_u = interp1(meshY(:,end),meshU(:,end),y_p,'linear');
            est_v = interp1(meshY(:,end),meshV(:,end),y_p,'linear');
            % region 3
        elseif y_p <= min(meshY(:)) && x_p > min(meshX(:)) && x_p < max(meshX(:))
            est_u = interp1(meshX(1,:),meshU(1,:),x_p,'linear');
            est_v = interp1(meshX(1,:),meshV(1,:),x_p,'linear');
        elseif y_p >= max(meshY(:)) && x_p > min(meshX(:)) && x_p < max(meshX(:))
            est_u = interp1(meshX(1,:),meshU(1,:),x_p,'linear');
            est_v = interp1(meshX(1,:),meshV(1,:),x_p,'linear');
            % region 4
        elseif x_p <= min(meshX(:)) && y_p <= min(meshY(:))
            est_u = meshU(1,1);
            est_v = meshV(1,1);
        elseif x_p >= max(meshX(:)) && y_p <= min(meshY(:))
            est_u = meshU(end,1);
            est_v = meshV(end,1);
        elseif x_p <= min(meshX(:)) && y_p >= max(meshY(:))
            est_u = meshU(1,end);
            est_v = meshV(1,end);
        elseif x_p >= max(meshX(:)) && y_p >= max(meshY(:))
            est_u = meshU(end,end);
            est_v = meshV(end,end);
        end

        % update the estimated centroids based on PIV
        x_p_0 = x_p - est_u;
        y_p_0 = y_p - est_v;

        bbox = Particles.BoundingBox(p_id,1:4);
        bbox_x1 = x_p_0 + 0.5;
        bbox_y1 = y_p_0 + 0.5;
        bbox_x2 = x_p_0 + bbox(1,3) - 0.5;
        bbox_y2 = y_p_0 + bbox(1,4) - 0.5;
        YN_bbox_mask = in_mask(bbox_x1,bbox_y1,static_mask) + ...
            in_mask(bbox_x1,bbox_y2,static_mask) + ...
            in_mask(bbox_x2,bbox_y1,static_mask) + ...
            in_mask(bbox_x2,bbox_y2,static_mask);

        if YN_bbox_mask < 4
            % a new flag of -7 indicates this particle may come from the
            % masked area
            A4_start(i,2) = -7;
        end

    end
    disp(i)
end

save('tracking_v2.mat','A4','A3','A2','A4_start');