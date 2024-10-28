function [maxy,maxx,corr,R_fft] = best_corr(x1d,x2d)

% by default settings, x1d is smaller than x2d
size_x1d = size(x1d);
% define a size for peak search 
dist = sqrt((size_x1d(1)-16)^2 + (size_x1d(2)-16)^2)/2;
size_x2d = size(x2d);
% calculate the excess size
d_x12d = (size_x2d - size_x1d + 1)/2;

% zero padding the x1d
x1d_pad = zeros(size(x2d));
x1d_pad(1:size_x1d(1),1:size_x1d(2)) = x1d;

R_fft = real(ifft2(conj(fft2(x1d_pad)).*fft2(x2d)));
R_fft = fftshift(fftshift(R_fft, 1), 2);
R_fft = R_fft/(std(x1d_pad(:))*std(x2d(:)))/numel(R_fft);
corr = max(abs(R_fft(:)));


if corr < 0.9

    % theoretically should be "== 1"
    % just make sure something can be output
    if numel(peaks2(R_fft,'MinPeakHeight',0.3,'MinPeakDistance',dist)) <= 1
        [maxy,maxx] = find(R_fft == max(R_fft(:)),1);
        % checked with PIV test: have to -1
        maxy = maxy - size_x2d(1)/2 - 1;
        maxx = maxx - size_x2d(2)/2 - 1;

    else
        [pks,locs_y,locs_x] = peaks2(R_fft,'MinPeakHeight',0.3,'MinPeakDistance',dist);
        pixels_from_center = zeros(length(pks),1);
        for i = 1:length(pks)
            pxl_y = abs(locs_y(i) - size_x2d(2)/2 - 1 - d_x12d(2));
            pxl_x = abs(locs_x(i) - size_x2d(1)/2 - 1 - d_x12d(1));
            pixels_from_center(i) = pxl_x + pxl_y;
        end
        pks_best = pixels_from_center == min(pixels_from_center);
        if sum(pks_best(:)) > 1
            corr = max(pks);
            [maxy,maxx] = find(R_fft == max(R_fft(:)),1);
            maxy = maxy - size_x2d(1)/2 - 1;
            maxx = maxx - size_x2d(2)/2 - 1;

        else
            corr = pks(pks_best);
            maxy = locs_y(pks_best) - size_x2d(1)/2 - 1;
            maxx = locs_x(pks_best) - size_x2d(2)/2 - 1;
        end

    end
    
else
    % corr >= 0.9
    % this is a step that just saves computational time
    [maxy,maxx] = find(R_fft == max(R_fft(:)),1);
    maxy = maxy - size_x2d(1)/2 - 1;
    maxx = maxx - size_x2d(2)/2 - 1;
end

end