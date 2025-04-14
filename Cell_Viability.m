function cell_Viability()
    try
        % Directly specify file paths
        live_path = 'live_path.jpg';  % Replace with your actual path
        dead_path = 'dead_path.jpg';  % Replace with your actual path
        
        % Read images and check dimensions
        live_img = imread(live_path);
        dead_img = imread(dead_path);
        
        % Validate image dimensions
        if ndims(live_img) < 2 || ndims(dead_img) < 2
            error('Images must have at least 2 dimensions');
        end
        
        % Convert to grayscale if RGB
        if size(live_img, 3) > 1
            live_img_gray = rgb2gray(live_img);
            live_g_norm = mat2gray(live_img(:,:,2));  % Normalized green channel
        else
            live_img_gray = live_img;
            live_g_norm = mat2gray(live_img_gray);
        end
        
        if size(dead_img, 3) > 1
            dead_img_gray = rgb2gray(dead_img);
            dead_r_norm = mat2gray(dead_img(:,:,1));  % Normalized red channel
        else
            dead_img_gray = dead_img;
            dead_r_norm = mat2gray(dead_img_gray);
        end
        
        % Ensure minimum image size for registration
        min_image_size = 50;
        if size(live_img_gray, 1) < min_image_size || size(live_img_gray, 2) < min_image_size || ...
           size(dead_img_gray, 1) < min_image_size || size(dead_img_gray, 2) < min_image_size
            error('Images are too small for analysis. Minimum size is %d pixels in each dimension.', min_image_size);
        end
        
        % Alternative drift correction if imregister fails
        try
            % Use normalized cross-correlation for image alignment
            [optimizer, metric] = imregconfig('multimodal');
            registered_dead_img_gray = imregister(dead_img_gray, live_img_gray, 'rigid', optimizer, metric);
            
            % Also register the normalized color channels
            registered_dead_r_norm = imregister(dead_r_norm, live_g_norm, 'rigid', optimizer, metric);
        catch ME
            % Fallback to manual alignment if registration fails
            warning('Automatic image registration failed. Using manual alignment method.');
            
            % Simple manual alignment using cross-correlation
            c = normxcorr2(dead_img_gray, live_img_gray);
            [max_val, imax] = max(abs(c(:)));
            [ypeak, xpeak] = ind2sub(size(c), imax);
            
            % Calculate offset
            offset_y = ypeak - size(live_img_gray, 1);
            offset_x = xpeak - size(live_img_gray, 2);
            
            % Simple translation-based alignment
            registered_dead_img_gray = imtranslate(dead_img_gray, [offset_x, offset_y]);
            registered_dead_r_norm = imtranslate(dead_r_norm, [offset_x, offset_y]);
        end
        
        % Configurable thresholds
        min_cell_area = 100; % Minimum area for cell detection
        min_box_size = 10;  % Minimum bounding box size (width or height)
        
        % Use adaptive thresholding with normalized image
        live_bw = imbinarize(live_img_gray, 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.4);
        dead_bw = imbinarize(registered_dead_img_gray, 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.4);
        
        % Clean up binary images (remove small objects)
        live_bw = bwareaopen(live_bw, min_cell_area);
        dead_bw = bwareaopen(dead_bw, min_cell_area);
        
        % Find connected components (cells)
        live_stats = regionprops(live_bw, 'Area', 'BoundingBox', 'Centroid');
        dead_stats = regionprops(dead_bw, 'Area', 'BoundingBox', 'Centroid');
        
        % Filter out small bounding boxes
        live_box_sizes = zeros(length(live_stats), 2);
        for i = 1:length(live_stats)
            live_box_sizes(i, :) = live_stats(i).BoundingBox(3:4);
        end
        
        dead_box_sizes = zeros(length(dead_stats), 2);
        for i = 1:length(dead_stats)
            dead_box_sizes(i, :) = dead_stats(i).BoundingBox(3:4);
        end
        
        % Create logical index for filtering
        live_valid_boxes = all(live_box_sizes >= min_box_size, 2);
        dead_valid_boxes = all(dead_box_sizes >= min_box_size, 2);
        
        % Filter stats based on valid boxes
        live_stats_filtered = live_stats(live_valid_boxes);
        dead_stats_filtered = dead_stats(dead_valid_boxes);
        
        % More robust overlap calculation
        overlap_bw = live_bw & dead_bw;
        overlap_stats = regionprops(overlap_bw, 'Area', 'BoundingBox', 'Centroid');
        
        % Filter overlap stats
        overlap_box_sizes = zeros(length(overlap_stats), 2);
        for i = 1:length(overlap_stats)
            overlap_box_sizes(i, :) = overlap_stats(i).BoundingBox(3:4);
        end
        overlap_valid_boxes = all(overlap_box_sizes >= min_box_size, 2);
        overlap_stats_filtered = overlap_stats(overlap_valid_boxes);
        
        % Calculate cell counts 
        live_cell_areas = [live_stats_filtered.Area];
        dead_cell_areas = [dead_stats_filtered.Area];
        overlap_cell_areas = [overlap_stats_filtered.Area];
        
        % More precise cell counting
        live_cells = sum(live_cell_areas >= min_cell_area & ~ismember(live_cell_areas, overlap_cell_areas));
        dead_cells = sum(dead_cell_areas >= min_cell_area);
        overlap_cells = sum(overlap_cell_areas >= min_cell_area);
        
        total_cells = live_cells + dead_cells;
        viability_percentage = (live_cells / total_cells) * 100;
        
        % Close any existing figures
        close all;
        
        % 1. Bounding Boxes on Original Color Images
        % Live Cells (Green Channel)
        figure('Name', 'Live Cells Bounding Boxes', 'NumberTitle', 'off');
        imshow(live_img);
        title('Live Cells with Bounding Boxes');
        hold on;
        for i = 1:numel(live_stats_filtered)
            rectangle('Position', live_stats_filtered(i).BoundingBox, 'EdgeColor', 'g', 'LineWidth', 2);
        end
        hold off;
        
        % Dead Cells (Red Channel)
        figure('Name', 'Dead Cells Bounding Boxes', 'NumberTitle', 'off');
        imshow(dead_img);
        title('Dead Cells with Bounding Boxes');
        hold on;
        for i = 1:numel(dead_stats_filtered)
            rectangle('Position', dead_stats_filtered(i).BoundingBox, 'EdgeColor', 'r', 'LineWidth', 2);
        end
        hold off;
        
        % Overlapping Cells
        figure('Name', 'Overlapping Cells', 'NumberTitle', 'off');
        imshow(dead_img);
        title('Overlapping Cells');
        hold on;
        for i = 1:numel(overlap_stats_filtered)
            rectangle('Position', overlap_stats_filtered(i).BoundingBox, 'EdgeColor', 'b', 'LineWidth', 2);
        end
        hold off;
        
        % 2. Histograms of Cell Areas
        figure('Name', 'Cell Area Distribution', 'NumberTitle', 'off');
        subplot(1,2,1);
        histogram(live_cell_areas, 'FaceColor', 'g', 'EdgeColor', 'k', 'BinWidth', 10);
        title('Live Cells Area Distribution');
        xlabel('Area (pixels)');
        ylabel('Frequency');
        
        subplot(1,2,2);
        histogram(dead_cell_areas, 'FaceColor', 'r', 'EdgeColor', 'k', 'BinWidth', 10);
        title('Dead Cells Area Distribution');
        xlabel('Area (pixels)');
        ylabel('Frequency');
        
        % 3. Binary Masks
        figure('Name', 'Binary Masks', 'NumberTitle', 'off');
        subplot(1,3,1);
        imshow(live_bw);
        title('Live Cells Mask');
        
        subplot(1,3,2);
        imshow(dead_bw);
        title('Dead Cells Mask');
        
        subplot(1,3,3);
        imshow(overlap_bw);
        title('Overlap Mask');
        
        % Display results
        fprintf('Minimum Cell Area Threshold: %d pixels\n', min_cell_area);
        fprintf('Minimum Bounding Box Size: %d pixels\n', min_box_size);
        fprintf('Total Live Cells: %d\n', live_cells);
        fprintf('Total Dead Cells: %d\n', dead_cells);
        fprintf('Overlapping Cells: %d\n', overlap_cells);
        fprintf('Cell Viability: %.2f%%\n', viability_percentage);
        
    catch ME
        % Comprehensive error handling
        fprintf('Error in cell viability analysis:\n');
        fprintf('Error Message: %s\n', ME.message);
        fprintf('Error Identifier: %s\n', ME.identifier);
        
        % Print stack trace for debugging
        for i = 1:length(ME.stack)
            fprintf('File: %s\n', ME.stack(i).file);
            fprintf('Name: %s\n', ME.stack(i).name);
            fprintf('Line: %d\n', ME.stack(i).line);
        end
    end
end