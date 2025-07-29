function [img]=filter_for_density(img,target_density)


method_none = 0;
method_setToZeroForDensity = 1;
method_tuneIntensity       = 2;
method_power               = 3;

density_method = method_none; 

if(density_method == method_power)
    % Parameters
    sigma = 2;             % Standard deviation for Gaussian smoothing
    threshold = 0.5;       % Initial threshold (adjust if necessary)
    power = 2;
    
    % Display the results
    figure;
    subplot(1, 3, 1)
    imshow(img)
    title('Original Image');

    is_desired_density = false;
    maxItersDensity = 100;
    idensity = 0;
    while(is_desired_density==false)

        img_smooth = imgaussfilt(img, sigma);
        threshold_val = mean(img);
        binary_img = img_smooth > threshold_val; % Convert to binary (black and white)
        current_density = sum(binary_img(:)) / numel(binary_img);
        
        fprintf(' %d -> density: %.2f%%\n',idensity, current_density * 100);
        imshow(img_smooth)

        is_desired_density= current_density<target_density;
    
        if(is_desired_density==false)
            img = img_smooth;
            imin = min(img(:));
            imax = max(img(:));
            img = (img-imin)/(imax-imin);
            
            img = img*2;
            img = img.^power;
    
            imax = max(img(:));
            img = img./imax;
        end
        
        idensity= idensity+1;
        if(idensity>maxItersDensity) 
            error('Not converged density filtering...')
        end
    end
    
    subplot(1, 3, 2)
    imshow(img_smooth)
    title('Smoothed Image');
    
    subplot(1, 3, 3); 
    imshow(binary_img);
    title(sprintf('Adjusted Density: %.2f%%', target_density * 100));

    img = img_smooth;

elseif(density_method == method_setToZeroForDensity)
    % Parameters
    sigma = 2;             % Standard deviation for Gaussian smoothing
    threshold = 0.55;       % Initial threshold (adjust if necessary)
    
    % Gaussian filter for smoothing
    img_smooth = imgaussfilt(img, sigma);
    
    % Threshold to binarize the image
    binary_img = img_smooth > threshold; % Convert to binary (black and white)
    
    % Compute the initial density
    current_density = sum(binary_img(:)) / numel(binary_img);
    fprintf('Initial density: %.2f%%\n', current_density * 100);
    
    % Adjust the density if needed
    maxItersDensity = 100;
    idensity = 0;
    if current_density > target_density
        % Calculate the number of pixels to retain
        num_pixels = numel(binary_img); 
        num_target_pixels = round(target_density * num_pixels); 
        
        % Find active (white) pixels
        active_indices = find(binary_img); 
        num_active_pixels = length(active_indices); 
        
        % Determine how many pixels to remove
        num_to_remove = num_active_pixels - num_target_pixels; 
        remove_indices = randperm(num_active_pixels, num_to_remove); 
        
        % Set selected pixels to black (0)
        binary_img(active_indices(remove_indices)) = 0;
    
        idensity= idensity+1;
        if(idensity>maxItersDensity) 
            error('Not converged density filtering...')
        end
    end
    
    % Display the results
    figure;
    subplot(1, 3, 1)
    imshow(img)
    title('Original Image');
    
    subplot(1, 3, 2)
    imshow(img_smooth)
    title('Smoothed Image');
    
    subplot(1, 3, 3); 
    imshow(binary_img);
    title(sprintf('Adjusted Density: %.2f%%', target_density * 100)); 

    img = binary_img;

elseif(density_method==method_tuneIntensity)
    % Calculate the current total intensity
    current_total_intensity = sum(img(:));
    
    % Scale the image intensities to match the target density
    scaling_factor = target_density / (current_total_intensity / numel(img));
    img_scaled = img * scaling_factor;
    
    % Cap values at 1 (to avoid intensities greater than the maximum)
    img_scaled(img_scaled > 1) = 1;
    
    % Display results
    figure;
    subplot(1, 2, 1); 
    imshow(img); 
    title('Original Image');
    subplot(1, 2, 2); 
    imshow(img_scaled); 
    title(sprintf('Scaled Density: %.2f%%', target_density * 100));

    img = img_scaled;
end