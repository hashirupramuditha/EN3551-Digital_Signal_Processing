%% Assignment 2 - Application of 2D-DCT for Image Compression

% Name        : A.A.H. Pramuditha
% Index No    : 200476P

clc;
clear;
close all;

%% 3.1) 2-D DCT - Application to Image Compression

% Load the files containing the image data
image1 = load('camera256.mat');
image2 = load('boat512.mat');
image3 = load('peppers512.mat');

% Load the extra image as a TIFF file and store it as a struct
extra_image = double(imread('lena256.tif'));
image4 = struct('data', extra_image);

% Store all the image data in a single struct object
images = struct();
images(1).Field1 = 'Image 1';
images(1).Field2 = image1;
images(2).Field1 = 'Image 2';
images(2).Field2 = image2;
images(3).Field1 = 'Image 3';
images(3).Field2 = image3;
images(4).Field1 = 'Image 4';
images(4).Field2 = image4;


% Define quality levels
quality_levels = [70, 25, 10];

% Define Standard Quality Matrix (Q_50)
Q_50 = [16, 11, 10, 16, 24, 40, 51, 61;
        12, 12, 14, 19, 16, 58, 60, 55;
        14, 13, 16, 24, 40, 57, 69, 56;
        14, 17, 22, 29, 51, 87, 80, 62;
        18, 22, 37, 56, 68, 109, 103, 77;
        24, 35, 55, 64, 81, 104, 113, 92;
        49, 64, 78, 87, 103, 121, 120, 101;
        72, 92, 95, 98, 112, 100, 103, 99;];

% Define the block size
block_size = 8;

% Image selection and compression procedure
for k = 1:4

    % Convert the struct object into a matrix
    image = struct2cell(images(k).Field2);
    data_matrix = cell2mat(image);
    
    % Display the dimensions of the matrix
    [rows, columns] = size(data_matrix);
    fprintf('Image Dimensions: %d rows x %d columns\n', rows, columns);
    
    % Calculate the number of blocks in each dimension
    x_blocks = size(data_matrix, 2) / block_size;
    y_blocks = size(data_matrix, 1) / block_size;
    
    figure();
    % Plot the original image
    % subplot(3, 4, 4*k - 3);
    subplot(2, 2, 1);
    imshow(data_matrix, []);
    title(['Original Image ', num2str(k)]);
    
    % Compress the image for a given quality level
    for n = 1:3
        
        % Calculate the required quality value
        quality_level = quality_levels(n);
        if (quality_level>50)
            J = (100-quality_level) / 50;
        else
            J = 50 / quality_level;
        end
        
        % Initialize an empty matrix to store the final image
        final_image = zeros(rows, columns);
        
        % Extract the blocks and do the operations
        for i = 1:y_blocks
            for j = 1:x_blocks

                % Extract the blocks with required size
                row_start = (i-1)*block_size + 1;
                row_end = i*block_size;
                column_start = (j-1)*block_size + 1;
                column_end = j*block_size;
                
                B = data_matrix(row_start:row_end, column_start:column_end);

                % Level-off the matrix by substracting 128 from each entry
                B_tilda = B - 128;
                
                % Apply 2-D DCT to the leveled-off matrix
                C = dct2(B_tilda);

                % Quantization
                Q = Q_50 * J;
                S = round(C ./ Q);

                % Decompress the S matrix
                R = Q .* S;

                % Apply 2-D inverse DCT  to matrix R
                E = idct2(R);

                % Correct the level-off operation done in compression stage
                decom_B = E + 128;
                
                % Reconstruct the total matrix for final image
                final_image(row_start:row_end, column_start:column_end) = decom_B;
            end
        end
        
        % Calculate the zero element percentage
        total_elements = numel(S);
        num_zeros = sum(S(:) == 0);
        percentage_zero = (num_zeros / total_elements) * 100;
        
        % Calculate Peak Signal to Noise Ratio (PSNR)
        squared_error_matrix = (final_image - data_matrix).^2;
        variance_sigma_e = mean(mean(squared_error_matrix));
        PSNR = 20 * log10(255 / (variance_sigma_e^0.5));
        
        % Plot the compressed versions of the image
        % subplot(3, 4, 4*k - 3 + n);
        subplot(2, 2, n+1);
        imshow(final_image, []);
        title(['Compressed, Quality Level: ', num2str(quality_level), ', zeros: ', num2str(percentage_zero), '%, PSNR: ', num2str(PSNR)]);
    end
end