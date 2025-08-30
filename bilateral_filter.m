function filtered_image = bilateral_filter(image, sigma_space, sigma_intensity)  
    % 确定核大小  
    kernel_size = 2 * ceil(sigma_space) + 1;  
    half_kernel_size = floor(kernel_size / 2);  
      
    % 获取图像尺寸  
    [rows, cols] = size(image);  
      
    % 初始化结果图像和权重图  
    result = zeros(rows, cols);  
      
    % 对每个像素应用双边滤波  
    for x = 1:rows  
        for y = 1:cols  
            Gspace_sum = 0;  
            pixel_sum = 0;  
            W_sum = 0;  
              
            % 遍历邻域内的像素  
            for dx = -half_kernel_size:half_kernel_size  
                for dy = -half_kernel_size:half_kernel_size  
                    nx = x + dx;  
                    ny = y + dy;  
                      
                    % 检查索引是否有效  
                    if nx > 0 && nx <= rows && ny > 0 && ny <= cols  
                        % 计算空间高斯权重  
                        Gspace = exp(-0.5 * ((dx^2 + dy^2) / (sigma_space^2)));  
                          
                        % 获取邻域像素值  
                        shifted_pixel = image(nx, ny);  
                          
                        % 计算强度差异  
                        intensity_diff = image(x, y) - shifted_pixel;  
                          
                        % 计算强度高斯权重  
                        Gintensity = exp(-0.5 * (intensity_diff^2 / (sigma_intensity^2)));  
                          
                        % 累加权重和像素值  
                        weight = Gspace * Gintensity;  
                        pixel_sum = pixel_sum + weight * shifted_pixel;  
                        W_sum = W_sum + weight;  
                    end  
                end  
            end  
              
            % 如果权重和大于0，则计算加权平均值  
            if W_sum > 0  
                result(x, y) = pixel_sum / W_sum;  
            else  
                % 防止除以零，使用原始像素值  
                result(x, y) = image(x, y);  
            end  
        end  
    end  
      
    % 返回结果图像  
    filtered_image = result;  
end