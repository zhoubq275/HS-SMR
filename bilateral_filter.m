function filtered_image = bilateral_filter(image, sigma_space, sigma_intensity)  
    % ȷ���˴�С  
    kernel_size = 2 * ceil(sigma_space) + 1;  
    half_kernel_size = floor(kernel_size / 2);  
      
    % ��ȡͼ��ߴ�  
    [rows, cols] = size(image);  
      
    % ��ʼ�����ͼ���Ȩ��ͼ  
    result = zeros(rows, cols);  
      
    % ��ÿ������Ӧ��˫���˲�  
    for x = 1:rows  
        for y = 1:cols  
            Gspace_sum = 0;  
            pixel_sum = 0;  
            W_sum = 0;  
              
            % ���������ڵ�����  
            for dx = -half_kernel_size:half_kernel_size  
                for dy = -half_kernel_size:half_kernel_size  
                    nx = x + dx;  
                    ny = y + dy;  
                      
                    % ��������Ƿ���Ч  
                    if nx > 0 && nx <= rows && ny > 0 && ny <= cols  
                        % ����ռ��˹Ȩ��  
                        Gspace = exp(-0.5 * ((dx^2 + dy^2) / (sigma_space^2)));  
                          
                        % ��ȡ��������ֵ  
                        shifted_pixel = image(nx, ny);  
                          
                        % ����ǿ�Ȳ���  
                        intensity_diff = image(x, y) - shifted_pixel;  
                          
                        % ����ǿ�ȸ�˹Ȩ��  
                        Gintensity = exp(-0.5 * (intensity_diff^2 / (sigma_intensity^2)));  
                          
                        % �ۼ�Ȩ�غ�����ֵ  
                        weight = Gspace * Gintensity;  
                        pixel_sum = pixel_sum + weight * shifted_pixel;  
                        W_sum = W_sum + weight;  
                    end  
                end  
            end  
              
            % ���Ȩ�غʹ���0��������Ȩƽ��ֵ  
            if W_sum > 0  
                result(x, y) = pixel_sum / W_sum;  
            else  
                % ��ֹ�����㣬ʹ��ԭʼ����ֵ  
                result(x, y) = image(x, y);  
            end  
        end  
    end  
      
    % ���ؽ��ͼ��  
    filtered_image = result;  
end