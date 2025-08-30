function [D_opt, E_opt, iter] = optimizeDisparity(Efun, D_init, IL, IR, edge_mask)
    % LMM 迭代优化主函数
    % 输入：
    %   D_init          - 初始视差图
    %   IL, IR      - 左右图像
    %   edge_mask   - 超像素边缘掩码
    % 输出：
    %   D_opt       - 优化后的视差图
    %   E_opt     - 最终能量值
    %   iter        - 实际迭代次数

    %% 参数设置
    max_iter       = 10000;
    epsilon        = 1e-5;
    alpha_init     = 0.01;
    alpha_decay    = 0.5;  
    rho            = 0.5;
    sigma          = 0.4;

    % 双边滤波参数
    sigma_space     = 0.7;
    sigma_intensity = 2.2;

    % 区域感知平滑权重
    lambda_min = 0.4;
    lambda_max = 0.6;
    lambda_matrix = lambda_min + (lambda_max - lambda_min) .* (1 - edge_mask);

    % 初始化变量
    D     = D_init;
    alpha = alpha_init;
    iter  = 0;

    %% 主迭代过程
    while iter < max_iter
        [E,grad] = buildEnergy(D, IL, IR, lambda_matrix);

        % 收敛判断
        grad_norm = norm(grad);
        if grad_norm < epsilon
            fprintf('[Iter %d] Gradient below threshold. Converged.\n', iter);
            break;
        end

        m = 0;
        success = false;

        while m < 5
            delta_D = -alpha * reshape(grad, size(D));
            D_new   = D + delta_D;
            E_new = buildEnergy(D_new, IL, IR, lambda_matrix);

            if E_new < E 
                % 应用双边滤波
                D_filtered = bilateral_filter(D_new, sigma_space, sigma_intensity);
                D_filtered;
                success = true;
                break;
            else
                alpha = alpha * alpha_decay; 
                m = m + 1;
            end
        end

        if ~success
            fprintf('[Iter %d] Line search failed to improve energy.\n', iter);
            break;
        end

        
        iter = iter + 1;
        fprintf('[Iter %d] Energy: %.6f | Grad norm: %.6f\n', iter, E_new, grad_norm);
    end

    % 更新视差图
    D_opt   = D;
    E_opt = E_new;
end