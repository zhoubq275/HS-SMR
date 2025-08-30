function [graph_nodes, graph_edges, edge_mask, superpixel_labels] = buildMRF(D, IL, numSuperpixels)
    % 基于初始视差图构建MRF图结构
    % 输入:
    %   D - 初始视差图
    %   IL - 左图
    %   numSuperpixels - 超像素数量
    % 输出:
    %   graph_nodes - 每个超像素块的结构体数组
    %   graph_edges - 超像素邻接边
    %   edge_mask - 边缘感知平滑权重图
    %   superpixel_labels - 超像素分割标签矩阵

    D_clean = D;
    D_clean(isnan(D)) = 0;  
    D_norm = mat2gray(D_clean);  

    [rows, cols] = size(D);
    [superpixel_labels, numLabels] = superpixels(D_norm, numSuperpixels, 'Compactness', 18);

    if size(IL, 3) == 3
        IL_gray = rgb2gray(IL);
    else
        IL_gray = IL;
    end
    [Gx, Gy] = imgradientxy(IL_gray, 'sobel');
    GradMag = sqrt(Gx.^2 + Gy.^2);

    % 构建图节点
    graph_nodes = struct('index', {}, 'pixels', {}, 'meanDisparity', {}, ...
                         'centroid', {}, 'localGradient', {});

    for k = 1:numLabels
        pixels = find(superpixel_labels == k);
        [y, x] = ind2sub([rows, cols], pixels);
        graph_nodes(k).index = k;
        graph_nodes(k).pixels = pixels;
        graph_nodes(k).meanDisparity = mean(D(pixels), 'omitnan');
        graph_nodes(k).centroid = [mean(x), mean(y)];
        graph_nodes(k).localGradient = mean(GradMag(pixels));
    end

    % 构建邻接边结构
    adjMatrix = sparse(numLabels, numLabels);
    neighbor_offsets = [0 1; 1 0; 0 -1; -1 0]; 
    for i = 1:rows
        for j = 1:cols
            label = superpixel_labels(i, j);
            for d = 1:size(neighbor_offsets, 1)
                ni = i + neighbor_offsets(d, 1);
                nj = j + neighbor_offsets(d, 2);
                if ni >= 1 && ni <= rows && nj >= 1 && nj <= cols
                    neighbor_label = superpixel_labels(ni, nj);
                    if label ~= neighbor_label
                        adjMatrix(label, neighbor_label) = 1;
                    end
                end
            end
        end
    end
    [r, c] = find(triu(adjMatrix));
    graph_edges = [r, c];

    % 生成边界掩膜图
    boundary_mask = boundarymask(superpixel_labels); 
    edge_mask = imgaussfilt(double(boundary_mask), 1.5);  
    edge_mask = mat2gray(edge_mask);  
end