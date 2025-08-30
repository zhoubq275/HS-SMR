function [graph_nodes, graph_edges, edge_mask, superpixel_labels] = buildMRF(D, IL, numSuperpixels)
    % ���ڳ�ʼ�Ӳ�ͼ����MRFͼ�ṹ
    % ����:
    %   D - ��ʼ�Ӳ�ͼ
    %   IL - ��ͼ
    %   numSuperpixels - ����������
    % ���:
    %   graph_nodes - ÿ�������ؿ�Ľṹ������
    %   graph_edges - �������ڽӱ�
    %   edge_mask - ��Ե��֪ƽ��Ȩ��ͼ
    %   superpixel_labels - �����طָ��ǩ����

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

    % ����ͼ�ڵ�
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

    % �����ڽӱ߽ṹ
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

    % ���ɱ߽���Ĥͼ
    boundary_mask = boundarymask(superpixel_labels); 
    edge_mask = imgaussfilt(double(boundary_mask), 1.5);  
    edge_mask = mat2gray(edge_mask);  
end