function [ship_map, centroids, areas, labels] = ...
    module6_post_processor(detection_map, min_area)
% MODULE6_POST_PROCESSOR
%   Applies morphological cleaning, connected component analysis,
%   size filtering, and centroid extraction (Sec. VI-E, Algorithm 1 step 15).
%
%   Steps:
%     1. Morphological opening (3x3) to remove isolated noise pixels
%     2. Connected component labelling
%     3. Size filter: discard components < min_area pixels
%     4. Centroid extraction

%% 1. Morphological opening
se = strel('square', 3);
ship_map = imopen(detection_map, se);

%% 2. Connected component labelling
[labels, n_comp] = bwlabel(ship_map, 8);   % 8-connectivity

%% 3. Size filter
props = regionprops(labels, 'Area', 'Centroid', 'BoundingBox', ...
                   'MajorAxisLength', 'MinorAxisLength', 'Orientation');
areas     = [props.Area];
keep_mask = areas >= min_area;

ship_map  = ismember(labels, find(keep_mask));
[labels, n_kept] = bwlabel(ship_map, 8);

%% 4. Centroid extraction
if n_kept > 0
    props2    = regionprops(labels, 'Area', 'Centroid');
    centroids = cat(1, props2.Centroid);   % [col, row] format
    areas     = [props2.Area]';
else
    centroids = zeros(0, 2);
    areas     = zeros(0, 1);
end

fprintf('  Post-processing: %d components -> %d ships (after size filter >= %d px)\n', ...
    n_comp, n_kept, min_area);

if n_kept > 0
    fprintf('  Ship centroids (col, row):\n');
    for i = 1:min(n_kept, 10)  % print first 10
        fprintf('    Ship %02d: (%.1f, %.1f) | area=%d px\n', ...
            i, centroids(i,1), centroids(i,2), areas(i));
    end
    if n_kept > 10
        fprintf('    ... and %d more.\n', n_kept-10);
    end
end
end
