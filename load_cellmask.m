function [mask,cellArea] = load_cellmask(data,designed_mask)
% load cellmask and calculate cell area

 mask = logical(imread([data.source 'Detection' filesep 'cellmask.tif']));
 px = data.pixelSize / data.M; 
    if ~isempty(designed_mask)
        mask2 = zeros(size(mask));
        for i = designed_mask(2):designed_mask(2) + designed_mask(3)
            for j = designed_mask(1):designed_mask(1) + designed_mask(4)
                mask2(i,j) = 1;
            end
        end
        mask2 = logical(mask2);
    else
        mask2 = true(size(mask));
    end
    mask = mask & mask2;
    cellArea = sum(mask(:)) * px^2 * 1e12; 
end

