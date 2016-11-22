function slices = getSlices(data)
    slices = zeros(size(data, 1)* size(data, 3), size(data, 2));
    for r=1:size(data, 2)
        currslice = data(:, r, :);
        slices(:, r) = currslice(:);
    end