function dataN = normalizeDataWithOrdering(data, runningOrder, params)

dataN = NormalizeData(permute(data, runningOrder), params);
[~, ic] = sort(runningOrder);
dataN = permute(dataN, ic);
