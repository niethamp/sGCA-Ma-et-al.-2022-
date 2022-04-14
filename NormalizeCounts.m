function [normCounts] = NormalizeCounts(counts)
%ratio of medians
% estimate pseudo-reference with geometric mean row by row
pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1);
% transform to common scale
normCounts = bsxfun(@rdivide,counts,sizeFactors);
end

