function xMAP = MAP_GMM(x, y, xEst, sigmaNoise, sigma, GMM)
xEst = im2col(xEst, [8, 8]);
xEstDC = mean(xEst);
xEst = bsxfun(@minus, xEst, xEstDC);
varMtx = sigma^2 * eye(64);
likelihoodMtx = zeros(GMM.ncomponents, size(xEst, 2));
for i = 1: GMM.ncomponents
    likelihoodMtx(i, :) = log(GMM.weights(i)) + loggausspdf2(xEst, GMM.mus(:, i), GMM.covs(:, :, i) + varMtx);
end
[~, ks] = max(likelihoodMtx);

gamma = 1/ (sigma^2);

for i = 1: GMM.ncomponents
    inds = find(ks == i);
    xEstNew(:, inds) = (eye(64) + gamma * GMM.covs(:, :, i)) \ ...
        (repmat(GMM.mus(:, i), 1, length(inds)) + gamma * GMM.covs(:, :, i) * xEst(:, inds));
end

xEstNew = bsxfun(@plus, xEstNew, xEstDC);
xEstNew = scol2im(xEstNew, 8, size(x, 1), size(x, 2), 'average');
varInvNoise = 1 / (sigmaNoise^2);
xMAP = varInvNoise / (varInvNoise + gamma) *  y + gamma / (varInvNoise + gamma) * xEstNew;
%psnr = cal_psnr(x, xMAP, 0, 0);
%ssim = cal_ssim(x*255, xMAP*255, 0, 0);
