function aGMM = EM_adaptation(GMM, xEst, sigma, rho)
%%-------------------------------------------------------------------------
% EM_adaptation: proposed EM adaptation algorithm
% Input:
%        GMM: generic GMM
%       xEst: new samples (it can be a clean image, pre-filtered image or example image)
%      sigma: remaining noise standard deviation in xEst
%        rho: relevance factor 
% Output:
%       aGMM: adapted GMM
%%-------------------------------------------------------------------------
xEst = im2col(xEst, [8, 8]);
xEstDC = mean(xEst);
xEst = bsxfun(@minus, xEst, xEstDC);
[d, N] = size(xEst);        % d: dimension of a patch; N: total number of patches
varMtx = sigma^2 * eye(64);
likelihoodMtx = zeros(GMM.ncomponents, size(xEst, 2));
for i = 1: GMM.ncomponents
    likelihoodMtx(i, :) = log(GMM.weights(i)) + loggausspdf2(xEst, GMM.mus(:, i), GMM.covs(:, :, i) + varMtx);
end
logsumexp_likelihoodMtx = logsumexp(likelihoodMtx, 1);
likelihoodMtx = bsxfun(@minus,likelihoodMtx,logsumexp_likelihoodMtx);
likelihoodMtx = exp(likelihoodMtx);
likelihoodMtx(likelihoodMtx == 0) = 10^-10;                              % assign a small value to avoid divide-by-zero
likelihoodMtx = bsxfun(@rdivide, likelihoodMtx, sum(likelihoodMtx, 1));  % posterior normalization
num = sum(likelihoodMtx, 2);                                             % effective number over GMM.ncomponents
likelihoodMtx = bsxfun(@rdivide, likelihoodMtx, sum(likelihoodMtx, 2));
weightsNew = num / N;                                                    % new weights
musNew = xEst * likelihoodMtx';                                          % new mus
alpha = num ./ (num + rho);                                              % weight for adaptation
weightsAdapted  = (alpha .* weightsNew) + ((1 - alpha) .* GMM.weights);
weightsAdapted = weightsAdapted / sum(weightsAdapted);
musAdapted = bsxfun(@times, musNew, alpha') + bsxfun(@times, GMM.mus, (1 - alpha)');

% % simplified covariance 
% tStart = tic;
%covsNew = zeros(1,1,N);
for i = 1: N
    covsNew(:, :, i) = xEst(:, i) * xEst(:, i)' - varMtx;
end
% disp('covsNew-1');
% disp(size(covsNew,1));
% disp('covsNew-2');
% disp(size(covsNew,2));
% disp('covsNew-3');
% disp(size(covsNew,3));
covsNew = reshape(covsNew, d * d, N);
covsNew = covsNew * likelihoodMtx';
covsNew = reshape(covsNew, [d, d, GMM.ncomponents]);                     % new covs
clear likelihoodMtx;
for k = 1: GMM.ncomponents
    covsAdapted(:, :, k) = alpha(k) * (covsNew(:, :, k)) + (1 - alpha(k)) * (GMM.covs(:, :, k) ...
        + GMM.mus(:, k) * GMM.mus(:, k)') - musAdapted(:, k) * musAdapted(:, k)';
end

% disp('covsAdapted-1');
% disp(size(covsAdapted,1));
% disp('covsAdapted-2');
% disp(size(covsAdapted,2));
% disp('covsAdapted-3');
% disp(size(covsAdapted,3));

% tElapse1 = toc(tStart);

% % unsimplified covariance 
% clear covsNew;
% tStart = tic;
% for k = 1: GMM.ncomponents
%     X_temp = bsxfun(@minus, xEst, musAdapted(:, k));
%     for i = 1: N
%         covsNew(:, :, i) = X_temp(:, i) * X_temp(:, i)' - varMtx;
%     end
%     covsNew = reshape(covsNew, d * d, N);
%     covsNew = covsNew * likelihoodMtx(k, :)';
%     covsNew = reshape(covsNew, [d, d]);
%     covsNew = topdm(covsNew);
%     covsAdapted(:, :, k) = alpha(k) * covsNew + (1 - alpha(k)) * (GMM.covs(:, :, k) + ...
%         (GMM.mus(:,k) - musAdapted(:,k)) * (GMM.mus(:,k) - musAdapted(:,k))');
% end
% tElapse2 = toc(tStart);
% fprintf('Time for simplified covariance is %.2f', tElapse1);
% fprintf('Time for unsimplified covariance is %.2f\n', tElapse2);
% fprintf('The time ratio is %.2f\n', tElapse2 / tElapse1);

for k = 1:GMM.ncomponents
    covsAdapted(:, :, k) = topdm(covsAdapted(:, :, k)); % ensure the covariance is positive semi-definite.
end

aGMM.ncomponents = GMM.ncomponents;
aGMM.mus = musAdapted;                                  % adapted: musAdapted; unadapted: GMM.mus
aGMM.covs = covsAdapted;                                % adapted: covsAdapted; unadapted: GMM.covs
aGMM.weights = weightsAdapted;                          % adapted: weightsAdapted; unadapted: GMM.weights
