function [phi,normaliser,pred] = hmmfilter(s,db,td,LIK)
%HMMFILTER Perform the forward sweep of the filtering.
%
%   This function is called by hmmgeolocate.m and likelihood.m
%
%   This function should not be called manually by the user.
%
%   Date: 7/8 - 2009, ver. 0.54
%   HMM geolocation toolbox, DTU Informatics and DTU Aqua

names = fieldnames(LIK);
names = names(~strcmp(names,'type'));
names = names(~strcmp(names,'mode'));
numnames = length(names);

% Initialise
[row,col] = size(db.depth);
icalc     = length(td.d24);
phi       = zeros(row,col,icalc); % Density function
pred      = phi;

% Posterior probability for initial position of fish
post = zeros(row,col); post(td.y0,td.x0) = 1; % Dirac delta

% Store probability distribution
phi(:,:,1) = post;
normaliser = ones(1,icalc-1);

par1.covmat = 2*s(1)*eye(2);
par2.covmat = 2*s(2)*eye(2);
kern1 = makekern2(par1);
kern2 = makekern2(par2);


% We know the initial distribution so iterations start at 2
for i=2:icalc 
    % Solve forward equations by convolution
    % Make time update, probability is spread out, prediction
    if     td.behav(i-1) == 1
        P = conv2(post,kern1,'same');
    elseif td.behav(i-1) == 2
        P = conv2(post,kern2,'same');
    end
    % Remove islands and normalise
    P(db.land) = 0;
    %P = normalise(P); % corrected on Oct 13th 2008 - this changes D
    %estimates to larger values because they are no longer "punished" for
    %spreading the probability on to land. Solution to this is to change
    %the program so that boundaries are taken more nicely care of (by doing
    %many small convolutions).
    pred(:,:,i) = P;
    % Use Bayes' theorem to find probability conditioned on depth 
    ltemp = ones(row,col);
    for j = 1:numnames
        ltemp = ltemp .* LIK.(names{j})(:,:,i-1); % for some reason it is faster to save LIK
    end
    post = ltemp .* P;
    %if sum(isnan(post(:))) ~= 0, error('NaN found in predicted distribution at day %i',i), end
    if sum(post(:)) == 0, error('Zero distribution at day %i, s = [%f,%f]\nIf you are estimating D try changing the bounds (in hmmgeolocate)',i,s(1),s(2)), end
    % Store likelihood function to be optimized later for D
    % The normalising constant corresponds to the conditional distribution
    % of the depth given the previous measurements
    [post,normaliser(i-1)] = normalise(post);
    % Store probability distribution
    phi(:,:,i) = post;
end