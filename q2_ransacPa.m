%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run('C:/VLFEATROOT/vlfeat-0.9.21/toolbox/vl_setup')
% Brittany Ally
% RANSAC Image Based Stitching 
% Affine Transform
%
% 1) PREPROCESSING: load both images, covert to single and greyscale
parlleft = imread("parliament-left.jpg");
parlright = imread("parliament-right.jpg");

%resizing
parlleft = imresize(parlleft, [2500 2500]);
parlright = imresize(parlright, [2500 2500]);

pleftc = im2single(parlleft); %colour single
prightc = im2single(parlright); %colour single

pleft = im2gray(pleftc);
pright = im2gray(prightc);

% 2) Detect keypoints and extract descriptors
[fleft,dleft] = vl_sift(pleft);
[fright,dright] = vl_sift(pright);

%%
%3) Match features - here, im using vl_ubcmatch
[matches, scores] = vl_ubcmatch(dleft, dright);
%%
[mat, sco] = vl_ubcmatch(dleft, dright, 70);

%% 
%4) Prune features
pruned = mat; %top matches

%% 
%5) Robust Transformation Estimation - RANSAC

N = 20; % # of iterations - set to small number just for simplification
best = 0; %update to best sample result w/ RANSAC (# of inliers)
bestaffest = [1 1 1; 1 1 1; 0 0 1]; %best estimate affine transform matrix

%holding for: pruned - outliers
fpruned = pruned;

% STEP 1: Select minimum point set (4 points for homography, 3 for affine)
for i = 1:N
    % STEP 2: Find Best Affine Transformation
    % take 3 random sample points
    [nrows,ncols] = size(sco);
    samples = randsample(ncols, 3); %or change to randperm
    
    v1 = samples(1,1);
    v2 = samples(2,1);
    v3 = samples(3,1);
    
    m1 = pruned(:,v1);
    m2 = pruned(:,v2);
    m3 = pruned(:,v3);
    
    %get x and y vals for leftp and rightp
    lx1 = fleft(1,m1(1,1));
    lx2 = fleft(1,m2(1,1));
    lx3 = fleft(1,m3(1,1));
    
    ly1 = fleft(2,m1(1,1));
    ly2 = fleft(2,m2(1,1));
    ly3 = fleft(2,m3(1,1));
    
    rx1 = fright(1,m1(2,1));
    rx2 = fright(1,m2(2,1));
    rx3 = fright(1,m3(2,1));
    
    ry1 = fright(2,m1(2,1));
    ry2 = fright(2,m2(2,1));
    ry3 = fright(2,m3(2,1));
    
    % calculate estimated affine transformation w/ matches on other img
    % using fit geotrans because it is fast and assignment doesn't say that
    % I can't - plus below I have commented code for manual derivation

    lpts = [lx1 ly1; lx2 ly2; lx3 ly3];
    rpts = [rx1 ry1; rx2 ry2; rx3 ry3];
    tform = fitgeotrans(rpts, lpts, 'affine');
%{    
    %or manually

    sumlx = lx1 + lx2 + lx3;
    sumly = ly1 + ly2 + ly3;
    
    sumrx = rx1 + rx2 + rx3;
    sumrx = ry1 + ry2 + ry3;
    
    rstuff = [sumrx sumry 0 0 1 0; 0 0 sumrx sumry 0 1];
    lstuff = [sumlx sumly]';
    
    Acol = rstuff \ lstuff;
    %convert to sq 3x3 for input
    Asq = [ Acol(1,1) Acol(2,1) Acol(3,1); Acol(4,1) Acol(5,1) Acol(6,1);
    0 0 1];
%}
    
    % STEP 3: Determine Inliers (set threshold, to count only pts within)
    threshold = 2; %5 pixels
    testr = imwarp(pright, tform);
    [~, prunedcol] = size(pruned);
    outs = [];
    totalin = 0;
    for j = 1:prunedcol
        m = pruned(:,j);
        
        %get x and y vals for leftp and rightp
        lx = fleft(1,m(1,1));
        ly = fleft(2,m(1,1));
        
        rx = fright(1,m(2,1));
        ry = fright(2,m(2,1));
        
        [tx,ty] = transformPointsForward(tform,rx,ry);
        pts = [lx ly; tx ty];
        thedist = pdist(pts);
        
        if thedist<threshold
            totalin = totalin + 1;
        else
            outs = [outs; j];
        end
    end
    
    if totalin > best
        best = totalin;
        bestaffest = tform;
        outs = fpruned;
    end
    totalin = 0;
    % STEP 4: repeat 
end

%%
%6) Compute optimal transformation
%use all inliers
[~, outsc] = size(outs);
for h = 1:outsc
    pruned(:,outs(h,1)) = [];
end
[fprunedrow, fprunedcol] = size(pruned);
bestaff = bestaffest;
leftpts = [];
rightpts = [];
for k = 1:fprunedcol
    m = fpruned(:,k);
   
    %get x and y vals for leftp and rightp
    lx = fleft(1,m(1,1));
    ly = fleft(2,m(1,1));
        
    rx = fright(1,m(2,1));
    ry = fright(2,m(2,1));
    
    lefttemp = [lx ly];
    righttemp = [rx ry];
    
    leftpts = [leftpts; lefttemp];
    rightpts = [rightpts; righttemp];
end
bestaff = fitgeotrans(rightpts, leftpts, 'affine');

%%
%7) Create panorama
finalrr = imwarp(prightc(:,:,1), bestaffest);
finalrg = imwarp(prightc(:,:,2), bestaffest);
finalrb = imwarp(prightc(:,:,3), bestaffest);
newimg = zeros(2800,4000,1,'single');
[heigth, length, ~] = size(finalrr);
[pxx,pyy] = transformPointsForward(bestaffest,1,1);
hoff = 2500 - heigth;
loff = 2500 - length;
ix1 = pxx+501-1400;
ix2 = pxx+3000-hoff-1400;
iy1 = 1+pyy+1250;
iy2 = 2500+pyy+1250-loff;

newimg(ix1:ix2, iy1:iy2, 1) = finalrr;
newimg(ix1:ix2, iy1:iy2, 2) = finalrg;
newimg(ix1:ix2, iy1:iy2, 3) = finalrb;

newimg(501:3000, 1:2500, 1:3) = pleftc;
imshow(newimg)


