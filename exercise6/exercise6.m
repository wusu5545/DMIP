%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostic Medical Image Processing (DMIP)
% WS 2014/15
% Exercise: Filtered Backprojection for fan beam geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fbp()
clear all;
close all;
close all hidden;
clc;

showRec = 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Initialise System Parameters + Load Sinogram  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projections have been created using the "CONRAD" software
% Here we set the parameters used in CONRAD and load the sinogram
% System parameters:
angleIncrement = 1; % angular increment
focalLength = 600; % source to detector distance
detLength = 200; % length of the detector
numberOfProjections = 134; % number of projections

% The scan was acquired as a short-scan, therefore we need to calculate the
% fan angle "halfFanAngle" and the angular short-scan range "maxBeta"
% TASK: Compute "halfFanAngle"
halfFanAngle = atan((detLength/2) / focalLength);
disp(['Half the fan angle: ' num2str(halfFanAngle*180/pi)]);
% TASK: Compute the short-scan range "maxBeta
maxBeta = pi + 2 * halfFanAngle;%180 degrees is not enough
disp(['Short-scan range: ' num2str(maxBeta*180/pi)]);

% From the short scan range we can now compute the individual angles for
% each projection
phis = linspace(0,maxBeta,numberOfProjections);

% TASK: Compute "detectorPositions"
% These are the distances from detector boundary to the detector cell
% centers, used to calculate ray angles (parkerWeighting, backprojection)
% (note: centered around 0, as central ray hits center of the detector)
detectorPositions = ((1:detLength) - (detLength/2 + 0.5))';

% Load the projections from the attached projection tif files
% Available: "Sinogram0.tif","Sinogram1.tif","Sinogram2.tif"
projs = double(imread('Sinogram2.tif'));

% Visualize the sinogram
if(showRec)
    figure(1);
    imagesc(projs);
    xlabel('Rotation angle');
    ylabel('Projections');
    axis image
    colormap gray;
    title('Full Sinogram');
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sinogram Weighting (Fan Beam) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% As we are now in the fan-beam case we need to apply different weights to
% the projection data, due to redundant rays and to compensate for different path lengths

% Get the cosine weights and multiply with the sinogram
cosW = cosineWeights(detectorPositions,focalLength,numberOfProjections);
projs_cos = projs.*cosW;

% Get the parker redundancy weights and multiply with the cosine weighted
% sinogram
parkerW = parkerWeights(detectorPositions,focalLength,phis,halfFanAngle);
projs_afterWeighting = projs_cos.*parkerW;

% Visualize weights and weighted sinogram
figure
subplot(2,1,1)
imagesc(cosW)
xlabel('Rotation angle');
ylabel('Detector position');
title('Cosine Weights');
subplot(2,1,2)
imagesc(parkerW)
xlabel('Rotation angle');
ylabel('Detector position');
title('Parker Weights');

figure
subplot(2,1,1)
imagesc(projs)
xlabel('Rotation angle');
ylabel('Detector position');
colormap gray;
title('Sinogram before weighting');
subplot(2,1,2)
imagesc(projs_afterWeighting)
xlabel('Rotation angle');
ylabel('Detector position');
colormap gray;
title('Sinogram after weighting');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtered Backprojection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define the output size of the reconstructed image
m=128;
n=128;
% With "scaling" we can increase/decrease the reconstruction resolution
scaling = 1;

% choose a reconstruction filter
fltr = 1;
%fltr = 0; % No filtering, no high-pass filter is applied!
%fltr = 1; % Ram-Lak
%fltr = 2; % Shepp-Logan

% Reconstructed image
ct = zeros(m*scaling,n*scaling);

% Half of the dimension
dimY2 = m/2;
dimX2 = n/2;

% Compute matrix of image pixel positions
% Note: Always sample in the space where you expect the result!
[I J] = meshgrid(linspace(1,m,m*scaling),linspace(1,n,n*scaling));

% Image pixels are centered around the origin
x = J(:)-(dimX2+0.5);
y = -(I(:)-(dimY2+0.5));

imInd = [x y];
nInd = size(imInd,1);

for phi = 1:numberOfProjections
    
    % Extract detector line "proj" and angle "rad" of this projection
    proj = projs_afterWeighting(:,phi);
    % "-maxBeta" rotates the output 
    % (Necessary as CONRAD angles apparently start at "-maxBeta")
    rad = phis(phi)-maxBeta;
    dimensions = [length(proj) length(proj)];
    
    % do the filtering on projections
    % Which filter should be used?
    if(fltr == 1)
        % Compute the filter step
        proj = RamLak(proj);
    elseif(fltr == 2)
        proj = SheppLogan(proj);
    else
        % No filter is applied
        proj = proj;
    end
    
    
    % Detector corner points (Central ray goes to center of detector)
    Po = [0; dimensions(2)/2];
    Pt = [0; -dimensions(2)/2];
    
    % Source point at oposite side
    S = [focalLength; 0];
    
    % Instead of rotating the image, the projection is rotated with the rotation matrix.
    % Define the rotation matrix
    R = [cos(rad), -sin(rad);
        sin(rad), cos(rad)];
    
    % Rotate the corner pixel positions of the 2-D image that we reconstruct.
    pPo = (R*Po);
    pPt = (R*Pt);
    pS = (R*S);
    
    if(showRec)
        figure(8);
        subplot(1,2,2);
        lim = ceil(sqrt(m^2+n^2)/100)*100;
        axis([-lim lim -lim lim]);
        plot(pPo(2),pPo(1),'--ro');
        %legend({'P_o'},'Location','EastOutside')
        hold on
        plot(pPt(2),pPt(1),'--bo');
        %legend({'P_o','P_t'},'Location','EastOutside')
        hold on
        plot([pPo(2) pPt(2)],[pPo(1) pPt(1)],'-g');
        hold on;
        %legend({'P_o','P_t','Detector'});
    end
    
    % Back projection:
    % Compute direction and normal of detector
    dirDet = pPt-pPo;
    dirDet = dirDet/norm(dirDet);
    
    % Hesse normal form: r*n0 - d = 0
    % d for Hesse normal form (distance of detector to origin - dot product)
    normalDet = [-dirDet(2);dirDet(1)];
    normalDet = normalDet/norm(normalDet);
    d = pPo'*normalDet;
    
    % intersect all rays with the detector
    % 1) get line equations for all rays
    % get ray directions (Source to pixel positions)
    StoPixels = (imInd-repmat(pS',nInd,1));
    StoPixelsN = StoPixels./repmat(sqrt(sum(StoPixels.^2,2)),1,2);
    % obtain the ray normals
    normalPix = StoPixelsN * [0,-1; 1,0]';
    % obtain line distances to world coordinate center
    dPix = pS'*normalPix';
    
    % TASK: 2) compute intersection points of ray lines and detector line
    pDet = zeros(nInd,2);
    for i=1:nInd
        pDet(i,:)=([normalDet';normalPix(i,:)]\[d;dPix(i)])';
    end
    
    % Distances between projected points and detector origin
    dis = pDet-repmat(pPo',nInd,1);
    % dis contains the vectors along the detector -> get eucledian lengths
    ts = sqrt((dis(:,1).*dis(:,1))+(dis(:,2).*dis(:,2)));
    
    % Get the interpolated value on the detector and assign to
    % corresponding pixel
    % use built in Matlab interpolation
    fbpv = interp1((1:size(proj,1))-0.5,proj,ts,'linear',0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% DISTANCE WEIGHTING %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply the distance weights necessary for fan-beam backprojection
    % Orthogonal project Source-to-Pixel vectors to central ray
    % 1) Transform pixel positions to radius and angle representation
    r = abs(imInd(:,1)+1i*imInd(:,2));
    alphas = pi/2 + angle(imInd(:,1)+1i*imInd(:,2));
    % TASK: 2) Compute the distance weights U for all pixels
    U = (focalLength + r .* sin (rad - alphas)) / focalLength;
    
    % Correct the pixel values with the distance weight
    fbpvWeighted = fbpv./(U.^2);
    
    % Accumulate FBPs
    FBP = reshape(fbpvWeighted,m*scaling,n*scaling);
    ct = ct + FBP;
    
    if(showRec)
        figure(8);
        subplot(1,2,1); imagesc(ct'); axis image; colormap gray; xlabel('x'); ylabel('y'); title(['Projection: ', num2str(phi), '/', num2str(numberOfProjections)]);
        drawnow;
        hold off
    end
end % loop projections

% adjust scale (needed because we have short-scan)
ct = ct*(maxBeta/numberOfProjections);
figure(8);
subplot(1,2,1); imagesc(ct'); axis image; colormap gray; xlabel('x'); ylabel('y'); title(['Projection: ', num2str(phi), '/', num2str(numberOfProjections)]);
drawnow;
hold off

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [cosW] = cosineWeights(detSize,nrProjections,fanAngle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cosW] = cosineWeights(detPositions,srd,nrProjections)
% Cosine weights for fan angle "fanAngle" and given detector Positions

sourceDetectorDistance = srd;

% TASK: Compute the cosine weights for ONE projection
cosW = sourceDetectorDistance ./ sqrt(sourceDetectorDistance.^2+detPositions.^2);

% Make sure "cosW" is a column vector and extend weights for all
% projections
cosW = squeeze(cosW);
cosW = repmat(cosW,1,nrProjections);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [parkerW] = parkerWeights(detSize,nrProjections,fanAngle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [parkerW] = parkerWeights(detPositions, srd, projAngles, halfFanAngle)
% Parker redundancy weights for fan angle "fanAngle" and given size.

% Calculate ray angles for each detector bin
sourceDetectorDistance = srd;
gammas = atan(detPositions/sourceDetectorDistance);

% extend projection angles and ray angles to the full sinogram
betas = projAngles;
gammas = repmat(gammas,1,length(projAngles));
betas = repmat(betas,length(detPositions),1);

% Compute parker weights
% initialize with ones: --> includes already the second requirement
parkerW = ones(length(detPositions),length(projAngles));
% TASK: Implement the first requirement
% 1st case: 0 <= beta <= (2*delta-2*gamma)
mask = (betas>=0 & betas<=(2*halfFanAngle-2*gammas));%mask is the condition index
parkerW(mask) = sin((pi/4)*betas(mask)./(halfFanAngle - gammas(mask))).^2;
% TASK: Implement the third requirement
% 3rd case: pi-2*gamma <= beta <= pi+2*halfFanAngle
mask = (betas>=pi-2*gammas & betas <= pi+2*halfFanAngle);
parkerW(mask) = sin((pi/4)*(pi+2*halfFanAngle - betas(mask))./(halfFanAngle + gammas(mask))).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [out] = SheppLogan(proj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out] = SheppLogan(proj)
% Shepp-Logan filtering with kernel size "width" is applied.

width = length(proj);
t = -floor((width-1)/2);

for tap=1:width
    shepp(tap) = -2/(pi^2*(4*t^2-1));
    t = t + 1;
end

out = real(ifft(fft([shepp.';zeros(width,1)]).*fft([proj;zeros(width,1)])));
out = out(width/2:width+width/2-1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [out] = RamLak(width)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out] = RamLak(proj)
% Ram-Lak filtering with kernel size "width" is applied.

width = length(proj);
t = -floor((width-1)/2);

for tap=1:width
    if(t==0)
        ramlak(tap) = 1/4;
    elseif(mod(t,2)==0)
        ramlak(tap) = 0;
    else
        ramlak(tap) = -1/(pi*t)^2;
    end
    t = t + 1;
end

out = real(ifft(fft([ramlak.';zeros(width,1)]).*fft([proj;zeros(width,1)])));
out = out(width/2:width+width/2-1);

end
