%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostic Medical Image Processing (DMIP)
% WS 2014/15
% Author: Marco Boegel, Yan Xia
% Exercise: Filtered backprojection for parallel beam geometry
% NOTE: Complete the '???' lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fbp()
clear all;
close all;
close all hidden;
clc;

showRec = 1;

im = phantom(64);

% Size of the input image
[m,n] = size(im);

% Normalize the image values
im = mat2gray(im);

figure(8);
subplot(1,3,1); imagesc(im); colormap gray; axis image; xlabel('x'); ylabel('y'); title('Original phantom image');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Forward Projection  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Extracting projections ...');
% Acquire the projection sequence
% Compute "numberOfProjections" 1-D projections along angle given the "startAngle" and the "angleIncrement".
% The result is a parallel projection of the image.
angleIncrement = 1;
startAngle = 0;
phi = startAngle;
numberOfProjections = 180;

% Create a cell array for projections
projs = cell(numberOfProjections, 1);

for i=1:numberOfProjections
    % Save the projection angles in an array
    phis(i) = phi;
    
    % Image is rotated by phi degrees counterclockwise, so the projection is lined up in the
    % columns & finally the columns are summed up
    rI = imrotate(im,phi,'bilinear','crop');
    
    if(showRec)
        figure(1);
        imagesc(rI);
        xlabel('x');
        ylabel('y');
        axis image
        colormap gray;
        title(['Angle: ', num2str(phi)]);
        drawnow;
    end
    
    % Sum up columnwise -> parallel beam
    projs{i} = sum(rI,1);
    
    % Compute the next rotation angle
    phi = phi + angleIncrement;
end

% Plot the sinogram of the image
sino = zeros(64, numberOfProjections); 
for i = 1:numberOfProjections
    sino(:,i) = projs{i};
end

figure(4);
imagesc(sino); colormap gray; axis image; xlabel('x'); ylabel('y'); title('Sinogram');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtered-Backprojection  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fltr = 0; % No filtering, no high-pass filter is applied!
%fltr = 1; % Ram-Lak
fltr = 2; % Shepp-Logan

% Reconstructed image
ct = zeros(size(im));

% Half of the dimension
dimY2 = m/2;
dimX2 = n/2;

% Compute matrix of image pixel positions
% Note: Always sample in the space where you expect the result!
[I J] = meshgrid(1:m,1:n);

x = J(:)-(dimX2+1);
y = -(I(:)-(dimY2+1));

if(showRec)
    figure(8);
    subplot(1,3,3);
    hold on;
    imagesc(x,y,im);
    hold on;
    daspect([1 1 1]);
end

imInd = [y x];
nInd = size(imInd,1);


for phi = 1:numberOfProjections
    disp('Processing projection angle: ')
    disp(phi)
    
    proj = projs{phi};
    ang = phis(phi);
    dimensions = [length(proj) length(proj)];
    
    % Detector border points
    % [-m/2,  -n/2]
    Po = [-dimensions(1)/2; -dimensions(2)/2];
    
    % [m/2, -n/2]
    Pt = [-dimensions(1)/2; dimensions(2)/2];
    
    % Which filter should be used?
    if(fltr == 1)
        fltm = RamLak(60);
        % Compute the filter step
        proj = conv(proj, fltm, 'same');
    elseif(fltr == 2)
        fltm = SheppLogan(60);
        % Compute the filter step
        proj = conv(proj, fltm, 'same');
    else
        % No filter is applied
        proj = proj;
    end
    
    rad = (pi/180)*ang;
    
    % Instead of rotating the image, the projection is rotated with the rotation matrix.
    % Define the rotation matrix
    R = [cos(rad), -sin(rad);
        sin(rad), cos(rad)];
    
    % Rotate the corner pixel positions of the 2-D image that we reconstruct.
    pPo = (R*Po);
    pPt = (R*Pt);
    
    if(showRec)
        figure(8);
        subplot(1,3,3);
        lim = ceil(sqrt(size(im,1)^2+size(im,2)^2)/100)*100;
        axis([-lim lim -lim lim]);
        plot(pPo(2),pPo(1),'--ro');
        hold on
        plot(pPt(2),pPt(1),'--bo');
        hold on
        plot([pPo(2) pPt(2)],[pPo(1) pPt(1)],'-g');
        hold on;
        legend({'P_o','P_t','Detector'});
    end
    
    % Back projection:
    % Compute direction and normal of detector
    dirDet = pPt-pPo;
    dirDet = dirDet/norm(dirDet);
    
    normalDet = [-dirDet(2);dirDet(1)];
    normalDet = normalDet/norm(normalDet);
    
    % Use Hesse normal form to calculate the distance of the image
    % points to the detector. Use it to project the points onto the
    % detector and to resolve the position on the detector line.
    
    % d for Hesse normal form
    d = pPo'*normalDet;
    % Distances between image points and detector
    dd = imInd*normalDet-d*ones(nInd,1);
    % Projection of image points onto detector
    pDet = imInd-repmat(normalDet',nInd,1).*repmat(dd,1,2);
    % Distances between projected points and detector origin
    dis = pDet-repmat(pPo',nInd,1);
    ts = sqrt((dis(:,1).*dis(:,1))+(dis(:,2).*dis(:,2)));
    
    % Calculate left and right neighbors on the detector line grid to
    % get affected detector cells
    li = floor(ts);
    li(li<1)=1;
    li(li>dimensions(2))=dimensions(2);
    ui = li+1;
    ui(ui>dimensions(2))=dimensions(2);
    
    % Calculate weights between neighbored cells
    ld = abs(ts-li);
    ud = abs(ts-ui);
    
    % Interpolate between detector cells and create the current FBP.
    fbpv = (ones(nInd,1)-ld).*proj(li)' + (ones(nInd,1)-ud).*proj(ui)';
    
    FBP = reshape(fbpv,m,n)';
    
    % Accumulate FBPs
    ct = ct + FBP;
    
    if(showRec)
        figure(8);
        subplot(1,3,2); imagesc(ct); axis image; colormap gray; xlabel('x'); ylabel('y'); title(['Projection: ', num2str(phi), '/', num2str(numberOfProjections)]);
        drawnow;
        hold off
    end
end % loop projections

ct = ct/numberOfProjections; % normalize according to projections
figure(8);
subplot(1,3,2); imagesc(ct); axis image; colormap gray; xlabel('x'); ylabel('y'); title(['Projection: ', num2str(phi), '/', num2str(numberOfProjections)]);
drawnow;
hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [shepp] = SheppLogan(width)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shepp] = SheppLogan(width)
% Shepp-Logan filtering with kernel size "width" is applied.
    t = -floor((width - 1)/2);
    for n=1:width
        shepp(n) = -(2/pi^2)*1/(4*t^2-1);
        t = t + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ramlak] = RamLak(width)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ramlak] = RamLak(width)
% Ram-Lak filtering with kernel size "width" is applied.
    t = -floor((width - 1)/2);
    for n=1:width
        if (t == 0)
            ramlak(n) = 1/4;
        elseif (mod(t,2) == 0)
            ramlak(n) = 0;
        else
            ramlak(n) = -1/(pi^2*t^2);
        end
        t = t+1;
    end
end