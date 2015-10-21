%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostic Medical Image Processing                      
% WS 2014/15                                           
% Exercise: rigid registration - MI
% NOTE: Complete the '???' lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function rigid_registration()
close all;
clear all;
clc;

   
     %%
    % initial information
    is = 256;
    paddings = 32;
    isp = is + paddings*2;
    
    startpositions = [5 5 5];
       
    % transformation function
    function t = transform(img, rot, transx, transy)
        % rotation
        t = imrotate(img, rot, 'crop'); % crop guarantees original size
        % if NaN appears
        t(isnan(t))=0;
        
        % translation
        is =size(t);
        [X,Y] = meshgrid(1:is(1), 1:is(2));
        t = interp2(X, Y, t, X+transx, Y+transy);
        % if NaN appears
        t(isnan(t))=0;
    end

% Gaussian filtering of the images to avoid local minima during optimization

    Image1 = imread('T1.png','png');
    Image1 =padarray(mat2gray(Image1), [paddings,paddings]);
    I1 = Image1;
    
    Image2 = imread('Proton.png','png');
    Image2 =padarray(mat2gray(Image2), [paddings,paddings]);
    I2 = Image2;
    
    size(Image1)
    nice = fspecial('gaussian', 10, 4);
    Image1 = imfilter(Image1, nice, 'same');
    Image2 = imfilter(Image2, nice, 'same');
    

    % visualization
    figure(2);
    colormap gray;
    subplot(2,3,1);
    imagesc(Image1);
    title('filtered brainT1 image');
    subplot(2,3,2);
    imagesc(Image2);
    title('filtered brainProton image');
    Image3 = transform(Image2, startpositions(1), startpositions(2), startpositions(3));
    subplot(2,3,3);
    imagesc(Image3);
    title('start image');
    
    values = [];
    
    % Perform the optimization using function fminsearch
    [x, fval] = fminsearch(@(x)sim(x,Image1,Image2), startpositions,optimset('MaxFunEvals',150));
     
    
     
    % define sim called by fminsearch
    function d = sim(position, fiximage, movingimage)
        
        % apply current transformation
        moved = transform(movingimage, position(1), position(2), position(3));
        
        % distance measure SSD
              diff = (fiximage - moved).^2;
              d = sum(sum(diff)) / is(1) /is(2);

        % distance measure mutual information  
            % calculate the difference between fiximage and moved
            diff = abs(fiximage - moved);
            % joint histogram of fiximage and moved
            hxy = jointH(im2uint8(fiximage), im2uint8(moved)); 
            % joint entropy of the joint histogram of fiximage and moved
            jointE = -sum(sum(hxy.*log2(hxy + (hxy == 0))));
            % mutual information formula
            d =  marginalE(hxy) + marginalE(hxy') - jointE;
            d = -d;
        
        values = [values d];
        
        % show images after each iteration step
        subplot(2,3,4);
        imagesc(moved);
        title('moving image');    
        subplot(2,3,5);
        imagesc(diff);
        title('difference image');
        subplot(2,3,6);
        plot(values);
        title('distance visualization');
        drawnow;     
    end


    hFig = figure(3);
    set(hFig, 'Position', [200 200 550 250])
    axis image
    colormap gray;
    subplot(1,2,1);
    I1 = imcrop(I1,[50 50 220 220]);
    imagesc(I1);
    title('brainT1 image');
    subplot(1,2,2);
    I2 = transform(I2, x(1), x(2), x(3));
    I2 = imcrop(I2,[50 50 220 220]);
    imagesc(I2);
    title('registered brainProton image');
    
    
        % joint histogram of 2 images
        function h = jointH(im1, im2)
        is = size(Image1);
        h = zeros(256);
        for i = 1:is(1)
            for j = 1:is(2)
                h(im1(i,j) + 1,im2(i,j) + 1) = h(im1(i,j) + 1,im2(i,j) + 1) + 1;
            end
        end
        h = h ./ is(1) /is(2);
        end
        
        % marginal entropy of one image
        function h = marginalE(hxy)
            marginal = sum(hxy);
            h=0;
            for k = 1 : 256
                if (marginal(k) ~= 0)
                    h = h + marginal(k) * log2(marginal(k));
                end
            end
            h = -h;
        end
  end