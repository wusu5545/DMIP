%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostic Medical Image Processing                      
% WS 2014/15                                           
% Exercise: rigid registration
% NOTE: Complete the '???' lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function rigid_registration()
close all;
clear all;
clc;

   %% transform a phantom image using imtransform (Rotate 45 degree and translate with t=(20, 1))
        % Load Shepp-logan phantom with size of 256
        image = phantom(256);
        figure(1);
        colormap bone;
        subplot(2,2,1);
        imagesc(image);
        title('original image');
        % Rotation angle
        a = pi/4;
        % Construct a 3x3 transformation matrx that contains the rotation matrix and translation vector
        tmat = [cos(a),-sin(a),0;sin(a),cos(a),0; 20,1,1];
        % Create spatial transformation structure (TFORM)
        tform = maketform('affine', tmat);
        transimage2 = imtransform(image, tform);
        subplot(2,2,2);
        imagesc(transimage2);
        title('transformed image 2');
        
     %%   
    % initial information
    is = 256;
    paddings = 32;
    isp = is + paddings*2;
    rotation = 45;
    transx = 2;
    transy = 1;
    startpositions = [10 10 10];
       
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

% Gaussian filtering of the image to avoid local minima during optimization and transform it using
% transform(), rotation, transx, and transy to obtain a moving image.
    Image1 =padarray(mat2gray(phantom(is)), [paddings,paddings]);
    nice = fspecial('gaussian', 10, 4);
    Image1 = imfilter(Image1,nice,'same'); 
    Image2 = transform(Image1,rotation,transx,transy)
    

    % visualization
    figure(2);
    colormap bone;
    subplot(2,3,1);
    imagesc(Image1);
    title('original image');
    subplot(2,3,2);
    imagesc(Image2);
    title('transformed image');
    Image3 = transform(Image2, startpositions(1), startpositions(2), startpositions(3));
    subplot(2,3,3);
    imagesc(Image3);
    title('start image');
    
    values = [];
    
    % Perform the optimization using function fminsearch 
    [x, fval] = fminsearch(@(x)sim(x,Image1,Image2),startpositions,optimset('MaxFunEvals',150));
    
   
   
    % define similarity function (SSD or MI)called by fminsearch
    function d = sim(position, fiximage, movingimage)
        
        % apply current transformation
        moved = transform(movingimage, position(1), position(2), position(3));
        
        % distance measure Sum of Squared Distances (SSD)
              diff = (moved - fiximage).^2;
              d = sum(sum(diff))/is(1)/is(2);

        % distance measure mutual information (MI)
            % calculate the difference between fiximage and moved
            %diff = abs(fiximage - moved);
            % joint histogram of fiximage and moved
            %???
            % joint entropy of the joint histogram of fiximage and moved
            %???
            % mutual information formula
            %???
        
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

        % joint histogram of 2 images
        function h = jointH(im1, im2)
        is = size(Image1);
        h = zeros(256);
        %???
        end
        
        % marginal entropy of one image
        function h = marginalE(hxy)
           marginal = sum(hxy);
           h=0;
           %???
        end
  end