%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostic Medical Image Processing (DMIP) 
% WS 2014/15
% Author: Marco Boegel & Yan Xia
% Exercise: RANSAC
% NOTE: Complete the '???' lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ransacline_ex
    clear all;
    close all;
    clc;
    
    % Point amount for model estimation (minimum number according to
    % RANSAC), probability for outliers, probability for correct model
    mn = 2;
    oP = 0.2;
    cP = 0.9999;
    
    % Data points to fit
    pts = [0   0;
           1   1;
           2   2;
           3   3;
           3.2 1.9;
           4   4;
           10  1.8];
       
    % Draw estimated line considering all points   
    figure(1);
    hold on;
    title('RANSAC for line fitting');
    xlabel('x');
    ylabel('y');
    
    x = [-min(pts(:,1))-5 max(pts(:,1))+5];
    aMdl = fitline(pts);
    am = aMdl(1);
    at = aMdl(2);
    
    ay = x*am+at;
    
    plot(x,ay);
    plot(pts(:,1)',pts(:,2)','x');
       
    % Estimate best model using RANSAC
    [eMdl mPts err] = commonransac(pts, mn, oP, cP, @fitline, @lineerror);
    
    % Draw estimated line and outliers
    em = eMdl(1);
    et = eMdl(2);
    
    ey = x*em+et;
    
    plot(x,ey,'r');
    plot(mPts(:,1)',mPts(:,2)','or');
    
    legend('Line using all points', 'Data ponts', 'RANSAC estimated line', 'Points used for model estimation', 'Location', 'NorthWest');
    
    err
    em
    et
    
    
end

function [bMdl mPts err] = commonransac(data, mn, oRelFr, cProb, mdlEstFct, errFct)
% RANSAC function to return best estimated model parameters
% 
%   [bMdl mPts err] = commonransac(data, mn, oProb, cProb, mdlEstFct,
%   errFct) returns the best estimated model parameters bMdl, the used
%   sample points mPts for building the model and the error err of the
%   estimated model. Used parameters are data as data points, mn as minimum
%   number of data points required to build the model, oRelFr as relative
%   frequency of outliers, cProb as probability to randomly choose the best
%   estimated model in one of the iterations, mdlEstFct as function for
%   model estimation and errFct as error function. data is assumed to
%   represent one data point in each row. mdlEstFct and errFct
%   require the following interfaces:
%   [mdl] = mdlEstFct(data)
%   [err] = errFct(mdl,data)
% 
    it = ceil(log(1-cProb)/log(1-(1-oRelFr)^mn));
    n = size(data,1);

    err = Inf;
    bMdl = [];
    mPts = [];
    
    for i=1:it
        rPerm = randperm(n);
        mData = data(rPerm(1:mn),:);
        
        cMdl = mdlEstFct(mData);
        lErr = errFct(cMdl,data);
        
        if(lErr < err)
            err = lErr;
            mPts = mData;
            bMdl = cMdl;
        end
    end
end

function [err] = lineerror(mdl,pts)
% Estimates the error for a line fitted through the data points pts using 
% model parameters mdl
    thr = 0.2;

    m = mdl(1);
    t = mdl(2);
    n = size(pts,1);
    
    x = [-min(pts(:,1)-5),max(pts(:,1))+5];
    y = m*x + t;
    cPt1 = [x(1);y(1)];
    cPt2 = [x(2);y(2)];
    r = cPt2 - cPt1;    %step 1
    
    nn = [-r(2);r(1)];
    nn = nn/norm(nn);   %step 2
    
    d = cPt1' * nn;     %step 3
    
    ds = abs(pts * nn - d);
    err = sum(ds > thr)/n;%step 4
end

function [mdl] = fitline(pts)
% Fits a line throug pts using least squares and returns line parameters
% [m t] in model vector mdl.

    % Build components for eq:
    %
    %         |m|                      
    %  M    * |t| =  ptsY
    %
    % and solve it to get the line eq for l
    
    x = pts(:,1);
    y = pts(:,2);
%     n = size(pts,1);
%     M = ones(n,2);
%     M(:,1) = x;
    M = [x,ones(length(x),1)]; 
    Mpinv = pinv(M);
    mdl = Mpinv * y;
    
end