%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostic Medical Image Processing (DMIP) 
% WS 2014/15
% Author: Marco Boegel, Yan Xia
% Exercise: Image Undistortion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data_scaling()

% Copy and paste your exercise 2 code here

% Find the proper scale via optimizing the condition number of the Gramian
% matrix

% Initial scales
x0 = [1 1];
% Start optimization
xs = fminsearch(@optFunction, x0); 
disp('Optimized scale');
xs
%Save old measurement matrix
AT = A;

% Build measurement matrix using optimized scales
for r = 1:NumCorresp
  c = 1;
  for i = 0:d
    for j = 0:(d-i)
      A(r,c) = (xs(1)*XU2vec(r))^i * (xs(2)*YU2vec(r))^j;
      c = c + 1;
    end
  end
end

% Compare condition numbers of scaled and unscaled measurement matrix
disp('cond(A^{T}A) optimized');
cond(A'*A)
disp('cond(A^{T}A) before');
cond(AT'*AT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimization function cond(A' * A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = optFunction(x)
    E = zeros( NumCorresp, NumKoeff );
    for r = 1:NumCorresp
        c = 1;
        for i = 0:d
            for j = 0:(d-i)
                E(r, c) = (x(1)*XU2vec(r))^i * (x(2)*YU2vec(r))^j;
                c = c + 1;
            end
        end
    end
    f = cond(E'*E);
end % optimization function

end