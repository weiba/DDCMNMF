function C = exp_colormap(colors, N)
% A simple function that returns a colormap, C, for visualizing input data.
% C is just a N x 3 matrix [R G B] describing the range 
% of color values.
%
% Example:
% >> C = exp_colormap('blue-yellow',64);
% >> colormap(C);
% 
% INPUT:
%
% colors: Options for it are: 'blue-yellow'; 'green-red'; 'yellow'; 
%                             'blue-white-red'; 'default'.
% N: The number of degrees of color to use.  The default is 64.
%
% OUTPUT:
%
% C: It returns a [N x 3] RGB colormap.
% 
% The colormaps returned range monotonically.
%  
  
if nargin < 2 
  N = 64;
end

if nargin < 1
  colors = 'green-red';
end

X = 0.5: -1/(N-1):-0.5;
X = abs(X).*2;

switch colors
case {'green-red'}
  R = [X(:, 1:N/2) zeros(1,N/2)];
  G = [zeros(1, N/2) X(:,(N/2 + 1):N)];
  B = zeros(1,N);

case {'blue-yellow'} 
  R = [zeros(1,N/2) X(:,(N/2 + 1):N)];
  B = [X(:,1:N/2) zeros(1,N/2)];
  G = [zeros(1,N/2) X(:,(N/2 + 1):N)];

case {'yellow'} 
  X = 0:1/(N - 1):1;
  R = X;
  B = zeros(1,N);
  G = X;

case {'blue-white-red'}
  R = [X(:,(N/2 + 1):N) ones(1,N/2) ];
  G = [X(:,(N/2 + 1):N) X(:, 1:N/2) ];
  B = [ones(1,N/2) X(:, 1:N/2) ];
  
case {'default'}
  C = colormap(summer);
  return;
  
otherwise
 error([colors ' is not a known option for coloring.']);
end

C = [R' G' B'];
end
