function OutputModule2TXT(Comodule,FeatureType,ResultsFile)
% Output the co-modules to a text file named [ResultsFile '.txt'], and then
% user will obtain a text file including (2*3) columns, the first
% three columns are the numbers of selected samples, N types of features. 
% And the next three columns are the index lists for the selected
% components.
%
% INPUT :
%
% Comodule: (K x (N+1)) cell array. Each row represents an identified
%            md-module.
% FeatureType: (1 x N) cell array storing the feature types. 
% ResultsFile: A string for the saved file name.
%
[K,N] = size(Comodule);
fid = fopen([ResultsFile '.txt'],'w+');
if(fid == -1)
    error('Fail in opening text file.')
end

OutStr = 'SampleIndex';
OutNum = '#Sample';
for i = 1:(N-1)
    OutStr = [OutStr ',' FeatureType{i} 'Index'];
    OutNum = [OutNum ',#' FeatureType{i}];
end
OutStr = [OutNum ',' OutStr];
OutStr = OutStr(1:end);
fprintf(fid,'%s\n',OutStr);
clear i OutStr OutNum

for k = 1:K
    OutStr = '';
    OutNum = '';
    for i = 1:N
        OutStr = [OutStr ','  ObtainStr(Comodule{k,i})];
        OutNum = [OutNum ',' num2str(length(Comodule{k,i}))];
    end
    OutStr = [OutNum  OutStr];
    OutStr = OutStr(2:length(OutStr));
    fprintf(fid,'%s\n',OutStr);
end
fclose(fid);
end

function str = ObtainStr(mat)
[m,n] = size(mat);
if (m == 1 && n == 1),
    str = ['[' mat2str(mat) ']'];
elseif(isempty(mat))
    str = '[]';
else
    str = mat2str(mat);
end
end
