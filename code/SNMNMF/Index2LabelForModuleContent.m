function Index2LabelForModuleContent(ModuleIndex,FeatureLabel,TypeName,ResultsFile)
%
% This function is used for converting the indexes into labels, and output 
% the lists into text files.
%
% INPUT :
%
% ModuleIndex: (K x 1) cell array, ModuleIndex{i} records the indexes of
%               one type of feature in one module. 
% FeatureLabel: (n x 1) cell array, all the feature names.
% TypeName: String for the component name, such as, 'Sample','Gene','miRNA'
% ResultsFile: The folder name to save these text files.
%

ListFile = [ResultsFile '/' TypeName 'Lists'];
if ~isdir(ListFile)
    mkdir(ListFile);
end

NCluster = size(ModuleIndex,1);

for k = 1:NCluster
    TxtName = [TypeName 'List_' num2str(k) '.txt'];
    ind = ModuleIndex{k,1};
    fid = fopen([ListFile '/' TxtName],'w+');
    if(fid ~= -1)
        len = length(ind);
        str = ['# The List concludes ' num2str(len) ' ' TypeName '(s).'];
        fprintf(fid,'%s\t\n',str);
        fprintf(fid,'%s\t\n',[TypeName 'Label']);
        for i = 1:len
            str = FeatureLabel{ind(i)};
            fprintf(fid,'%s\t\n',str);
        end
        fclose(fid);
        clear len str ind TxtName fid
    else
        error(['Something wrong with writing' TypeName 'list.'])
    end
end
clear k ListFile
end

