function writeJsonAtomic(fileName,value)
%WRITEJSONATOMIC Replace a JSON file only after the new content is written.

fileName = char(fileName);
temporary = [fileName '.tmp'];
text = jsonencode(tdgl.io.jsonSafe(value),'PrettyPrint',true);
fileId = fopen(temporary,'w','n','UTF-8');
if fileId < 0
    error('tdgl:io:CannotWriteJson','Cannot open %s for writing.',temporary);
end
cleanup = onCleanup(@() fclose(fileId));
fprintf(fileId,'%s\n',text);
clear cleanup;
[success,message] = movefile(temporary,fileName,'f');
if ~success
    error('tdgl:io:CannotCommitJson','Cannot replace %s: %s',fileName,message);
end
end
