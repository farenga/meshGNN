rootpath = fileparts(which('main.m'));
FolderName = [rootpath,'/output_files/images'];   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = get(FigHandle, 'Name');
    saveas(FigHandle, fullfile(FolderName, [FigName,'.png']));
end
