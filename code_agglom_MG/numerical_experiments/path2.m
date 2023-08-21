function P = path2(folder)
rootpath = fileparts(which('main.m'));
P = [rootpath,'/output_files/',folder,'/'];
end