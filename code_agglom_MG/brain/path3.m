function P = path3(folder)
rootpath = fileparts(which('main.m'));
P = [rootpath,'/brain/',folder,'/'];
end