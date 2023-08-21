%==========================================================================
%                      Build mesh for time test
%==========================================================================
%clc;
clear all;
close all;

%range_Ne = [8 16 32 64 128 256 512 1024 2048 4096];
range_Ne = [8 32];

index = 1;
for Ne = range_Ne
    
    X = ['  Ne = ',num2str(Ne),'.'];
    disp(X);
    
    Tau{index} = generate_mesh(@MbbDomain,Ne,100);
    neighbour{index} = neighbours(Tau{index});
    
    [error] = Plot_PolyMesh(Tau{index},index,'k',neighbour{index},'k');
    
    index = index + 1;
end

%filename = 'Tau_region.mat';
%save(filename,'Tau','neighbour');
