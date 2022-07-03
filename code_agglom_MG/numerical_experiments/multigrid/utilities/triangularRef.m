function [Tau] = triangularRef(Tau_coarse)

addpath(genpath('/Users/mac/Desktop/Multilevel_NonNested_QuadratureFree_MOX30/mesh2d-master'));

nodes = Tau_coarse.coord;

edges = [];
for ie = 1:Tau_coarse.ne
    edges = [edges; [Tau_coarse.edges_phys{ie}, Tau_coarse.tag(ie)*ones(size(Tau_coarse.edges_phys{ie},1),1)] ];
end

part{1} = edges;
for k = Tau_coarse.tag
    part{k+1} = [find(edges(:,3)==k)];
end

hmax = +0.045 ;
[vlfs,tlfs, hlfs] = lfshfn2(nodes,edges,part) ;
hlfs = min(hmax,hlfs) ;
[slfs] = idxtri2(vlfs,tlfs) ;


[vert,etri,tria,tnum] = refine2(nodes,edges,part,[],hfun, ...
                         vlfs,tlfs,slfs,hlfs);