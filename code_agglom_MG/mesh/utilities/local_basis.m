function vert = local_basis(mesh)
    dim = manifold(mesh);
    vert = vertices(mesh);
    %%%
    if mesh.dim > dim
        [~, score] = pca(vert);
        vert = score(:,1:dim);
    end
    %%%
%     if mesh.dim > dim
%         vert = vert-vert(1,:);
%         coeff = rand(mesh.elem_num(1),dim);
%         basis = vert'*coeff;
%         while rank(basis) < dim
%             coeff = randn(mesh.elem_num(1),dim);
%             basis = vert'*coeff;
%         end
%         vert = vert*basis;
%     end
end