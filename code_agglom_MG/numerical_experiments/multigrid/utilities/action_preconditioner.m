%============================================================
% compute the action of the preconditioner P
% on the rhb b in the linear system:
% P A x= P b
% with P=sum_{i=0...NS}{Pi}
%============================================================

function [z] = action_preconditioner(b, local_solvers, coarse_solver, type_preconditioner, A)


NSUB = coarse_solver.NSUB;
ndof = length(b);

z = sparse(ndof,1);


switch type_preconditioner
    
    case 'ADDITIVE'
        
        %----------------------------------------------------------------
        % ADDITIVE PRECONDITIONER
        %----------------------------------------------------------------
        
        % COARSE SOLVER
        A0 = coarse_solver.A0;
        R0_T = coarse_solver.R0_T;
        b0 = (R0_T')*b;
        z0 = A0\b0;
        z = z + R0_T*z0;
        
        % LOCAL SOLVER
        for i = 1:NSUB
            Ai = local_solvers(i).Ai;
            index = local_solvers(i).index_dof;
            bi = b(index);
            zi = Ai\bi;
            z(index) = z(index) + zi;
        end
        
    case 'MULTIPLICATIVE'
        %----------------------------------------------------------------
        % MULTIPLICATIVE PRECONDITIONER
        %----------------------------------------------------------------
        
        % COARSE SOLVER
        A0 = coarse_solver.A0;
        R0_T = coarse_solver.R0_T;
        z = z +R0_T*(A0\((R0_T')*b));
        
        % LOCAL SOLVER
        for i = 1:NSUB
            Ai = local_solvers(i).Ai;
            index = local_solvers(i).index_dof;
            ndof_i = length(local_solvers(i).index_dof);
            Ri_T = local_solvers(i).Ri_T;
            z = z + Ri_T*(Ai\((Ri_T')*(b-A*z)));
        end
        
        
    case 'SYMMETRIZED'
        %----------------------------------------------------------------
        % SYMMETRIZED PRECONDITIONER
        %----------------------------------------------------------------
        
        % COARSE SOLVER
        A0 = coarse_solver.A0;
        R0_T = coarse_solver.R0_T;
        z = z +R0_T*(A0\((R0_T')*b));
        
        % LOCAL SOLVER
        for i = 1:NSUB
            Ai = local_solvers(i).Ai;
            index = local_solvers(i).index_dof;
            ndof_i = length(local_solvers(i).index_dof);
            Ri_T = local_solvers(i).Ri_T;
            z = z + Ri_T*(Ai\((Ri_T')*(b-A*z)));
        end
        
        % SYMMETRIZED LOCAL SOLVER
        for i = NSUB:-1:1
            Ai_trans = (local_solvers(i).Ai)';
            index = local_solvers(i).index_dof;
            ndof_i = length(local_solvers(i).index_dof);
            Ri_T = local_solvers(i).Ri_T;
            z = z + Ri_T*(Ai_trans\((Ri_T')*(b-A*z)));
        end
        % SYMMETRIZED COARSE SOLVER
        A0 = coarse_solver.A0;
        R0_T = coarse_solver.R0_T;
        z = z +R0_T*(A0\((R0_T')*b));
        disp('CHECKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK')
        
        
end
