%======================================================================
% Here we build the coarse solver
%======================================================================

function [coarse_solver] = MY_COARSE_SOLVER(femregion,A,R0_T,Nelem)


%disp('-------- BEGIN COMPUTING COARSE SOLVER --------');

%===========================================================
% OUTPUT
%===========================================================
coarse_solver = struct('NSUB',Nelem,...
                     'A0',A,...
                     'R0_T',R0_T,...
                     'coarse_femregion',femregion);

%disp('-------- END COMPUTING COARSE SOLVER -------- ');
end
