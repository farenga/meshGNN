clc;
clear all;
close all;

% Model problem:
%   - \delta u = f    on \Omega = (0,1)^2
%   u = 0     on \partial \Omega
%
% let f to be the forcing term such that uex = (x.^2-x).*(y.^2-y).*exp(x)
% is the exact solution of the problem
%
% we solve using SIPG method, with Voronoi mesh, and verify the estimate
% ||uh - uex||_{L2} \le C h^{p+1} ||uex||_{L2},
% where uh is the numerical solution given by the SIPG method, and h is the
% diameter of the mesh.
%
% Here the stiffness matrix is computed using the Sukumar method for the
% exact integration of homogeneous function.


range_deg_approx = [1 2 3 4];
range_elements = [10 40 160 640 2560];

for p = range_deg_approx
    
    disp(['Degree of approximation p = ',num2str(p)]);
    EE = [];
    Dati = dati;
    Dati.fem = p;
    Dati.nqn = 2*Dati.fem + 1;
    
    for N = range_elements
        
        disp(['  Number of elements: ',num2str(N)]);
        
        % Polygonal SIPG method, assemble matrix
        [region] = generate_mesh(@MbbDomain,N,100); % generate a Voronoi mesh of N elements on \Omega
        [neighbour] = neighbours(region);  % struct describing the 
        [femregion] = create_dof(Dati,region);
        
        % get the coefficient lists of the Legendre 1D polynomials up to degree deg, 
        % and their derivatives
        [P,dP] = get_coeff_Legendre1D(p);
        
        % Compute matrix using Sukumar method for homogeneous functions
        [Matrices]= SukumarMatrix2D(femregion,neighbour,P,dP,Dati);
        
        % here we compute the solution and evaluate the error 
        u_h = Matrices.A\Matrices.f;
        [E_L2]= compute_errors(Dati,femregion,u_h);
        EE = [EE E_L2];
        disp(['  => ||uex - uh||_{L2} = ',num2str(EE(end))]);
        
    end
    
    ErrorL2{p} = EE;
    
    % reference slope using the degree of approximation p
    e4slope1(1) = EE(1);
    for i = 2:length(range_elements)
        e4slope1(i) = e4slope1(i-1)/(2^(p+1));
    end
    e4slope{p} = e4slope1;
    
end


% PLOT ERRORS
for p = range_deg_approx
    
    loglog(range_elements,ErrorL2{p},'b*-','linewidth',2 );
    grid on;
    hold on;
    loglog(range_elements,e4slope{p},'r.-','linewidth',2 );
    labx = xlabel('log (Number elements)');
    laby = ylabel('log ||uex - uh||_{L2}');
    titolo = title('Error Norm L2');
    
    string4legend = ['Error L2 with p = ',num2str(p)];
    AX = legend(string4legend,'Reference slope','Location','SouthWest');
    %set(labx,'FontSize',20);
    %set(laby,'FontSize',20);
    %set(titolo,'FontSize',20);
    
end


