function rho = get_diffusion_coefficient(femregion,Dati,rho_1,rho_2)

rho = zeros(femregion.ne,1);

switch Dati.permeability
    
    case 'EvenAndOdd'
        
        for ie = 1:femregion.ne
            if( mod(ie,2) == 0 )
                rho(ie) = rho_1;
            else
                rho(ie) = rho_2;
            end
        end
        
    case 'BlackAzure'
        
        rho_Azure = rho_1;
        rho_Black = rho_2;
        for ie = 1:femregion.ne
            
            BBox = femregion.BBox(ie,:);
            x1B = BBox(1); x2B = BBox(2);
            y1B = BBox(3); y2B = BBox(4);
            
            if (x1B <= 0.5 && y1B < 0.5) || (x2B > 0.5 && y2B > 0.5)
                % Black
                rho(ie) = rho_Azure;
            else
                % White
                rho(ie) = rho_Black;
            end
            
        end
        
    case 'PolygonalChessboard'
        
        ChessRegion = Dati.DiffusionCoefficientRegion;
        ChessData = Dati;
        ChessData.permeability = 'EvenAndOdd';
        ChessRegion.rho = get_diffusion_coefficient(ChessRegion,ChessData,rho_1,rho_2);
        for ie = 1:femregion.ne
            
            CoordElement = femregion.coords_element{ie};
            CenterMass(1) = sum(CoordElement(:,1))/size(CoordElement,1); CenterMass(2) = sum(CoordElement(:,2))/size(CoordElement,1);
            for il = 1:ChessRegion.ne
                
                ChessCoordElement = ChessRegion.coords_element{il};
                if inpolygon(CenterMass(1),CenterMass(2),ChessCoordElement(:,1),ChessCoordElement(:,2))
                    rho(ie) = ChessRegion.rho(il);
                end
                
            end
            
        end
        
end