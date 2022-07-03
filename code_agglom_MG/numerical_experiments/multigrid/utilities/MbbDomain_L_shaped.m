%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = MbbDomain_L_shaped(Demand,Arg)
BdBox = [0 2 0 2];
BdBoxA = [0 2 0 2];
BdBoxB = [1 2 1 2];
switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBoxA,BdBoxB);
    case('BC');    x = BndryCnds(Arg,BdBox,BdBoxA,BdBoxB);
    case('BdBox'); x = BdBox;
    case('BdBoxes'); x{1} = BdBoxA; x{2} = BdBoxB;
end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBoxA,BdBoxB)
dA = dRectangle(P,BdBoxA(1),BdBoxA(2),BdBoxA(3),BdBoxA(4));
dB = dRectangle(P,BdBoxB(1),BdBoxB(2),BdBoxB(3),BdBoxB(4));
Dist = dDiff(dA,dB);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,BdBox,BdBoxA,BdBoxB)
eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
LeftEdgeNodes = find(abs(Node(:,1)-BdBox(1))<eps);
LeftUpperNode = find(abs(Node(:,1)-BdBox(1))<eps & ...
    abs(Node(:,2)-BdBox(4))<eps);
RigthBottomNode = find(abs(Node(:,1)-BdBox(2))<eps & ...
    abs(Node(:,2)-BdBox(3))<eps);
FixedNodes = [LeftEdgeNodes; RigthBottomNode];
Supp = zeros(length(FixedNodes),3);
Supp(:,1)=FixedNodes; Supp(1:end-1,2)=1; Supp(end,3)=1;
Load = [LeftUpperNode,0,-1];
x = {Supp,Load};
%-------------------------------------------------------------------------%