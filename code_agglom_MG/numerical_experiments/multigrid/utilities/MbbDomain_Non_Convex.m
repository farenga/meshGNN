%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = MbbDomain_Non_Convex(Demand,Arg)
BdBox = [-1 1 0 2];
BdBoxA = [-1 1 0 2];
BdBoxA = [0 1 0 1];
BdBoxB = [-0.5 0.5 0 0.5];
BdBoxC = [-0.5 0.5 1.5 2];
switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBoxA,BdBoxB);
    case('BC');    x = BndryCnds(Arg,BdBox,BdBoxA,BdBoxB);
    case('BdBox'); x = BdBox;
end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBoxA,BdBoxB)
dA = dRectangle(P,BdBoxA(1),BdBoxA(2),BdBoxA(3),BdBoxA(4));
%dB = dRectangle(P,BdBoxB(1),BdBoxB(2),BdBoxB(3),BdBoxB(4));
dB = dCircle(P,0,0,0.5);
%dC = dLine(P,0,0,1,0);
dC = dCircle(P,0,2,0.5);
%Dist = dIntersect(dC,ddiff(dA,dB));
%Dist = dDiff(dDiff(dA,dB),dC);
Dist = dIntersect(dA,dDiff(dA,dB));
d1 = dLine(P,0,0,1,0);
d2 = dLine(P,1,0,1.2,1);
d3 = dLine(P,1.2,1,0.8,2);
d4 = dLine(P,0.8,2,0,1);
d5 = dLine(P,0,1,0,0);
Dist = dIntersect(d1,dIntersect(d2,dIntersect(d3,dIntersect(d4,d5))));
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