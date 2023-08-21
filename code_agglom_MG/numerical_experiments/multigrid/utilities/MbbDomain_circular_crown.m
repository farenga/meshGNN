%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = MbbDomain_circular_crown(Demand,Arg)
  %BdBox = [0 1 0 1];
  circx=0.5; circr=2;
  BdBox = [circx-circr,circx+circr,-circr,circr];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg,BdBox);
    case('BdBox'); x = BdBox;
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  circx=0; circr=2;
  d1 = dcircle(P,circx,0,circr);
  d2 = dcircle(P,circx,0,circr/3);
  %d2 = dRectangle(P,-0.5,0.5,-0.5,0.5);
  Dist = ddiff(d1,d2);
  %Dist = d1;
  %Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
  %Dist = dCircle(P,(BdBox(1)+BdBox(2))/2,(BdBox(3)+BdBox(4))/2, (BdBox(4)-BdBox(3))/2);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,BdBox)
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