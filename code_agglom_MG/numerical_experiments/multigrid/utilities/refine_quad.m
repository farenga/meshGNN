function [pn, tn]= refine_quad(p, t)

xnod = p(1,:);
ynod = p(2,:);

nodes = t;

nele=length(nodes(1,:));
nno=length(xnod);
nodesn=[];
bnodn=[];
sum=ones(1,4)/4.;
check=sparse(1,1,1,nno,nno,10*nno);
check(1,1)=0;

for iel=1:nele
    
    %     new nodes
    
    iv=nodes(1:end-1,iel);
    xnod=[xnod, sum*xnod(iv)'];
    ynod=[ynod, sum*ynod(iv)'];
    nno=nno+1; nodmid=nno;
    if check(iv(1),iv(2))==0
        nno=nno+1;
        xnod=[xnod, (xnod(iv(1))+xnod(iv(2)))/2.];
        ynod=[ynod, (ynod(iv(1))+ynod(iv(2)))/2.];
        check(iv(1),iv(2))=nno;check(iv(2),iv(1))=nno;
    end
    
    if check(iv(3),iv(2))==0
        nno=nno+1;
        xnod=[xnod, (xnod(iv(3))+xnod(iv(2)))/2.];
        ynod=[ynod, (ynod(iv(3))+ynod(iv(2)))/2.];
        check(iv(3),iv(2))=nno;check(iv(2),iv(3))=nno;
    end
    
    if check(iv(3),iv(4))==0
        nno=nno+1;
        xnod=[xnod, (xnod(iv(3))+xnod(iv(4)))/2.];
        ynod=[ynod, (ynod(iv(3))+ynod(iv(4)))/2.];
        check(iv(3),iv(4))=nno;check(iv(4),iv(3))=nno;
    end
    if check(iv(1),iv(4))==0
        nno=nno+1;
        xnod=[xnod, (xnod(iv(1))+xnod(iv(4)))/2.];
        ynod=[ynod, (ynod(iv(1))+ynod(iv(4)))/2.];
        check(iv(1),iv(4))=nno;check(iv(4),iv(1))=nno;
    end
    ivn=[iv(1) check(iv(1),iv(2)) nodmid check(iv(1),iv(4))];
    nodesn=[nodesn, ivn'];
    ivn=[check(iv(1),iv(2)) iv(2) check(iv(2),iv(3)) nodmid];
    nodesn=[nodesn, ivn'];
    ivn=[nodmid check(iv(3),iv(2)) iv(3) check(iv(3),iv(4))];
    nodesn=[nodesn, ivn'];
    ivn=[check(iv(1),iv(4)) nodmid check(iv(3),iv(4)) iv(4)];
    nodesn=[nodesn, ivn'];
end

pn(1,:)=xnod; pn(2,:)=ynod;

tn = [nodesn; ones(1,length(nodesn))];

