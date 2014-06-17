function [A,b] = vert2con(V)
% VERT2CON - convert a set of points to the set of inequality constraints
%            which most tightly contain the points; i.e., create
%            constraints to bound the convex hull of the given points
%
% [A,b,Aeq,Beq] = vert2con(V)
%
% V = a set of points, each ROW of which is one point
% A,b = a set of constraints such that A*x <= b defines
%       the region of space enclosing the convex hull of
%       the given points
%
% For n dimensions:
% V = p x n matrix (p vertices, n dimensions)
% A = m x n matrix (m constraints, n dimensions)
% b = m x 1 vector (m constraints)
%
% NOTES: (1) In higher dimensions, duplicate constraints can
%            appear. This program detects duplicates at up to 6
%            digits of precision, then returns the unique constraints.
%        (2) See companion function CON2VERT.
%        (3) ver 1.0: initial version, June 2005.
%        (4) ver 1.1: enhanced redundancy checks, July 2005
%        (5) Written by Michael Kleder
% 
% EXAMPLE:
%
% V=rand(20,2)*6-2;
% [A,b]=vert2con(V)
% figure('renderer','zbuffer')
% hold on
% plot(V(:,1),V(:,2),'r.')
% [x,y]=ndgrid(-3:.01:5);
% p=[x(:) y(:)]';
% p=(A*p <= repmat(b,[1 length(p)]));
% p = double(all(p));
% p=reshape(p,size(x));
% h=pcolor(x,y,p);
% set(h,'edgecolor','none')
% set(h,'zdata',get(h,'zdata')-1) % keep in back
% axis equal
% set(gca,'color','none')
% title('A*x <= b  (1=True, 0=False)')
% colorbar

% do manual conversion for 2-D
if size(V,1) == 2 && size(V,2) == 2
    p1 = V(1,:);
    p2 = V(2,:);
    slope = (p1(2) - p2(2))/(p1(1) - p2(1));
    if abs(slope) == inf
        A = [1 0];
        b = [p1(1)];
    elseif abs(slope) == 0
        A = [0 1];
        b = [p1(2)];
    else
        A = [-slope 1];
        b = [p1(2) - slope*p1(1)];
    end
else
    k = convhulln(V);
    c = mean(V(unique(k),:));
    V=V-repmat(c,[size(V,1) 1]);
    A  = NaN*zeros(size(k,1),size(V,2));
    rc=0;
    for ix = 1:size(k,1)
        F = V(k(ix,:),:);
        if rank(F,1e-5) == size(F,1)
            rc=rc+1;
            A(rc,:)=F\ones(size(F,1),1);
        end
    end
    A=A(1:rc,:);
    b=ones(size(A,1),1);
    b=b+A*c';
    % eliminate dumplicate constraints:
    [null,I]=unique(num2str([A b],6),'rows'); % note precision here
    A=A(I,:); % rounding is NOT done for actual returned results
    b=b(I);
end
return
