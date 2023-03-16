clear pts
G = graph;
G = G.addnode(1);
G.Nodes.x = [0 0 0 0];

for i = 1:size(Polytopes,2)
    %     if length(Polytopes{i})~=0
    %         [A,b] = vert2lcon(Polytopes{i});
    %
    %         pts{i} = [A b];
    %     else
    %         pts{i} = [];
    %     end
    pts{i} = Polytopes{i};
    G = G.addnode(1);
    G.Nodes.x(end,:) = PolyCenter{i};
end
G = rmnode(G,1);


s_dim = size(pts{1},2)-1;
opts = optimset('Display','off');
hold on

PolyInt = {};
ind = 1;

for i = 1:size(pts,2)
    for j = (i+1):size(pts,2)
        if (length(Polytopes{i})~=0 && length(Polytopes{j})~=0)
            touching = false;
%             [x,fval,ef] = linprog(zeros(1,s_dim),[pts{i}(:,1:s_dim); pts{j}(:,1:s_dim)],[pts{i}(:,end); pts{j}(:,end)]+0.01,[],[],[],[],[],opts);

objective=[0 0];
if overlap
     IN=struct('obj',objective,'A',[pts{i}(:,1:s_dim); pts{j}(:,1:s_dim)],'B',[pts{i}(:,end); pts{j}(:,end)]);
else
        IN=struct('obj',objective,'A',[pts{i}(:,1:s_dim); pts{j}(:,1:s_dim)],'B',[pts{i}(:,end); pts{j}(:,end)]+0.01);
end
        OUT = cddmex('solve_lp',IN);
        if OUT.how==1
                touching = true;
            end
            if touching == true
                tmp = [pts{i}; pts{j}];
                %                 V = lcon2vert(tmp(:,1:end-1),tmp(:,end)+0.01);
                %                 if (size(V,1)==1)
                %                     break
                %                 end
                PolyInt{ind} = tmp;
                ind = ind+1;
                G = G.addedge(i,j,norm(PolyCenter{i} - PolyCenter{j})); % weight is euclidean distance of centers
                %                 G = G.addedge(j,i,norm(PolyCenter{i} - PolyCenter{j})); % weight is euclidean distance of centers
                %                 sprintf("%i is touching %i", i, j)
                %                 plot([PolyCenter{i}(1) PolyCenter{j}(1)], [PolyCenter{i}(2) PolyCenter{j}(2)],'b','linewidth',1)

            end
        end
    end
end

dual_G = digraph;
dual_G = dual_G.addnode(1);
dual_G.Nodes.P = {0};
dual_G.Nodes.C = [0 0 0 0];
% Dualize the graph
for j = 1:size(G.Edges,1)
    dual_G = dual_G.addnode(1);
    dual_G.Nodes.P{end} = PolyInt{j};
    if overlap
V = lcon2vert(PolyInt{j}(:,1:end-1),PolyInt{j}(:,end));
    else
    V = lcon2vert(PolyInt{j}(:,1:end-1),PolyInt{j}(:,end)+0.01);
    end
    dual_G.Nodes.C(end,:) = mean(V);
end
dual_G = rmnode(dual_G,1);
for i = 1:size(G.Edges,1)
    start_N = G.Edges.EndNodes(i,1);
    end_N = G.Edges.EndNodes(i,2);
    for j = (i+1):size(G.Edges,1)
        if G.Edges.EndNodes(j,1) == start_N || G.Edges.EndNodes(j,2) == start_N || G.Edges.EndNodes(j,1) == end_N ||  G.Edges.EndNodes(j,2) == end_N
            dual_G = dual_G.addedge(i,j);
            dual_G = dual_G.addedge(j,i);
            %             plot([dual_G.Nodes.C(i,1) dual_G.Nodes.C(j,1)], [dual_G.Nodes.C(i,2) dual_G.Nodes.C(j,2)],'r','linewidth',1)
        end
    end
end


for i = 1:ind-1
    if overlap
V = lcon2vert(PolyInt{i}(:,1:end-1), PolyInt{i}(:,end));
    else
    V = lcon2vert(PolyInt{i}(:,1:end-1), PolyInt{i}(:,end)+0.01);
    end
    c_h = convhull(V(:,1:2));
    c_h = V(c_h,:);
    patch(c_h(:,1),c_h(:,2),'k','facealpha',0.01)
    %     scatter(dual_G.Nodes.C(i,1), dual_G.Nodes.C(i,2), 100,'r','filled')
    %     line([V(1,1) V(2,1)],[V(1,2) V(2,2)],'b','facealpha',0.05)
end

if overlap
    buffer = 0.0;
    PC = cell2mat(PolyCenter)';
    [~,SI] = min(vecnorm(IC'-PC'));
    [~,EI] = min(vecnorm(EC'-PC'));
    SI = [];
    EI = [];
    for i=1:size(PolyInt,2)
        if PolyInt{i}(:,1:end-1)*IC' <= PolyInt{i}(:,end)
            SI = [SI i];
        end
        if PolyInt{i}(:,1:end-1)*EC' <= PolyInt{i}(:,end)
            EI = [EI i];
        end
    end

else
    buffer = 0.01;
    for i=1:size(Polytopes,2)
        if Polytopes{i}(:,1:end-1)*IC' <= Polytopes{i}(:,end)
            SI = i;
        end
        if Polytopes{i}(:,1:end-1)*EC' <= Polytopes{i}(:,end)
            EI = i;
        end
    end

    ind_to = find(G.Edges.EndNodes(:,1) == SI);
    ind_from = find(G.Edges.EndNodes(:,2) == SI);
    % SI_to = dual_G.Edges.EndNodes(ind_to,2);
    % SI_from = dual_G.Edges.EndNodes(ind_from,1);
    SI = [ind_to; ind_from];
    ind_to = find(G.Edges.EndNodes(:,1) == EI);
    ind_from = find(G.Edges.EndNodes(:,2) == EI);
    % EI_to = dual_G.Edges.EndNodes(ind_to,2);
    % EI_from = dual_G.Edges.EndNodes(ind_from,1);
    EI = [ind_to; ind_from];
end
drawnow

