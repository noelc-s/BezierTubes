
if overlap
    buffer = 0.0;
else
    buffer = 0.01;
end

fileID = fopen('poly.txt','w');
fprintf(fileID, "polyhedra = (\n");


str = "";
for i = 1:size(PolyInt,2)
    V = lcon2vert(PolyInt{i}(:,1:end-1),PolyInt{i}(:,end)+buffer);
    str1 = "Polyhedron.from_vertices((";
    for j = 1:size(V,1)
        str1 = append(str1,sprintf("[%i,%i,%i,%i],",V(j,1),V(j,2),V(j,3),V(j,4)));
    end
    str1 = append(str1,")),");
    str = append(str,str1);
end
fprintf(fileID, str);
fprintf(fileID, "\n)");


e_str = "'s': (";
for i = 1:length(SI)
    e_str = e_str + "'p"+num2str(SI(i)-1)+"',";
end
e_str = e_str + "),";
target_added=false;
for i = 1:dual_G.Edges.EndNodes(end,1)
    ind = find(dual_G.Edges.EndNodes(:,1)==i);
    str1 = "'p"+num2str(i-1)+"': (";
    for j = 1:length(ind)
        str1 = str1+"'p"+num2str(dual_G.Edges.EndNodes(ind(j),2)-1)+"',";
    end
    if any(EI==i)
        target_added = true;
        str1=str1+"'t',";
    end
    str1 = str1+")";
    e_str = e_str+str1+", ";
end

if ~target_added
    e_str = e_str+"'p"+num2str(EI-1)+"': ('t',)";
end
fprintf(fileID, "\nedges = {\n");
fprintf(fileID, e_str);
fprintf(fileID, "\n}");


vel_scaling = 0.01;
Cost = [];
for m = 1:4
I_m = zeros(2,8);
I_m(1,(m-1)*2+1) = 1;
I_m(2,(m-1)*2+2) = 1;
Cost = [Cost;[I_m*kron(eye(4),eye(2))'; I_m*kron(H,eye(2))']*kron(D_nT,eye(2)) - 0.999*[zeros(4) eye(4)]]; % weird condition so B'*B has nonnegative eigenvalues...
end
Cost(:,[3 4 7 8]) = Cost(:,[3 4 7 8])*vel_scaling;
fprintf(fileID, "\n\nH = np.array([");
for i= 1:size(Cost,1)
    fprintf(fileID,"[")
fprintf(fileID, "%f,",Cost(i,:));
fprintf(fileID,"],\n")
end
fprintf(fileID,"])")


fclose(fileID);