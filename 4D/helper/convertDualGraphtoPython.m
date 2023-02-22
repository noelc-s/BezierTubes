
buffer = 0.01;

fileID = fopen('poly.txt','w');
fprintf(fileID, "polyhedra = (\n");


str = "";
for i = 1:size(PolyInt,2)
    V = lcon2vert(PolyInt{i}(:,1:end-1),PolyInt{i}(:,end)+buffer);
    str1 = "Polyhedron.from_vertices((";
    for j = 1:size(V,1)
        str1 = append(str1,sprintf("[%i,%i],",V(j,1),V(j,2)));
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


fclose(fileID);