function [I, tree] = relabel(tree)

w=[];
for i=1:length(tree)
     w(:,i)=tree{i}.clustering;
end;

for i=2:length(tree)-1
     h=w(:,i);
     [~,b]=sort(h);
     w=w(b,:);
 end;
I=w(:,1);

for i=2:length(tree)-1
    w(:,i)=1+[0;cumsum(abs(diff(w(:,i)))>0)];
end;

[~, b]=sort(I);
w=w(b,:);


for i=1:length(tree)-1
    
    tree{i}.clustering=w(:,i)';
    foldersize=[];
    superfolder=[];
    for j=1:tree{i}.folder_count
        foldersize=[foldersize sum(w(:,i)==j)];
        superfolder=[superfolder w(find(w(:,i)==j, 1 ),i+1)]; 
    end;
    tree{i}.folder_sizes=foldersize;
    tree{i}.super_folders=superfolder;
end;