function savex(x,filename)
if nargin<2
    filename=sprintf('./Data/%s.dat',inputname(1));
end
f=fopen(filename,'w');
fprintf(f,'%s\n',inputname(1));
for i=1:numel(x)
    fprintf(f,'%8.5f\n',x(i));
end
fclose(f);
end