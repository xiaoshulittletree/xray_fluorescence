for i=1:81 
    content =fileread(['BMXStarget_',num2str(i),'.mac']);
    content=regexprep(content,'250 eV','1 keV');
    content=regexprep(content,'/run/setCut  0.001 mm', '/run/setCut  10 mm');
    fid=fopen(['BMXStarget_',num2str(i),'.mac'],'w');
    fwrite(fid, content);
    fclose(fid);
    clear content
end
    
