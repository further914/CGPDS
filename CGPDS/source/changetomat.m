function changetomat
options.dirName='../missa/gray/'; options.from=0; options.to=149;
width=360;
height=288;
Y=preprocessVideo('missa',width,height,'ras',options);

playMov(height, width,0.1, Y);
end