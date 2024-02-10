%
clear all
xyz=dlmread('TorusSurfDirK1K220171201.xyz');
plot3( xyz(:,1),xyz(:,2),xyz(:,3),'r.');
hold on
axis equal
len = 0.3;
Pos=xyz(:,1:3);
Dir1=xyz(:,4:6);
DirA = Pos - len*Dir1;
DirB = Pos+ len*Dir1;
Lx=[DirA(:,1), DirB(:,1)]';
Ly=[DirA(:,2), DirB(:,2)]';
Lz=[DirA(:,3), DirB(:,3)]';
plot3(Lx,Ly,Lz,'b-')