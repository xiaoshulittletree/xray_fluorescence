
folder='/home/shuzhang/GEANT4/XrayFluo/'

TRM=zeros(81,150,36);
edges=0.001*10.^((0:150)*5/150);
electronenergy=10*10.^((0:80)*4/80); %keV
thetaedges=(0:5:180)/180*pi;
sr=2*pi*(cos(thetaedges(1:36))-cos(thetaedges(2:37)));
for i=1:81
    filename=[folder,num2str(electronenergy(i),'%.0f'),'keV_nt_gamma.csv']
    [EventID,Energy,Theta,Phi,Rho]=importfile(filename);
     
    anglebin=discretize(Theta,thetaedges);
    for j=1:36
        TRM(i,:,j)=histcounts(Energy(anglebin==j),edges)/EventID(end)/sr(j);
    end
end
