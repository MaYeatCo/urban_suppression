%setting the parameters
N=3;
D3=2.5;%1.9928608732968183 - this separates between super (below) and sub (above)
D2=D3;
ubar=0.5;
alpha=0.05;
D1=(1/32).*D3.^(-2).*ubar.^(-1).*(16.*alpha.*D3+(-8).*alpha.*D3.^2+ ...
  2.*alpha.*ubar+(-16).*D3.*ubar+(-1).*alpha.*D3.*ubar+(-2).* ...
  ubar.^2+D3.*ubar.^2)+(1/32).*(D3.^(-4).*ubar.^(-2).*(512.* ...
  alpha.^2.*D3.^2+(-256).*alpha.^2.*D3.^3+64.*alpha.^2.*D3.*ubar+( ...
  -1024).*alpha.*D3.^2.*ubar+(-64).*alpha.^2.*D3.^2.*ubar+256.* ...
  alpha.*D3.^3.*ubar+16.*alpha.^2.*D3.^3.*ubar+4.*alpha.^2.*ubar.^2+ ...
  (-128).*alpha.*D3.*ubar.^2+(-4).*alpha.^2.*D3.*ubar.^2+512.* ...
  D3.^2.*ubar.^2+96.*alpha.*D3.^2.*ubar.^2+alpha.^2.*D3.^2.*ubar.^2+ ...
  (-16).*alpha.*D3.^3.*ubar.^2+(-8).*alpha.*ubar.^3+64.*D3.*ubar.^3+ ...
  8.*alpha.*D3.*ubar.^3+(-32).*D3.^2.*ubar.^3+(-2).*alpha.*D3.^2.* ...
  ubar.^3+4.*ubar.^4+(-4).*D3.*ubar.^4+D3.^2.*ubar.^4)).^(1/2);
beta = (-1).*D1.^(-1).*ubar.^(-2).*((-2).*alpha.^2+(-1).*alpha.^2.*D2+4.* ...
  alpha.*ubar+alpha.*D2.*ubar+(-8).*alpha.*D1.*D2.*ubar+(-2).* ...
  ubar.^2+8.*D1.*D2.*ubar.^2);

vbar= 1-alpha/ubar;
wbar = ubar*(beta/(ubar-alpha)-1);

K=(1/2).*D1.^(-1).*D2.^(-1).*ubar.^(-1).*((-1).*alpha+ubar).^(-1).*( ...
  alpha.^2.*(2+D2)+(-1).*alpha.*(4+D2).*ubar+(2+(-1).*beta.*D1).* ...
  ubar.^2);
chi=D3.*(alpha+beta+(-1).*ubar).^(-1).*(alpha+D1.*D2.*K.^2+(-1).*ubar+ ...
  (-1).*alpha.*beta.*((-1).*alpha+ubar).^(-1));

%these are coefficients for IC
Q=(1/2).*D2.^(-1).*ubar.^(-2).*((-1).*alpha+ubar).^(-1).*((-1).* ...
  alpha.^2.*((-2)+D2)+alpha.*((-4)+D2).*ubar+(2+(-1).*beta.*D1).* ...
  ubar.^2);

R=(1/4).*D1.^(-1).*D2.^(-1).*ubar.^(-2).*((-1).*alpha+ubar).^(-3).*( ...
  alpha.^4.*(2+D2).^2+(-2).*alpha.^3.*(8+6.*D2+D2.^2).*ubar+ ...
  alpha.^2.*(24+2.*beta.*D1.*((-2)+D2)+4.*(3+alpha.*D1).*D2+D2.^2).* ...
  ubar.^2+(-2).*alpha.*(beta.*D1.*((-4)+D2)+2.*(4+D2+3.*alpha.*D1.* ...
  D2)).*ubar.^3+(4+(-4).*beta.*D1+beta.^2.*D1.^2+12.*alpha.*D1.*D2) ...
  .*ubar.^4+(-4).*D1.*D2.*ubar.^5);

%setting time and space
Tmax=5000;
wk=sqrt(K);
a = 0; 
b = 4*pi/wk;  
L = b-a;
eps = -0.01;
P0=0.08;
lowerLeft  = [a   ,a   ];
lowerRight = [b , a  ];
upperRight = [b , b];
upperLeft =  [a , b];


%Setting the geometry
S = [3,4 lowerLeft(1), lowerRight(1), upperRight(1), upperLeft(1), ...
         lowerLeft(2), lowerRight(2), upperRight(2), upperLeft(2)];                     
gdm = S';
ns = 'S';
sf = 'S';
g = decsg(gdm,ns,sf');

%parameters for PDE
d=[1;1;1];
Pf=[N, alpha, beta ];
Pc=[N, D1, D2,D3,chi,eps];
model = createpde(N);
geometryFromEdges(model,g);
generateMesh(model);

%BC for PDE
applyBoundaryCondition(model,"neumann","Edge",1:4);
specifyCoefficients(model,"m",0,"d",d,"c",@(location,state)ccoeffunction(location,state,Pc),"a",0,"f",@(location,state)fcoeffunction(location,state,Pf));

%Initial condition function
icFcn = @(region) [ubar + P0*(cos(region.x*wk)+cos(region.y*wk));
                   vbar + P0*(cos(region.x*wk)+cos(region.y*wk))*Q;
                   wbar + P0*(cos(region.x*wk)+cos(region.y*wk))*R];
setInitialConditions(model,icFcn);

tlist = linspace(0,Tmax,Tmax+1);

%solve the PDE and measure the time
tic
results = solvepde(model,tlist);
toc

%mesh needed for plotting and the correct rms of ODE part
X=model.Mesh.Nodes

%we focus on A equation
u1=reshape(results.NodalSolution(:,1,:),[size(results.NodalSolution(:,1,:),1),size(results.NodalSolution(:,1,:),3)]);
rmsu1=rms(u1-ubar);

%if you need to run longer, you can set repeats
rpts=0;
if rpts>0
    tic
for j=1:rpts
setInitialConditions(model,results);
tlist=linspace(Tmax*j+1,Tmax*(j+1),Tmax);
results=solvepde(model,tlist);
u1=reshape(results.NodalSolution(:,1,:),[size(results.NodalSolution(:,1,:),1),size(results.NodalSolution(:,1,:),3)]);
rmsu1=[rmsu1,rms(u1-ubar)];
end
    toc
end

%plot rms
figure()
plot(linspace(0,Tmax*(rpts+1),Tmax*(rpts+1)+1),rmsu1)

%reaction terms of PDE
function f = fcoeffunction(location,state,P)
N=P(1); alpha=P(2); beta=P(3); 
nr = length(location.x); 
f = zeros(N,nr); 
f(1,:) = -state.u(1,:)+state.u(1,:).*state.u(2,:)+alpha;
f(2,:) = -state.u(1,:).*state.u(2,:)-state.u(3,:).*state.u(2,:)+beta;
end

%diffusion and advection terms of PDE
function cmatrix = ccoeffunction(location,state,P)
N=P(1); D1=P(2); D2=P(3); D3=P(4); chi=P(5); eps=P(6);
n1 = 2*N;
nr = numel(location.x);
cmatrix = zeros(n1^2,nr);
cmatrix(1,:)=D1*ones(1,nr);
cmatrix(4,:)=D1*ones(1,nr);
cmatrix(5,:)=-2*state.u(2,:)./state.u(1,:);
cmatrix(8,:)=-2*state.u(2,:)./state.u(1,:);
cmatrix(9,:)=-(chi-eps)*state.u(3,:)./state.u(1,:);
cmatrix(12,:)=-(chi-eps)*state.u(3,:)./state.u(1,:);
cmatrix(2*n1+5,:)=D2*ones(1,nr);
cmatrix(2*n1+8,:)=D2*ones(1,nr);
cmatrix(n1^2-3,:)=D3*ones(1,nr);
cmatrix(n1^2,:)=D3*ones(1,nr);
end