%parameters
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

%these are used for IC
Q=(1/2).*D2.^(-1).*ubar.^(-2).*((-1).*alpha+ubar).^(-1).*((-1).* ...
  alpha.^2.*((-2)+D2)+alpha.*((-4)+D2).*ubar+(2+(-1).*beta.*D1).* ...
  ubar.^2);

R=(1/4).*D1.^(-1).*D2.^(-1).*ubar.^(-2).*((-1).*alpha+ubar).^(-3).*( ...
  alpha.^4.*(2+D2).^2+(-2).*alpha.^3.*(8+6.*D2+D2.^2).*ubar+ ...
  alpha.^2.*(24+2.*beta.*D1.*((-2)+D2)+4.*(3+alpha.*D1).*D2+D2.^2).* ...
  ubar.^2+(-2).*alpha.*(beta.*D1.*((-4)+D2)+2.*(4+D2+3.*alpha.*D1.* ...
  D2)).*ubar.^3+(4+(-4).*beta.*D1+beta.^2.*D1.^2+12.*alpha.*D1.*D2) ...
  .*ubar.^4+(-4).*D1.*D2.*ubar.^5);

%epsilon and initial magnitude
eps = -0.01;
P0=0.08;

%various intesities for hotspot arrests
kappa1=5;
kappa2=1;

%Defining time and space
Tmax=4000; %length of pre-suppression
Tdamp=500; %length of suppression
Trem=4000; %length of post-suppression

wk=sqrt(K);
a = 0; 
b = 4*pi/wk;  
L = b-a;

lowerLeft  = [a   ,a   ];
lowerRight = [b , a  ];
upperRight = [b , b];
upperLeft =  [a , b];
% Setting geometry
S = [3,4 lowerLeft(1), lowerRight(1), upperRight(1), upperLeft(1), ...
         lowerLeft(2), lowerRight(2), upperRight(2), upperLeft(2)];                     
gdm = S';
ns = 'S';
sf = 'S';
g = decsg(gdm,ns,sf');
d=[1;1;1];
Pf=[N, alpha, beta ];
Pc=[N, D1, D2,D3,chi,eps];

model = createpde(N);
geometryFromEdges(model,g);
generateMesh(model);%,"Hmax",0.1);

applyBoundaryCondition(model,"neumann","Edge",1:4);
specifyCoefficients(model,"m",0,"d",d,"c",@(location,state)ccoeffunction(location,state,Pc),"a",0,"f",@(location,state)fcoeffunction(location,state,Pf));
%initial conditions
icFcn = @(region) [ubar + P0*(cos(region.x*wk)+cos(region.y*wk));
                   vbar + P0*(cos(region.x*wk)+cos(region.y*wk))*Q;
                   wbar + P0*(cos(region.x*wk)+cos(region.y*wk))*R];
setInitialConditions(model,icFcn);

%solve for pre-suppression
tlist1 = linspace(0,Tmax,Tmax+1);
results1 = solvepde(model,tlist1);

u1=reshape(results1.NodalSolution(:,1,:),[size(results1.NodalSolution(:,1,:),1),size(results1.NodalSolution(:,1,:),3)]);
v1=reshape(results1.NodalSolution(:,2,:),[size(results1.NodalSolution(:,2,:),1),size(results1.NodalSolution(:,2,:),3)]);
w1=reshape(results1.NodalSolution(:,3,:),[size(results1.NodalSolution(:,3,:),1),size(results1.NodalSolution(:,3,:),3)]);

%saving mesh
X=model.Mesh.Nodes(1,:);
Y=model.Mesh.Nodes(2,:);

%ic for the next part
utail=u1(:,end);

%setting and solving suppression part
Pf2=[N, alpha, beta,kappa1,kappa2];
specifyCoefficients(model,"m",0,"d",d,"c",@(location,state)ccoeffunction(location,state,Pc),"a",0,"f",@(location,state)fcoeffunction2(location,state,Pf2,utail,X,Y));
setInitialConditions(model,results1)
tlist2=linspace(Tmax+1,Tmax+Tdamp,Tdamp);
results2 = solvepde(model,tlist2);

u2=reshape(results2.NodalSolution(:,1,:),[size(results2.NodalSolution(:,1,:),1),size(results2.NodalSolution(:,1,:),3)]);
v2=reshape(results2.NodalSolution(:,2,:),[size(results2.NodalSolution(:,2,:),1),size(results2.NodalSolution(:,2,:),3)]);
w2=reshape(results2.NodalSolution(:,3,:),[size(results2.NodalSolution(:,3,:),1),size(results2.NodalSolution(:,3,:),3)]);

%setting and solving post-suppression part
specifyCoefficients(model,"m",0,"d",d,"c",@(location,state)ccoeffunction(location,state,Pc),"a",0,"f",@(location,state)fcoeffunction(location,state,Pf));
setInitialConditions(model,results2)
tlist3=linspace(Tmax+Tdamp+1,Tmax+Tdamp+Trem,Trem);
results3 = solvepde(model,tlist3);

u3=reshape(results3.NodalSolution(:,1,:),[size(results3.NodalSolution(:,1,:),1),size(results3.NodalSolution(:,1,:),3)]);
v3=reshape(results3.NodalSolution(:,2,:),[size(results3.NodalSolution(:,2,:),1),size(results3.NodalSolution(:,2,:),3)]);
w3=reshape(results3.NodalSolution(:,3,:),[size(results3.NodalSolution(:,3,:),1),size(results3.NodalSolution(:,3,:),3)]);

%plotting rms with highlighted/colored backgrounds
figure()
hold on
r1=rectangle('Position',[tlist1(1),0,tlist1(end)-tlist1(1),0.45],'FaceColor',[0.8500, 0.3250, 0.0980],'FaceAlpha',0.1,'EdgeColor','k',...
    'LineWidth',0.1);
r2=rectangle('Position',[tlist2(1),0,tlist2(end)-tlist2(1),0.45],'FaceColor',[0.4940, 0.1840, 0.5560],'FaceAlpha',0.1,'EdgeColor','k',...
    'LineWidth',0.1);
r3=rectangle('Position',[tlist3(1),0,tlist3(end)-tlist3(1),0.45],'FaceColor',[0, 0.4470, 0.7410],'FaceAlpha',0.1,'EdgeColor','k',...
    'LineWidth',0.1);
plot(tlist1,rms(u1-mean(u1(:,end))),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2)
plot(tlist2,rms(u2-mean(u2(:,end))),'Color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
plot(tlist3,rms(u3-mean(u3(:,end))),'Color',[0, 0.4470, 0.7410],'LineWidth',2)
xlim([tlist1(1),tlist3(end)])
ylabel('A_{amp}(t)')
xlabel('t')
hold off

%plotting the end of pre-suppression
figure()
pdeplot(results1.Mesh,XYData=u1(:,end),ColorBar='on')
colormap('jet')
xlim([a,b])
xlabel('x')
ylim([a,b])
ylabel('y')
pbaspect([1 1 1])

%plotting the beginning of suppression
figure()
pdeplot(results2.Mesh,XYData=u2(:,3),ColorBar='on')
colormap('jet')
xlim([a,b])
xlabel('x')
ylim([a,b])
ylabel('y')
pbaspect([1 1 1])

%plotting the middle of post-suppression
figure()
pdeplot(results3.Mesh,XYData=u3(:,5500),ColorBar='on')
colormap('jet')
xlim([a,b])
xlabel('x')
ylim([a,b])
ylabel('y')
clim([0.49,0.51])
pbaspect([1 1 1])

%reaction terms for post- and pre-suppresion - identical to the original
function f = fcoeffunction(location,state,P)
N=P(1); alpha=P(2); beta=P(3); 
nr = length(location.x); 
f = zeros(N,nr); 
f(1,:) = -state.u(1,:)+state.u(1,:).*state.u(2,:)+alpha;
f(2,:) = -state.u(1,:).*state.u(2,:)-state.u(3,:).*state.u(2,:)+beta;
end

%reaction terms for suppression
function f = fcoeffunction2(location,state,Pf2,utail,X,Y)
N=Pf2(1); alpha=Pf2(2); beta=Pf2(3); kappa1=Pf2(4); kappa2=Pf2(5);
nr = length(location.x); 
f = zeros(N,nr); 
us=griddata(X, Y, utail, location.x, location.y);
As=quantile(utail,0.5,'all'); %Acutoff is defined as median/quantile
damp=kappa1*(us>As)+kappa2*(us<=As);
f(1,:) = -state.u(1,:)+state.u(1,:).*state.u(2,:)+alpha;
f(2,:) = -state.u(1,:).*state.u(2,:)-damp.*state.u(3,:).*state.u(2,:)+beta; %change suppression function here
end

%diffusion and advection terms are identical to the original
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

