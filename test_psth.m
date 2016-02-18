
% run setup.m first to compile the models in stan

%% Simulate PSTHs with COM-Poisson noise...
%   Following Fig 6 from Stevenson, JCNS, 2016
%   generates simulated variably-dispersed tuning curve data

ntrials = 50;
n = 0:1000; % range to calculate pdf over
S = []; ey=[]; vy=[];

% varying mean and fano factor
% mu = linspace(5,1,100);
% nu = linspace(2,0.25,100);

% % fixed mean, varying fano factor
mu = ones(1,100)*5;
mu(31:70)=5.255;
mu(71:end)=5;
nu  = ones(1,100);
nu(31:70)=2;
nu(71:end)=1;

for i=1:ntrials
    S(i,:)=com_rnd(mu.^nu,nu);
end
for j=1:length(mu)
    pdf = com_pdf(n, mu(j).^nu(j), nu(j));
    pdf = pdf/sum(pdf);
    ey(j) = n*pdf';
    vy(j) = n.^2*pdf' - ey(j)^2;
end

figure(1); clf
subplot(2,1,1)
bar(mean(S),1)
hold on
errorbar(mean(S),std(S),'.')
plot(ey,'r')
hold off
axis tight
box off; set(gca,'TickDir','out')
ylabel('Spike Count')
subplot(2,1,2)
stairs(vy./ey,'r')
hold on
plot(var(S)./mean(S),'o')
hold off
box off; set(gca,'TickDir','out')
xlabel('Time')
ylabel('Fano Factor')

%%

t = reshape(repmat(linspace(0,1,size(S,2)),size(S,1),1),[],1);
y = S(:);

mfit = cell(0);
mfit{1} = fitModels(y,t,false,model(2),'cubicbspline',1,1);
mfit{2} = fitModels_v2(y,theta,model(2),true,1:20,1);
mfit{3} = fitModels_v2(y,theta,model(3),true,1:20,1);
mfit{4} = fitModels_v2(y,theta,model(4),true,1:20,2:10);

