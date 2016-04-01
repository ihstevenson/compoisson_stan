
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

%% MAP Estimation...

t = reshape(repmat(linspace(0,1,size(S,2)),size(S,1),1),[],1);
y = S(:);

figure(2); clf
mfit = cell(0);
for m = 1:3
    mfit{neuron,m} = model_optimize(y,t,false,model(m),'cubicbspline',5,9,[100 100],false);
    subplot(2,length(model),m)
    bar(linspace(0,1,size(S,2)),mean(S),1)
    hold on
    errorbar(linspace(0,1,size(S,2)),mean(S),std(S),'.')
    plot(mfit{neuron,m}.x_fine,mfit{neuron,m}.ey_fine,'r')
    hold off
    title(strrep(mfit{neuron,m}.model.model_name,'_',' '));
    axis tight
    ylabel('Spike Count')
    subplot(2,length(model),length(model)+m)
    plot(linspace(0,1,size(S,2)),var(S)./mean(S),'ko');
    hold on
    stairs(linspace(0,1,size(S,2)),vy./ey,'b')
    plot(mfit{neuron,m}.x_fine,mfit{neuron,m}.vy_fine./mfit{neuron,m}.ey_fine,'r')        
    hold off
    ylabel('Fano Factor')
    drawnow
end

