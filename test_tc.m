
% run setup.m first to compile the models in stan

%% Simulate tuning curves with COM-Poisson noise...
%   Following Fig 3 from Stevenson, JCNS, 2016
%   generates simulated under-dispersed, over-dispersed, and
%   variably-dispersed tuning curve data

rep = 16;
n = 0:1000; % range to calculate pdf over

tcAll=cell(1); ey=[]; vy=[];
x0 = linspace(0,2*pi-2*pi/16,16);
Xrbf = getBasis('fourier',x0',2);

b = [2.5 1 .1 0 0;...
    2.5 0.5 -0.5 0 0;...
    2 0.1 -1 0 0];

nu = [(3+sin(x0)*1);...
    (0.4-cos(x0)*0.05-sin(x0)*-0.05);...
    cos(x0)/2+1];

for neuron=1:size(b,1)
    tcAll{neuron} = [];
    for i=1:length(x0)
        mu = exp(b(neuron,:)*Xrbf(i,:)');
        nun = nu(neuron,i);
        samp = com_rnd(mu.^nun, nun, rep);
        tcAll{neuron} = [tcAll{neuron}; [repmat(x0(i),rep,1) samp]];
        
        pdf = com_pdf(n, mu.^nun, nun);
        pdf = pdf/sum(pdf);
        ey(neuron,i) = n*pdf';
        vy(neuron,i) = n.^2*pdf' - ey(neuron,i)^2;
        ff(neuron,i) = var(samp)/mean(samp);
    end
    
    figure(neuron); clf;
    subplot(2,1,1)
    hold on
    plot(tcAll{neuron}(:,1)+randn(size(tcAll{neuron}(:,1)))/20,tcAll{neuron}(:,2),'o')
    plot(x0,ey(neuron,:),'k')
    axis tight
    ylabel('Spike Count')
    subplot(2,1,2)
    hold on
    plot(x0,vy(neuron,:)./ey(neuron,:),'k')
    plot(x0,ff(neuron,:),'o')
    axis tight
    ylim([0 3])
    xlabel('Stimulus Direction')
    ylabel('Fano Factor')
end

figure(4);
clf
plot(ey',vy')
line(xlim,xlim)
box off; set(gca,'TickDir','out')
xlabel('Mean Spike Count')
ylabel('Spike Count Variance');

%% MAP Estimation...

for neuron = 1
    theta = tcAll{neuron}(:,1);
    y = tcAll{neuron}(:,2);
    figure(neuron); clf
    mfit = cell(0);
    for m = 1:3
        mfit{neuron,m} = model_optimize(y,theta,false,model(m),'fourier',2,0,[100 1000],false);
        subplot(2,length(model),m)
        plot(theta,y,'bo')
        hold on
        plot(mfit{neuron,m}.x_fine,mfit{neuron,m}.ey_fine,'r')
        plot(x0,ey(neuron,:),'k');
        hold off
        title(strrep(mfit{neuron,m}.model.model_name,'_',' '));
        ylabel('Spike Count')
        subplot(2,length(model),length(model)+m)
        hold on
        plot(mfit{neuron,m}.x_fine,mfit{neuron,m}.vy_fine./mfit{neuron,m}.ey_fine,'r')
        plot(x0,vy(neuron,:)./ey(neuron,:),'k');
        plot(x0,ff(neuron,:),'o')
        hold off
        ylabel('Fano Factor')
        drawnow
    end
end

%% Bayesian Estimation...

for neuron = 1
    theta = tcAll{neuron}(:,1);
    y = tcAll{neuron}(:,2);
    figure(neuron); clf
    mfit = cell(0);
    for m = 1
        mfit{neuron,m} = model_sample(y,theta,false,model(m),'fourier',2,0,[100 1000],false);
        subplot(2,length(model),m)
        plot(theta,y,'bo')
        hold on
        plot(mfit{neuron,m}.x_fine,mfit{neuron,m}.ey_fine,'r')
        plot(x0,ey(neuron,:),'k');
        hold off
        title(strrep(mfit{neuron,m}.model.model_name,'_',' '));
        ylabel('Spike Count')
        subplot(2,length(model),length(model)+m)
        hold on
        plot(mfit{neuron,m}.x_fine,mfit{neuron,m}.vy_fine./mfit{neuron,m}.ey_fine,'r')
        plot(x0,vy(neuron,:)./ey(neuron,:),'k');
        plot(x0,ff(neuron,:),'o')
        hold off
        ylabel('Fano Factor')
        drawnow
    end
end

