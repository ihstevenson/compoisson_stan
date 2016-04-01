% model_optimize
% given     y       [trials x neurons] spike counts
%           theta   [trials x 1] stimulus/mov parameter
%           getSuffStats boolean whether to use sufficient statistics
%           model   compiled stan model
%           basis   string for basis name
%           nvec    number of bases for mean/lambda
%           mvec    number of basis for r/nu
%           sigma   s.d. of normal priors
%           verbose boolean, generates summary figure

function fout = model_optimize(y,theta,getSuffStats,model,basis,nvec,mvec,sigma,verbose)

% Compute sufficient statistics for speed...
if getSuffStats
    utheta = unique(theta);
    for neuron=1:size(y,2)
        for i=1:length(utheta)
            uy(i,neuron)   = sum(y(theta==utheta(i),neuron));
            uylf(i,neuron) = sum(gammaln(y(theta==utheta(i),neuron)+1));
        end
    end
else
    utheta = theta;
    uy = y;
    uylf = gammaln(y+1);
end

data.N = size(uy,1);
data.smax = max(1000,ceil(max(uy(:))*2));
data.lgm = gammaln(1:data.smax);

if ~isempty(sigma)
    data.sigma_b = sigma(1);
    data.sigma_c = sigma(2);
else
    data.sigma_b = 100;
    data.sigma_c = 10;
end

fout=[];
c=1;
for neuron=1:size(y,2)
    for m=1:length(mvec)
        for n=1:length(nvec)
            fprintf('Neuron %02i, k=%02i, l=%02i...\n',neuron,nvec(n),mvec(m))
            
            data.X=getBasis(basis,utheta,nvec(n));
            data.G=getBasis(basis,utheta,mvec(m));
            data.P = size(data.X,2);
            data.Q = size(data.G,2);
            
            data.y = uy(:,neuron);
            data.lyf = uylf(:,neuron);
            
            d = RData(data);
            d.type('X') = 'matrix';
            d.type('G') = 'matrix';
            
            fit = model.optimizing('data',d);
            addlistener(fit,'exit',@exitHandler);
            fit.block();
            fit.print();
            
            if fit.exit_value==0
                fout(neuron,c).model.model_name = model.model_name;
                fout(neuron,c).model.basis = basis;
                fout(neuron,c).n = nvec(n);
                fout(neuron,c).m = mvec(m);
                
                fout(neuron,c).x_fine = linspace(min(theta),max(theta),256)';
                Xrbf0 = getBasis(basis,fout(neuron,c).x_fine,nvec(n));
                Grbf0 = getBasis(basis,fout(neuron,c).x_fine,mvec(m));
                [fout(neuron,c).res,fout(neuron,c).ey_fine,fout(neuron,c).vy_fine,fullp] = getFitRes(model,fit,data,{Xrbf0,Grbf0});
                
                fout(neuron,c).x = unique(theta);
                Xrbf0 = getBasis(basis,fout(neuron,c).x,nvec(n));
                Grbf0 = getBasis(basis,fout(neuron,c).x,mvec(m));
                [~,fout(neuron,c).ey,fout(neuron,c).vy,fout(neuron,c).fullp] = getFitRes(model,fit,data,{Xrbf0,Grbf0});
            end
            
            if verbose
                try
                    figure(1); clf
                    subplot(2,1,1)
                    plot(theta,y,'bo')
                    hold on
                    plot(fout(neuron,c).x_fine,fout(neuron,c).ey_fine,'r')
                    hold off
                    title(strrep(fout(neuron,c).model.model_name,'_',' '));
                    ylabel('Spike Count')
                    subplot(2,1,2)
                    hold on
                    plot(fout(neuron,c).x_fine,fout(neuron,c).vy_fine./fout(neuron,c).ey_fine,'r')
                    hold off
                    ylabel('Fano Factor')
                    drawnow
                end
            end
            c=c+1;
        end
    end
end
end

