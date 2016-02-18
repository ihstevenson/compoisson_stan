
% Given samples, return log-posterior, expected value of y, variance of y,
% and full posterior predictive

function [params,ey,vy,fullp] = getFitRes(model,fit,data,Xrbf0,maxsamps)

if nargin<5, maxsamps=Inf; end

try
    if strcmp(fit.model.command{2},'sample ')
        isPermuted=true;
    elseif strcmp(fit.model.command{2},'optimize ')
        isPermuted=false;
    else
        error('Unknown method STAN command');
    end
    params = fit.extract('permuted',false);
catch
    isPermuted=true;
    params=[];
end

switch model.model_name
    case 'poisson_regression'
        if isstruct(fit)
            b = fit.b;
            lp = fit.lp__;
        else
            b = fit.extract('permuted',isPermuted).b;
            lp = fit.extract('permuted',isPermuted).lp__;
        end
        if min(size(b))==1, b=b'; end
        
        if iscell(Xrbf0), Xrbf0=Xrbf0{1}; end
        ey=[]; vy=[]; fullp=zeros(data.smax+1,size(Xrbf0,1));
        for i=1:min(length(lp),maxsamps)
            ey(i,:) = exp(Xrbf0*b(i,:)');
            vy(i,:) = ey(i,:);
            for j=1:size(Xrbf0,1)
                fullp(:,j) = fullp(:,j)+poisspdf(0:data.smax,ey(i,j))';
            end
        end
    case 'negbino_regression'
        if isstruct(fit)
            b = fit.b;
            c = fit.c;
            lp = fit.lp__;
        else
            b = fit.extract('permuted',isPermuted).b;
            c = fit.extract('permuted',isPermuted).c;
            lp = fit.extract('permuted',isPermuted).lp__;
        end
        if min(size(b))==1, b=b'; c=c'; end
        
        if iscell(Xrbf0),
            XrbfG0=Xrbf0{2};
            Xrbf0=Xrbf0{1};
        else
            XrbfG0=Xrbf0;
        end
        ey=[]; vy=[]; fullp=zeros(data.smax+1,size(Xrbf0,1));
        for i=1:min(length(lp),maxsamps)
            ey(i,:) = exp(Xrbf0*b(i,:)');
            phi(i,:) = exp(XrbfG0*c(i,:)');
            vy(i,:) = ey(i,:) + ey(i,:).^2./phi(i,:);
            for j=1:size(Xrbf0,1)
                p = phi(i,j)./(ey(i,j)+phi(i,j));
                fullp(:,j) = fullp(:,j)+nbinpdf(0:data.smax,phi(i,j),p)';
            end
        end
    case {'compoisson_regression'}
        if isstruct(fit)
            b = fit.b;
            c = fit.c;
            lp = fit.lp__;
        else
            b = fit.extract('permuted',isPermuted).b;
            c = fit.extract('permuted',isPermuted).c;
            lp = fit.extract('permuted',isPermuted).lp__;
        end
        if min(size(b))==1, b=b'; c=c'; end
        
        if iscell(Xrbf0),
            XrbfG0=Xrbf0{2};
            Xrbf0=Xrbf0{1};
        else
            XrbfG0=Xrbf0;
        end
        ey=[]; vy=[];  fullp=zeros(data.smax+1,size(Xrbf0,1));
        for i=1:min(length(lp),maxsamps)
            lambda = exp(Xrbf0*b(i,:)');
            nu = exp(XrbfG0*c(i,:)');
            for j=1:length(lambda)
                pdf = com_pdf(0:data.smax, lambda(j), nu(j));
                pdf = pdf/sum(pdf);
                if any(~isfinite(pdf)), pdf(~isfinite(pdf))=0; end
                ey(i,j) = [0:data.smax]*pdf';
                vy(i,j) = [0:data.smax].^2*pdf' - ey(i,j)^2;
                fullp(:,j) = fullp(:,j)+pdf';
            end
        end
    case {'gencount_regression'}
        if isstruct(fit)
            b = fit.b;
            c = fit.c;
            lp = fit.lp__;
        else
            b = fit.extract('permuted',isPermuted).b;
            c = fit.extract('permuted',isPermuted).c;
            lp = fit.extract('permuted',isPermuted).lp__;
        end
        if min(size(b))==1, b=b'; c=c'; end
        
        if iscell(Xrbf0),
            Xrbf0=Xrbf0{1};
        else
            XrbfG0=Xrbf0;
        end
        ey=[]; vy=[];  fullp=zeros(data.smax+1,size(Xrbf0,1));
        for i=1:min(length(lp),maxsamps)
            theta = Xrbf0*b(i,:)';
            for j=1:length(theta)
                pdf = gc_pdf(0:data.smax, theta(j), c');
                pdf = pdf/sum(pdf);
                if any(~isfinite(pdf)), pdf(~isfinite(pdf))=0; end
                ey(i,j) = [0:data.smax]*pdf';
                vy(i,j) = [0:data.smax].^2*pdf' - ey(i,j)^2;
                fullp(:,j) = fullp(:,j)+pdf';
            end
        end
    otherwise
        error('Model type not found...')
end
fullp = bsxfun(@rdivide,fullp,sum(fullp));