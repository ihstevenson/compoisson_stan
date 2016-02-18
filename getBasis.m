function X = getBasis(basis,x,n)

switch lower(basis)
    case 'cubicbspline_circ'
        X = getCubicBSplineBasis(x,n,true);
    case 'cubicbspline'
        X = getCubicBSplineBasis(x,n,false);
    case 'fourier'
        X = getFourierBasis(x,n);
    case 'polynomial'
        X = getPolynomialBasis(x,n);
    otherwise
        error('Unknown basis');
end

function b = getCubicBSplineBasis(x,nknots,isCirc)

if ~isCirc & nknots>1
    nknots=nknots-1; % for consistency across circ/non-circ
    weights = ones(nknots+1,1);
    knots = linspace(-2/nknots,1+2/nknots,nknots+5);
    s = fastBSpline(knots,weights);
    b = s.getBasis(x);
else
    knots = linspace(-2*2*pi/nknots,2*pi+2*2*pi/nknots,nknots+5);
    weights = ones(nknots+1,1);
    s = fastBSpline(knots,weights);
    b = zeros(length(x),nknots+1);
    for k=-4:4
        xk = x+k*2*pi;
        b = b+s.getBasis(xk);
    end
    b = b(:,1:end-1);
end

if nknots>1
    b = [b(:,1)*0+1 b];
end

function b = getFourierBasis(x,n)

b=[];
for i=1:n
    b(:,2*i-1) = sin(x*i);
    b(:,2*i)   = cos(x*i);
end
b = [x*0+1 b];

function b = getPolynomialBasis(x,n)

b=[];
for i=1:n
    b(:,i)   = x.^i;
end
b = [x*0+1 b];