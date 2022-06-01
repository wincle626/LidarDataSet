%% Copyright @ Heriot-Watt University UDRC WP 2.1
% Author: Yun Wu
function x = pgd(para, AtA, Atb, lambda, beta, gamma)

MAX_ITER = para.MAX_ITER;
% ABSTOL   = 1e-4;
x = zeros(size(AtA, 2),1);
% f = @(u) 0.5*sum_square(A*u-b);

for k = 1:MAX_ITER
%     while 1
        grad_x = AtA*x - Atb;
        z = prox_l1(x - lambda*grad_x, lambda*gamma);
%         if f(z) <= f(x) + grad_x'*(z - x) + (1/(2*lambda))*sum_square(z - x)
%             break;
%         end
        lambda = beta*lambda;
%     end
    x = z;

%     h.prox_optval(k) = objective(A, b, gamma, x, x);
%     if k > 1 && abs(h.prox_optval(k) - h.prox_optval(k-1)) < ABSTOL
%         break;
%     end
end

end

function p = objective(A, b, lambda, x, z)
    p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) );
end

function y = sum_square( x, dim )
    narginchk(1,2);
    y = x .* x;
    if nargin == 2
        y = sum( y, dim );
    else
        y = sum( y );
    end
end

function z = prox_l1(x, kappa)
    z = max( 0, x - kappa ) - max( 0, -x - kappa );
end