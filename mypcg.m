function [u,flag,relres,it]=mypcg(A,f,tol,maxit,M,u0)
% MYPCG solves linear equations with conjugate gradients.
%   [u,flag,relres,it]=mypcg(A,f,tol,maxit,M,u0)
%   solves Au=f in maximal maxit iterations and within
%   tol, the relative tolerance in the energy-norm defined by A.
%   M is optional and is used as a preconditioner, default M=I_n.
%   u0 is the optional starting vector, default u0=0.
%   A and M can be matrices or function handles.
%   flag 0 converged to desired tolerance in iter iterations
%   flag 1 not converged

if nargin<2,                    fprintf('Error not enough input arguments.\n'); end;
if nargin<3 || isempty(tol),    tol=1e-6; end;
if nargin<4 || isempty(maxit),  maxit=100; end;
if nargin<5,                    M=@(x)speye(size(A))*x; end;
if nargin<6,                    u0=0*f; end;

if isempty(M), M=@(x)speye(size(A))*x; end;
if ~(isa(M,'function_handle')), M=@(x) M*x;end;
if ~(isa(A,'function_handle')), A=@(x) A*x;end;

if length(u0)==length(f)
    u=u0;
else
    u=0*f;
end

d=f-A(u);
p=M(d);
it=0;
dMd=d(:)'*p(:);
dMd0=dMd;
for it=1:maxit
    Ap    = A(p);
    alpha = dMd/(p(:)'*Ap(:));
    u     = u + alpha*p;
%     res = f - A(u);
%     if res'*res < tol^2, flag=0;break, end
    dMd_old = dMd;
    d     = d - alpha*Ap;
    Md    = M(d);
    dMd   = d(:)'*Md(:);
    beta  = dMd/dMd_old;
    p     = Md + beta*p;  
    fprintf('    it %d relres %e tol %e\n',it,sqrt(dMd/dMd0),tol);
    if dMd/dMd0 < tol^2, flag=0; break,end
end
relres=(dMd/dMd0)/abs(dMd/dMd0)*sqrt(abs(dMd/dMd0));
if it==maxit, flag=1;end