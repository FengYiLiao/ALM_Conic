clc; clear; close all
yalmip('clear')
%ALM for Maxcut
%Package requirement: Yalmip and Mosek

datapath = '';
savepath = 'results/';
name     = {'n20r5MC'};
idx        = 1;





load([datapath,name{idx},'.mat']);

% n   = 20;
% At  = At(1:n,:);
% b   = b(1:n);
% C   = reshape(c,100,100);
% C   = C(1:n,1:n);
% K.s = n;
% c   = reshape(C,[],1);
% At = At_sdp;
% b  = b_sdp;
% c  = c_sdp;
% K  = K_sdp;
% save('n20m20dr10.mat','At','b','c','K','Optimal');%


%At = At';



[info,mosektime]=SolveMosek(At,b,c,K);
% TrueCost  = c.'*xx;
Optimal.Cost = info.sol.itr.pobjval;
% XX = [M,M;M,M];



m        = height(At);
n        = sqrt(width(At));
r        = 1; %augmetned term

x        = sdpvar(n); %Yalmip variable
u        = sdpvar(n);


Xtrue                    = zeros(n);
[IndSym,IndDiag,IndOffDiag,ShrinkIndDiag,ShrinkIndOffDiag,IndOffDiagCounter] = SymmetricIndices(n,false);
Xtrue(IndSym)            = info.sol.itr.barx;
Xtrue(IndOffDiagCounter) = Xtrue(IndOffDiag);

ytrue                    = info.sol.itr.y;
Ztrue                    = reshape(c -At.'*ytrue,n,n);


yk       = zeros(m,1);
z        = c - At.'*yk;
zk       = reshape(z,n,n);
C        = reshape(z,n,n);


Max_iter     = 10;
cost         = [];
Dcost        = [];
xstar        = [];
Affinefeasi  = [];
DAffinefeasi = [];
Conefeasi    = [];
DualGap      = [];
Dist        = [];
DDist        = [];
normb        = norm(b);
normC        = norm(C,'fro');
normXtrue    = norm(Xtrue,'fro');
normytrue    = norm(ytrue,'fro');
normZtrue    = norm(Ztrue,'fro');




for k = 1:Max_iter
    %lk        = trace(reshape(c,n,n)*x) + trace(x*u) + y.'*(b-At*vec(x)) + r/2*trace(u*u) +r/2*(At*vec(x)-b).'*(At*vec(x)-b);
    %lk        = trace(reshape(c,n,n)*x) + y.'*(b-At*vec(x)) +r/2*(At*vec(x)-b).'*(At*vec(x)-b);
    %lk        = trace(C*x) - trace(zk*(x-u))  + r/2*trace((x-u)*(x-u)) ;
    bAx       = b-At*vec(x);
    xu        = x-u; %trace(C*x)
    %lk        =  c(:).'*x(:) - zk(:).'*xu(:) + r/2*xu(:).'*xu(:) + yk.'*bAx +r/2*(bAx).'*(bAx);
    %lk        =  c(:).'*x(:)  + yk.'*bAx +r/2*(bAx).'*(bAx);
    lk        =  c(:).'*x(:) + 1/(2*r)*norm(r*x-zk-u,'fro')^2 + 1/(2*r)*norm(yk+r*bAx)^2;
    %lk        =  c(:).'*x(:) + 1/(2*r)*trace((r*x-zk-u).'*(r*x-zk-u)) + 1/(2*r)*(yk+r*bAx).'*(yk+r*bAx);


    %cons      = [At*vec(x)==b,u>=0];
    cons      = [u>=0];
    %cons      = [x>=0];

    ops = sdpsettings('solver','mosek');
    %ops.cdcs.relTol = 1e-5;
    %optimize(cons,lk,ops);
    optimize(cons,lk);

    xk        = value(x);
    cost(k)   = c.'*x(:);
    xstar{k}  = xk;
    
    %dual update
    tempzk    = zk - r*xk;
    [V,D]     = eig(tempzk); 
    D(D<=0)   = 0;
    zk        = V*D*V';
    zk        = (zk+zk')/2;

    yk        = yk+r*value(bAx);
    
    Ztest     = reshape(c - At.'*yk,n,n);
    
    %residual
    Dcost        = [Dcost,b.'*yk]; 
    Affinefeasi  = [Affinefeasi,norm(b-At*vec(xk))/(1+normb)];
    DAffinefeasi = [DAffinefeasi,norm(Ztest-zk)/(1+normC)];
    %DualGap      = [DualGap,abs(cost(k) - b'*yk)];
    [V,D]        = eig(xk);
    Conefeasi    = [Conefeasi,norm(D(D<0))/(1+norm(xk))];
    Dist         = [Dist,norm(xk-Xtrue,'fro')/(1+normXtrue)];
    DDist        = [DDist,sqrt(norm(yk-ytrue)^2+norm(zk-Ztrue,'fro')^2)/(1+normytrue+normZtrue)];
end


% semilogy(abs(cost-Optimal.Cost)/abs(Optimal.Cost));
% hold on 

% semilogy(Conefeasi);
%semilogy(Dist);
%hold on
%semilogy(Affinefeasi);
%semilogy(DDist);
%xlim([0,5]);
%semilogy(abs(Optimal.Cost-Dcost)/abs(Optimal.Cost));
%semilogy(DAffinefeasi);
%semilogy(DualGap);
Out.PCostgap     = abs(cost-Optimal.Cost)/abs(Optimal.Cost);
Out.DCostgap     = abs(Dcost-Optimal.Cost)/abs(Optimal.Cost);

%
Out.Conefeasi    = Conefeasi;
Out.PCost        = cost;
Out.DCost        = Dcost;
Out.Affinefeasi  = Affinefeasi;
Out.DAffinefeasi = DAffinefeasi;
Out.Dist         = Dist;
Out.DDist        = DDist;
Out.OptimalCost  = info.sol.itr.pobjval;


semilogy(Out.PCostgap);
%semilogy(Out.Dist);
hold on
%semilogy(Out.DDist);
%semilogy(Out.DCostgap);
semilogy(Out.Conefeasi);

save([savepath,name{idx},'_result'],'Out');
