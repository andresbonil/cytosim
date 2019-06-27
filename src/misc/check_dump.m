function check_dump(path)

% This is used to explore Cytosim's linear system in matlab
% - load the matrices and vector from Cytosim's dump
% - plot convergence pattern of BICGstab, with and without preconditionning
%
% F. Nedelec, 16 Oct. 2014, March 2018, June 2018, 26 Jan 2019

if nargin < 1
    path = 'dump';
end

%% Loading

if isdir(path)

    cwd = pwd;
    cd(path);
    
    ord = load('ord.txt');
    time_step = load('stp.txt');
    fprintf(1, 'system of size %i with time_step %f\n', ord, time_step);
    
    obj = fread(fopen('obj.bin'), ord, 'double');
    drg = fread(fopen('drg.bin'), ord, 'double');
    sys = fread(fopen('sys.bin'), [ord, ord], 'double');
    ela = fread(fopen('ela.bin'), [ord, ord], 'double');  %elasticity matrix
    mob = fread(fopen('mob.bin'), [ord, ord], 'double');  %projection matrix
    con = fread(fopen('con.bin'), [ord, ord], 'double');  %preconditionner
    pts = fread(fopen('pts.bin'), ord, 'double');
    rhs = fread(fopen('rhs.bin'), ord, 'double');
    sol = fread(fopen('sol.bin'), ord, 'double');
    
    cd(cwd);
else
    error(['cannot find dump directory ',path]);
end

%% Check matrix

show(abs(sys)); set(gcf, 'name', 'System matrix');
%show(abs(mob)); set(gcf, 'name','Projection matrix');
%show(abs(ela)); set(gcf, 'name','Elasticity matrix');

if ( 1 )
    mat = eye(ord) - time_step * mob * ela;
    err0 = max(max(abs(mat-sys)));
    fprintf(1, 'norm8(system matrix - reconstituted matrix) : %e\n', err0);
    if ( err0 > 1e-8 )
        show(abs(mat));
        set(gcf, 'name', 'Reconstituted matrix');
    end
end
if ( 0 )
    figure('Position', [50 50 1000 1000]);
    plot(reshape(mat,1,ord*ord), reshape(sys,1,ord*ord), '.')
    xl = xlim;
    ylim(xl);
    xlabel('Reconstituted matrix');
    ylabel('cytosim matrix');
end
if ( 0 )
    figure; hold on;
    plot(abs(sys), '^b');
    plot(abs(mat), 'vr');
end
if ( 1 )
    figure('name', 'System matrix structure');
    spy(mat)
    drawnow;
end

%% check iterative solver
    
    function y = mfun1(x)
        y = con * x;
    end

    function y = mfun2(x, mode)
        if ( strcmp(mode, 'notransp') )
            y = con * x;
        else
            y = con' * x;
        end
    end

if 1
    tol = 0.0001;
    maxit = ord;
    
    sss = sparse(sys);
    
    % without preconditionning:
    [x0,fl0,rr0,itr,rv0] = bicgstab(sss, rhs, tol, maxit);
    fprintf(1, 'BCGS        converved after %6.1f vecmuls %f\n', 2*itr, rr0);
    
    figure('Position', [100 300 1400 600]);
    subplot(1,2,1);
    plot(x0, sol, 'k.');
    xlabel('matlab solution');
    ylabel('cytosim solution');
    xl = xlim;
    ylim(xl);
    
    subplot(1,2,2);
    semilogy(rv0/rv0(1),'b:', 'Linewidth', 2);
    
    xlabel('Number of M*V operations');
    ylabel('Relative residual');
    title('Solver convergence');
    hold on;
    
    %% Try different preconditionners
    
    iCON = inv(con);
    
    % with cytosim's preconditionner:
    [x0,fl0,rr0,itr,rv0] = bicgstab(sss, rhs, tol, maxit, @mfun1);
    semilogy(rv0/rv0(1),'b-', 'Linewidth', 2);
    fprintf(1, 'BCGS-P      converved after %6.1f vecmuls %f\n', 2*itr, rr0);
    
    for i = 3:8
        RS = 2^i;
        [x0,fl0,rr0,itr,rv0] = gmres(sss, rhs, RS, tol, maxit);
        semilogy(rv0/rv0(1),'k:', 'Linewidth', 2);
        fprintf(1, 'GMRES   %03i converved after %6.1f vecmuls %f\n', RS, (itr(1)-1)*RS+itr(2), rr0);
    end
    for i = 2:8
        RS = 2^i;
        [x0,fl0,rr0,itr,rv0] = gmres(sss, rhs, RS, tol, maxit, @mfun1);
        semilogy(rv0/rv0(1),'k-', 'Linewidth', 2);
        fprintf(1, 'GMRES-P %03i converved after %6.1f vecmuls %f\n', RS, (itr(1)-1)*RS+itr(2), rr0);
    end
    
    % preconditionner = incomplete LU factorization
    [L, U] = ilu(sss);
    [x0,fl0,rr0,itr,rv0] = bicgstab(sss, rhs, tol, maxit, L, U);
    semilogy(rv0/rv0(1),'b--', 'Linewidth', 2);
    fprintf(1, 'BCGS-LU      converved after %6.1f vecmuls %f\n', 2*itr, rr0);
   
    for i = 2:7
        RS = 2^i;
        [x0,fl0,rr0,itr,rv0] = gmres(sss, rhs, RS, tol, maxit, L, U);
        semilogy(rv0/rv0(1),'k--', 'Linewidth', 2);
        fprintf(1, 'GMRES-LU %03i converved after %6.1f vecmuls %f\n', RS, (itr(1)-1)*RS+itr(2), rr0);
    end
    
    if ( 0 )
        % preconditionner = Symmetric successive over-relaxation
        L = tril(sss);
        U = triu(sss);
        V = diag(sss);
        D = diag(V);
        M = (D+L)*diag(1./V)*(D+U);
        [x0,fl0,rr0,itr,rv0] = bicgstab(sss, rhs, tol, maxit, M);
        semilogy(rv0/rv0(1),'b--', 'Linewidth', 2);
        fprintf(1, 'BCGS-SSOR     converved after %6.1f vecmuls %f\n', 2*itr, rr0);
        
        for i = 2:7
            RS = 2^i;
            [x0,fl0,rr0,itr,rv0] = gmres(sss, rhs, RS, tol, maxit, M);
            semilogy(rv0/rv0(1),'k--', 'Linewidth', 2);
            fprintf(1, 'GMRES-SSOR %03i converved after %6.1f vecmuls %f\n', RS, (itr(1)-1)*RS+itr(2), rr0);
        end
    end
    
    if ( 0 )
        % QMR method
        [x0,fl0,rr0,itr,rv0] = qmr(sss, rhs, tol, maxit);
        semilogy(rv0/rv0(1),'m:', 'Linewidth', 2);
        fprintf(1, 'QMR         converved after %6.1f vecmuls %f\n', 2*itr, rr0);
        
        [x0,fl0,rr0,itr,rv0] = qmr(sss, rhs, tol, maxit, @mfun2);
        semilogy(rv0/rv0(1),'m-', 'Linewidth', 2);
        fprintf(1, 'QMR-P       converved after %6.1f vecmuls %f\n', 2*itr, rr0);
    end
    
    % IDRS method
    OPT.smoothing = 1;
    for i = 3:8
        RS = 2^i;
        [x0,fl0,rr0,itr,rv0] = idrs(sss, rhs, RS, tol, maxit, [], [], [], OPT);
        semilogy(rv0/rv0(1),'r:', 'Linewidth', 2);
        fprintf(1, 'IDRS   %03i  converved after %6.1f vecmuls %f\n', RS, itr, rr0);
    end
    for i = 2:7
        RS = 2^i;
        [x0,fl0,rr0,itr,rv0] = idrs(sss, rhs, RS, tol, maxit, iCON, [], [], OPT);
        semilogy(rv0/rv0(1),'r--', 'Linewidth', 2);
        fprintf(1, 'IDRS-P %03i  converved after %6.1f vecmuls %f\n', RS, itr, rr0);
    end
end

%% check different preconditionners

if 1 
    
end

 
end