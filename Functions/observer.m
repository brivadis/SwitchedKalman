function dz = observer(t, z, zpast, DeltaT, n, theta, gamma, A, B, b, C)

x = z(1:n);
xhat = z((n+1):2*n);
xhatpast = zpast((n+1):2*n, 1);

phase = z(end);
phasepast = zpast(end);

u = cont(phase, xhat);
Au = A + u*B;

upast = cont(phasepast, xhatpast);
Aupast = A + upast*B;

Sinv = vect2matrix(z((2*n+1):(2*n+n^2)));
S = vect2matrix(z((2*n+n^2+1):(2*n+2*n^2)));

Res = vect2matrix(z((2*n+2*n^2+1):(2*n+3*n^2)));
Gram = vect2matrix(z((2*n+3*n^2+1):(2*n+4*n^2)));

dx = Au*x + u*b;
dxhat = Au*xhat + u*b - Sinv*(C')*C*(xhat-x);

dSinv = Sinv*Au' + Au*Sinv + theta*Sinv - Sinv*(C')*C*Sinv;
dS = -Au'*S - S*Au - theta*S + C'*C;

if t>DeltaT
    dRes = -Res*(Au+gamma/2*eye(n))+(Aupast+gamma/2*eye(n))*Res;
    dGram = -(Au+gamma/2*eye(n))'*Gram - Gram*(Au+gamma/2*eye(n)) + C'*C - (C*Res)'*C*Res;
else
    dRes = -Res*(Au+gamma/2*eye(n));
    dGram = -(Au+gamma/2*eye(n))'*Gram - Gram*(Au+gamma/2*eye(n)) + C'*C;
end

% dRes = -Res*Au;
% 
% if t>DeltaT
%     dGram = -Au'*Gram - Gram*Au - (-Aupast'*Gram - Gram*Aupast);
% else
%     dGram = -Au'*Gram - Gram*Au + C'*C;
% end

dz = [dx; dxhat; matrix2vect(dSinv); matrix2vect(dS); matrix2vect(dRes); matrix2vect(dGram); 0];

end