function[H]=BFGS(yk,dk,H0)
H=(eye(3)-(dk*yk')/(dk'*yk))*H0*(eye(3)-(yk*dk')/(dk'*yk));
H=H+(dk*dk')/(dk'*yk);