function  Jdot = BodyJacobianTimeDerivative( model, q, qd, body_num, Xend, Xend_dot)

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q)
    [q,qd] = confVecToCell(model,q,qd);
end

if nargin == 4
    Xend = eye( 6 );
end
if nargin <= 5
    Xend_dot = 0*Xend;
end

Jdot = q{1}(1)*0 + zeros( size(Xend,1), model.NV );

v = {};
for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
  end
  Sd{i} = crm(v{i})*S{i};
end

X  = Xend;
X2 = Xend_dot -Xend * crm( v{body_num} );

j = body_num;
while j > 0
    jj = model.vinds{j};
    Jdot(:,jj) = X * Sd{j} + X2 * S{j};
    if model.parent(j) > 0
        X  = X  * Xup{j};
        X2 = X2 * Xup{j};
    end
    j = model.parent(j);
end