function  J = BodyJacobian( model, q, body_num, Xend)

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q)
    [q] = confVecToCell(model,q);
end

J = zeros( size(Xend,1), model.NV );

X = Xend;
j = body_num;
while j > 0
    [ XJ, S ] = jcalc( model.jtype{j}, q{j} );
    jj = model.vinds(j);
    J(:, jj) = X * S;
    if model.parent(j) > 0
        X = X *  XJ * model.Xtree{ j };
    end
    j = model.parent(j);
end
