function  J = BodyJacobian( model, q, body_num, Xend)

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q)
    [q] = confVecToCell(model,q);
end

J = q{1}(1)*0 + zeros( size(Xend,1), model.NV );

X = Xend;
j = body_num;
while j > 0
    [ Xup, S ] = model.joint{j}.kinematics( model.Xtree{j} , q{j} );
    jj = model.vinds{j};
    J(:, jj) = X * S;
    if model.parent(j) > 0
        X = X *  Xup;
    end
    j = model.parent(j);
end
