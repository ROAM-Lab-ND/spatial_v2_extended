classdef RBD_model
    %RBD_Model 
    
    properties
        parent 
        NB % TODO - global change to NJ?
        NL
        NV
        NQ
        N_RB % number of rigid bodies
        I 
        Il
        Ir
        I_RB
        gravity
        appearance
        Xtree
        XtreeKin
        XtreeKinRot
        joint
        qinds
        vinds
        ancestors
        ancestor_vinds
        subtree_vinds
        successor_vinds
        param_inds
    end
    
    methods
        function a_grav = getGravity(model) 
            
            if length( model.gravity) > 0
              g = model.gravity;
            else
              g = [0;0;-9.81];
            end

            if size(model.Xtree{1},1) == 3		% is model planar?
              a_grav = [0;g(1);g(2)];
            else
              a_grav = [0;0;0;g(1);g(2);g(3)];
            end
        end
        
        %% Utilities
        function var_parent = getParentVariable(model, i, vars, default)

            if nargin == 3
                default = zeros(6,1);
            end

            p = model.parent;

            if p(i) == 0
                var_parent = default;
            else
                var_parent = vars{ p(i) };
            end
        end
        
        %%
        function model = postProcessModel(model)
            qi = 0;
            vi = 0;
            model.N_RB = 0;
            
            model.qinds = {};
            model.vinds = {};
            model.ancestors = {};
            model.ancestor_vinds = {};
            model.subtree_vinds = {};
            model.param_inds = {};

            % Connectivity processing
            for i = 1:model.NB
                nqi = model.joint{i}.nq;
                nvi = model.joint{i}.nv;
                
                model.qinds{i} = qi+1 : qi+nqi;
                model.vinds{i} = vi+1 : vi+nvi;

                if model.parent(i) == 0
                    model.ancestors{i} = [];
                    model.ancestor_vinds{i} = [];
                else
                    model.ancestors{i} = model.ancestors{ model.parent(i)};
                    model.ancestor_vinds{i} = model.ancestor_vinds{ model.parent(i)};
                end
                model.ancestors{i} = [model.ancestors{i}  i];
                model.ancestor_vinds{i} = [model.ancestor_vinds{i} model.vinds{i}];

                qi = qi + nqi;
                vi = vi + nvi;
                model.subtree_vinds{i} = [];
                model.successor_vinds{i} = [];
            end

            % more connectivity processing
            for i = model.NB:-1:1
                ii = model.vinds{i};
                model.subtree_vinds{i} = [ii model.subtree_vinds{i}];
                if model.parent(i) > 0
                    p = model.parent(i);
                    model.subtree_vinds{p} = [model.subtree_vinds{i} model.subtree_vinds{p}];
                    model.successor_vinds{p} = [ii model.successor_vinds{i} model.successor_vinds{p}];
                end
            end
            
            % Xtree initialization when previous group has more than one
            % body
            for i = 1:model.NB
               if model.parent(i) > 0
                  p = model.parent(i);
                  
                  if size(model.Xtree{i},2) < model.joint{p}.bodies*6
                     newXtree = zeros( size(model.Xtree{i},1) , model.joint{p}.bodies*6 );
                     
                     bod = model.joint{p}.output_body;
                     inds = 6*(bod-1)+1 : 6*bod;
                     
                     newXtree(:,inds) = model.Xtree{i};
                     model.Xtree{i} = newXtree;
                  end
               end
               
               model.param_inds{i} = (model.N_RB*10+1): ((model.N_RB+model.joint{i}.bodies)*10);
               model.N_RB = model.N_RB + model.joint{i}.bodies;
            end
            model.NV = vi;
            model.NQ = qi;
        end
       
        %%
        function q = normalizeConfVec(model,q)

            if ~isfield(model,'nq')
                model = model.postProcessModel();
            end

            if ~iscell(q)
                [q] = model.confVecToCell(q);
            end

            for i = 1:model.NB
                switch class( model.joint{i})
                    case {'floatingBaseJoint','sphericalJoint'}
                        q{i}(1:4) = q{i}(1:4)/norm( q{i}(1:4) );
                    otherwise
                        q{i} = q{i};
                end
            end
            q = cell2mat(q);
        end
        
        %%
        function [q_cell, v_cell, vd_cell, vd2_cell] = confVecToCell(model,q,v,vd, vd2)
            qj = 0;
            vj = 0;

            q_cell = cell(model.NB,size(q,2));
            if nargin > 2
                v_cell = cell(model.NB,size(v,2));
            end
            if nargin > 3
                vd_cell = cell(model.NB,size(vd,2));
            end
            if nargin > 4
                vd2_cell = cell(model.NB,size(vd2,2));
            end

            for i = 1:model.NB

                q_cell{i} = q(qj+1 : qj+model.joint{i}.nq,: );
                if nargin > 2
                    v_cell{i} = v( vj+1 : vj+model.joint{i}.nv,: );
                end
                if nargin > 3
                    vd_cell{i} = vd( vj+1 : vj+model.joint{i}.nv,: );
                end
                if nargin > 4
                    vd2_cell{i} = vd2( vj+1 : vj+model.joint{i}.nv,: );
                end

                qj = qj+model.joint{i}.nq;
                vj = vj+model.joint{i}.nv;
            end
        end
        
    end
end


