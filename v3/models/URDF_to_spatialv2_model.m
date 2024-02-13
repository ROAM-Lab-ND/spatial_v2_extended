function [model,  robot] = URDF_to_spatialv2_model(file, addRandomInertia)
    robot = importrobot(file);
    model = struct();
    model.NB = robot.NumBodies;
    if nargin == 1
        addRandomInertia = 0;
    end
    
    % Loop over bodies to populate spatial_v2 model structure
    for i =1:robot.NumBodies
        
        body = robot.Bodies{i};
        parentName = body.Parent.Name;
        p = findString(robot.BodyNames, parentName); % Index of parent
        model.parent(i) = p; % spatial_v2 parent array
        parent = body.Parent;
        joint = body.Joint;

        %% Parse frame attachment information 

        % Frames attached to body i are:
        %   i+ [frame after joint i <- main frame used for spatial_v2 computations]
        %   i  [URDF Frame i at CoM -- not used in spatial_v2]
        %   for each child j, we attach frame j- to body i [just before joint j]
        % 
        % For spatial_v2 we need:
        %   Xtree{i} = {}^{i-} X_{p+} where p is the parent

        T_p_iminus = joint.JointToParentTransform;
        if p == 0 
            % There is no joint preceeding the base.
            T_pplus_p = eye(4);
        else
            T_pplus_p = parent.Joint.ChildToJointTransform;
        end
        T_pplus_iminus = T_pplus_p * T_p_iminus;
        
        % Check for +/-1, 0 and overwrite to exact values where sensible
        rounded = round(T_pplus_iminus);
        flag = abs( rounded-T_pplus_iminus) < 1e-5;
        T_pplus_iminus( flag) = rounded(flag);


        % X_iminus_pplus 
        model.Xtree{i} = AdjointRepresentation( inv( T_pplus_iminus ) );
        model.Xtree_sym{i} =  model.Xtree{i}; % Needed for RPNA

        %% Parse Joint information

        % Get joint type
        if strcmp(joint.Type,'revolute')
            jtype = 'R';
        elseif strcmp(joint.Type,'prismatic')
            jtype = 'P';
        else
            % TODO: Implemented 'fixed' connections
            assert(1==0,'Only revolute and prismatic supported.')
        end
        
        % Get joint axis
        if all( joint.JointAxis == [1 0 0])
            axis = 'x';
        elseif all( joint.JointAxis == [0 1 0])
            axis = 'y';
        elseif all( joint.JointAxis == [0 0 1])
            axis = 'z';
        elseif all( joint.JointAxis == [-1 0 0])
            axis = 'x-';
        elseif all( joint.JointAxis == [0 -1 0])
            axis = 'y-';
        elseif all( joint.JointAxis == [0 0 -1])
            axis = 'z-';
        else
            % TODO: support non-axially aligned joints by redefining joint
            % frames
            assert(1==0,'Only axially aligned joints currently supported');
        end
        jtype = [jtype axis];
        model.jtype{i}= jtype;
   
        %% Parse inertial properties in URDF
        m = body.Mass;
        % If no mass, then assign random inertial properties.
        if m == 0 && addRandomInertia==1
            fprintf('Body %s has no mass listed in URDF:\n   adding random inertial properties\n', body.Name);
            m = rand();
            body.Mass = m;
            body.CenterOfMass = rand(1,3);
            I3 = randomPositiveDefinite(3);
            I3 = diag([1 2 3]);
            body.Inertia = [I3(1,1) I3(2,2) I3(3,3) I3(2,3) I3(1,3) I3(1,2)];
        end

        c = body.CenterOfMass;
        Ivec = body.Inertia;
        Ixx = Ivec(1); Iyy =Ivec(2); Izz = Ivec(3);
        Iyz = Ivec(4); Ixz =Ivec(5); Ixy = Ivec(6);
        I_3D = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];
        I_i = [I_3D m*skew(c) ; -m*skew(c) m*eye(3)];
    
        % Transform inertia to correct frame:
        %    URDF has a link frame at the link COM that is separate from the 
        %    frame after each joint. spatial_v2 does not. So, for each body i, 
        %     we need to convert from [URDF] link i frame data to frame i+ data.
        T_iplus_i = joint.ChildToJointTransform;
        X_i_iplus = AdjointRepresentation( inv(T_iplus_i) );
        model.I{i} = X_i_iplus'*I_i*X_i_iplus;
    end
end

% Find a string in a list of strings. Return 0 if not found.
function index = findString(list, str)
    index = 0;
    for i = 1:length(list)
        if strcmp(list{i},str)
            index = i;
            break
        end
    end
end

% Convert from homogeneous transformation to spatial transformation
% i.e., from SE(3) homogeneous transform to its Adjoint matrix.
function X = AdjointRepresentation(T)
    R = T(1:3,1:3);
    p = T(1:3,4);
    X = [R zeros(3); skew(p)*R R];
end