function [model robot] = UDRF_to_spatialv2_model(file)
    robot = importrobot(file);
    model = struct();
    model.NB = robot.NumBodies;
    
    for i =1:robot.NumBodies
        body = robot.Bodies{i};
        parentName = body.Parent.Name;
        p = findString(robot.BodyNames, parentName);
        model.parent(i) = p;
        parent = body.Parent;
         
        joint = body.Joint;

        T_pi_iminus = joint.JointToParentTransform;
        if p == 0 
            T_piplus_pi = eye(4);
        else
            T_piplus_pi = parent.Joint.ChildToJointTransform;
        end
        T_piplus_iminus = T_piplus_pi * T_pi_iminus;
        
        rounded = round(T_piplus_iminus);
        
        flag = abs( rounded-T_piplus_iminus) < 1e-5;
        T_piplus_iminus( flag) = rounded(flag);




        % X_iminus_piplus 
        model.Xtree{i} = AdjointRepresentation( inv( T_piplus_iminus ) );
        model.Xtree_sym{i} =  model.Xtree{i};

        % Get joint type
        if strcmp(joint.Type,'revolute')
            jtype = 'R';
        elseif strcmp(joint.Type,'prismatic')
            jtype = 'P';
        else
            assert(1==0,'Only revolute and prismatic supported.')
        end
        
        % Get joint axis
        if all( joint.JointAxis == [1 0 0])
            axis = 'x';
        elseif all( joint.JointAxis == [0 1 0])
            axis = 'y';
        elseif all( joint.JointAxis == [0 0 1])
            axis = 'z';
        else
            assert(1==0,'Only axially aligned joints currently supported');
        end
        jtype = [jtype axis];
        model.jtype{i}= jtype;
        
        % Get inertia in correct frame

        m = body.Mass;
        if m == 0
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
    

        T_iplus_i = joint.ChildToJointTransform;
        X_i_iplus = AdjointRepresentation( inv(T_iplus_i) );

        model.I{i} = X_i_iplus'*I_i*X_i_iplus;
    end
end

function index = findString(list, str)
    index = 0;
    for i = 1:length(list)
        if strcmp(list{i},str)
            index = i;
            break
        end
    end
end

function X = AdjointRepresentation(T)
    R = T(1:3,1:3);
    p = T(1:3,4);
    X = [R zeros(3); skew(p)*R R];
end