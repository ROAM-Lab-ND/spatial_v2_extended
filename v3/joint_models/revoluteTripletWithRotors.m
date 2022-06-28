classdef revoluteTripletWithRotors
    %revoluteJointWithRotor 
    
    properties
        nv = 3
        nq = 3
        jointAxis
        rotorAxis
        
        gearRatio
        beltRatio
        
        XtreeInternal
        output_body = 6
        bodies = 6
    end
    
    methods        
        function [Xup, S, Sd, v] = kinematics(obj, Xtree, q, qdot, vp)
            % Assuming that the generalized coordaintes are the realtive
            % angles. All the challenge comes from rotors.
            
            [XJ1 , S1 ] = jcalc( ['R' obj.jointAxis{1}] , q(1) );
            [XJ2 , S2 ] = jcalc( ['R' obj.jointAxis{2}] , q(2) );
            [XJ3 , S3 ] = jcalc( ['R' obj.jointAxis{3}] , q(3) );
           
            gr1 = obj.gearRatio{1};
            gr2 = obj.gearRatio{2};
            gr3 = obj.gearRatio{3};
            b1  = obj.beltRatio{1};
            b2  = obj.beltRatio{2};
            b3  = obj.beltRatio{3};
            
            n1 = gr1*b1;
            n2 = gr2*b2;
            n3 = gr3*b3;
            
            qr1 = n1*q(1);
            qr2 = n2*q(2) + gr2*b1*q(1);
            qr3 = n3*q(3) - gr3*b1*q(1) - gr3*b2*q(2);

            [XR1 , Sr1 ] = jcalc( ['R' obj.rotorAxis{1}], qr1 );
            [XR2 , Sr2 ] = jcalc( ['R' obj.rotorAxis{2}], qr2 );
            [XR3 , Sr3 ] = jcalc( ['R' obj.rotorAxis{3}], qr3 );

            X1p  = XJ1 * Xtree(1:6,:);
            Xr1p = XR1 * Xtree(7:12,:);
            X21  = XJ2 * obj.XtreeInternal(1:6,:);
            X2p  = X21 * X1p;
            Xr2p = XR2 * Xtree(13:18,:);
            Xr3p = XR3 * Xtree(19:24,:);
            X32  = XJ3 * obj.XtreeInternal(7:12,:);
            X3p  = X32 * X2p;

            Xup = [X1p ; Xr1p ; X2p; Xr2p ; Xr3p ; X3p];

            
            z = zeros(6,1);
            S  = [S1           z          z; 
                  Sr1*n1       z          z;
                  X21*S1       S2         z;
                  Sr2*gr2*b1   Sr2*n2     z;
                  Sr3*gr3*b1   Sr3*gr3*b2  Sr3*n3;
                  X32*X21*S1   X32*S2     S3
                  ];
              
            if nargout > 2
                v  = Xup*vp + S*qdot;
                v1 = v(1:6);
                vr1= v(7:12);
                v2 = v(13:18);
                vr2= v(19:24);
                vr3= v(25:30);
                v3 = v(31:36);
                
                S1d = crm(v1)*S1;
                S2d = crm(v2)*S2;
                S3d = crm(v3)*S3;
                Sr1d = crm(vr1)*Sr1;
                Sr2d = crm(vr2)*Sr2;
                Sr3d = crm(vr3)*Sr3;

                Sd  = [S1d           z          z; 
                       Sr1d*n1       z          z;
                       X21*S1d       S2d         z;
                       Sr2d*gr2*b1   Sr2d*n2     z;
                       Sr3d*gr3*b1   Sr3d*gr3*b2  Sr3d*n3;
                       X32*X21*S1d   X32*S2d     S3d
                      ];
            end
        end
        
    end
end

