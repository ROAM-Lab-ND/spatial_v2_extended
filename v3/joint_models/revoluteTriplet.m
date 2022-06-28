classdef revoluteTriplet
    %revoluteJointWithRotor 
    
    properties
        nv = 3
        nq = 3
        axis
        XtreeInternal
        output_body = 3
        bodies = 3
    end
    
    methods        
        function [Xup, S, Sd, v] = kinematics(obj, Xtree, q, qdot, vp)
            
            
            [XJ1, S1 ] = jcalc( ['R' obj.axis{1}], q(1) );
            [XJ2, S2 ] = jcalc( ['R' obj.axis{2}], q(2) );
            [XJ3, S3 ] = jcalc( ['R' obj.axis{3}], q(3) );
            
            X1p = XJ1*Xtree;
            X21 = XJ2*obj.XtreeInternal(1:6,:);
            X2p = X21*X1p;
            X32 = XJ3*obj.XtreeInternal(7:12,:);
            X3p = X32*X2p;
            
            Xup = [X1p;X2p;X3p];
            
            S  = [
                  S1         zeros(6,1) zeros(6,1);
                  X21*S1     S2         zeros(6,1);
                  X32*X21*S1 X32*S2     S3
                ];
            
            if nargout > 2
                v  = Xup*vp + S*qdot;
                v1 = v(1:6);
                v2 = v(7:12);
                v3 = v(13:18);
                
                S1d = crm(v1)*S1;
                S2d = crm(v2)*S2;
                S3d = crm(v3)*S3;
                
                Sd = [
                      S1d         zeros(6,1)  zeros(6,1);
                      X21*S1d     S2d         zeros(6,1);
                      X32*X21*S1d X32*S2d     S3d
                     ];
            end
        end
        
    end
end

