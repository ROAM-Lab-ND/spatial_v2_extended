printHeader('CMM'); 
UnitTest_CMM

printHeader('Coriolis'); 
UnitTest_Coriolis

printHeader('Derivatives'); 
UnitTest_Derivatives

printHeader('Derivatives (Floating)'); 
UnitTest_Derivatives_FB

printHeader('Main Dynamics'); 
UnitTest_Dynamics

printHeader('Orientation'); 
UnitTest_Orientation

printHeader('Orientation Rates'); 
UnitTest_OrientationRates

function printHeader(st)
    fprintf('************************************\n');
    fprintf('%s\n',st);
    fprintf('************************************\n');
end