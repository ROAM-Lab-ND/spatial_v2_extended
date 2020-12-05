function [angs1, angs2] = rotToFixedAngles(R,seq)

    [angs1,angs2] = rotToEulerAngles(R, seq([3 2 1]));
    angs1 = angs1([3 2 1]);
    angs2 = angs2([3 2 1]);

end

