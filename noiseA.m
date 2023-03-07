function NA = noiseA(Es,bps, Eb_N0_lin,r)

Eb = Es / (bps*r);
NP = Eb ./ Eb_N0_lin;
NA = sqrt(NP);