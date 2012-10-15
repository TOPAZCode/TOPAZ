(* {DipAlpha=00000, DipAlpha=00002}, dl=dilepton, fh=fully hadronic *)

  Vidl = {-0.27611456 10^+01,-0.19007229 10^+02};
errVidl = {0.83551649 10^-03,0.57564967 10^-02};

Redl = {0.52802892 10^+00,0.16778684 10^+02};
errRedl = {0.21616606 10^-03,0.80858866 10^-02};

sumdl = Redl+Vidl;
errdl = Sqrt[errVidl^2+errRedl^2];


Vifh = {-0.56581680 10^+02,-0.13174657 10^+04};
errVifh = {0.17139635 10^-01,0.39881709 10^+00};

Refh = {0.83259169 10^+01,0.12700363 10^+04};
errRefh = {0.29177137 10^-01,0.51060443 10^+00}

sumfh = Refh+Vifh;
errfh = Sqrt[errVifh^2+errRefh^2];

(******************************************************)

  resLOnoDK = 0.10366486 10^4; 
errLOnoDK = 0.17903682 10^+01;

BdNLOtop = (1.337450-1.465331)/1.465331;

preddl = 2 resLOnoDK /81 * BdNLOtop;


BdNLOwjj = 2/3 * 0.109517038/Pi;

predfh = 2 resLOnoDK ( 4/9 BdNLOtop + 2/3 BdNLOwjj);
