


zero4={
{ 0, 0, 0, 0},
{ 0, 0, 0, 0},
{ 0, 0, 0, 0},
{ 0, 0, 0, 0}};


(* 4D Dirac represenation *)
gamma4[0]= {
{ 1, 0, 0, 0},
{ 0, 1, 0, 0},
{ 0, 0,-1, 0},
{ 0, 0, 0,-1}};

gamma4[1]= {
{ 0, 0, 0, 1},
{ 0, 0, 1, 0},
{ 0,-1, 0, 0},
{-1, 0, 0, 0}};

gamma4[2]= {
{ 0, 0, 0,-I},
{ 0, 0, I, 0},
{ 0, I, 0, 0},
{-I, 0, 0, 0}};

gamma4[3]= {
{ 0, 0, 1, 0},
{ 0, 0, 0,-1},
{-1, 0, 0, 0},
{ 0, 1, 0, 0}};

gamma4[5]= {
{ 0, 0, 1, 0},
{ 0, 0, 0, 1},
{ 1, 0, 0, 0},
{ 0, 1, 0, 0}};


gamma4Hat = gamma4[5];
gamma4Five= gamma4[5];

testGamma5 = I*gamma4[0].gamma4[1].gamma4[2].gamma4[3] // Factor;
% // MatrixForm
test = testGamma5 - gamma4[5] // Factor
% // MatrixForm


test= gamma4[0].gamma4Five + gamma4Five.gamma4[0] // Factor;
% // MatrixForm
test= gamma4[1].gamma4Five + gamma4Five.gamma4[1] // Factor;
% // MatrixForm
test= gamma4[2].gamma4Five + gamma4Five.gamma4[2] // Factor;
% // MatrixForm
test= gamma4[3].gamma4Five + gamma4Five.gamma4[3] // Factor;
% // MatrixForm




omega4Minus = 1/2( IdentityMatrix[4] - gamma4[5]);
omega4Plus  = 1/2( IdentityMatrix[4] + gamma4[5]);


(* checking projector properties *)
test = omega4Minus + omega4PlusD // Factor;
% // MatrixForm

test = omega4Minus.omega4Minus - omega4Minus// Factor;
% // MatrixForm

test = omega4Plus.omega4Plus - omega4Plus// Factor;
% // MatrixForm

test = omega4Plus.omega4Minus// Factor;
% // MatrixForm

test = omega4Minus.omega4Plus// Factor;
% // MatrixForm




spi4 = {sp1,sp2,sp3,sp4};

spb2 = spi4.(v1*gamma4[0] - v2*gamma4[1] - v3*gamma4[2] - v4*gamma4[3])  // Simplify;
% // MatrixForm
% // FortranForm

spi2 = (v1*gamma4[0] - v2*gamma4[1] - v3*gamma4[2] - v4*gamma4[3]).spi4  // Simplify;
% // MatrixForm
% // FortranForm


spi4Left = omega4Minus.spi4
% // MatrixForm
% // FortranForm

spi4Right = omega4Plus.spi4
% // MatrixForm
% // FortranForm














(* 6D Dirac representation (8x8) *)

gamma6[0] = { {gamma4[0],zero4},{zero4,gamma4[0]} };
% // MatrixForm
gamma6[0] = ArrayFlatten[gamma6[0]];
% // MatrixForm

gamma6[1] = { {gamma4[1],zero4},{zero4,gamma4[1]} };
% // MatrixForm
gamma6[1] = ArrayFlatten[gamma6[1]];
% // MatrixForm

gamma6[2] = { {gamma4[2],zero4},{zero4,gamma4[2]} };
% // MatrixForm
gamma6[2] = ArrayFlatten[gamma6[2]];
% // MatrixForm

gamma6[3] = { {gamma4[3],zero4},{zero4,gamma4[3]} };
% // MatrixForm
gamma6[3] = ArrayFlatten[gamma6[3]];
% // MatrixForm

gamma6[4] = { {zero4,gamma4[5]},{-gamma4[5],zero4} };
% // MatrixForm
gamma6[4] = ArrayFlatten[gamma6[4]];
% // MatrixForm

gamma6[5] = { {zero4,I*gamma4[5]},{I*gamma4[5],zero4} };
% // MatrixForm
gamma6[5] = ArrayFlatten[gamma6[5]];
% // MatrixForm


(* this is 'GammaHat' in D=6 *)
gamma6[6] = - gamma6[0].gamma6[1].gamma6[2].gamma6[3].gamma6[4].gamma6[5]  // Factor;
% // MatrixForm
gamma6Hat = gamma6[6];

(* this is 'Gamma5' in D=6 *)
gamma6Five= { {gamma4Five,zero4},{zero4,gamma4Five} };
% // MatrixForm
gamma6Five = ArrayFlatten[gamma6Five];
% // MatrixForm




test= gamma6[0].gamma6Five + gamma6Five.gamma6[0] // Factor;
% // MatrixForm
test= gamma6[1].gamma6Five + gamma6Five.gamma6[1] // Factor;
% // MatrixForm
test= gamma6[2].gamma6Five + gamma6Five.gamma6[2] // Factor;
% // MatrixForm
test= gamma6[3].gamma6Five + gamma6Five.gamma6[3] // Factor;
% // MatrixForm
test= gamma6[4].gamma6Five - gamma6Five.gamma6[4] // Factor;
% // MatrixForm
test= gamma6[5].gamma6Five - gamma6Five.gamma6[5] // Factor;
% // MatrixForm




omega6Minus = 1/2( IdentityMatrix[8] - gamma6Five);
% // MatrixForm

omega6Plus  = 1/2( IdentityMatrix[8] + gamma6Five);
% // MatrixForm

test = omega6Minus + omega6Plus // Factor;
% // MatrixForm

test = omega6Minus.omega6Minus - omega6Minus// Factor;
% // MatrixForm

test = omega6Plus.omega6Plus - omega6Plus// Factor;
% // MatrixForm

test = omega6Plus.omega6Minus// Factor;
% // MatrixForm

test = omega6Minus.omega6Plus// Factor;
% // MatrixForm



spi6 = {sp1,sp2,sp3,sp4,sp5,sp6,sp7,sp8};

spb2 = spi6.(v1*gamma6[0] - v2*gamma6[1] - v3*gamma6[2] - v4*gamma6[3] - v5*gamma6[4]- v6*gamma6[5])  // Simplify;
% // MatrixForm
% // FortranForm

spi2 = (v1*gamma6[0] - v2*gamma6[1] - v3*gamma6[2] - v4*gamma6[3] - v5*gamma6[4]- v6*gamma6[5]).spi6  // Simplify;
% // MatrixForm
% // FortranForm


spi6Left = omega6Minus.spi6  // Factor;
% // MatrixForm
% // FortranForm

spi6Right = omega6Plus.spi6  // Factor;
% // MatrixForm
% // FortranForm






(* 8D Dirac representation (16x16) *)

zero6={
{ 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0}};
% // MatrixForm


gamma8[0] = { {gamma6[0],zero6},{zero6,gamma6[0]} };
% // MatrixForm
gamma8[0] = ArrayFlatten[gamma8[0]];
% // MatrixForm


gamma8[1] = { {gamma6[1],zero6},{zero6,gamma6[1]} };
% // MatrixForm
gamma8[1] = ArrayFlatten[gamma8[1]];
% // MatrixForm


gamma8[2] = { {gamma6[2],zero6},{zero6,gamma6[2]} };
% // MatrixForm
gamma8[2] = ArrayFlatten[gamma8[2]];
% // MatrixForm


gamma8[3] = { {gamma6[3],zero6},{zero6,gamma6[3]} };
% // MatrixForm
gamma8[3] = ArrayFlatten[gamma8[3]];
% // MatrixForm


gamma8[4] = { {gamma6[4],zero6},{zero6,gamma6[4]} };
% // MatrixForm
gamma8[4] = ArrayFlatten[gamma8[4]];
% // MatrixForm


gamma8[5] = { {gamma6[5],zero6},{zero6,gamma6[5]} };
% // MatrixForm
gamma8[5] = ArrayFlatten[gamma8[5]];
% // MatrixForm


gamma8[6] = { {zero6,gamma6Hat},{-gamma6Hat,zero6} };
% // MatrixForm
gamma8[6] = ArrayFlatten[gamma8[6]];
% // MatrixForm


gamma8[7] = { {zero6,I*gamma6Hat},{I*gamma6Hat,zero6} };
% // MatrixForm
gamma8[7] = ArrayFlatten[gamma8[7]];
% // MatrixForm


(* this is GammaHat in D=8 *)
gamma8[8] = -I*gamma8[0].gamma8[1].gamma8[2].gamma8[3].gamma8[4].gamma8[5].gamma8[6].gamma8[7]  // Factor;
% // MatrixForm

(* this is 'Gamma5' in D=8 *)
gamma8Five= { {gamma6Five,zero6},{zero6,gamma6Five} };
% // MatrixForm
gamma8Five = ArrayFlatten[gamma8Five];
% // MatrixForm



test= gamma8[0].gamma8Five + gamma8Five.gamma8[0] // Factor;
% // MatrixForm
test= gamma8[1].gamma8Five + gamma8Five.gamma8[1] // Factor;
% // MatrixForm
test= gamma8[2].gamma8Five + gamma8Five.gamma8[2] // Factor;
% // MatrixForm
test= gamma8[3].gamma8Five + gamma8Five.gamma8[3] // Factor;
% // MatrixForm
test= gamma8[4].gamma8Five - gamma8Five.gamma8[4] // Factor;
% // MatrixForm
test= gamma8[5].gamma8Five - gamma8Five.gamma8[5] // Factor;
% // MatrixForm
test= gamma8[6].gamma8Five - gamma8Five.gamma8[6] // Factor;
% // MatrixForm
test= gamma8[7].gamma8Five - gamma8Five.gamma8[7] // Factor;
% // MatrixForm




omega8Minus = 1/2( IdentityMatrix[16] - gamma8Five);
% // MatrixForm
omega8Plus  = 1/2( IdentityMatrix[16] + gamma8Five);
% // MatrixForm

test = omega8Minus + omega8Plus // Factor;
% // MatrixForm

test = omega8Minus.omega8Minus - omega8Minus// Factor;
% // MatrixForm

test = omega8Plus.omega8Plus - omega8Plus// Factor;
% // MatrixForm

test = omega8Plus.omega8Minus// Factor;
% // MatrixForm

test = omega8Minus.omega8Plus// Factor;
% // MatrixForm




spi8 = {sp1,sp2,sp3,sp4,sp5,sp6,sp7,sp8,sp9,sp10,sp11,sp12,sp13,sp14,sp15,sp16};

spb2 = spi8.(v1*gamma8[0] - v2*gamma8[1] - v3*gamma8[2] - v4*gamma8[3] - v5*gamma8[4] - v6*gamma8[5] - v7*gamma8[6]- v8*gamma8[7])  // Simplify;
% // MatrixForm
% // FortranForm

spi2 = (v1*gamma8[0] - v2*gamma8[1] - v3*gamma8[2] - v4*gamma8[3] - v5*gamma8[4] - v6*gamma8[5] - v7*gamma8[6]- v8*gamma8[7]).spi8  // Simplify;
% // MatrixForm
% // FortranForm


spi8Left = omega8Minus.spi8
% // MatrixForm
% // FortranForm

spi8Right = omega8Plus.spi8
% // MatrixForm
% // FortranForm










(* ---------------------------------------------------- *)



(* standard 4D Weyl representation, see e.g. HELAS or Wikipedia... *)
gammaWstd[0]= {
{ 0, 0, 1, 0},
{ 0, 0, 0, 1},
{ 1, 0, 0, 0},
{ 0, 1, 0, 0}};

gammaWstd[1]= {
{ 0, 0, 0, 1},
{ 0, 0, 1, 0},
{ 0,-1, 0, 0},
{-1, 0, 0, 0}};

gammaWstd[2]= {
{ 0, 0, 0,-I},
{ 0, 0, I, 0},
{ 0, I, 0, 0},
{-I, 0, 0, 0}};

gammaWstd[3]= {
{ 0, 0, 1, 0},
{ 0, 0, 0,-1},
{-1, 0, 0, 0},
{ 0, 1, 0, 0}};

gammaWstd[5]= {
{-1, 0, 0, 0},
{ 0,-1, 0, 0},
{ 0, 0, 1, 0},
{ 0, 0, 0, 1}};


(*  conversion to K.Ellis' Weyl repr.  *)
gammaW[0]=+gammaWstd[0];
gammaW[1]=-gammaWstd[1];
gammaW[2]=-gammaWstd[2];
gammaW[3]=-gammaWstd[3];
gammaW[5]=-gammaWstd[5];

G= {
{ 1, 0, 0, 0},
{ 0,-1, 0, 1},
{ 0, 0, -1,0},
{ 0, 0, 0,-1}};



(* unitary matrix to convert Dirac to K.Ellis' Weyl repr.  *)
U= 1/Sqrt[2]*{
{ 1, 0, 1, 0},
{ 0, 1, 0, 1},
{ 1, 0, -1, 0},
{ 0, 1, 0, -1}};
UT= Transpose[U];
(* check unitarity *)
check = U.UT
check = U.U
check = UT.UT
check = Det[U]


(* check the transformation *)
check = U.gamma4[0].UT - gammaW[0]
check = U.gamma4[1].UT - gammaW[1]
check = U.gamma4[2].UT - gammaW[2]
check = U.gamma4[3].UT - gammaW[3]
check = U.gamma4[5].UT - gammaW[5]




(* this is a vertex ubar.gamma^mu.v *)
ubar={ub1,ub2,ub3,ub4};
v={v1,v2,v3,v4};
ubar.gammaW[0].v
ubar.gammaW[1].v
ubar.gammaW[2].v
ubar.gammaW[3].v


ubar.gammaWstd[0].v  // FortranForm
ubar.gammaWstd[1].v  // FortranForm
ubar.gammaWstd[2].v  // FortranForm
ubar.gammaWstd[3].v  // FortranForm


ubar.gamma4[0].v
ubar.gamma4[1].v
ubar.gamma4[2].v
ubar.gamma4[3].v




(* check Clifford relation *)
Do[ Do[ Print[ mu,nu, "   ", gammaW[mu].gammaW[nu] + gammaW[nu].gammaW[mu] ] ,{nu,0,3}]  ,{mu,0,3}]


omegaMinusW = 1/2( IdentityMatrix[4] - gammaW[5]);
omegaPlusW  = 1/2( IdentityMatrix[4] + gammaW[5]);

ubar = {sp1,sp2,sp3,sp4};

ChiPL = ubar.omegaPlusW  // Simplify
% // FortranForm;

ChiML = ubar.omegaMinusW  // Simplify
% // FortranForm;

ChiPR = omegaPlusW.ubar  // Simplify
% // FortranForm;

ChiMR = omegaMinusW.ubar  // Simplify
% // FortranForm;

