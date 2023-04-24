Sets
i               /1*2/
*i1_76-16-4    i2_420-46-2
j               /1*102/
l               /1*3/
;
Alias(j,k,kk)
Alias(i,ii)


$call gdxxrw.exe R116R143a.xlsx par=psigma rng=psigma!a1:cY3 rdim=1 cdim=1  
Parameter psigma(i,j);
$GDXIN R116R143a.gdx
$LOAD psigma
$GDXIN
option decimals = 8;
display psigma;

$call gdxxrw.exe R116R143a.xlsx par=psigma_pure rng=psigma!a6:cY8 rdim=1 cdim=1  
Parameter psigma_pure(i,j);
$GDXIN R116R143a.gdx
$LOAD psigma_pure
$GDXIN
option decimals = 8;
display psigma_pure;

$call gdxxrw.exe R116R143a.xlsx par=square_sum_rps1_rps2 rng=square_sum_rps1_rps2!A1:CY103 rdim=1 cdim=1
Parameter square_sum_rps1_rps2(j,k);
$gdxin R116R143a.gdx
$load square_sum_rps1_rps2
$gdxin
display square_sum_rps1_rps2

$call gdxxrw.exe R116R143a.xlsx par=square_rps1_rps2 rng=square_rps1_rps2!A1:CY103 rdim=1 cdim=1
Parameter square_rps1_rps2(j,k);
$gdxin R116R143a.gdx
$load square_rps1_rps2
$gdxin
display square_rps1_rps2

$call gdxxrw.exe R116R143a.xlsx par=gammaold rng=gammaold!A1:b102 rdim=1 
Parameter gammaold(j);
$gdxin R116R143a.gdx
$load gammaold
$gdxin
display gammaold;

$call gdxxrw.exe R116R143a.xlsx par=chb rng=chb!A1:CY103 rdim=1 cdim=1
Parameter chb(j,k);
$gdxin R116R143a.gdx
$load chb
$gdxin
display chb;

Parameters
*provide perimental data for initialization
Tini                       /273.3/
Pini                       /1.8243/
x1                         /0.8754/
y1                         /0.8966/
*pure material parameters
Tc(i)                   /'1' 293.03,'2' 345.86/
Pc(i)                   /'1' 3.048,'2' 3.761/
omega(i)                /'1' 0.2566,'2' 0.2615/
Area(i)                 /'1' 119.2712848,'2' 96.78839793/
Vol(i)                  /'1' 108.7691499,'2' 81.96211169/
disp(i)                 /'1' 68.624425,'2' 78.04/
*COSMO parameters 
q0                      /79.53/
r0                      /66.69/
z_coordination          /10/
w                       /0.27027/
A_ES                    /6525.69/
B_ES                    /148590000/
c_OT_OT                 /932.31/
a_eff                   /7.25/
*other parameters
R                       /8.31433/
*kij parameters
xkij(i)                 /'1' 0.5,'2' 0.5/
Pkij                    /0.1/
kij                     /0.1600/
;
Table Twu(i,l)
           1                     2                  3
1   0.506824734000000   0.826828382000000   1.224304913000000
2   0.244999902000000   0.849065209000000   2.129778311000000
;
Variable
x(i)
y(i)
T
P
* the combinatorial part of ln(gamma_i)
qq(i)
rr(i)
ll(i)
phioverxl(i)
phioverxv(i)
thetaoverphil(i)
thetaoverphiv(i)

lngammacombl(i)

lngammacombv(i)
*the dispersion part of ln(gamma_i)
Adisp
lngamma_displ(i)
lngamma_dispv(i)
*the residual part of ln(gamma_i)
C_ES
Delta_W(j,k)
psigma_sl(j)
psigma_sv(j)
lngamma_sl(k)
lngamma_sv(k)
lngamma(i,j)
lngamma_residl(i)
lngamma_residv(i)
* cosmo gamma
gamma_cosmol(i)
gamma_cosmov(i)
*WS_COSMO_PR
alpha(i)
a_c(i)
a_i(i)
b_i(i)
b_i_a_i(i)
*kij
b_aRT_ij(i,ii)
g_E_resl
g_E_resv
bl
bv
al
av
AAl
AAv
BBl
BBv
CF_l
CF_v
QQQQl
QQQQv
DDl
DDv
dQQl(ii)
dQQv(i)
dDDl(i)
dDDv(i)
dbl(i)
dbv(i)
dal(i)
dav(i)
ln_phi_0l(i)
ln_phi_1l
ln_phi_2l(i)
ln_phi_3l
ln_phi_0v(i)
ln_phi_1v
ln_phi_2v(i)
ln_phi_3v
ln_phil(i)
ln_phiv(i)

obj
;
x.fx('1')=x1;
x.fx('2')=1-x1;
y.l('1')=y1;
y.l('2')=1-y1;
y.lo(i)=0;
T.fx=Tini;
P.l=Pini;
P.lo=0.5*Pini;
P.up=5*Pini;
*kij.l=0;
* the combinatorial part of ln(gamma_i)
qq.l(i)=Area(i)/q0;
rr.l(i)=Vol(i)/r0;
ll.l(i)=z_coordination/2*(rr.l(i)-qq.l(i))-(rr.l(i)-1);
phioverxl.l(i)=rr.l(i)/(x.l('1')*rr.l('1')+x.l('2')*rr.l('2'));
phioverxv.l(i)=rr.l(i)/(y.l('1')*rr.l('1')+y.l('2')*rr.l('2'));
thetaoverphil.l(i)=qq.l(i)/rr.l(i)/(x.l('1')*qq.l('1')+x.l('2')*qq.l('2'))*(x.l('1')*rr.l('1')+x.l('2')*rr.l('2'));
thetaoverphiv.l(i)=qq.l(i)/rr.l(i)/(y.l('1')*qq.l('1')+y.l('2')*qq.l('2'))*(y.l('1')*rr.l('1')+y.l('2')*rr.l('2'));
lngammacombl.l(i)=log(phioverxl.l(i))+z_coordination/2*qq.l(i)*log(thetaoverphil.l(i))+ll.l(i)-phioverxl.l(i)*(x.l('1')*ll.l('1')+x.l('2')*ll.l('2'));
lngammacombv.l(i)=log(phioverxv.l(i))+z_coordination/2*qq.l(i)*log(thetaoverphiv.l(i))+ll.l(i)-phioverxv.l(i)*(y.l('1')*ll.l('1')+y.l('2')*ll.l('2'));

*the dispersion part of ln(gamma_i)
Adisp.l=w*(0.5*sum(i,disp(i))-sqrt(disp('1')*disp('2')));
lngamma_displ.l('1')=Adisp.l*x.l('2')*x.l('2');
lngamma_displ.l('2')=Adisp.l*x.l('1')*x.l('1');
lngamma_dispv.l('1')=Adisp.l*y.l('2')*y.l('2');
lngamma_dispv.l('2')=Adisp.l*y.l('1')*y.l('1');
*the residual part of ln(gamma_i)
C_ES.l=A_ES+B_ES/T.l/T.l;
Delta_W.l(j,k)=-(C_ES.l*square_sum_rps1_rps2(j,k)-c_OT_OT*(chb(j,k)*square_rps1_rps2(j,k)))/R/4184/T.l;
psigma_sl.l(j)=sum(i,psigma_pure(i,j)*(x.l(i)*Area(i)))/sum(i,x.l(i)*Area(i));
psigma_sv.l(j)=sum(i,psigma_pure(i,j)*(y.l(i)*Area(i)))/sum(i,y.l(i)*Area(i));
lngamma_sl.l(k)=(1/(sum(j,(exp(Delta_W.l(j,k)))*(psigma_sl.l(j)*gammaold(j))))+gammaold(k))/2;
lngamma_sv.l(k)=(1/(sum(j,(exp(Delta_W.l(j,k)))*(psigma_sv.l(j)*gammaold(j))))+gammaold(k))/2;

lngamma.l(i,k)=(1/(sum(j,(exp(Delta_W.l(j,k)))*(psigma_pure(i,j)*gammaold(j))))+gammaold(k))/2;

lngamma_residl.l(i)=Area(i)/a_eff*sum(j,psigma_pure(i,j)*(lngamma_sl.l(j)-lngamma.l(i,j)));
lngamma_residv.l(i)=Area(i)/a_eff*sum(j,psigma_pure(i,j)*(lngamma_sv.l(j)-lngamma.l(i,j)));
* cosmo gamma
gamma_cosmol.l(i)=lngammacombl.l(i)+lngamma_displ.l(i)+lngamma_residl.l(i);
gamma_cosmov.l(i)=lngammacombv.l(i)+lngamma_dispv.l(i)+lngamma_residv.l(i);
*WS_COSMO_PR
alpha.l(i)=(T.l/Tc(i))**(Twu(i,'3')*(Twu(i,'2')-1))*exp(Twu(i,'1')*(1-(T.l/Tc(i))**(Twu(i,'3')*Twu(i,'2'))));
a_c.l(i)=0.45724*R*R*(Tc(i)**2)/Pc(i);
a_i.l(i)=alpha.l(i)*a_c.l(i);
b_i.l(i)=0.0778*R*Tc(i)/Pc(i);
b_i_a_i.l(i)=b_i.l(i)-a_i.l(i)/R/T.l;
*kij.l=0;
b_aRT_ij.l('1','1')=2*b_i_a_i.l('1')*1*0.5;
b_aRT_ij.l('2','1')=(b_i_a_i.l('1')+b_i_a_i.l('2'))*0.5*(1-kij);
b_aRT_ij.l('1','2')=(b_i_a_i.l('1')+b_i_a_i.l('2'))*0.5*(1-kij);
b_aRT_ij.l('2','2')=2*b_i_a_i.l('2')*1*0.5;
g_E_resl.l=R*T.l*sum(i,x.l(i)*gamma_cosmol.l(i));
g_E_resv.l=R*T.l*sum(i,y.l(i)*gamma_cosmov.l(i));

*
bl.l=sum(ii,(sum(i,x.l(i)*b_aRT_ij.l(i,ii)))*x.l(ii))/(1-sum(i,x.l(i)*(a_i.l(i)/b_i.l(i)/R/T.l))+g_E_resl.l/R/T.l/0.62323);
bv.l=sum(ii,(sum(i,y.l(i)*b_aRT_ij.l(i,ii)))*y.l(ii))/(1-sum(i,y.l(i)*(a_i.l(i)/b_i.l(i)/R/T.l))+g_E_resv.l/R/T.l/0.62323);
al.l=bl.l*sum(i,x.l(i)*(a_i.l(i)/b_i.l(i)))+g_E_resl.l/(-0.62323);
av.l=bv.l*sum(i,y.l(i)*(a_i.l(i)/b_i.l(i)))+g_E_resv.l/(-0.62323);
AAl.l=al.l*P.l/R/R/T.l/T.l;
AAv.l=av.l*P.l/R/R/T.l/T.l;
BBl.l=bl.l*P.l/R/T.l;
BBv.l=bv.l*P.l/R/T.l;

QQQQl.l=(x.l('1')*b_aRT_ij.l('1','1')+x.l('2')*b_aRT_ij.l('2','1'))*x.l('1')+(x.l('1')*b_aRT_ij.l('1','2')+x.l('2')*b_aRT_ij.l('2','2'))*x.l('2');
QQQQv.l=(y.l('1')*b_aRT_ij.l('1','1')+y.l('2')*b_aRT_ij.l('2','1'))*y.l('1')+(y.l('1')*b_aRT_ij.l('1','2')+y.l('2')*b_aRT_ij.l('2','2'))*y.l('2');
DDl.l=sum(i,x.l(i)*(a_i.l(i)/b_i.l(i)/R/T.l))-g_E_resl.l/R/T.l/0.62323;
DDv.l=sum(i,y.l(i)*(a_i.l(i)/b_i.l(i)/R/T.l))-g_E_resv.l/R/T.l/0.62323;
dQQl.l(i)=2*sum(ii,b_aRT_ij.l(i,ii)*x.l(i));
dQQv.l(i)=2*sum(ii,b_aRT_ij.l(i,ii)*y.l(i));
dDDl.l(i)=a_i.l(i)/b_i.l(i)/R/T.l-gamma_cosmol.l(i)/(-0.62323);
dDDv.l(i)=a_i.l(i)/b_i.l(i)/R/T.l-gamma_cosmov.l(i)/(-0.62323);
dbl.l(i)=1/(1-DDl.l)*dQQl.l(i)-QQQQl.l/(1-DDl.l)/(1-DDl.l)*(1-dDDl.l(i));
dbv.l(i)=1/(1-DDv.l)*dQQv.l(i)-QQQQv.l/(1-DDv.l)/(1-DDv.l)*(1-dDDv.l(i));
dal.l(i)=R*T.l*DDl.l*dbl.l(i)+R*T.l*bl.l*dDDl.l(i);
dav.l(i)=R*T.l*DDv.l*dbv.l(i)+R*T.l*bv.l*dDDv.l(i);
CF_v.l=0.7;
CF_l.l=1.01*BBl.l;
CF_l.lo=BBl.l+0.0000001;
CF_v.lo=BBv.l+0.0000001;
ln_phi_0l.l(i)=(dbl.l(i)/bl.l)*(CF_l.l-1);
ln_phi_1l.l=log(CF_l.l-BBl.l);
ln_phi_2l.l(i)=dal.l(i)/al.l-dbl.l(i)/bl.l;
ln_phi_3l.l=AAl.l/BBl.l/2.82483*log((2*CF_l.l-0.82843*BBl.l)/(2*CF_l.l+4.82843*BBl.l));
ln_phi_0v.l(i)=(dbv.l(i)/bv.l)*(CF_v.l-1);
ln_phi_1v.l=log(CF_v.l-BBv.l);
ln_phi_2v.l(i)=dav.l(i)/av.l-dbv.l(i)/bv.l;
ln_phi_3v.l=AAv.l/BBv.l/2.82483*log((2*CF_v.l-0.82843*BBv.l)/(2*CF_v.l+4.82843*BBv.l));
ln_phil.l(i)=ln_phi_0l.l(i)-ln_phi_1l.l+ln_phi_2l.l(i)*ln_phi_3l.l;
ln_phiv.l(i)=ln_phi_0v.l(i)-ln_phi_1v.l+ln_phi_2v.l(i)*ln_phi_3v.l;

Equations e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e18,e19,e20,e21,e22,e24,e25,e26,e27,e28,e29,e30,e31,e32,e33,e34,e35,
e36,e37,e38,e39,e40,e41,e42,e43,e44,e45,e46,e47,e48,e49,e50,e51,e52,e53,e54,e55,e56,e57,e571,e58,e581,e59,e60,e61,e62,e63,e64,e65,e66,e67,e68,e69,e70,
e71,e72,e73,e74,e75,e76,e;
*e75,e76,e77,e78,e79,e80,e81,e82,e83,e84,e85,e86,e87,e88,e89,e90,e91,e92,e93,e94,e95,e96,e97,e98,e99,e100,e101,e102,e103,e104,
*e105,e106,e107,e108,e109,e110,e111,e112,e113,e114,e115,e116,e117,e118,e119,
*the combinatorial part of ln(gamma_i)
* the combinatorial part of ln(gamma_i)
e1(i)..             qq(i)=e=Area(i)/q0;
e2(i)..             rr(i)=e=Vol(i)/r0;
e3(i)..             ll(i)=e=z_coordination/2*(rr(i)-qq(i))-(rr(i)-1);
e4(i)..             phioverxl(i)=e=rr(i)/(x('1')*rr('1')+x('2')*rr('2'));
e5(i)..             phioverxv(i)=e=rr(i)/(y('1')*rr('1')+y('2')*rr('2'));
e6(i)..             thetaoverphil(i)=e=qq(i)/rr(i)/(x('1')*qq('1')+x('2')*qq('2'))*(x('1')*rr('1')+x('2')*rr('2'));
e7(i)..             thetaoverphiv(i)=e=qq(i)/rr(i)/(y('1')*qq('1')+y('2')*qq('2'))*(y('1')*rr('1')+y('2')*rr('2'));
e8(i)..             lngammacombl(i)=e=log(phioverxl(i))+z_coordination/2*qq(i)*log(thetaoverphil(i))+ll(i)-phioverxl(i)*(x('1')*ll('1')+x('2')*ll('2'));
e9(i)..             lngammacombv(i)=e=log(phioverxv(i))+z_coordination/2*qq(i)*log(thetaoverphiv(i))+ll(i)-phioverxv(i)*(y('1')*ll('1')+y('2')*ll('2'));
  
*the dispersion part of ln(gamma_i)
e10..               Adisp=e=w*(0.5*sum(i,disp(i))-sqrt(disp('1')*disp('2')));
e11('1')..          lngamma_displ('1')=e=Adisp*x('2')*x('2');
e12('2')..          lngamma_displ('2')=e=Adisp*x('1')*x('1');
e13('1')..          lngamma_dispv('1')=e=Adisp*y('2')*y('2');
e14('2')..          lngamma_dispv('2')=e=Adisp*y('1')*y('1');
*the residual part of ln(gamma_i)
e15..               C_ES=e=A_ES+B_ES/T/T;
e16(j,k)..          Delta_W(j,k)=e=-(C_ES*square_sum_rps1_rps2(j,k)-c_OT_OT*(chb(j,k)*square_rps1_rps2(j,k)))/R/4184/T;
e18(j)..            psigma_sl(j)=e=sum(i,psigma_pure(i,j)*(x(i)*Area(i)))/sum(i,x(i)*Area(i));
e19(j)..            psigma_sv(j)=e=sum(i,psigma_pure(i,j)*(y(i)*Area(i)))/sum(i,y(i)*Area(i));  
e20(k)..            lngamma_sl(k)=e=(1/(sum(j,exp(Delta_W(j,k))*(psigma_sl(j)*gammaold(j))))+gammaold(k))/2;
e21(k)..            lngamma_sv(k)=e=(1/(sum(j,exp(Delta_W(j,k))*(psigma_sv(j)*gammaold(j))))+gammaold(k))/2;

e22(i,k)..          lngamma(i,k)=e=(1/(sum(j,exp(Delta_W(j,k))*(psigma_pure(i,j)*gammaold(j))))+gammaold(k))/2;

e24(i)..            lngamma_residl(i)=e=Area(i)/a_eff*sum(j,psigma_pure(i,j)*(lngamma_sl(j)-lngamma(i,j)));
e25(i)..            lngamma_residv(i)=e=Area(i)/a_eff*sum(j,psigma_pure(i,j)*(lngamma_sv(j)-lngamma(i,j)));
* cosmo gamma
e26(i)..            gamma_cosmol(i)=e=lngammacombl(i)+lngamma_displ(i)+lngamma_residl(i);
e27(i)..            gamma_cosmov(i)=e=lngammacombv(i)+lngamma_dispv(i)+lngamma_residv(i);
e28(i)..            alpha(i)=e=(T/Tc(i))**(Twu(i,'3')*(Twu(i,'2')-1))*exp(Twu(i,'1')*(1-(T/Tc(i))**(Twu(i,'3')*Twu(i,'2'))));
e29(i)..            a_c(i)=e=0.45724*R*R*(Tc(i)**2)/Pc(i);
e30(i)..            a_i(i)=e=alpha(i)*a_c(i);
e31(i)..            b_i(i)=e=0.0778*R*Tc(i)/Pc(i);
e32(i)..            b_i_a_i(i)=e=b_i(i)-a_i(i)/R/T;
e33(i,ii)..         b_aRT_ij('1','1')=e=2*b_i_a_i('1')*1*0.5;
e34(i,ii)..         b_aRT_ij('2','1')=e=(b_i_a_i('1')+b_i_a_i('2'))*0.5*(1-kij);
e35(i,ii)..         b_aRT_ij('1','2')=e=(b_i_a_i('1')+b_i_a_i('2'))*0.5*(1-kij);
e36(i,ii)..         b_aRT_ij('2','2')=e=2*b_i_a_i('2')*1*0.5;
e37..               g_E_resl=e=R*T*sum(i,x(i)*gamma_cosmol(i));
e38..               g_E_resv=e=R*T*sum(i,y(i)*gamma_cosmov(i));
e39..               bl=e=sum(ii,(sum(i,x(i)*b_aRT_ij(i,ii)))*x(ii))/(1-sum(i,x(i)*(a_i(i)/b_i(i)/R/T))+g_E_resl/R/T/0.62323);
e40..               bv=e=sum(ii,(sum(i,y(i)*b_aRT_ij(i,ii)))*y(ii))/(1-sum(i,y(i)*(a_i(i)/b_i(i)/R/T))+g_E_resv/R/T/0.62323);
e41..               al=e=bl*sum(i,x(i)*(a_i(i)/b_i(i)))+g_E_resl/(-0.62323);
e42..               av=e=bv*sum(i,x(i)*(a_i(i)/b_i(i)))+g_E_resv/(-0.62323);
e43..               AAl=e=al*P/R/R/T/T;
e44..               AAv=e=av*P/R/R/T/T;
e45..               BBl=e=bl*P/R/T;
e46..               BBv=e=bv*P/R/T;
e47..               0=e=CF_l*CF_l*CF_l-(1-BBl)*CF_l*CF_l-(2*BBl-AAl+3*BBl*BBl)*CF_l-(AAl*BBl-BBl*BBl-BBl*BBl*BBl);
e48..               0=l=3*CF_l*CF_l-2*(1-BBl)*CF_l-(2*BBl-AAl+3*BBl*BBl);                   
e49..               0=g=6*CF_l-2*(1-BBl);
e50..               0=e=CF_v*CF_v*CF_v-(1-BBv)*CF_v*CF_v-(2*BBv-AAv+3*BBv*BBv)*CF_v-(AAv*BBv-BBv*BBv-BBv*BBv*BBv);
e51..               0=l=3*CF_v*CF_v-2*(1-BBv)*CF_v-(2*BBv-AAv+3*BBv*BBv);
e52..               0=l=6*CF_v-2*(1-BBv);
e53..               QQQQl=e=(x('1')*b_aRT_ij('1','1')+x('2')*b_aRT_ij('2','1'))*x('1')+(x('1')*b_aRT_ij('1','2')+x('2')*b_aRT_ij('2','2'))*x('2');
e54..               QQQQv=e=(y('1')*b_aRT_ij('1','1')+y('2')*b_aRT_ij('2','1'))*y('1')+(y('1')*b_aRT_ij('1','2')+y('2')*b_aRT_ij('2','2'))*y('2');
e55..               DDl=e=sum(i,x(i)*(a_i(i)/b_i(i)/R/T))-g_E_resl/R/T/0.62323;
e56..               DDv=e=sum(i,y(i)*(a_i(i)/b_i(i)/R/T))-g_E_resv/R/T/0.62323;
e57(i)..            dQQl('1')=e=2*(b_aRT_ij('1','1')*x('1')+b_aRT_ij('1','2')*x('2'));
e571(i)..           dQQl('2')=e=2*(b_aRT_ij('2','1')*x('1')+b_aRT_ij('2','2')*x('2'));
e58(i)..            dQQv('1')=e=2*(b_aRT_ij('1','1')*y('1')+b_aRT_ij('1','2')*y('2'));
e581(i)..           dQQv('2')=e=2*(b_aRT_ij('2','1')*y('1')+b_aRT_ij('2','2')*y('2'));
e59(i)..            dDDl(i)=e=a_i(i)/b_i(i)/R/T-gamma_cosmol(i)/(-0.62323);
e60(i)..            dDDv(i)=e=a_i(i)/b_i(i)/R/T-gamma_cosmov(i)/(-0.62323);
e61(i)..            dbl(i)=e=1/(1-DDl)*dQQl(i)-QQQQl/(1-DDl)/(1-DDl)*(1-dDDl(i));
e62(i)..            dbv(i)=e=1/(1-DDv)*dQQv(i)-QQQQv/(1-DDv)/(1-DDv)*(1-dDDv(i));
e63(i)..            dal(i)=e=R*T*DDl*dbl(i)+R*T*bl*dDDl(i);
e64(i)..            dav(i)=e=R*T*DDv*dbv(i)+R*T*bv*dDDv(i);
e65(i)..            ln_phi_0l(i)=e=(dbl(i)/bl)*(CF_l-1);
e66..               ln_phi_1l=e=log(CF_l-BBl);
e67(i)..            ln_phi_2l(i)=e=dal(i)/al-dbl(i)/bl;
e68..               ln_phi_3l=e=AAl/BBl/2.82483*log((2*CF_l-0.82843*BBl)/(2*CF_l+4.82843*BBl));
e69(i)..            ln_phi_0v(i)=e=(dbv(i)/bv)*(CF_v-1);
e70..               ln_phi_1v=e=log(CF_v-BBv);
e71(i)..            ln_phi_2v(i)=e=dav(i)/av-dbv(i)/bv;
e72..               ln_phi_3v=e=AAv/BBv/2.82483*log((2*CF_v-0.82843*BBv)/(2*CF_v+4.82843*BBv));
e73(i)..            ln_phil(i)=e=ln_phi_0l(i)-ln_phi_1l+ln_phi_2l(i)*ln_phi_3l;
e74(i)..            ln_phiv(i)=e=ln_phi_0v(i)-ln_phi_1v+ln_phi_2v(i)*ln_phi_3v;

e75..               sum(i,y(i))=e=1;
e76(i)..            x(i)*exp(ln_phil(i))=e=y(i)*exp(ln_phiv(i));



e..                 obj=e=0.6171-y('1');


Model vlecosmo /all/;

*option NLP = conopt;
option MINLP = ANTIGONE;
option sys12 = 1;
option optcr  = 0.000001;
option iterlim= 1000000;
option reslim = 1000000;
option DOMLIM=999999;

SOLVE vlecosmo using NLP minimizing obj;

option decimals = 8;
display qq.l,rr.l,ll.l,phioverxl.l,phioverxv.l,thetaoverphil.l,thetaoverphiv.l,lngammacombl.l,lngammacombv.l,
Adisp.l,lngamma_displ.l,lngamma_dispv.l,C_ES.l,Delta_W.l,psigma_sl.l,psigma_sv.l,lngamma_sl.l,lngamma_sv.l,
lngamma.l,lngamma_residl.l,lngamma_residv.l,gamma_cosmol.l,gamma_cosmov.l,alpha.l,a_c.l,a_i.l,b_i.l,b_i_a_i.l,
b_aRT_ij.l,g_E_resl.l,g_E_resv.l,bl.l,bv.l,al.l,av.l,AAl.l,AAv.l,BBl.l,BBv.l,CF_l.l,CF_v.l,QQQQl.l,QQQQv.l,
DDl.l,DDv.l,dQQl.l,dQQv.l,dDDl.l,dDDv.l,dbl.l,dbv.l,dal.l,dav.l,ln_phi_0l.l,ln_phi_1l.l,ln_phi_2l.l,ln_phi_3l.l,
ln_phi_0v.l,ln_phi_1v.l,ln_phi_2v.l,ln_phi_3v.l,ln_phil.l,ln_phiv.l,x.l,y.l;
