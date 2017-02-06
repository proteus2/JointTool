function statics = div(constant, statics,alphar,rho)

% Input is (elm,alphar,rho)

frdof = constant.str.frdof;

V = sqrt(2/rho);

cd ..
cd('Kernel')

statics = structure(constant,statics,0);
cd('Panelaero')
statics = panaero(constant,statics,V,alphar,rho);
cd ../..
cd('Postprocessor')

ev = eig(statics.str.Ks(frdof,frdof),statics.aero.Ka(frdof,frdof));
[ev_s ev_i] = sort(ev,'ascend');

qdiv = min(ev_s);
Vdiv = sqrt(2*qdiv/rho);