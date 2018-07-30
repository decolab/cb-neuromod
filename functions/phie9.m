function result=phie9(x)
global wgaine wgaini Receptor;
g=0.16;
I=125.;
c=310;
y=(c.*x-I).*(1+Receptor*wgaine);
 if y~=0
  result = y./(1-exp(-g.*y));
 else
  result=0;
 end
end
