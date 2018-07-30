function result=phii9(x)
global wgaine wgaini Receptor;
g=0.087;
I=177.;
c=615; 
y=(c.*x-I).*(1+Receptor*wgaini);
 if y~=0
  result = y./(1-exp(-g.*y));
 else
  result=0;
 end
end
