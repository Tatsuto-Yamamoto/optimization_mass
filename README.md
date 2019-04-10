# optimization_mass

<mass_function>
パラメーターは7個で、a,b,c,d,f,g,h
x は (a-b) を意味する
u, d, t, b, c, s, ele, muon, tau
8種類の粒子について、

m_u/m_t = ups = ((x)*Exp[(x)/3]*Sinh[c/3 - c*d] + (x)*Sinh[c*d] - 
      c*Exp[(x)/3]*Cosh[c/3 - c*d] + 
      c*Cosh[c*d])/((x)*Exp[(x)/3] Sinh[c - c*d] + (-x)*
       Sinh[2*c/3 - c*d] - c*Exp[(x)/3] Cosh[c - c*d] + 
      c*Cosh[2*c/3 - c*d])/(2.2/173000) // Simplify

      → PDFのp22,式(60)
      /(2.2/173000)は、実験値で得られたm_u/m_tである。
      
      charms = ((x)*Exp[(x)/3]*Sinh[2 c/3 - c*d] + (-x)*Sinh[c/3 - c*d] - 
     c*Exp[(x)/3]*Cosh[2 c/3 - c*d] + 
     c*Cosh[c/3 - c*d])/((x)*Exp[(x)/3] Sinh[c - c*d] + (-x)*
      Sinh[2*c/3 - c*d] - c*Exp[(x)/3] Cosh[c - c*d] + 
     c*Cosh[2*c/3 - c*d])/(1275/173000)
     
     downs = ((a - b)^2 - c^2)/((a - f)^2 - 
     c^2) Sqrt[(f*(-Exp[-2*b/3] + 1))/(b*(-Exp[-2*f/3] + 
        1))] ((a - f)*Exp[(a - f)/3]*Sinh[c/3 - c*d] + (a - f)*
       Sinh[c*d] - c*Exp[(a - f)/3]*Cosh[c/3 - c*d] + 
      c*Cosh[-c*d])/((a - b)*Exp[(a - b)/3] Sinh[c - c*d] + (-a + b)*
       Sinh[2*c/3 - c*d] - c*Exp[(a - b)/3] Cosh[c - c*d] + 
      c*Cosh[2*c/3 - c*d])/(4.7/173000)
      
      stranges = ((a - b)^2 - c^2)/((a - f)^2 - 
     c^2) Sqrt[(f*(-Exp[-2*b/3] + 1))/(b*(-Exp[-2*f/3] + 
        1))] ((a - f)*Exp[(a - f)/3]*Sinh[2 c/3 - c*d] + (-a + f)*
       Sinh[c/3 - c*d] - c*Exp[(a - f)/3]*Cosh[2 c/3 - c*d] + 
      c*Cosh[c/3 - c*d])/((a - b)*
       Exp[(a - b)/3] Sinh[c - c*d] + (-a + b)*Sinh[2*c/3 - c*d] - 
      c*Exp[(a - b)/3] Cosh[c - c*d] + c*Cosh[2*c/3 - c*d])/(95/
     173000)
     
     botoms = ((a - b)^2 - c^2)/((a - f)^2 - 
     c^2) Sqrt[(f*(-Exp[-2*b/3] + 1))/(b*(-Exp[-2*f/3] + 
        1))] ((a - f)*Exp[(a - f)/3]*Sinh[c - c*d] + (-a + f)*
       Sinh[2 c/3 - c*d] - c*Exp[(a - f)/3]*Cosh[c - c*d] + 
      c*Cosh[2 c/3 - c*d])/((a - b)*
       Exp[(a - b)/3] Sinh[c - c*d] + (-a + b)*Sinh[2*c/3 - c*d] - 
      c*Exp[(a - b)/3] Cosh[c - c*d] + c*Cosh[2*c/3 - c*d])/(41.8/
     1730)
     
     electrons = ((a - b)^2 - c^2)/((g - h)^2 - c^2)*
  Sqrt[(g*(Exp[2*a/3] - 1))/(a*(Exp[2*g/3] - 1))]*
  Sqrt[(h*(-Exp[-2*b/3] + 1))/(b*(-Exp[-2*h/3] + 1))] ((g - h)*
       Exp[(g - h)/3]*Sinh[c/3 - c*d] + (g - h)*Sinh[c*d] - 
      c*Exp[(g - h)/3]*Cosh[c/3 - c*d] + 
      c*Cosh[-c*d])/((a - b)*Exp[(a - b)/3] Sinh[c - c*d] + (-a + b)*
       Sinh[2*c/3 - c*d] - c*Exp[(a - b)/3] Cosh[c - c*d] + 
      c*Cosh[2*c/3 - c*d])/(0.510999/173000)
      
      muons = ((a - b)^2 - c^2)/((g - h)^2 - c^2)*
  Sqrt[(g*(Exp[2*a/3] - 1))/(a*(Exp[2*g/3] - 1))]*
  Sqrt[(h*(-Exp[-2*b/3] + 1))/(b*(-Exp[-2*h/3] + 1))] ((g - h)*
       Exp[(g - h)/3]*Sinh[2*c/3 - c*d] + (-g + h)*Sinh[c/3 - c*d] - 
      c*Exp[(g - h)/3]*Cosh[2*c/3 - c*d] + 
      c*Cosh[c/3 - c*d])/((a - b)*
       Exp[(a - b)/3] Sinh[c - c*d] + (-a + b)*Sinh[2*c/3 - c*d] - 
      c*Exp[(a - b)/3] Cosh[c - c*d] + c*Cosh[2*c/3 - c*d])/(105.65/
     173000)
     
     tauons = ((a - b)^2 - c^2)/((g - h)^2 - c^2)*
  Sqrt[(g*(Exp[2*a/3] - 1))/(a*(Exp[2*g/3] - 1))]*
  Sqrt[(h*(-Exp[-2*b/3] + 1))/(b*(-Exp[-2*h/3] + 1))] ((g - h)*
       Exp[(g - h)/3]*Sinh[c - c*d] + (-g + h)*Sinh[2*c/3 - c*d] - 
      c*Exp[(g - h)/3]*Cosh[c - c*d] + 
      c*Cosh[2*c/3 - c*d])/((a - b)*
       Exp[(a - b)/3] Sinh[c - c*d] + (-a + b)*Sinh[2*c/3 - c*d] - 
      c*Exp[(a - b)/3] Cosh[c - c*d] + c*Cosh[2*c/3 - c*d])/(1776.86/
     173000)

