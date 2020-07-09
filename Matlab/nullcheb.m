   % CHEBYSHEV LOWPASS FILTER with SPECIFIED NULL

   % filter length
   N  = 31;
   N2 = N - 2;

   % set band-edges and weighting
   wp = 0.26*pi;
   ws = 0.34*pi;
   wo = 0.30*pi;
   Kp = 1;
   Ks = 4;

   % set null
   wn = 0.59*pi;

   L  = 1000;
   w  = [0:L]*pi/L;
   W  = Kp*(w<=wp) + Ks*(w>=ws);
   D  = (w<=wo);

   A1 = 2*(cos(w)-cos(wn));
   W2 = W.*abs(A1);
   D2 = D./A1;

   h2 = fircheb(N2,D2,W2);
   h1 = [1 -2*cos(wn) 1];
   h  = conv(h2,h1);
