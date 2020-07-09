   % CHEBYSHEV NOTCH FILTER 

   % filter length
   N  = 51;
   N2 = N - 4;

   % set band-edges 
   w0 = 0.55*pi;
   w1 = 0.65*pi;

   % set notch frequency
   wn = 0.60*pi;

   L  = 1000;
   w  = [0:L]*pi/L;
   W  = (w<=w0) + (w>=w1);
   D  = ones(size(w));

   A1 = 4*(cos(w)-cos(wn)).^2;
   W2 = W.*abs(A1);
   D2 = D./A1;

   h1 = conv([1 -2*cos(wn) 1],[1 -2*cos(wn) 1]);
   h2 = fircheb(N2,D2,W2);
   h  = conv(h2,h1);
