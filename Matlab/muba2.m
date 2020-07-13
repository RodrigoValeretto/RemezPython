   N   = 105;
   ws1 = 0.20*pi;
   wp1 = 0.25*pi;
   wp2 = 0.63*pi;
   ws2 = 0.68*pi;
   Ks1 = 10;
   Kp  = 1;
   Ks2 = 1;

   wo = (ws1+wp1)/2;
   w1 = (wp1+ws2)/2;

   L = 4000;
   w = [0:L]*pi/L;

   W = Ks1*(w<=ws1) + Kp*((w>=wp1)&(w<=wp2)) + Ks2*(w>=ws2);
   D = (wp1<=w)&(w<=wp2);
   h = remez(N,D,W);
   n = [1:N];
   stem(n,h)
