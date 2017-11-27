(* Created with the Wolfram Language : www.wolfram.com *)
{{Derivative[0, 2][a][x, z] + 
   (Qzz[x, z]*(a[x, z] + (-1 + z)*(Derivative[0, 1][a][x, z] - 
         z^2*Qxz[x, z]*Derivative[1, 0][a][x, z]))*
      (-Derivative[0, 1][Qtt][x, z] + z^2*Qxz[x, z]*Derivative[1, 0][Qtt][x, 
         z]) + (Qtt[x, z]*Qzz[x, z]*(a[x, z] + 
        (-1 + z)*(Derivative[0, 1][a][x, z] - z^2*Qxz[x, z]*
           Derivative[1, 0][a][x, z]))*(Derivative[0, 1][Qxx][x, z] - 
        z^2*Qxz[x, z]*Derivative[1, 0][Qxx][x, z]))/Qxx[x, z] - 
     2*Qtt[x, z]*Qzz[x, z]*Derivative[0, 1][a][x, z]*
      (-2 + (-1 + z)*z^2*Derivative[1, 0][Qxz][x, z]) + 
     (Qtt[x, z]*Qzz[x, z]*(a[x, z] + (-1 + z)*(Derivative[0, 1][a][x, z] - 
          z^2*Qxz[x, z]*Derivative[1, 0][a][x, z]))*
       (Derivative[0, 1][Qyy][x, z] - z^2*Qxz[x, z]*Derivative[1, 0][Qyy][x, 
          z]))/Qyy[x, z] + Qtt[x, z]*(a[x, z] + 
       (-1 + z)*(Derivative[0, 1][a][x, z] - z^2*Qxz[x, z]*
          Derivative[1, 0][a][x, z]))*(-Derivative[0, 1][Qzz][x, z] + 
       z^2*Qxz[x, z]*Derivative[1, 0][Qzz][x, z]) - 
     2*z*Qtt[x, z]*Qzz[x, z]*(z*((-1 + z)*Derivative[0, 1][Qxz][x, z]*
          Derivative[1, 0][a][x, z] + a[x, z]*Derivative[1, 0][Qxz][x, z]) - 
       2*Qxz[x, z]*(Derivative[1, 0][a][x, z]*(1 - 2*z + 
           (-1 + z)*z^3*Derivative[1, 0][Qxz][x, z]) - 
         (-1 + z)*z*Derivative[1, 1][a][x, z]) - (-1 + z)*z^3*Qxz[x, z]^2*
        Derivative[2, 0][a][x, z]) + 
     (Qzz[x, z]*(Qtt[x, z]*Qyy[x, z]*Qzz[x, z]*Derivative[1, 0][a][x, z]*
         Derivative[1, 0][Qxx][x, z] + Qxx[x, z]*
         (-(Qtt[x, z]*Qzz[x, z]*Derivative[1, 0][a][x, z]*
            Derivative[1, 0][Qyy][x, z]) + Qyy[x, z]*
           (-(Qtt[x, z]*Derivative[1, 0][a][x, z]*Derivative[1, 0][Qzz][x, 
               z]) + Qzz[x, z]*(Derivative[1, 0][a][x, z]*Derivative[1, 0][
                 Qtt][x, z] - 2*Qtt[x, z]*Derivative[2, 0][a][x, z])))))/
      (P[z]*Qxx[x, z]^2*Qyy[x, z]))/(2*(-1 + z)*Qtt[x, z]*Qzz[x, z]), 
  Derivative[0, 2][Qtt][x, z] + 
   (Qzz[x, z]*(((-2 + z)^2*P[z]^2)/z^2 + (2*(-1 + z)*z^2*a[x, z]^2*P[z])/
       Qzz[x, z] + (P[z]^2*Qtt[x, z]*(2*(2 - 3*z + z^2)*Qyy[x, z]*Qzz[x, z] + 
         Qxx[x, z]*(2*(2 - 3*z + z^2)*Qzz[x, z] + z^2*Qyy[x, z]*
            (-1 + 2*(-1 + z)^3*z^3*Qxz[x, z]^2*Derivative[1][P][z]))))/
       (z^2*Qxx[x, z]*Qyy[x, z]*Qzz[x, z]) - 
      (2*(-1 + z)^2*P[z]^2*Derivative[0, 1][Qtt][x, z])/(z*Qxx[x, z]) - 
      (2*(-1 + z)^2*P[z]^2*Derivative[0, 1][Qtt][x, z])/(z*Qyy[x, z]) + 
      ((-1 + z)*(-2 + 3*z)*P[z]^2*Derivative[0, 1][Qtt][x, z])/
       (z*Qzz[x, z]) + (2*(-1 + z)^2*z^2*P[z]^3*Qxz[x, z]^2*
        (-((-2 + z)*Qtt[x, z]) + (-1 + z)*z*Derivative[0, 1][Qtt][x, z]))/
       Qzz[x, z] - (4*(-1 + z)^2*z^2*a[x, z]*P[z]*
        (-Derivative[0, 1][a][x, z] + z^2*Qxz[x, z]*Derivative[1, 0][a][x, 
           z]))/Qzz[x, z] - ((-1 + z)^2*(Qtt[x, z]*Qxx[x, z]*
          Derivative[1][P][z]^2 + Qzz[x, z]*
          (-(Qxx[x, z]*Derivative[1][P][z]^2) + 
           2*z^2*Derivative[1, 0][a][x, z]^2)))/(Qxx[x, z]*Qzz[x, z]) - 
      ((-1 + z)*P[z]^2*((-2 + z)*Qzz[x, z]*Derivative[0, 1][Qtt][x, z] + 
         2*(-1 + z)*z*(Derivative[0, 1][Qtt][x, z] - z^2*Qxz[x, z]*
             Derivative[1, 0][Qtt][x, z])^2))/(z*Qtt[x, z]*Qzz[x, z]) + 
      (2*(-1 + z)^2*z*P[z]^2*Qxz[x, z]*(-2*Derivative[1, 0][Qtt][x, z] - 
         2*z*Derivative[1, 1][Qtt][x, z] + z^3*Qxz[x, z]*
          Derivative[2, 0][Qtt][x, z]))/Qzz[x, z] + 
      ((-1 + z)*P[z]*(2*Qtt[x, z]^2*(-((-1 + z)*z*Qyy[x, z]*Qzz[x, z]*
             Derivative[1][P][z]) + Qxx[x, z]*
            (-((-1 + z)*z*Qzz[x, z]*Derivative[1][P][z]) + 
             Qyy[x, z]*(6*Qzz[x, z] + z^2*(Derivative[1][P][z] + 
                 (-1 + z)*Derivative[2][P][z])))) + z^2*Qyy[x, z]*Qzz[x, z]*
          ((-1 + z)*Qxx[x, z]*Derivative[1][P][z]*Derivative[0, 1][Qtt][x, 
             z] + 2*Derivative[1, 0][Qtt][x, z]^2) + z*Qtt[x, z]*Qyy[x, z]*
          (Qxx[x, z]*(-2*(-2 + z)*Qzz[x, z]*Derivative[1][P][z] + 
             (-1 + z)*z*(2*(-1 + z)*z^2*Derivative[0, 1][a][x, z]^2 + 
               Derivative[1][P][z]*Derivative[0, 1][Qtt][x, z] - 4*(-1 + z)*
                z^4*Qxz[x, z]*Derivative[0, 1][a][x, z]*Derivative[1, 0][a][
                 x, z] + 2*(-1 + z)*z^6*Qxz[x, z]^2*Derivative[1, 0][a][x, z]^
                 2)) - 2*z*Qzz[x, z]*Derivative[2, 0][Qtt][x, z])))/
       (z^2*Qtt[x, z]*Qxx[x, z]*Qyy[x, z]*Qzz[x, z])))/(2*(-1 + z)^2*P[z]^2), 
  Derivative[0, 2][Qzz][x, z] + 
   (2*(-1 + z)*z^4*a[x, z]^2*P[z]*Qtt[x, z]*Qxx[x, z]^2*Qyy[x, z]^2*
      Qzz[x, z]^2 - 4*(-1 + z)^2*z^4*a[x, z]*P[z]*Qtt[x, z]*Qxx[x, z]^2*
      Qyy[x, z]^2*Qzz[x, z]^2*(-Derivative[0, 1][a][x, z] + 
       z^2*Qxz[x, z]*Derivative[1, 0][a][x, z]) + (-1 + z)^2*z^2*Qtt[x, z]*
      Qxx[x, z]*Qyy[x, z]^2*Qzz[x, z]^2*
      (Qtt[x, z]*Qxx[x, z]*Derivative[1][P][z]^2 + 
       Qzz[x, z]*(-(Qxx[x, z]*Derivative[1][P][z]^2) + 
         2*z^2*Derivative[1, 0][a][x, z]^2)) + 2*(-1 + z)^2*z^4*P[z]^3*
      Qtt[x, z]^2*Qxx[x, z]^2*Qxz[x, z]*Qyy[x, z]^2*Qzz[x, z]*
      (4*(-1 + z)*z*Qzz[x, z]*Derivative[0, 1][Qxz][x, z] + 
       Qxz[x, z]*(-((-1 + z)*z*Derivative[0, 1][Qzz][x, z]) + 
         Qzz[x, z]*(-8 + 11*z - 4*(-1 + z)*z^3*Derivative[1, 0][Qxz][x, 
             z])) + 2*(-1 + z)*z^3*Qxz[x, z]^2*Derivative[1, 0][Qzz][x, z]) + 
     (-1 + z)*P[z]*Qxx[x, z]*Qyy[x, z]*Qzz[x, z]*
      (2*(-1 + z)*z^2*Qxx[x, z]*Qyy[x, z]*Qzz[x, z]^2*Derivative[1][P][z]*
        (-Derivative[0, 1][Qtt][x, z] + z^2*Qxz[x, z]*Derivative[1, 0][Qtt][
           x, z]) + z^2*Qtt[x, z]*Qxx[x, z]*Qyy[x, z]*Qzz[x, z]*
        (2*Qzz[x, z]*(Derivative[1][P][z] + (-1 + z)*Derivative[2][P][z]) + 
         (-1 + z)*(2*(-1 + z)*z^2*Derivative[0, 1][a][x, z]^2 - 
           4*(-1 + z)*z^4*Qxz[x, z]*Derivative[0, 1][a][x, z]*
            Derivative[1, 0][a][x, z] + 2*(-1 + z)*z^6*Qxz[x, z]^2*
            Derivative[1, 0][a][x, z]^2 + Derivative[1][P][z]*
            (2*Derivative[0, 1][Qtt][x, z] + Derivative[0, 1][Qzz][x, z] - 
             2*z^2*Qxz[x, z]*Derivative[1, 0][Qtt][x, z]))) + 
       Qtt[x, z]^2*(Qxx[x, z]*(-2*(-1 + z)*z*Qzz[x, z]^2*Derivative[1][P][
             z] + Qyy[x, z]*(12*Qzz[x, z]^2 - 2*(-2 + z)*z*Qzz[x, z]*
              Derivative[1][P][z] + (-1 + z)*z^2*Derivative[1][P][z]*
              Derivative[0, 1][Qzz][x, z])) + 2*z*Qyy[x, z]*
          (-((-1 + z)*Qzz[x, z]^2*Derivative[1][P][z]) + 
           z*Derivative[1, 0][Qzz][x, z]^2 - z*Qzz[x, z]*
            Derivative[2, 0][Qzz][x, z]))) + 
     P[z]^2*(z*Qtt[x, z]*Qxx[x, z]^2*Qyy[x, z]^2*Qzz[x, z]^2*
        (-(z*Qzz[x, z]) - (2 - 3*z + z^2)*(2*Derivative[0, 1][Qtt][x, z] + 
           Derivative[0, 1][Qzz][x, z] - 2*z^2*Qxz[x, z]*
            Derivative[1, 0][Qtt][x, z])) + (-1 + z)*z*Qxx[x, z]^2*
        Qyy[x, z]^2*Qzz[x, z]^2*(Derivative[0, 1][Qtt][x, z] - 
         z^2*Qxz[x, z]*Derivative[1, 0][Qtt][x, z])*(2*(-2 + z)*Qzz[x, z] + 
         (-1 + z)*z*(Derivative[0, 1][Qtt][x, z] - z^2*Qxz[x, z]*
            Derivative[1, 0][Qtt][x, z])) + Qtt[x, z]^2*
        ((-1 + z)^2*z*Qyy[x, z]^2*Qzz[x, z]^2*(Derivative[0, 1][Qxx][x, z] - 
           z^2*Qxz[x, z]*Derivative[1, 0][Qxx][x, z])*(4*Qzz[x, z] + 
           z*Derivative[0, 1][Qxx][x, z] - z^3*Qxz[x, z]*
            Derivative[1, 0][Qxx][x, z]) + 2*(-1 + z)*z*Qxx[x, z]*Qyy[x, z]^2*
          Qzz[x, z]^2*(-Qzz[x, z] - (-1 + z)*(Derivative[0, 1][Qzz][x, z] + 
             2*Derivative[0, 1][Qxx][x, z]*(1 + z^3*Derivative[1, 0][Qxz][x, 
                 z]) - 2*z^2*Qxz[x, z]*Derivative[1, 0][Qxx][x, z]*
              (1 + z^3*Derivative[1, 0][Qxz][x, z]))) + 
         Qxx[x, z]^2*((-1 + z)^2*z*Qzz[x, z]^2*(Derivative[0, 1][Qyy][x, z] - 
             z^2*Qxz[x, z]*Derivative[1, 0][Qyy][x, z])*(4*Qzz[x, z] + 
             z*Derivative[0, 1][Qyy][x, z] - z^3*Qxz[x, z]*Derivative[1, 0][
                Qyy][x, z]) + 2*(-1 + z)*z*Qyy[x, z]*Qzz[x, z]^2*
            (-Qzz[x, z] - (-1 + z)*(2*Derivative[0, 1][Qyy][x, z] + 
               Derivative[0, 1][Qzz][x, z] - 2*z^2*Qxz[x, z]*
                Derivative[1, 0][Qyy][x, z])) + Qyy[x, z]^2*
            (Qzz[x, z]^2*(12 - 20*z + 9*z^2 + 6*(-1 + z)^3*z^5*Qxz[x, z]^2*
                Derivative[1][P][z] + 8*(-1 + z)^2*z^3*Derivative[1, 0][Qxz][
                 x, z] + 4*(-1 + z)^2*z^6*Derivative[1, 0][Qxz][x, z]^2) - 
             3*(-1 + z)^2*z^2*(Derivative[0, 1][Qzz][x, z] - z^2*Qxz[x, z]*
                 Derivative[1, 0][Qzz][x, z])^2 + (-1 + z)*z*Qzz[x, z]*
              ((-2 + 3*z)*Derivative[0, 1][Qzz][x, z] + 2*(-1 + z)*z^2*
                Qxz[x, z]*(-2*Derivative[1, 0][Qzz][x, z] - 
                 2*z*Derivative[1, 1][Qzz][x, z] + z^3*Qxz[x, z]*
                  Derivative[2, 0][Qzz][x, z])))))))/
    (2*(-1 + z)^2*z^2*P[z]^2*Qtt[x, z]^2*Qxx[x, z]^2*Qyy[x, z]^2*Qzz[x, z]), 
  Derivative[0, 2][Qxx][x, z] + 
   (Qzz[x, z]*(4*(-1 + z)*P[z] - (2*z^4*a[x, z]^2*Qxx[x, z])/
       (Qtt[x, z]*Qzz[x, z]) - (2*(-1 + z)*z*P[z]*Derivative[0, 1][Qxx][x, 
         z])/Qyy[x, z] + (z*(-2 + 3*z)*P[z]*Derivative[0, 1][Qxx][x, z])/
       Qzz[x, z] - ((-2 + z)*P[z]*(-2*Qxx[x, z] + 
         z*Derivative[0, 1][Qxx][x, z]))/Qtt[x, z] - 
      (2*(-1 + z)*z^4*Derivative[1, 0][a][x, z]^2)/(P[z]*Qtt[x, z]) + 
      (4*(-1 + z)*z^4*a[x, z]*Qxx[x, z]*(-Derivative[0, 1][a][x, z] + 
         z^2*Qxz[x, z]*Derivative[1, 0][a][x, z]))/(Qtt[x, z]*Qzz[x, z]) + 
      (2*(-2 + z)*z^3*P[z]*Qxx[x, z]*Qxz[x, z]*Derivative[1, 0][Qtt][x, z])/
       Qtt[x, z]^2 - (2*(-1 + z)*z*P[z]*
        (Qzz[x, z]*(Derivative[0, 1][Qxx][x, z] - 2*z^2*Qxz[x, z]*
            Derivative[1, 0][Qxx][x, z]) + 
         z*(Derivative[0, 1][Qxx][x, z] - z^2*Qxz[x, z]*Derivative[1, 0][Qxx][
              x, z])^2))/(Qxx[x, z]*Qzz[x, z]) + 
      (2*(-1 + z)^2*z^4*P[z]^2*Qxz[x, z]^2*
        (z*Qzz[x, z]*Derivative[0, 1][Qxx][x, z] - 2*Qxx[x, z]*
          (Qzz[x, z]*(1 - 2*z^3*Derivative[1, 0][Qxz][x, z]) + 
           z^3*Qxz[x, z]*Derivative[1, 0][Qzz][x, z])))/Qzz[x, z]^2 - 
      (2*P[z]*Qxx[x, z]*(-2*(-1 + z)*Qyy[x, z]*Qzz[x, z]^2 - 
         2*(-1 + z)*z^3*Qxz[x, z]*Qzz[x, z]^2*Derivative[1, 0][Qyy][x, z] + 
         z*Qyy[x, z]^2*(Qzz[x, z]*(1 + 4*(-1 + z)*z^2*Derivative[1, 0][Qxz][
               x, z] + 2*(-1 + z)*z^5*Derivative[1, 0][Qxz][x, z]^2) + 
           z^2*(2*(-1 + z)*z*Derivative[0, 1][Qxz][x, z] + 
             Qxz[x, z]*(-2 + 3*z - 2*(-1 + z)*z^3*Derivative[1, 0][Qxz][x, 
                 z]))*Derivative[1, 0][Qzz][x, z])))/
       (Qyy[x, z]^2*Qzz[x, z]^2) + (2*(-1 + z)*z^3*P[z]*Qxz[x, z]*
        (-2*Derivative[1, 0][Qxx][x, z] - 2*z*Derivative[1, 1][Qxx][x, z] + 
         z^3*Qxz[x, z]*Derivative[2, 0][Qxx][x, z]))/Qzz[x, z] + 
      (-((-1 + z)*z*Qtt[x, z]*Qxx[x, z]^2*Qyy[x, z]^2*Qzz[x, z]*
          (-(z*Qzz[x, z]*Derivative[1][P][z]*Derivative[0, 1][Qxx][x, z]) + 
           2*Qxx[x, z]*(Qzz[x, z]*Derivative[1][P][z] + (-1 + z)*z^3*
              (Derivative[0, 1][a][x, z] - z^2*Qxz[x, z]*Derivative[1, 0][a][
                  x, z])^2))) - z^2*Qxx[x, z]^2*Qyy[x, z]^2*Qzz[x, z]^2*
         Derivative[1, 0][Qtt][x, z]*(2*(-1 + z)*z^2*Qxx[x, z]*Qxz[x, z]*
           Derivative[1][P][z] + Derivative[1, 0][Qtt][x, z]) + 
        Qtt[x, z]^2*(3*z^2*Qyy[x, z]^2*Qzz[x, z]^2*
           Derivative[1, 0][Qxx][x, z]^2 + 2*Qxx[x, z]^3*Qyy[x, z]^2*
           (6*Qzz[x, z]^2 - (-1 + z)*z*Qzz[x, z]*Derivative[1][P][z] - 
            (-1 + z)*z^4*Qxz[x, z]*Derivative[1][P][z]*Derivative[1, 0][Qzz][
              x, z]) + z^2*Qxx[x, z]^2*(-(Qzz[x, z]^2*Derivative[1, 0][Qyy][
                x, z]^2) + Qyy[x, z]^2*((-1 + z)*Qzz[x, z]*Derivative[1][P][
                z]*Derivative[0, 1][Qxx][x, z] - Derivative[1, 0][Qzz][x, 
                z]^2)) - 2*z^2*Qxx[x, z]*Qyy[x, z]^2*Qzz[x, z]^2*
           Derivative[2, 0][Qxx][x, z]))/(Qtt[x, z]^2*Qxx[x, z]^2*Qyy[x, z]^2*
        Qzz[x, z]^2)))/(2*(-1 + z)*z^2*P[z]), Derivative[0, 2][Qyy][x, z] + 
   (Qzz[x, z]*(4*(-1 + z)*P[z] - (2*z*P[z]*Qyy[x, z])/Qzz[x, z] - 
      (2*z^4*a[x, z]^2*Qyy[x, z])/(Qtt[x, z]*Qzz[x, z]) + 
      (z*(-2 + 3*z)*P[z]*Derivative[0, 1][Qyy][x, z])/Qzz[x, z] + 
      ((-2 + z)*P[z]*(2*Qyy[x, z] - z*Derivative[0, 1][Qyy][x, z]))/
       Qtt[x, z] - (2*(-1 + z)*P[z]*(-2*Qyy[x, z] + 
         z*Derivative[0, 1][Qyy][x, z]))/Qxx[x, z] + 
      (2*(-1 + z)^2*z^4*P[z]^2*Qxz[x, z]^2*(-2*Qyy[x, z] + 
         z*Derivative[0, 1][Qyy][x, z]))/Qzz[x, z] + 
      (2*(-1 + z)*z^4*Qyy[x, z]*Derivative[1, 0][a][x, z]^2)/
       (P[z]*Qtt[x, z]*Qxx[x, z]) + (4*(-1 + z)*z^4*a[x, z]*Qyy[x, z]*
        (-Derivative[0, 1][a][x, z] + z^2*Qxz[x, z]*Derivative[1, 0][a][x, 
           z]))/(Qtt[x, z]*Qzz[x, z]) - 
      (2*(-1 + z)*z*P[z]*(Qzz[x, z]*Derivative[0, 1][Qyy][x, z] + 
         z*(Derivative[0, 1][Qyy][x, z] - z^2*Qxz[x, z]*Derivative[1, 0][Qyy][
              x, z])^2))/(Qyy[x, z]*Qzz[x, z]) + 
      (2*(-1 + z)*z^3*P[z]*Qxz[x, z]*(-2*Derivative[1, 0][Qyy][x, z] - 
         2*z*Derivative[1, 1][Qyy][x, z] + z^3*Qxz[x, z]*
          Derivative[2, 0][Qyy][x, z]))/Qzz[x, z] + 
      (-((-1 + z)*z*Qxx[x, z]*Qyy[x, z]*(-(z*Qzz[x, z]*Derivative[1][P][z]*
             Derivative[0, 1][Qyy][x, z]) + 2*Qyy[x, z]*
            (Qzz[x, z]*Derivative[1][P][z] + (-1 + z)*z^3*
              (Derivative[0, 1][a][x, z] - z^2*Qxz[x, z]*Derivative[1, 0][a][
                  x, z])^2))) + Qtt[x, z]*(Qxx[x, z]*Qyy[x, z]*
           (2*Qyy[x, z]*(6*Qzz[x, z] - (-1 + z)*z*Derivative[1][P][z]) + 
            (-1 + z)*z^2*Derivative[1][P][z]*Derivative[0, 1][Qyy][x, z]) + 
          2*z^2*Qzz[x, z]*(Derivative[1, 0][Qyy][x, z]^2 - 
            Qyy[x, z]*Derivative[2, 0][Qyy][x, z])))/(Qtt[x, z]*Qxx[x, z]*
        Qyy[x, z]*Qzz[x, z])))/(2*(-1 + z)*z^2*P[z]), 
  Derivative[0, 2][Qxz][x, z] + 
   ((-1 + z)*z*Qxx[x, z]^2*Qyy[x, z]^2*Qzz[x, z]^2*
      (-4*z^2*a[x, z]*Qtt[x, z]*Derivative[1, 0][a][x, z] + 
       Qzz[x, z]*Derivative[1][P][z]*Derivative[1, 0][Qtt][x, z] + 
       Qtt[x, z]*(-4*(-1 + z)*z^2*Derivative[0, 1][a][x, z]*
          Derivative[1, 0][a][x, z] + 4*(-1 + z)*z^4*Qxz[x, z]*
          Derivative[1, 0][a][x, z]^2 - Derivative[1][P][z]*
          Derivative[1, 0][Qtt][x, z])) + 2*(-1 + z)^2*z^5*P[z]^3*Qtt[x, z]^2*
      Qxx[x, z]^3*Qxz[x, z]^2*Qyy[x, z]^2*
      (3*(-1 + z)*z*Qzz[x, z]*Derivative[0, 1][Qxz][x, z] + 
       Qxz[x, z]*(-((-1 + z)*z*Derivative[0, 1][Qzz][x, z]) + 
         Qzz[x, z]*(-7 + 9*z - 2*(-1 + z)*z^3*Derivative[1, 0][Qxz][x, z])) + 
       (-1 + z)*z^3*Qxz[x, z]^2*Derivative[1, 0][Qzz][x, z]) + 
     (-1 + z)*z*P[z]^2*Qxx[x, z]*(Qtt[x, z]*Qxx[x, z]^2*Qyy[x, z]^2*
        Qzz[x, z]^2*((6 - 4*z)*Qxz[x, z] - (-2 + z)*z*Derivative[0, 1][Qxz][
           x, z]) + (-2 + z)*z*Qxx[x, z]^2*Qxz[x, z]*Qyy[x, z]^2*Qzz[x, z]^2*
        (Derivative[0, 1][Qtt][x, z] - z^2*Qxz[x, z]*Derivative[1, 0][Qtt][x, 
           z]) + Qtt[x, z]^2*(-2*(-1 + z)*z*Qxz[x, z]*Qyy[x, z]^2*Qzz[x, z]^2*
          (-Derivative[0, 1][Qxx][x, z] + z^2*Qxz[x, z]*Derivative[1, 0][Qxx][
             x, z]) + 2*Qxx[x, z]*Qyy[x, z]^2*Qzz[x, z]*
          (-((-1 + z)*z*Qzz[x, z]*Derivative[0, 1][Qxz][x, z]) + 
           Qxz[x, z]*Qzz[x, z]*(3 - 4*z - 2*(-1 + z)*z^3*Derivative[1, 0][
                Qxz][x, z]) + (-1 + z)*z^3*Qxz[x, z]^2*Derivative[1, 0][Qzz][
             x, z]) + Qxx[x, z]^2*(4*(-1 + z)^2*z^5*Qxz[x, z]^3*Qyy[x, z]^2*
            Qzz[x, z]*Derivative[1][P][z] + z*Qyy[x, z]*Derivative[0, 1][Qxz][
             x, z]*(-2*(-1 + z)*Qzz[x, z]^2 + Qyy[x, z]*((-10 + 13*z)*
                Qzz[x, z] - 2*(-1 + z)*z*Derivative[0, 1][Qzz][x, z])) + 
           Qxz[x, z]*(2*(3 - 4*z)*Qyy[x, z]*Qzz[x, z]^2 + 2*(-1 + z)*z*
              Qzz[x, z]^2*Derivative[0, 1][Qyy][x, z] + Qyy[x, z]^2*
              (z*(Derivative[0, 1][Qzz][x, z]*(2 - 3*z + 2*(-1 + z)*z^3*
                    Derivative[1, 0][Qxz][x, z]) + 2*(-1 + z)*z^3*
                  Derivative[0, 1][Qxz][x, z]*Derivative[1, 0][Qzz][x, z]) - 
               2*Qzz[x, z]*(3 - 6*z + z^3*(-6 + 7*z)*Derivative[1, 0][Qxz][x, 
                   z] + 2*(-1 + z)*z^4*Derivative[1, 1][Qxz][x, z]))) + 
           z^3*Qxz[x, z]^2*(-2*(-1 + z)*Qzz[x, z]^2*Derivative[1, 0][Qyy][x, 
               z] + Qyy[x, z]^2*(-2 + 3*z - 2*(-1 + z)*z^3*Derivative[1, 0][
                  Qxz][x, z])*Derivative[1, 0][Qzz][x, z] + 
             2*(-1 + z)*z^3*Qyy[x, z]^2*Qzz[x, z]*Derivative[2, 0][Qxz][x, 
               z])))) + P[z]*(Qtt[x, z]*Qxx[x, z]^2*Qyy[x, z]^2*Qzz[x, z]^2*
        ((-1 + z)*z^2*Qxx[x, z]*(Qxz[x, z]*((-2 + 4*z)*Derivative[1][P][z] + 
             (-1 + z)*z*Derivative[2][P][z]) + (-1 + z)*z*Derivative[1][P][z]*
            Derivative[0, 1][Qxz][x, z]) + (-2 + z)*Derivative[1, 0][Qtt][x, 
           z]) + Qxx[x, z]^2*Qyy[x, z]^2*Qzz[x, z]^2*
        ((-1 + z)^2*z^3*Qxx[x, z]*Qxz[x, z]*Derivative[1][P][z]*
          (-Derivative[0, 1][Qtt][x, z] + z^2*Qxz[x, z]*Derivative[1, 0][Qtt][
             x, z]) + Derivative[1, 0][Qtt][x, z]*(-((-2 + z)*Qzz[x, z]) + 
           (-1 + z)*z*(-Derivative[0, 1][Qtt][x, z] + z^2*Qxz[x, z]*
              Derivative[1, 0][Qtt][x, z]))) + (-1 + z)*Qtt[x, z]^2*
        (Qyy[x, z]^2*Qzz[x, z]^2*Derivative[1, 0][Qxx][x, z]*
          (-2*Qzz[x, z] + z*Derivative[0, 1][Qxx][x, z] - 
           z^3*Qxz[x, z]*Derivative[1, 0][Qxx][x, z]) + 
         z^2*Qxx[x, z]^3*Qyy[x, z]*(3*(-1 + z)*z*Qyy[x, z]*Qzz[x, z]*
            Derivative[1][P][z]*Derivative[0, 1][Qxz][x, z] + 
           Qxz[x, z]*(-2*(-1 + z)*Qzz[x, z]^2*Derivative[1][P][z] + 
             Qyy[x, z]*(-((-1 + z)*z*Derivative[1][P][z]*Derivative[0, 1][
                   Qzz][x, z]) + Qzz[x, z]*((-1 + z)*z*Derivative[2][P][z] + 
                 Derivative[1][P][z]*(-6 + 8*z - 2*(-1 + z)*z^3*
                    Derivative[1, 0][Qxz][x, z])))) + (-1 + z)*z^3*
            Qxz[x, z]^2*Qyy[x, z]*Derivative[1][P][z]*Derivative[1, 0][Qzz][
             x, z]) - 2*Qxx[x, z]*Qyy[x, z]^2*Qzz[x, z]*
          (Qzz[x, z]*Derivative[1, 0][Qxx][x, z]*
            (1 + z^3*Derivative[1, 0][Qxz][x, z]) + 
           z*(Derivative[0, 1][Qxx][x, z] - z^2*Qxz[x, z]*Derivative[1, 0][
                Qxx][x, z])*Derivative[1, 0][Qzz][x, z]) + 
         Qxx[x, z]^2*(2*Qyy[x, z]*Qzz[x, z]^2*Derivative[1, 0][Qyy][x, z] - 
           Qzz[x, z]^2*(2*Qzz[x, z] + z*Derivative[0, 1][Qyy][x, z])*
            Derivative[1, 0][Qyy][x, z] + z^2*Qxz[x, z]*
            (z*Qzz[x, z]^2*Derivative[1, 0][Qyy][x, z]^2 - 
             Qyy[x, z]^2*(2*(-1 + z)*Qzz[x, z]^2*Derivative[1][P][z] + z*
                Derivative[1, 0][Qzz][x, z]^2)) + Qyy[x, z]^2*
            (z*Derivative[0, 1][Qzz][x, z]*Derivative[1, 0][Qzz][x, z] + 
             4*Qzz[x, z]*(1 + z^3*Derivative[1, 0][Qxz][x, z])*
              Derivative[1, 0][Qzz][x, z] - 2*z^3*Qzz[x, z]^2*
              Derivative[2, 0][Qxz][x, z])))))/(2*(-1 + z)^2*z^3*P[z]^2*
     Qtt[x, z]^2*Qxx[x, z]^3*Qyy[x, z]^2*Qzz[x, z])}, 
 (Qxx[x, z]^2*Qyy[x, z]^2*Qzz[x, z]^2*
    ((-1 + z)*z^4*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*
      (4 - 2*z^4*\[Mu]1^2 + z^3*(2 + \[Mu]1^2))^2*Qxx[x, z]^2*Qxz[x, z]^2*
      Qzz[x, z] + 4*z^2*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*Qzz[x, z]*
      Derivative[1, 0][Qtt][x, z]^2 + 2*Qxx[x, z]*
      ((4 - 2*z^4*\[Mu]1^2 + z^3*(2 + \[Mu]1^2))^2*Qzz[x, z]^2 - 
       2*z*(-8 + 2*z^8*\[Mu]1^4 + 2*z^3*(2 + \[Mu]1^2) - 
         3*z^7*\[Mu]1^2*(2 + \[Mu]1^2) + z^6*(2 + \[Mu]1^2)^2)*Qzz[x, z]*
        Derivative[0, 1][Qtt][x, z] + (-1 + z)^2*z^2*
        (2 + 2*z + 2*z^2 - z^3*\[Mu]1^2)^2*(Derivative[0, 1][Qtt][x, z] - 
          z^2*Qxz[x, z]*Derivative[1, 0][Qtt][x, z])^2)) - 
   2*Qtt[x, z]*Qxx[x, z]*Qyy[x, z]*Qzz[x, z]*
    (4*z^2*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*Qyy[x, z]*Qzz[x, z]^2*
      Derivative[1, 0][Qtt][x, z]*Derivative[1, 0][Qxx][x, z] + 
     (-1 + z)*z^4*(8 + 8*z + 8*z^2 + 2*z^7*\[Mu]1^4 - 2*z^3*(-2 + \[Mu]1^2) - 
       2*z^4*(-2 + \[Mu]1^2) - 2*z^5*(-2 + \[Mu]1^2) - 
       z^6*\[Mu]1^2*(6 + \[Mu]1^2))*Qxx[x, z]^3*Qxz[x, z]*Qzz[x, z]*
      (z^4*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))^2*Qxz[x, z]^3*Qyy[x, z] + 
       2*z*(-2 - z^4*\[Mu]1^2 + z^3*(2 + \[Mu]1^2))*Qyy[x, z]*
        Derivative[0, 1][Qxz][x, z] + Qxz[x, z]*
        (2*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*Qzz[x, z] + 
         Qyy[x, z]*(-4 - 6*z^4*\[Mu]1^2 + 5*z^3*(2 + \[Mu]1^2) + 
           2*(-1 + z)*z^3*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*
            Derivative[1, 0][Qxz][x, z]))) - 
     2*(-1 + z)*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*Qxx[x, z]*Qzz[x, z]*
      (2*z^2*Qzz[x, z]*Derivative[1, 0][Qtt][x, z]*Derivative[1, 0][Qyy][x, 
         z] + Qyy[x, z]*((8 - 4*z^4*\[Mu]1^2 + 2*z^3*(2 + \[Mu]1^2))*
          Qzz[x, z]^2 + z*Qzz[x, z]*(2*(2 + z^4*\[Mu]1^2 - 
             z^3*(2 + \[Mu]1^2))*Derivative[0, 1][Qtt][x, z] + 
           (-4 + 2*z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*
            (-Derivative[0, 1][Qxx][x, z] + 2*z^2*Qxz[x, z]*Derivative[1, 0][
                Qxx][x, z])) + z^2*((-1 + z)*(-2 - 2*z - 2*z^2 + 
             z^3*\[Mu]1^2)*Derivative[0, 1][Qtt][x, z]*
            (Derivative[0, 1][Qxx][x, z] - z^2*Qxz[x, z]*Derivative[1, 0][
                Qxx][x, z]) + Derivative[1, 0][Qtt][x, z]*
            (z^2*(-2 - z^4*\[Mu]1^2 + z^3*(2 + \[Mu]1^2))*Qxz[x, z]*
              Derivative[0, 1][Qxx][x, z] + z^4*(2 + z^4*\[Mu]1^2 - z^3*
                (2 + \[Mu]1^2))*Qxz[x, z]^2*Derivative[1, 0][Qxx][x, z] + 
             2*Derivative[1, 0][Qzz][x, z])))) + 
     2*Qxx[x, z]^2*(-((-1 + z)*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*Qzz[x, z]*
         ((8 - 4*z^4*\[Mu]1^2 + 2*z^3*(2 + \[Mu]1^2))*Qzz[x, z]^2 + 
          z*Qzz[x, z]*(2*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*
             Derivative[0, 1][Qtt][x, z] + (4 - 2*z^4*\[Mu]1^2 + 
              z^3*(2 + \[Mu]1^2))*Derivative[0, 1][Qyy][x, z]) + 
          (-1 + z)*z^2*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*
           (-Derivative[0, 1][Qtt][x, z] + z^2*Qxz[x, z]*
             Derivative[1, 0][Qtt][x, z])*(-Derivative[0, 1][Qyy][x, z] + 
            z^2*Qxz[x, z]*Derivative[1, 0][Qyy][x, z]))) + 
       Qyy[x, z]*((-4 + 2*z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*Qzz[x, z]^2*
          (-12 + 6*z^3 + 3*z^3*\[Mu]1^2 - 2*z^4*\[Mu]1^2 + 
           2*z^4*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))^2*Qxz[x, z]^2 + 
           (-4*z^3 - 2*z^7*\[Mu]1^2 + 2*z^6*(2 + \[Mu]1^2))*
            Derivative[1, 0][Qxz][x, z]) + 
         z^2*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))^2*
          (Derivative[0, 1][Qtt][x, z] - z^2*Qxz[x, z]*Derivative[1, 0][Qtt][
             x, z])*(Derivative[0, 1][Qzz][x, z] - z^2*Qxz[x, z]*
            Derivative[1, 0][Qzz][x, z]) - (-1 + z)*z*(-2 - 2*z - 2*z^2 + 
           z^3*\[Mu]1^2)*Qzz[x, z]*((-4 + 2*z^4*\[Mu]1^2 - 
             z^3*(2 + \[Mu]1^2))*Derivative[0, 1][Qzz][x, z] + 
           Derivative[0, 1][Qtt][x, z]*(-12 + 6*z^3 + 3*z^3*\[Mu]1^2 - 
             2*z^4*\[Mu]1^2 + z^4*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))^2*
              Qxz[x, z]^2 + (-4*z^3 - 2*z^7*\[Mu]1^2 + 2*z^6*(2 + \[Mu]1^2))*
              Derivative[1, 0][Qxz][x, z]) + 
           2*z^2*(z*(-2 - z^4*\[Mu]1^2 + z^3*(2 + \[Mu]1^2))*Derivative[0, 1][
                Qxz][x, z]*Derivative[1, 0][Qtt][x, z] + 
             Qxz[x, z]*(Derivative[1, 0][Qtt][x, z]*(4 - 2*z^4*\[Mu]1^2 + 
                 z^3*(2 + \[Mu]1^2) + 2*(-1 + z)*z^3*(-2 - 2*z - 2*z^2 + 
                   z^3*\[Mu]1^2)*Derivative[1, 0][Qxz][x, z]) + 
               (4 - 2*z^4*\[Mu]1^2 + z^3*(2 + \[Mu]1^2))*Derivative[1, 0][
                  Qzz][x, z])))))) + Qtt[x, z]^2*
    (4*z^2*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*Qyy[x, z]^2*Qzz[x, z]^3*
      Derivative[1, 0][Qxx][x, z]^2 + (-1 + z)*z^4*(-2 - 2*z - 2*z^2 + 
       z^3*\[Mu]1^2)*Qxx[x, z]^4*Qzz[x, z]*
      (z^4*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))^2*Qxz[x, z]^3*Qyy[x, z] + 
        2*z*(-2 - z^4*\[Mu]1^2 + z^3*(2 + \[Mu]1^2))*Qyy[x, z]*
         Derivative[0, 1][Qxz][x, z] + Qxz[x, z]*
         (2*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*Qzz[x, z] + 
          Qyy[x, z]*(-4 - 6*z^4*\[Mu]1^2 + 5*z^3*(2 + \[Mu]1^2) + 
            2*(-1 + z)*z^3*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*
             Derivative[1, 0][Qxz][x, z])))^2 + 
     2*(-1 + z)*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*Qxx[x, z]*Qyy[x, z]*
      Qzz[x, z]^2*(-4*z^2*Qzz[x, z]*Derivative[1, 0][Qxx][x, z]*
        Derivative[1, 0][Qyy][x, z] + Qyy[x, z]*
        (4*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*Qzz[x, z]^2 - 
         4*(-1 + z)*z*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*Qzz[x, z]*
          (-Derivative[0, 1][Qxx][x, z] + 2*z^2*Qxz[x, z]*
            Derivative[1, 0][Qxx][x, z]) + 
         z^2*((2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*
            Derivative[0, 1][Qxx][x, z]^2 + 2*z^2*(-2 - z^4*\[Mu]1^2 + 
             z^3*(2 + \[Mu]1^2))*Qxz[x, z]*Derivative[0, 1][Qxx][x, z]*
            Derivative[1, 0][Qxx][x, z] + Derivative[1, 0][Qxx][x, z]*
            (z^4*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*Qxz[x, z]^2*
              Derivative[1, 0][Qxx][x, z] - 4*Derivative[1, 0][Qzz][x, 
               z])))) + 4*(-1 + z)*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*
      Qxx[x, z]^2*Qzz[x, z]*(z^2*Qzz[x, z]^2*Derivative[1, 0][Qyy][x, z]^2 + 
       Qyy[x, z]*Qzz[x, z]*(4*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*
          Qzz[x, z]^2 - 2*(-1 + z)*z*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*
          Qzz[x, z]*(-Derivative[0, 1][Qxx][x, z] - Derivative[0, 1][Qyy][x, 
            z] + 2*z^2*Qxz[x, z]*Derivative[1, 0][Qxx][x, z]) + 
         z^2*(z^2*(-2 - z^4*\[Mu]1^2 + z^3*(2 + \[Mu]1^2))*Qxz[x, z]*
            Derivative[0, 1][Qyy][x, z]*Derivative[1, 0][Qxx][x, z] + 
           z^4*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*Qxz[x, z]^2*
            Derivative[1, 0][Qxx][x, z]*Derivative[1, 0][Qyy][x, z] + 
           (-1 + z)*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*Derivative[0, 1][Qxx][
             x, z]*(Derivative[0, 1][Qyy][x, z] - z^2*Qxz[x, z]*
              Derivative[1, 0][Qyy][x, z]) + 2*Derivative[1, 0][Qyy][x, z]*
            Derivative[1, 0][Qzz][x, z])) + Qyy[x, z]^2*
        (Qzz[x, z]^2*(-24 + 12*z^3 + 6*z^3*\[Mu]1^2 - 4*z^4*\[Mu]1^2 + 
           3*z^4*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))^2*Qxz[x, z]^2 + 
           (-8*z^3 - 4*z^7*\[Mu]1^2 + 4*z^6*(2 + \[Mu]1^2))*
            Derivative[1, 0][Qxz][x, z]) + 
         z^2*(z^2*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*Qxz[x, z]*
            Derivative[0, 1][Qzz][x, z]*Derivative[1, 0][Qxx][x, z] + 
           z^4*(-2 - z^4*\[Mu]1^2 + z^3*(2 + \[Mu]1^2))*Qxz[x, z]^2*
            Derivative[1, 0][Qxx][x, z]*Derivative[1, 0][Qzz][x, z] + 
           Derivative[1, 0][Qzz][x, z]^2 + (-1 + z)*(-2 - 2*z - 2*z^2 + 
             z^3*\[Mu]1^2)*Derivative[0, 1][Qxx][x, z]*
            (-Derivative[0, 1][Qzz][x, z] + z^2*Qxz[x, z]*Derivative[1, 0][
                Qzz][x, z])) + z*Qzz[x, z]*(Derivative[0, 1][Qxx][x, z]*
            (-12 + 6*z^3 + 3*z^3*\[Mu]1^2 - 2*z^4*\[Mu]1^2 + 
             z^4*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))^2*Qxz[x, z]^2 + 
             (-4*z^3 - 2*z^7*\[Mu]1^2 + 2*z^6*(2 + \[Mu]1^2))*
              Derivative[1, 0][Qxz][x, z]) - 2*(-1 + z)*(-2 - 2*z - 2*z^2 + 
             z^3*\[Mu]1^2)*(Derivative[0, 1][Qzz][x, z] + 
             z^2*(z^4*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*Qxz[x, z]^3*
                Derivative[1, 0][Qxx][x, z] - z*Derivative[0, 1][Qxz][x, z]*
                Derivative[1, 0][Qxx][x, z] - 2*Qxz[x, z]*
                (2*Derivative[1, 0][Qxx][x, z] + Derivative[1, 0][Qzz][x, 
                  z])))))) + 2*Qxx[x, z]^3*
      ((2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))^2*Qzz[x, z]^2*
        (4*Qzz[x, z]^2 + 4*z*Qzz[x, z]*Derivative[0, 1][Qyy][x, z] + 
         (z*Derivative[0, 1][Qyy][x, z] - z^3*Qxz[x, z]*Derivative[1, 0][Qyy][
             x, z])^2) + Qyy[x, z]^2*
        (Qzz[x, z]^2*(-16*z^4*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))^3*
            Qxz[x, z]^2 + 3*z^8*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))^4*
            Qxz[x, z]^4 - 4*z^5*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))^3*
            Qxz[x, z]*Derivative[0, 1][Qxz][x, z] + 
           (12 + 2*z^4*\[Mu]1^2 - 3*z^3*(2 + \[Mu]1^2) + 2*(-1 + z)*z^3*
              (-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*Derivative[1, 0][Qxz][x, z])^
            2) + (-1 + z)^2*z^2*(2 + 2*z + 2*z^2 - z^3*\[Mu]1^2)^2*
          (Derivative[0, 1][Qzz][x, z] - z^2*Qxz[x, z]*Derivative[1, 0][Qzz][
              x, z])^2 + 2*(-1 + z)*z*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*
          Qzz[x, z]*(-(Derivative[0, 1][Qzz][x, z]*(-12 + 6*z^3 + 
              3*z^3*\[Mu]1^2 - 2*z^4*\[Mu]1^2 + z^4*(2 + z^4*\[Mu]1^2 - 
                 z^3*(2 + \[Mu]1^2))^2*Qxz[x, z]^2 + (-4*z^3 - 
                2*z^7*\[Mu]1^2 + 2*z^6*(2 + \[Mu]1^2))*Derivative[1, 0][Qxz][
                x, z])) + 2*(-1 + z)*z^2*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*
            (-4*Qxz[x, z] + (-1 + z)*z^4*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*
              Qxz[x, z]^3 - z*Derivative[0, 1][Qxz][x, z])*
            Derivative[1, 0][Qzz][x, z])) + 2*(-1 + z)*(-2 - 2*z - 2*z^2 + 
         z^3*\[Mu]1^2)*Qyy[x, z]*Qzz[x, z]*
        (2*Qzz[x, z]^2*(-12 + 6*z^3 + 3*z^3*\[Mu]1^2 - 2*z^4*\[Mu]1^2 + 
           2*z^4*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))^2*Qxz[x, z]^2 + 
           (-4*z^3 - 2*z^7*\[Mu]1^2 + 2*z^6*(2 + \[Mu]1^2))*
            Derivative[1, 0][Qxz][x, z]) - (-1 + z)*z^2*(-2 - 2*z - 2*z^2 + 
           z^3*\[Mu]1^2)*(-Derivative[0, 1][Qyy][x, z] + z^2*Qxz[x, z]*
            Derivative[1, 0][Qyy][x, z])*(-Derivative[0, 1][Qzz][x, z] + 
           z^2*Qxz[x, z]*Derivative[1, 0][Qzz][x, z]) + 
         z*Qzz[x, z]*(Derivative[0, 1][Qyy][x, z]*(-12 + 6*z^3 + 
             3*z^3*\[Mu]1^2 - 2*z^4*\[Mu]1^2 + z^4*(2 + z^4*\[Mu]1^2 - 
                z^3*(2 + \[Mu]1^2))^2*Qxz[x, z]^2 + 
             (-4*z^3 - 2*z^7*\[Mu]1^2 + 2*z^6*(2 + \[Mu]1^2))*
              Derivative[1, 0][Qxz][x, z]) + 
           2*((-2 - z^4*\[Mu]1^2 + z^3*(2 + \[Mu]1^2))*Derivative[0, 1][Qzz][
               x, z] + z^2*(z*(-2 - z^4*\[Mu]1^2 + z^3*(2 + \[Mu]1^2))*
                Derivative[0, 1][Qxz][x, z]*Derivative[1, 0][Qyy][x, z] + 
               Qxz[x, z]*((4 - 2*z^4*\[Mu]1^2 + z^3*(2 + \[Mu]1^2) + 
                   2*(-1 + z)*z^3*(-2 - 2*z - 2*z^2 + z^3*\[Mu]1^2)*
                    Derivative[1, 0][Qxz][x, z])*Derivative[1, 0][Qyy][x, 
                   z] + 2*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*
                  Derivative[1, 0][Qzz][x, z]))))))))/
  (16*(2 + z^4*\[Mu]1^2 - z^3*(2 + \[Mu]1^2))*Qtt[x, z]^2*Qxx[x, z]^3*
   Qyy[x, z]^2*Qzz[x, z]^3), {Derivative[0, 1][a][x, y] - 
   (8*a[x, y]^3*Qxx[x, y]^3*Qyy[x, y]^2 - 2*Qtt[x, y]*Qxx[x, y]*Qyy[x, y]*
      (-2*(-6 + \[Mu]1^2)*Qxx[x, y]^2*Qxz[x, y]*Qyy[x, y]*
        Derivative[1, 0][a][x, y] - Qtt[x, y]*Qyy[x, y]*
        Derivative[1, 0][a][x, y]*Derivative[1, 0][Qxx][x, y] + 
       Qtt[x, y]*Qxx[x, y]*(Derivative[1, 0][a][x, y]*Derivative[1, 0][Qyy][
           x, y] + 2*Qyy[x, y]*Derivative[2, 0][a][x, y])) + 
     a[x, y]*(2*Qxx[x, y]^2*Qyy[x, y]^2*Derivative[1, 0][Qtt][x, y]*
        (-2*(-6 + \[Mu]1^2)*Qxx[x, y]*Qxz[x, y] + 
         3*Derivative[1, 0][Qtt][x, y]) + Qtt[x, y]*Qxx[x, y]^2*Qyy[x, y]*
        (Qxx[x, y]*(2*Qyy[x, y]*(-6*(-2 + \[Mu]1^2) + (-6 + \[Mu]1^2)*
              Derivative[1, 0][Qxz][x, y]) + (-6 + \[Mu]1^2)*Qxz[x, y]*
            Derivative[1, 0][Qyy][x, y]) + Qyy[x, y]*
          ((-6 + \[Mu]1^2)*Qxz[x, y]*Derivative[1, 0][Qxx][x, y] - 
           4*Derivative[2, 0][Qtt][x, y])) + Qtt[x, y]^2*
        (2*(-6 + \[Mu]1^2)*Qxx[x, y]^3*Qyy[x, y] - 3*Qyy[x, y]^2*
          Derivative[1, 0][Qxx][x, y]^2 + 2*Qxx[x, y]*Qyy[x, y]^2*
          Derivative[2, 0][Qxx][x, y] + Qxx[x, y]^2*
          (2*(-6 + \[Mu]1^2)*Qyy[x, y]^2 - Derivative[1, 0][Qyy][x, y]^2 + 
           2*Qyy[x, y]*Derivative[2, 0][Qyy][x, y]))))/
    (4*(-6 + \[Mu]1^2)*Qtt[x, y]*Qxx[x, y]^3*Qyy[x, y]^2), 
  -Qtt[x, y] + Qzz[x, y], Derivative[0, 1][Qtt][x, y] - 
   (2*Qtt[x, y]^2*((-6 + \[Mu]1^2)*Qyy[x, y] + 
       Qxx[x, y]*(-6 + \[Mu]1^2 + 12*Qyy[x, y])) + 
     Qyy[x, y]*(4*a[x, y]^2*Qxx[x, y] - (-6 + \[Mu]1^2)*Qxx[x, y]*
        Derivative[0, 1][Qzz][x, y] + 4*Derivative[1, 0][Qtt][x, y]^2) - 
     4*Qtt[x, y]*Qyy[x, y]*(2*\[Mu]1^2*Qxx[x, y] + Derivative[2, 0][Qtt][x, 
        y]))/((-6 + \[Mu]1^2)*Qxx[x, y]*Qyy[x, y]), 
  Derivative[0, 1][Qxx][x, y] - (-2*a[x, y]^2*Qxx[x, y]^3*Qyy[x, y]^2 + 
     2*(-6 + \[Mu]1^2)*Qtt[x, y]*Qxx[x, y]^3*Qyy[x, y]^2 + 
     2*Qxx[x, y]^2*Qyy[x, y]^2*((-6 + \[Mu]1^2)*Qxx[x, y]*Qxz[x, y] - 
       Derivative[1, 0][Qtt][x, y])*Derivative[1, 0][Qtt][x, y] + 
     Qtt[x, y]^2*(12*Qxx[x, y]^3*Qyy[x, y]^2 + 3*Qyy[x, y]^2*
        Derivative[1, 0][Qxx][x, y]^2 - Qxx[x, y]^2*
        Derivative[1, 0][Qyy][x, y]^2 - 2*Qxx[x, y]*Qyy[x, y]^2*
        Derivative[2, 0][Qxx][x, y]))/((-6 + \[Mu]1^2)*Qtt[x, y]*Qxx[x, y]^2*
     Qyy[x, y]^2), Derivative[0, 1][Qyy][x, y] - 
   (2*(-(a[x, y]^2*Qxx[x, y]*Qyy[x, y]^2) + 
      Qtt[x, y]*((-6 + \[Mu]1^2 + 6*Qtt[x, y])*Qxx[x, y]*Qyy[x, y]^2 + 
        Qtt[x, y]*(Derivative[1, 0][Qyy][x, y]^2 - Qyy[x, y]*
           Derivative[2, 0][Qyy][x, y]))))/((-6 + \[Mu]1^2)*Qtt[x, y]*
     Qxx[x, y]*Qyy[x, y]), Derivative[0, 1][Qxz][x, y] - 
   (Qxx[x, y]^3*Qyy[x, y]^3*(2*a[x, y]^2*Qxx[x, y]*
        ((-6 + \[Mu]1^2)*Qxx[x, y]*Qxz[x, y] - 2*Derivative[1, 0][Qtt][x, 
           y]) - Derivative[1, 0][Qtt][x, y]*((-6 + \[Mu]1^2)^2*Qxx[x, y]^2*
          Qxz[x, y]^2 + 4*Derivative[1, 0][Qtt][x, y]^2 - 
         2*(-6 + \[Mu]1^2)*Qxx[x, y]*(2*Derivative[0, 1][Qzz][x, y] - 
           Qxz[x, y]*Derivative[1, 0][Qtt][x, y]))) + 
     Qtt[x, y]*Qxx[x, y]^2*Qyy[x, y]^2*
      (-2*Qyy[x, y]*Derivative[1, 0][Qtt][x, y]^2*Derivative[1, 0][Qxx][x, 
         y] + (-6 + \[Mu]1^2)*Qxx[x, y]^3*Qxz[x, y]*Qyy[x, y]*
        (36 - 14*\[Mu]1^2 + (-6 + \[Mu]1^2)*Derivative[1, 0][Qxz][x, y]) - 
       2*Qxx[x, y]*Qyy[x, y]*(a[x, y]^2*Derivative[1, 0][Qxx][x, y] - 
         2*Derivative[1, 0][Qtt][x, y]*((-6 + \[Mu]1^2)*Qxz[x, y]*
            Derivative[1, 0][Qxx][x, y] + 2*Derivative[2, 0][Qtt][x, y])) + 
       2*Qxx[x, y]^2*(4*a[x, y]*Qyy[x, y]*Derivative[1, 0][a][x, y] + 
         a[x, y]^2*Derivative[1, 0][Qyy][x, y] + Qyy[x, y]*
          (2*Derivative[1, 0][Qtt][x, y]*(4*\[Mu]1^2 + (-6 + \[Mu]1^2)*
              Derivative[1, 0][Qxz][x, y]) - (-6 + \[Mu]1^2)*Qxz[x, y]*
            Derivative[2, 0][Qtt][x, y]))) + Qtt[x, y]^2*Qxx[x, y]*Qyy[x, y]*
      (2*(-6 + \[Mu]1^2)*Qxx[x, y]^4*Qxz[x, y]*Qyy[x, y]*
        (-6 + \[Mu]1^2 + 6*Qyy[x, y]) - 6*Qyy[x, y]^2*Derivative[1, 0][Qtt][
         x, y]*Derivative[1, 0][Qxx][x, y]^2 + Qxx[x, y]^2*
        (-2*(-6 + \[Mu]1^2)*Qyy[x, y]^2*(2*Derivative[1, 0][Qtt][x, y] + 
           Derivative[1, 0][Qxx][x, y]*Derivative[1, 0][Qxz][x, y]) + 
         2*Derivative[1, 0][Qtt][x, y]*Derivative[1, 0][Qyy][x, y]^2) + 
       Qxx[x, y]*Qyy[x, y]^2*(-((-6 + \[Mu]1^2)*Qxz[x, y]*
           Derivative[1, 0][Qxx][x, y]^2) + 4*Derivative[1, 0][Qtt][x, y]*
          Derivative[2, 0][Qxx][x, y]) + Qxx[x, y]^3*
        ((-6 + \[Mu]1^2)*Qxz[x, y]*(2*(-6 + \[Mu]1^2)*Qyy[x, y]^2 + 
           Derivative[1, 0][Qyy][x, y]^2) - 2*Qyy[x, y]*
          (2*(-6 + \[Mu]1^2 + 18*Qyy[x, y])*Derivative[1, 0][Qtt][x, y] + 
           (-6 + \[Mu]1^2)*Qyy[x, y]*Derivative[2, 0][Qxz][x, y]))) - 
     Qtt[x, y]^3*(-3*Qyy[x, y]^3*Derivative[1, 0][Qxx][x, y]^3 + 
       2*Qxx[x, y]^4*Qyy[x, y]*(-6 + \[Mu]1^2 + 6*Qyy[x, y])*
        Derivative[1, 0][Qyy][x, y] + Qxx[x, y]^2*Qyy[x, y]*
        Derivative[1, 0][Qxx][x, y]*(2*(-6 + \[Mu]1^2)*Qyy[x, y]^2 + 
         Derivative[1, 0][Qyy][x, y]^2) + 2*Qxx[x, y]*Qyy[x, y]^3*
        Derivative[1, 0][Qxx][x, y]*Derivative[2, 0][Qxx][x, y] - 
       2*Qxx[x, y]^3*(6*Qyy[x, y]^3*Derivative[1, 0][Qxx][x, y] - 
         Derivative[1, 0][Qyy][x, y]^3 + Qyy[x, y]*Derivative[1, 0][Qyy][x, 
           y]*Derivative[2, 0][Qyy][x, y])))/(2*(-6 + \[Mu]1^2)^2*Qtt[x, y]*
     Qxx[x, y]^5*Qyy[x, y]^3)}}
