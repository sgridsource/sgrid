Comparison with arXiv:0802.1249v2
---------------------------------
My modes differ from the $h^{lm}$ in Eq. (9.3) and the 
$\hat{H}^{lm}$ in Eq. (9.4) in arXiv:0802.1249v2.

Note: Eq. (9.3) of arXiv:0802.1249v2 gives:
$$
h^{lm} = (2 \mu / D) (m \omega)^{2/3}
         \sqrt{16\pi/5)} e^{-im\psi} \hat{H}^{lm}
$$


My modes are then:
$$
h_{WT}^{lm} = h^{lm} m/[r (m \omega)^{2/3}]
            = [2 m_1 m_2 / (r D)]  \sqrt{16\pi/5)} e^{-im\psi} \hat{H}^{lm}
$$
The factor $m/[r (m \omega)^{2/3}]$ is one at leading PN order, but not in
my code.

I checked that my l=m=2 mode agrees.

BUT arXiv:0802.1249v2 has a constant l=2,m=0 mode at leading order! I don't
have that. This is a memory term that depends on the past trajectory...

