import numpy as np
import math

def damppole(r,rorder,alphai,alphak,dmpi,dmpk,dmpik,pentyp):
    eps = 0.001
    diff = abs(alphai-alphak)
    dampi = alphai * r
    dampk = alphak * r
    expi = math.exp(-dampi)
    expk = math.exp(-dampk)
    # core-valence charge penetration damping for Gordon f1
    if pentyp == "GORDON1":
        dampi2 = dampi * dampi
        dampi3 = dampi * dampi2
        dampi4 = dampi2 * dampi2
        dampi5 = dampi2 * dampi3
        dmpi[0] = 1. - (1. + 0.5*dampi)*expi
        dmpi[2] = 1. - (1. + dampi + 0.5*dampi2)*expi
        dmpi[4] = 1. - (1. + dampi + 0.5*dampi2 \
                     + dampi3/6.)*expi
        dmpi[6] = 1. - (1. + dampi + 0.5*dampi2 \
                     + dampi3/6. + dampi4/30.)*expi
        dmpi[8] = 1. - (1. + dampi + 0.5*dampi2 \
                     + dampi3/6. + 4.*dampi4/105. \
                     + dampi5/210.)*expi
        if diff < eps:
            dmpk[0] = dmpi[0]
            dmpk[2] = dmpi[2]
            dmpk[4] = dmpi[4]
            dmpk[6] = dmpi[6]
            dmpk[8] = dmpi[8]
        else:
            dampk2 = dampk * dampk
            dampk3 = dampk * dampk2
            dampk4 = dampk2 * dampk2
            dampk5 = dampk2 * dampk3
            dmpk[0] = 1. - (1. + 0.5*dampk)*expk
            dmpk[2] = 1. - (1. + dampk + 0.5*dampk2)*expk
            dmpk[4] = 1. - (1. + dampk + 0.5*dampk2 \
                         + dampk3/6.)*expk
            dmpk[6] = 1. - (1. + dampk + 0.5*dampk2 \
                         + dampk3/6. + dampk4/30.)*expk
            dmpk[8] = 1. - (1. + dampk + 0.5*dampk2 \
                         + dampk3/6. + 4.*dampk4/105. \
                         + dampk5/210.)*expk
        # valence-valence charge penetration damping for Gordon f1
        if diff < eps:
            dampi6 = dampi3 * dampi3
            dampi7 = dampi3 * dampi4
            dmpik[0] = 1. - (1. + 11.*dampi/16. \
                          + 3.*dampi2/16. \
                          + dampi3/48.)*expi
            dmpik[2] = 1. - (1. + dampi + 0.5*dampi2 \
                          + 7.*dampi3/48. \
                          + dampi4/48.)*expi
            dmpik[4] = 1. - (1. + dampi + 0.5*dampi2 \
                          + dampi3/6. + dampi4/24. \
                          + dampi5/144.)*expi
            dmpik[6] = 1. - (1. + dampi + 0.5*dampi2 \
                          + dampi3/6. + dampi4/24. \
                          + dampi5/120. + dampi6/720.)*expi
            dmpik[8] = 1. - (1. + dampi + 0.5*dampi2 \
                          + dampi3/6. + dampi4/24. \
                          + dampi5/120. + dampi6/720. \
                          + dampi7/5040.)*expi
            if rorder >= 11:
                dampi8 = dampi4 * dampi4
                dmpik[10] = 1. - (1. + dampi + 0.5*dampi2 \
                               + dampi3/6. + dampi4/24. \
                               + dampi5/120. + dampi6/720. \
                               + dampi7/5040. + dampi8/45360.)*expi
        else:
            alphai2 = alphai * alphai
            alphak2 = alphak * alphak
            termi = alphak2 / (alphak2-alphai2)
            termk = alphai2 / (alphai2-alphak2)
            termi2 = termi * termi
            termk2 = termk * termk
            dmpik[0] = 1. - termi2*(1. + 2.*termk \
                          + 0.5*dampi)*expi \
                       - termk2*(1. + 2.*termi \
                            + 0.5*dampk)*expk
            dmpik[2] = 1. - termi2*(1.+dampi+0.5*dampi2)*expi \
                          - termk2*(1.+dampk+0.5*dampk2)*expk \
                          - 2.*termi2*termk*(1.+dampi)*expi \
                          - 2.*termk2*termi*(1.+dampk)*expk
            dmpik[4] = 1. - termi2*(1. + dampi + 0.5*dampi2 \
                          + dampi3/6.)*expi \
                       - termk2*(1. + dampk + 0.5*dampk2 \
                            + dampk3/6.)*expk \
                       - 2.*termi2*termk \
                            *(1. + dampi + dampi2/3.)*expi \
                       - 2.*termk2*termi \
                            *(1. + dampk + dampk2/3.)*expk
            dmpik[6] = 1. - termi2*(1. + dampi + 0.5*dampi2 \
                          + dampi3/6. + dampi4/30.)*expi \
                       - termk2*(1. + dampk + 0.5*dampk2 \
                            + dampk3/6. + dampk4/30.)*expk \
                       - 2.*termi2*termk*(1. + dampi \
                            + 2.*dampi2/5. + dampi3/15.)*expi \
                       - 2.*termk2*termi*(1. + dampk \
                            + 2.*dampk2/5. + dampk3/15.)*expk
            dmpik[8] = 1. - termi2*(1. + dampi + 0.5*dampi2 \
                          + dampi3/6. + 4.*dampi4/105. \
                          + dampi5/210.)*expi \
                       - termk2*(1. + dampk + 0.5*dampk2 \
                            + dampk3/6. + 4.*dampk4/105. \
                            + dampk5/210.)*expk \
                       - 2.*termi2*termk*(1. + dampi \
                            + 3.*dampi2/7. \
                            + 2.*dampi3/21. \
                            + dampi4/105.)*expi \
                       - 2.*termk2*termi*(1. + dampk \
                            + 3.*dampk2/7. \
                            + 2.*dampk3/21. \
                            + dampk4/105.)*expk
            if rorder >= 11:
                dampi6 = dampi3 * dampi3
                dampk6 = dampk3 * dampk3
                dmpik[10] = 1. - termi2*(1. + dampi \
                               + 0.5*dampi2 + dampi3/6. \
                               + 5.*dampi4/126. \
                               + 2.*dampi5/315. \
                               + dampi6/1890.)*expi \
                            - termk2*(1. + dampk \
                                 + 0.5*dampk2 + dampk3/6. \
                                 + 5.*dampk4/126. \
                                 + 2.*dampk5/315. \
                                 + dampk6/1890.)*expk \
                            - 2.*termi2*termk*(1. + dampi \
                                 + 4.*dampi2/9. + dampi3/9. \
                                 + dampi4/63. + dampi5/945.)*expi \
                            - 2.*termk2*termi*(1. + dampk \
                                 + 4.*dampk2/9. + dampk3/9. \
                                 + dampk4/63. + dampk5/945.)*expk 
    # core-valence charge penetration damping for Gordon f2
    elif pentyp == "GORDON2":
        dampi2 = dampi * dampi
        dampi3 = dampi * dampi2
        dmpi[0] = 1. - expi
        dmpi[2] = 1. - (1. + dampi)*expi
        dmpi[4] = 1. - (1. + dampi + dampi2/3.)*expi
        dmpi[6] = 1. - (1. + dampi + 0.4*dampi2 \
                     + dampi3/15.)*expi
        if diff < eps:
            dmpk[0] = dmpi[0]
            dmpk[2] = dmpi[2]
            dmpk[4] = dmpi[4]
            dmpk[6] = dmpi[6]
        else:
            dampk2 = dampk * dampk
            dampk3 = dampk * dampk2
            dmpk[0] = 1. - expk
            dmpk[2] = 1. - (1. + dampk)*expk
            dmpk[4] = 1. - (1. + dampk + dampk2/3.)*expk
            dmpk[6] = 1. - (1. + dampk + 0.4*dampk2 \
                          + dampk3/15.)*expk
        # valence-valence charge penetration damping for Gordon f2
        dampi4 = dampi2 * dampi2
        dampi5 = dampi2 * dampi3
        if diff < eps:
            dampi6 = dampi3 * dampi3
            dmpik[0] = 1. - (1. + 0.5*dampi)*expi
            dmpik[2] = 1. - (1. + dampi + 0.5*dampi2)*expi
            dmpik[4] = 1. - (1. + dampi + 0.5*dampi2 \
                          + dampi3/6.)*expi
            dmpik[6] = 1. - (1. + dampi + 0.5*dampi2 \
                          + dampi3/6. + dampi4/30.)*expi
            dmpik[8] = 1. - (1. + dampi + 0.5*dampi2 \
                          + dampi3/6. + 4.*dampi4/105. \
                          + dampi5/210.)*expi
            if rorder >= 11:
                dmpik[10] = 1. - (1. + dampi + 0.5*dampi2 \
                               + dampi3/6. + 5.*dampi4/126. \
                               + 2.*dampi5/315. \
                               + dampi6/1890.)*expi
        else:
            dampk4 = dampk2 * dampk2
            dampk5 = dampk2 * dampk3
            alphai2 = alphai * alphai
            alphak2 = alphak * alphak
            termi = alphak2 / (alphak2-alphai2)
            termk = alphai2 / (alphai2-alphak2)
            dmpik[0] = 1. - termi*expi - termk*expk
            dmpik[2] = 1. - termi*(1. + dampi)*expi \
                          - termk*(1. + dampk)*expk
            dmpik[4] = 1. - termi*(1. + dampi + dampi2/3.)*expi \
                          - termk*(1. + dampk + dampk2/3.)*expk
            dmpik[6] = 1. - termi*(1. + dampi + 0.4*dampi2 \
                          + dampi3/15.)*expi \
                       - termk*(1. + dampk + 0.4*dampk2 \
                            + dampk3/15.)*expk
            dmpik[8] = 1. - termi*(1. + dampi + 3.*dampi2/7. \
                          + 2.*dampi3/21. + dampi4/105.)*expi \
                       - termk*(1. + dampk + 3.*dampk2/7. \
                            + 2.*dampk3/21. + dampk4/105.)*expk
            if rorder >= 11:
                dmpik[10] = 1. - termi*(1. + dampi \
                               + 4.*dampi2/9. + dampi3/9. \
                               + dampi4/63. + dampi5/945.)*expi \
                            - termk*(1. + dampk \
                                 + 4.*dampk2/9. + dampk3/9. \
                                 + dampk4/63. + dampk5/945.)*expk
