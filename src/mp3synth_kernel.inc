// Ο πυρήνας της τράπεζας σύνθεσης, γραμμένος ΜΙΑ φορά και μεταγλωττισμένος μία φορά για κάθε δρόμο (σκέτη C++, SSE2, AVX2+FMA,
// NEON): το mp3synth.cpp ορίζει πριν από κάθε #include τι είναι «διάνυσμα» και πώς γίνονται οι πράξεις του. Έτσι οι δρόμοι δεν
// μπορούν να διαφέρουν στα μαθηματικά, μόνο στο πόσες χρονικές θέσεις ή πόσα δείγματα χωρά ένα διάνυσμα.
//
// Τι πρέπει να έχει ορίσει όποιος το συμπεριλαμβάνει:
//   SYN_NAME(x)          το όνομα x με την κατάληξη του δρόμου
//   SYN_ATTR             ό,τι χρειάζεται ο μεταγλωττιστής για να βγάλει εντολές αυτού του δρόμου (ή τίποτα)
//   SYN_LANES            πόσοι float σε ένα διάνυσμα: 1, 4 ή 8
//   SYN_V                ο τύπος του διανύσματος
//   SYN_ADD SYN_SUB SYN_MUL           a + b · a − b · a·b
//   SYN_MADD(a, b, c)    a·b + c      SYN_NMADD(a, b, c)   c − a·b      (μία εντολή όπου υπάρχει FMA)
//   SYN_ZERO()           όλο μηδέν
//   SYN_TW(p)            SYN_LANES φορές ο ίδιος float, διαβασμένος από τον πίνακα p των σταθερών (στοιχισμένος)
//   SYN_LOAD(p)          SYN_LANES float από ΔΙΚΗ ΜΑΣ μνήμη (πίνακες, SynthState: στοιχισμένη)
//   SYN_STORE(p, v)      SYN_LANES float σε δική μας μνήμη (στοιχισμένη)
//   SYN_STOREU(p, v)     SYN_LANES float σε μνήμη χωρίς στοίχιση (η έξοδος)
//   SYN_LOADU(p)         … και ανάγνωση χωρίς στοίχιση
//   SYN_REV(v)           το διάνυσμα με τις θέσεις του ανάποδα
//   SYN_GATHER4(in, o, x0, x1, x2, x3)     η αναστροφή 4×4 της εισόδου του DCT, βλ. gather παρακάτω
//   SYN_SCATTER4(row, o, q0, q1, q2, q3)   η αναστροφή 4×4 της εξόδου του, βλ. scatter παρακάτω
//   SYN_KEEP(x0, y0, x1, y1)               προαιρετικό: «αυτές οι τέσσερις τιμές του παραθύρου να μείνουν σε καταχωρητές»
//   SYN_WANT_DCT         1 αν θέλει τη συνάρτηση dct (gather, core, scatter μαζί για μία ομάδα)· ο AVX2 τα καλεί χωριστά
//   SYN_WANT_WINDOW      1 αν θέλει και το βήμα του παραθύρου (ο βοηθητικός δρόμος «AVX2 σε μισό πλάτος» θέλει μόνο τον DCT)

// ---------------------------------------------------------------------------------------------------------------------
// ΜΕΡΟΣ 1. Οι 32 καινούριοι αριθμοί κάθε χρονικής θέσης: C[n] = Σ_k S[k]·cos((2k+1)·n·π/64), n = 0 … 31 (DCT-II 32 σημείων).
//
// Κάθε «μεταβλητή» εδώ είναι ένα διάνυσμα με την ΙΔΙΑ ποσότητα για SYN_LANES διαδοχικές χρονικές θέσεις. Γι' αυτό ο κώδικας
// διαβάζεται σαν να ήταν για μία θέση, δεν έχει καμία ανακατάταξη μέσα σε διάνυσμα, και κάθε πράξη του δουλεύει για 8 (ή 4)
// θέσεις μαζί. Η παραγοντοποίηση (κάθε βήμα αποδεικνύεται στο mp3synth.cpp, «Ο DCT»):
//
//   DCT-II N σημείων   =  N/2 αθροίσματα και N/2 διαφορές  →  DCT-II N/2 (οι ζυγές έξοδοι)  +  DCT-IV N/2 (οι μονές)
//   DCT-IV M σημείων   =  M/2 περιστροφές  →  μιγαδικός DFT M/2 σημείων  →  M/2 περιστροφές
//
// Όλα τα βήματα είναι αθροίσματα, διαφορές και περιστροφές (ορθογώνιοι πίνακες): κανένας συντελεστής πάνω από 1, καμία
// αφαίρεση μεγάλων αριθμών που σχεδόν αλληλοαναιρούνται. Γι' αυτό το σφάλμα στρογγύλευσης μένει στα λίγα 1e-7.

#define SYN_INL SYN_INLINE SYN_ATTR
#ifndef SYN_KEEP
#define SYN_KEEP(x0, y0, x1, y1)
#endif
// S = A + B, D = A − B. Διαβάζει και τα δύο πριν γράψει, ώστε τα S, D να μπορούν να είναι τα ίδια τα A, B.
#define SYN_BF(S, D, A, B)                                                                                              \
    {                                                                                                                   \
        const SYN_V bfA = (A), bfB = (B);                                                                               \
        (S) = SYN_ADD(bfA, bfB);                                                                                        \
        (D) = SYN_SUB(bfA, bfB);                                                                                        \
    }

// Η περιστροφή ΠΡΙΝ από τον DFT: (a − j·b)·(c + j·s) = (a·c + b·s) + j·(a·s − b·c), με c = cos φ, s = sin φ στο tw[0], tw[1].
SYN_INL void SYN_NAME(pre)(const float (*tw)[8], SYN_V a, SYN_V b, SYN_V &re, SYN_V &im)
{
    const SYN_V c = SYN_TW(tw[0]), s = SYN_TW(tw[1]);
    re = SYN_MADD(a, c, SYN_MUL(b, s));
    im = SYN_NMADD(b, c, SYN_MUL(a, s));
}

// Η περιστροφή ΜΕΤΑ τον DFT: (r + j·m)·(c + j·s) = (r·c − m·s) + j·(r·s + m·c).
SYN_INL void SYN_NAME(post)(const float (*tw)[8], SYN_V r, SYN_V m, SYN_V &re, SYN_V &im)
{
    const SYN_V c = SYN_TW(tw[0]), s = SYN_TW(tw[1]);
    re = SYN_NMADD(m, s, SYN_MUL(r, c));
    im = SYN_MADD(r, s, SYN_MUL(m, c));
}

// Η ίδια για γωνία π/4, όπου c = s = √½: (r·√½ − m·√½) + j·(r·√½ + m·√½). Το r·√½ υπολογίζεται μία φορά: τρεις πράξεις αντί για
// τέσσερις όπου υπάρχει FMA.
SYN_INL void SYN_NAME(post45)(const SynthConsts &k, SYN_V r, SYN_V m, SYN_V &re, SYN_V &im)
{
    const SYN_V h = SYN_TW(k.half), rh = SYN_MUL(r, h);
    re = SYN_NMADD(m, h, rh);
    im = SYN_MADD(m, h, rh);
}

// Μιγαδικός DFT 4 σημείων με ΘΕΤΙΚΟ εκθέτη: U[p] = Σ_n u[n]·(+j)^(n·p). Με a0 = u0 + u2, a1 = u0 − u2, b0 = u1 + u3,
// b1 = u1 − u3 (οι δύο DFT 2 σημείων των ζυγών και των μονών n):
//   U[0] = a0 + b0     U[2] = a0 − b0     U[1] = a1 + j·b1     U[3] = a1 − j·b1,     όπου j·(x + j·y) = −y + j·x.
SYN_INL void SYN_NAME(dft4)(SYN_V u0r, SYN_V u0i, SYN_V u1r, SYN_V u1i, SYN_V u2r, SYN_V u2i, SYN_V u3r, SYN_V u3i, SYN_V *Ur, SYN_V *Ui)
{
    SYN_V a0r, a0i, a1r, a1i, b0r, b0i, b1r, b1i;
    SYN_BF(a0r, a1r, u0r, u2r)
    SYN_BF(a0i, a1i, u0i, u2i)
    SYN_BF(b0r, b1r, u1r, u3r)
    SYN_BF(b0i, b1i, u1i, u3i)
    SYN_BF(Ur[0], Ur[2], a0r, b0r)
    SYN_BF(Ui[0], Ui[2], a0i, b0i)
    SYN_BF(Ur[3], Ur[1], a1r, b1i)   // πραγματικό μέρος: a1r ∓ b1i  (το «+» πάει στο U[3])
    SYN_BF(Ui[1], Ui[3], a1i, b1r)   // φανταστικό μέρος: a1i ± b1r
}

// DCT-IV M σημείων: X[m] = Σ_k d[k]·cos((2k+1)(2m+1)π/(4M)). Για όλα τα M ο ίδιος τύπος (απόδειξη στο mp3synth.cpp):
//   u[n] = (d[2n] − j·d[M−1−2n])·e^(+jπ(4n+1)/(4M))        n = 0 … M/2 − 1
//   U[p] = Σ_n u[n]·e^(+2πj·n·p/(M/2))                      (DFT M/2 σημείων, θετικός εκθέτης)
//   z[p] = U[p]·e^(+jπp/M)
//   X[2p] = Re z[p]        X[M−1−2p] = Im z[p]
// Η έξοδος X[m] γράφεται στο X[step·m]: έτσι κάθε DCT-IV βάζει τις τιμές του κατευθείαν στις θέσεις C[n] που του ανήκουν.

// M = 2: ένα μόνο u, κανένας DFT, z[0] = u[0].
SYN_INL void SYN_NAME(dct4x2)(const SynthConsts &k, SYN_V d0, SYN_V d1, SYN_V *X, int step)
{
    SYN_NAME(pre)(k.pre2, d0, d1, X[0], X[step]);
}

// M = 4: U[0] = u0 + u1, U[1] = u0 − u1 · z[0] = U[0], z[1] = U[1]·e^(jπ/4).
SYN_INL void SYN_NAME(dct4x4)(const SynthConsts &k, const SYN_V *d, SYN_V *X, int step)
{
    SYN_V u0r, u0i, u1r, u1i, U0r, U0i, U1r, U1i;
    SYN_NAME(pre)(k.pre4 + 0, d[0], d[3], u0r, u0i);
    SYN_NAME(pre)(k.pre4 + 2, d[2], d[1], u1r, u1i);
    SYN_BF(U0r, U1r, u0r, u1r)
    SYN_BF(U0i, U1i, u0i, u1i)
    X[0] = U0r;
    X[3 * step] = U0i;
    SYN_NAME(post45)(k, U1r, U1i, X[2 * step], X[1 * step]);
}

// M = 8: DFT 4 σημείων · z[p] = U[p]·e^(jπp/8), p = 0 … 3 (το p = 0 δεν θέλει πράξη, το p = 2 είναι η γωνία π/4).
SYN_INL void SYN_NAME(dct4x8)(const SynthConsts &k, const SYN_V *d, SYN_V *X, int step)
{
    SYN_V ur[4], ui[4], Ur[4], Ui[4];
    SYN_NAME(pre)(k.pre8 + 0, d[0], d[7], ur[0], ui[0]);
    SYN_NAME(pre)(k.pre8 + 2, d[2], d[5], ur[1], ui[1]);
    SYN_NAME(pre)(k.pre8 + 4, d[4], d[3], ur[2], ui[2]);
    SYN_NAME(pre)(k.pre8 + 6, d[6], d[1], ur[3], ui[3]);
    SYN_NAME(dft4)(ur[0], ui[0], ur[1], ui[1], ur[2], ui[2], ur[3], ui[3], Ur, Ui);
    X[0] = Ur[0];
    X[7 * step] = Ui[0];
    SYN_NAME(post)(k.post8 + 2, Ur[1], Ui[1], X[2 * step], X[5 * step]);
    SYN_NAME(post45)(k, Ur[2], Ui[2], X[4 * step], X[3 * step]);
    SYN_NAME(post)(k.post8 + 6, Ur[3], Ui[3], X[6 * step], X[1 * step]);
}

// M = 16: DFT 8 σημείων από δύο DFT 4 σημείων, E των ζυγών n και O των μονών:
//   U[p] = E[p] + w^p·O[p],   U[p+4] = E[p] − w^p·O[p],   w = e^(+2πj/8) = (1 + j)·√½,   p = 0 … 3
//   w^0 = 1     w^1·(x + jy) = ((x − y) + j(x + y))·√½     w^2·(x + jy) = −y + jx     w^3·(x + jy) = (−(x + y) + j(x − y))·√½
// και μετά z[p] = U[p]·e^(jπp/16), p = 0 … 7 (p = 0 χωρίς πράξη, p = 4 η γωνία π/4).
SYN_INL void SYN_NAME(dct4x16)(const SynthConsts &k, const SYN_V *d, SYN_V *X, int step)
{
    SYN_V ur[8], ui[8], Er[4], Ei[4], Or[4], Oi[4], Ur[8], Ui[8];
    SYN_NAME(pre)(k.pre16 + 0, d[0], d[15], ur[0], ui[0]);
    SYN_NAME(pre)(k.pre16 + 2, d[2], d[13], ur[1], ui[1]);
    SYN_NAME(pre)(k.pre16 + 4, d[4], d[11], ur[2], ui[2]);
    SYN_NAME(pre)(k.pre16 + 6, d[6], d[9], ur[3], ui[3]);
    SYN_NAME(pre)(k.pre16 + 8, d[8], d[7], ur[4], ui[4]);
    SYN_NAME(pre)(k.pre16 + 10, d[10], d[5], ur[5], ui[5]);
    SYN_NAME(pre)(k.pre16 + 12, d[12], d[3], ur[6], ui[6]);
    SYN_NAME(pre)(k.pre16 + 14, d[14], d[1], ur[7], ui[7]);
    SYN_NAME(dft4)(ur[0], ui[0], ur[2], ui[2], ur[4], ui[4], ur[6], ui[6], Er, Ei);
    SYN_NAME(dft4)(ur[1], ui[1], ur[3], ui[3], ur[5], ui[5], ur[7], ui[7], Or, Oi);
    const SYN_V h = SYN_TW(k.half);
    // p = 0
    SYN_BF(Ur[0], Ur[4], Er[0], Or[0])
    SYN_BF(Ui[0], Ui[4], Ei[0], Oi[0])
    // p = 1: w·O = ((Or − Oi) + j(Or + Oi))·√½
    {
        const SYN_V wr = SYN_SUB(Or[1], Oi[1]), wi = SYN_ADD(Or[1], Oi[1]);
        Ur[1] = SYN_MADD(wr, h, Er[1]);
        Ur[5] = SYN_NMADD(wr, h, Er[1]);
        Ui[1] = SYN_MADD(wi, h, Ei[1]);
        Ui[5] = SYN_NMADD(wi, h, Ei[1]);
    }
    // p = 2: w²·O = −Oi + j·Or
    SYN_BF(Ur[6], Ur[2], Er[2], Oi[2])
    SYN_BF(Ui[2], Ui[6], Ei[2], Or[2])
    // p = 3: w³·O = (−(Or + Oi) + j(Or − Oi))·√½
    {
        const SYN_V wr = SYN_ADD(Or[3], Oi[3]), wi = SYN_SUB(Or[3], Oi[3]);
        Ur[3] = SYN_NMADD(wr, h, Er[3]);
        Ur[7] = SYN_MADD(wr, h, Er[3]);
        Ui[3] = SYN_MADD(wi, h, Ei[3]);
        Ui[7] = SYN_NMADD(wi, h, Ei[3]);
    }
    X[0] = Ur[0];
    X[15 * step] = Ui[0];
    SYN_NAME(post)(k.post16 + 2, Ur[1], Ui[1], X[2 * step], X[13 * step]);
    SYN_NAME(post)(k.post16 + 4, Ur[2], Ui[2], X[4 * step], X[11 * step]);
    SYN_NAME(post)(k.post16 + 6, Ur[3], Ui[3], X[6 * step], X[9 * step]);
    SYN_NAME(post45)(k, Ur[4], Ui[4], X[8 * step], X[7 * step]);
    SYN_NAME(post)(k.post16 + 10, Ur[5], Ui[5], X[10 * step], X[5 * step]);
    SYN_NAME(post)(k.post16 + 12, Ur[6], Ui[6], X[12 * step], X[3 * step]);
    SYN_NAME(post)(k.post16 + 14, Ur[7], Ui[7], X[14 * step], X[1 * step]);
}

// Ο DCT για SYN_LANES διαδοχικές χρονικές θέσεις, σε τρία κομμάτια που ο κάθε δρόμος τα βάζει στη σειρά που τον συμφέρει (βλ.
// mp3synth.cpp, «Η διάταξη της δουλειάς»):
//   gather   από την `in` (η πρώτη θέση της ομάδας, 32 υποζώνες η καθεμία) στα x[0 … 31], μία υποζώνη ανά διάνυσμα
//   core     x[0 … 31] → c[0 … 31] = C[0 … 31]   (χαλάει τα x)
//   scatter  από τα c στις σειρές της ιστορίας `row` (η σειρά της πρώτης θέσης της ομάδας) και στο `mid`
//
// SYN_GATHER4(in, o, x0 … x3): οι υποζώνες o … o+3 όλων των θέσεων της ομάδας, μία υποζώνη ανά διάνυσμα (αναστροφή 4×4).
// SYN_SCATTER4(row, o, q0 … q3): το αντίστροφο, τα τέσσερα διανύσματα γίνονται οι τιμές o … o+3 της σειράς κάθε θέσης.
SYN_INL void SYN_NAME(gather)(const float *in, SYN_V *x)
{
    SYN_GATHER4(in, 0, x[0], x[1], x[2], x[3])
    SYN_GATHER4(in, 4, x[4], x[5], x[6], x[7])
    SYN_GATHER4(in, 8, x[8], x[9], x[10], x[11])
    SYN_GATHER4(in, 12, x[12], x[13], x[14], x[15])
    SYN_GATHER4(in, 16, x[16], x[17], x[18], x[19])
    SYN_GATHER4(in, 20, x[20], x[21], x[22], x[23])
    SYN_GATHER4(in, 24, x[24], x[25], x[26], x[27])
    SYN_GATHER4(in, 28, x[28], x[29], x[30], x[31])
}

SYN_INL void SYN_NAME(core)(const SynthConsts &k, SYN_V *x, SYN_V *c)
{
    SYN_V d1[16], d2[8], d3[4];
    // N = 32: x[k] ← x[k] + x[31−k] (για τις ζυγές C), d1[k] = x[k] − x[31−k] (για τις μονές). C[2m+1] = DCT-IV16(d1)[m].
    SYN_BF(x[0], d1[0], x[0], x[31])
    SYN_BF(x[1], d1[1], x[1], x[30])
    SYN_BF(x[2], d1[2], x[2], x[29])
    SYN_BF(x[3], d1[3], x[3], x[28])
    SYN_BF(x[4], d1[4], x[4], x[27])
    SYN_BF(x[5], d1[5], x[5], x[26])
    SYN_BF(x[6], d1[6], x[6], x[25])
    SYN_BF(x[7], d1[7], x[7], x[24])
    SYN_BF(x[8], d1[8], x[8], x[23])
    SYN_BF(x[9], d1[9], x[9], x[22])
    SYN_BF(x[10], d1[10], x[10], x[21])
    SYN_BF(x[11], d1[11], x[11], x[20])
    SYN_BF(x[12], d1[12], x[12], x[19])
    SYN_BF(x[13], d1[13], x[13], x[18])
    SYN_BF(x[14], d1[14], x[14], x[17])
    SYN_BF(x[15], d1[15], x[15], x[16])
    SYN_NAME(dct4x16)(k, d1, c + 1, 2);

    // N = 16 πάνω στα 16 αθροίσματα: C[4m+2] = DCT-IV8(d2)[m].
    SYN_BF(x[0], d2[0], x[0], x[15])
    SYN_BF(x[1], d2[1], x[1], x[14])
    SYN_BF(x[2], d2[2], x[2], x[13])
    SYN_BF(x[3], d2[3], x[3], x[12])
    SYN_BF(x[4], d2[4], x[4], x[11])
    SYN_BF(x[5], d2[5], x[5], x[10])
    SYN_BF(x[6], d2[6], x[6], x[9])
    SYN_BF(x[7], d2[7], x[7], x[8])
    SYN_NAME(dct4x8)(k, d2, c + 2, 4);

    // N = 8: C[8m+4] = DCT-IV4(d3)[m].
    SYN_BF(x[0], d3[0], x[0], x[7])
    SYN_BF(x[1], d3[1], x[1], x[6])
    SYN_BF(x[2], d3[2], x[2], x[5])
    SYN_BF(x[3], d3[3], x[3], x[4])
    SYN_NAME(dct4x4)(k, d3, c + 4, 8);

    // N = 4: C[8] και C[24] = DCT-IV2 των δύο διαφορών.
    {
        SYN_V e0, e1;
        SYN_BF(x[0], e0, x[0], x[3])
        SYN_BF(x[1], e1, x[1], x[2])
        SYN_NAME(dct4x2)(k, e0, e1, c + 8, 16);
    }
    // N = 2: C[0] = x0 + x1, C[16] = (x0 − x1)·cos(π/4).
    c[0] = SYN_ADD(x[0], x[1]);
    c[16] = SYN_MUL(SYN_SUB(x[0], x[1]), SYN_TW(k.half));
}

// Η σειρά μιας θέσης (βλ. mp3synth.cpp, «Το παράθυρο»):   row[m] = C[17 + m], m = 0 … 14 · row[15] = 0 · row[16 + m] = C[15 − m],
// m = 0 … 15 · και mid = C[16].
SYN_INL void SYN_NAME(scatter)(const SYN_V *c, float *row, float *mid)
{
    SYN_SCATTER4(row, 0, c[17], c[18], c[19], c[20])
    SYN_SCATTER4(row, 4, c[21], c[22], c[23], c[24])
    SYN_SCATTER4(row, 8, c[25], c[26], c[27], c[28])
    SYN_SCATTER4(row, 12, c[29], c[30], c[31], SYN_ZERO())
    SYN_SCATTER4(row, 16, c[15], c[14], c[13], c[12])
    SYN_SCATTER4(row, 20, c[11], c[10], c[9], c[8])
    SYN_SCATTER4(row, 24, c[7], c[6], c[5], c[4])
    SYN_SCATTER4(row, 28, c[3], c[2], c[1], c[0])
    SYN_STOREU(mid, c[16]);
}

#if SYN_WANT_DCT
// Τα τρία μαζί, για μία ομάδα.
SYN_ATTR static void SYN_NAME(dct)(const SynthConsts &k, const float *in, float *row, float *mid)
{
    SYN_V x[32], c[32];
    SYN_NAME(gather)(in, x);
    SYN_NAME(core)(k, x, c);
    SYN_NAME(scatter)(c, row, mid);
}
#endif

#if SYN_WANT_WINDOW
// ---------------------------------------------------------------------------------------------------------------------
// ΜΕΡΟΣ 2. Το παράθυρο: 16 όροι για καθένα από τα 32 δείγματα εξόδου, για τις 18 θέσεις του granule.
//
// `old`: οι 18 σειρές του προηγούμενου granule · `cur`: οι 18 αυτού · `midOld`, `midCur`: τα C[16] τους · `out`: 576 δείγματα.
// Η θέση s χρειάζεται τις 16 σειρές s − 15 … s. Όσες έχουν αριθμό ≥ 0 είναι στο `cur`, οι αρνητικές είναι οι τελευταίες του `old`
// (η −1 είναι η σειρά 17 του). SYN_ROW(q) είναι η σειρά s − 15 + q, q = 0 … 16 (η 16 για τη θέση s + 1).
// Για τη θέση s, ο όρος t (0 = η ίδια η θέση, 15 = η πιο παλιά) είναι η σειρά s + 15 − t. Με G[t][m] = D[32t + m + 1] και
// j = m + 1 = 1 … 16 (απόδειξη στο mp3synth.cpp, «Το παράθυρο»):
//
//   out[j]      = Σ_{t ζυγό} G[t][m]·E_t[m]  −  Σ_{t μονό} G[t][m]·O_t[m]              E = το πρώτο μισό της σειράς (C[17+m])
//   out[32 − j] = Σ_t G[15−t][m]·X_t[m],   X_t = E_t για t ζυγό, O_t για t μονό        O = το δεύτερο μισό (C[15−m])
//
// Οι όροι t = i και t = 15 − i (i = 0 … 7) χρειάζονται τα ΙΔΙΑ τέσσερα διανύσματα, a = G[i], b = G[15−i], x = X_i, y = X_(15−i):
//   out[j]      παίρνει  ±(a·x − b·y)      («+» για i ζυγό)
//   out[32 − j] παίρνει    b·x + a·y
// άρα τέσσερα γινόμενα από τέσσερα διανύσματα. Οι τέσσερις αθροιστές (la, lb, ua, ub, με out[j] = la − lb και out[32−j] = ua + ub)
// είναι ανεξάρτητοι μεταξύ τους, ώστε ο επεξεργαστής να τρέχει τέσσερις κοντές αλυσίδες μαζί αντί για μία τέσσερις φορές πιο
// μακριά. Η σειρά i = 0 … 7 πάει από τις άκρες του παραθύρου (μικροί συντελεστές) προς το κέντρο (μεγάλοι): οι μικροί όροι
// αθροίζονται πρώτοι, που είναι η σειρά με το μικρότερο σφάλμα στρογγύλευσης.
// Δύο διαδοχικές θέσεις (s και s + 1) σε κάθε γύρο: οι σειρές του παραθύρου a, b είναι οι ίδιες για τις δύο, άρα διαβάζονται μία
// φορά για οκτώ γινόμενα, και οι αλυσίδες που τρέχουν μαζί είναι οκτώ. Η σειρά της θέσης s + 1 είναι 32 τιμές πιο πέρα.
#define SYN_PAIR(i, XOFF, YOFF, ACC)                                                                                     \
    {                                                                                                                   \
        const SYN_V a = SYN_LOAD(k.g[i] + B), b = SYN_LOAD(k.g[15 - i] + B);                                            \
        SYN_V x0 = SYN_LOAD(SYN_ROW(15 - i) + XOFF + B), y0 = SYN_LOAD(SYN_ROW(i) + YOFF + B);                          \
        SYN_V x1 = SYN_LOAD(SYN_ROW(16 - i) + XOFF + B), y1 = SYN_LOAD(SYN_ROW(i + 1) + YOFF + B);                      \
        SYN_KEEP(x0, y0, x1, y1)                                                                                        \
        la0 = ACC(a, x0, la0);                                                                                          \
        lb0 = ACC(b, y0, lb0);                                                                                          \
        ua0 = SYN_MADD(b, x0, ua0);                                                                                     \
        ub0 = SYN_MADD(a, y0, ub0);                                                                                     \
        la1 = ACC(a, x1, la1);                                                                                          \
        lb1 = ACC(b, y1, lb1);                                                                                          \
        ua1 = SYN_MADD(b, x1, ua1);                                                                                     \
        ub1 = SYN_MADD(a, y1, ub1);                                                                                     \
    }
#define SYN_PAIR_EVEN(i) SYN_PAIR(i, 0, 16, SYN_MADD)    // i ζυγό: x = E της σειράς του t = i, y = O της σειράς του t = 15 − i
#define SYN_PAIR_ODD(i) SYN_PAIR(i, 16, 0, SYN_NMADD)    // i μονό: x = O, y = E, και τα δύο γινόμενα του out[j] με «−»

// Το σώμα για δύο θέσεις (s, s + 1).
//   m = B … B + SYN_LANES − 1:  out[1 + m] = la − lb,  out[31 − m] = ua + ub (γι' αυτό γράφεται ανάποδα). Το out[16] γράφεται δύο
//   φορές (m = 15 από τις δύο μεριές) με την ίδια ακριβώς τιμή: εκεί E = 0 και G[15−t] = −G[t], άρα οι δύο τύποι έχουν τα ίδια
//   γινόμενα με την ίδια σειρά. Στο i = 0 οι αθροιστές αρχίζουν από το πρώτο γινόμενο, όχι από το μηδέν.
#define SYN_TWO_BODY \
    for (int B = 0; B < 16; B += SYN_LANES) { \
        SYN_V la0, lb0, ua0, ub0, la1, lb1, ua1, ub1; \
        { \
            const SYN_V a = SYN_LOAD(k.g[0] + B), b = SYN_LOAD(k.g[15] + B); \
            const SYN_V x0 = SYN_LOAD(SYN_ROW(15) + B), y0 = SYN_LOAD(SYN_ROW(0) + 16 + B); \
            const SYN_V x1 = SYN_LOAD(SYN_ROW(16) + B), y1 = SYN_LOAD(SYN_ROW(1) + 16 + B); \
            la0 = SYN_MUL(a, x0); \
            lb0 = SYN_MUL(b, y0); \
            ua0 = SYN_MUL(b, x0); \
            ub0 = SYN_MUL(a, y0); \
            la1 = SYN_MUL(a, x1); \
            lb1 = SYN_MUL(b, y1); \
            ua1 = SYN_MUL(b, x1); \
            ub1 = SYN_MUL(a, y1); \
        } \
        SYN_PAIR_ODD(1) \
        SYN_PAIR_EVEN(2) \
        SYN_PAIR_ODD(3) \
        SYN_PAIR_EVEN(4) \
        SYN_PAIR_ODD(5) \
        SYN_PAIR_EVEN(6) \
        SYN_PAIR_ODD(7) \
        SYN_STOREU(o + 1 + B, SYN_SUB(la0, lb0)); \
        SYN_STOREU(o + 32 - SYN_LANES - B, SYN_REV(SYN_ADD(ua0, ub0))); \
        SYN_STOREU(o + 32 + 1 + B, SYN_SUB(la1, lb1)); \
        SYN_STOREU(o + 32 + 32 - SYN_LANES - B, SYN_REV(SYN_ADD(ua1, ub1))); \
    }

// Οι θέσεις s και s + 1 όταν οι σειρές τους είναι μοιρασμένες στα δύο μισά: `ro` + 32q είναι η σειρά q αν βρίσκεται στο παλιό
// μισό, `rc` + 32q αν βρίσκεται στο νέο, και στο παλιό είναι όσες έχουν q < q0 (q0 = 15 − s). Η επιλογή γίνεται την ώρα που
// τρέχει, ΕΝΑ αντίγραφο του κώδικα για όλα τα s (βλ. mp3synth.cpp, «Τι μετρήθηκε»). 64 δείγματα στο `o`, εκτός από τα o[0] και
// o[32] (τα γράφει ο καλών).
#define SYN_ROW(q) (((q) >= q0 ? rc : ro) + 32 * (q))
SYN_INL void SYN_NAME(windowTwoSplit)(const SynthConsts &k, const float *ro, const float *rc, const int q0, float *o)
{
    SYN_TWO_BODY
}
#undef SYN_ROW

// Οι ίδιες δύο θέσεις όταν οι 17 σειρές τους είναι συνεχόμενες στη μνήμη από το `r` και πέρα (βλ. window).
#define SYN_ROW(q) (r + 32 * (q))
SYN_INL void SYN_NAME(windowTwoLinear)(const SynthConsts &k, const float *r, float *o)
{
    SYN_TWO_BODY
}
#undef SYN_ROW

SYN_ATTR static void SYN_NAME(window)(const SynthConsts &k, const float *old, const float *cur, const float *mid, float *out)
{
    // out[0] = Σ_t (±D[32t])·C_t[16] («+» για t ζυγό): 16 όροι πάνω στα mid, για SYN_LANES θέσεις μαζί (κατά μήκος του χρόνου,
    // όπως ο DCT). `mid`: τα 33 C[16] στη σειρά του χρόνου (15 του προηγούμενου granule, 18 αυτού), βλ. SYN_MID_BEGIN.
    // Το g0 έχει τους συντελεστές γυρισμένους στη σειρά της μνήμης. Η τελευταία ομάδα πιάνει τις τελευταίες SYN_LANES θέσεις,
    // όποιες κι αν είναι, και ξαναϋπολογίζει όσες χρειαστεί.
    float first[18];
    for (int s = 0;; s += SYN_LANES) {
        if (s > 18 - SYN_LANES) s = 18 - SYN_LANES;
        SYN_V acc = SYN_MUL(SYN_TW(k.g0[0]), SYN_LOADU(mid + s));
        for (int u = 1; u < 16; ++u) acc = SYN_MADD(SYN_TW(k.g0[u]), SYN_LOADU(mid + s + u), acc);
        SYN_STOREU(first + s, acc);
        if (s == 18 - SYN_LANES) break;
    }
    // Κάθε δεύτερη κλήση το νέο μισό βρίσκεται στη μνήμη αμέσως μετά το παλιό: οι 33 σειρές είναι τότε συνεχόμενες. Στις άλλες
    // κλήσεις κάθε ανάγνωση διαλέγει μισό (windowTwoSplit).
    if (old + 18 * 32 == cur) {
        for (int s = 0; s < 18; s += 2) {
            SYN_NAME(windowTwoLinear)(k, cur + 32 * (s - 15), out + 32 * s);
            out[32 * s] = first[s];
            out[32 * s + 32] = first[s + 1];
        }
        return;
    }
    const float *const base = old < cur ? old : cur;   // τα δύο μισά είναι ο ίδιος πίνακας (SynthState::rows)
    const std::ptrdiff_t oo = old - base, co = cur - base;   // ptrdiff_t, όχι long: στα Windows το long είναι 32 bit (MSVC C4244)
    for (int s = 0; s < 18; s += 2) {
        SYN_NAME(windowTwoSplit)(k, base + (oo + 32 * (s + 3)), base + (co + 32 * (s - 15)), 15 - s, out + 32 * s);
        out[32 * s] = first[s];
        out[32 * s + 32] = first[s + 1];
    }
}
#undef SYN_PAIR
#undef SYN_PAIR_EVEN
#undef SYN_PAIR_ODD
#endif   // SYN_WANT_WINDOW

#undef SYN_INL
#undef SYN_BF
#undef SYN_KEEP

// Ό,τι όρισε ο δρόμος για αυτή τη μεταγλώττιση σβήνεται εδώ, ώστε ο επόμενος να αρχίσει από καθαρό χαρτί.
#undef SYN_NAME
#undef SYN_ATTR
#undef SYN_LANES
#undef SYN_V
#undef SYN_ADD
#undef SYN_SUB
#undef SYN_MUL
#undef SYN_MADD
#undef SYN_NMADD
#undef SYN_ZERO
#undef SYN_TW
#undef SYN_LOAD
#undef SYN_LOADU
#undef SYN_STORE
#undef SYN_STOREU
#undef SYN_REV
#undef SYN_SHUF
#undef SYN_GATHER4
#undef SYN_SCATTER4
#undef SYN_WANT_WINDOW
#undef SYN_WANT_DCT
