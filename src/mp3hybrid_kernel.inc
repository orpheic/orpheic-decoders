// Ο πυρήνας του υβριδικού σταδίου, γραμμένος ΜΙΑ φορά πάνω σε έναν μικρό «τύπο διανύσματος». Το mp3hybrid.cpp τον
// συμπεριλαμβάνει μία φορά ανά δρόμο επεξεργαστή, μέσα σε δικό του namespace, αφού ορίσει:
//   HK_LANES            πόσες υποζώνες χωρά ένα διάνυσμα (1, 4 ή 8)
//   V                   ο τύπος (float, __m128, __m256, float32x4_t)
//   HK_FN               τι πρέπει να γράφει κάθε συνάρτηση για να επιτρέπονται οι εντολές του δρόμου (target(...)) ή τίποτα
//   vld, vst            ανάγνωση και γραφή HK_LANES συνεχόμενων float ΧΩΡΙΣ υπόθεση στοίχισης (δεδομένα)
//   vk                  ανάγνωση σταθεράς από τους πίνακες (στοιχισμένοι στα 32 byte)
//   vadd, vsub, vmul    πράξεις ανά θέση
//   vmadd(a, b, c)      a·b + c        vmsub(a, b, c)   a·b − c        vnmadd(a, b, c)   c − a·b
//   vzero()             μηδενικό διάνυσμα
//   HK_ENTER()          ό,τι χρειάζεται στην είσοδο των δύο μεγάλων συναρτήσεων (στοίχιση της στοίβας στον AVX)
//   HK_LEAVE()          ό,τι χρειάζεται πριν από την επιστροφή (vzeroupper στον AVX)
// Οι αριθμοί «Βήμα 1 … 4» στα σχόλια είναι της παραγωγής στην κορυφή του mp3hybrid.cpp.
//
// Επειδή κάθε «αριθμός» εδώ είναι ένα διάνυσμα με μία θέση ανά υποζώνη, δεν υπάρχει πουθενά ανακάτεμα θέσεων: μόνο
// αναγνώσεις, πράξεις ανά θέση, γραφές. Γι' αυτό ο ίδιος κώδικας βγάζει και τους τέσσερις δρόμους.
//
// ΓΙΑΤΙ ΟΛΑ ΕΙΝΑΙ ΓΡΑΜΜΕΝΑ ΑΝΑΠΤΥΓΜΕΝΑ, με ονόματα X0 … X17 αντί για πίνακες και βρόχους: με -O2 ο GCC δεν ξετυλίγει τους
// μικρούς βρόχους και κρατά τους πίνακες στη στοίβα (μετρημένο στην πρώτη μορφή του πυρήνα: 95 ns ανά μακρύ μπλοκ αντί 67).
// Έτσι κάθε τιμή είναι ένας καταχωρητής σε κάθε μεταγλωττιστή, χωρίς να βασιζόμαστε σε διακόπτες.

#if defined(MP3_HYBRID_MUTANT) && ((MP3_HYBRID_MUTANT == 16 && HK_LANES == 8) || (MP3_HYBRID_MUTANT == 17 && HK_LANES == 4) || (MP3_HYBRID_MUTANT == 18 && HK_LANES == 1))
#define HK_MUTANT_PATH 1   // σκόπιμο σφάλμα ΜΟΝΟ σε αυτόν τον δρόμο: δείχνει ότι η δοκιμή κρίνει κάθε δρόμο χωριστά
#else
#define HK_MUTANT_PATH 0
#endif

// Ο DCT-II 9 σημείων (Βήμα 3), c_l = Σ p_k·cos(l·(2k+1)·10°), σε δύο μισά που δεν έχουν τίποτα κοινό: μετά το δίπλωμα
// k ↔ 8 − k τα ζυγά l βλέπουν μόνο τα αθροίσματα s_k = p_k + p_(8−k) (και s4 = p4) και τα περιττά μόνο τις διαφορές
// d_k = p_k − p_(8−k). Ολόκληρος: 16 πολλαπλασιασμοί και 32 προσθέσεις (μαζί με τις 8 του διπλώματος)· όπου υπάρχει FMA
// οι 14 προσθέσεις είναι η ίδια εντολή με έναν πολλαπλασιασμό: 34 εντολές.

// Το ζυγό μισό: από τα s0 … s4 τα c0, c2, c4, c6, c8.
HK_FN HK_INLINE void dct9even(V s0, V s1, V s2, V s3, V s4, V &c0, V &c2, V &c4, V &c6, V &c8)
{
    const Tables &t = kT;
    // l = 0 και l = 6: συντελεστές μόνο 1 και ½ (cos 60°, cos 180°, cos 300°, cos 420°).
    const V a = vadd(vadd(s0, s2), s3), b = vadd(s1, s4);
    c0 = vadd(a, b);
    c6 = vmsub(a, vk(t.half), b);
    // l = 2, 4, 8: ο όρος του s1 είναι ±½, του s4 ±1· οι τρεις άλλοι συντελεστές είναι cos 20°, 40°, 80° σε κύκλο, και
    // επειδή cos 20° = cos 40° + cos 80° γράφονται με δύο μόνο διαφορές u, v.
    const V e = vmsub(s1, vk(t.half), s4);
    const V u = vsub(s0, s2), v = vsub(s3, s2);
    c2 = vnmadd(v, vk(t.c40), vmadd(u, vk(t.c20), e));   //  e + c20·u − c40·v
    c4 = vmadd(v, vk(t.c80), vmsub(u, vk(t.c40), e));    // −e + c40·u + c80·v
    c8 = vnmadd(v, vk(t.c20), vmsub(u, vk(t.c80), e));   // −e + c80·u − c20·v
}

// Το περιττό μισό: από τα d0 … d3 τα c1, c3, c5, c7.
HK_FN HK_INLINE void dct9odd(V d0, V d1, V d2, V d3, V &c1, V &c3, V &c5, V &c7)
{
    const Tables &t = kT;
    // l = 3: cos 30°, 90°, 150°, 210° = c30·(1, 0, −1, −1).
#if HK_MUTANT_PATH
    c3 = vmul(vsub(vsub(d0, d2), d3), vk(t.c10));
#else
    c3 = vmul(vsub(vsub(d0, d2), d3), vk(t.c30));
#endif
    // l = 1, 5, 7: ο όρος του d1 είναι ±cos 30°· οι άλλοι cos 10°, 50°, 70° σε κύκλο, με cos 10° = cos 50° + cos 70°.
    const V g = vmul(d1, vk(t.c30));
    const V up = vadd(d0, d3), vp = vsub(d3, d2);
    c1 = vnmadd(vp, vk(t.c50), vmadd(up, vk(t.c10), g));   //  g + c10·u' − c50·v'
    c5 = vmadd(vp, vk(t.c70), vmsub(up, vk(t.c50), g));    // −g + c50·u' + c70·v'
    c7 = vnmadd(vp, vk(t.c10), vmsub(up, vk(t.c70), g));   // −g + c70·u' − c10·v'
}

// Οι δύο θέσεις της εξόδου που παίρνουν την τιμή C_j (j = 9 … 17): σειρές j − 9 και 26 − j. Παράθυρο, πρόσημο και η
// πρόσθεση της παλιάς επικάλυψης σε μία εντολή ανά θέση (Βήμα 4).
HK_FN HK_INLINE void longOut(V c, int j, const LongWin &w, int par, const float *HK_RESTRICT ov, float *HK_RESTRICT out)
{
    vst(out + 32 * (j - 9), vmadd(c, vk(w.a[j - 9] + par), vld(ov + 32 * (j - 9))));
    vst(out + 32 * (26 - j), vmadd(c, vk(w.b[j - 9] + par), vld(ov + 32 * (26 - j))));
}
// Οι δύο θέσεις της νέας επικάλυψης που παίρνουν την τιμή C_j (j = 0 … 8): σειρές 8 − j και 9 + j.
HK_FN HK_INLINE void longOv(V c, int j, const LongWin &w, int par, float *HK_RESTRICT ov)
{
    vst(ov + 32 * (8 - j), vmul(c, vk(w.c[j] + par)));
    vst(ov + 32 * (9 + j), vmul(c, vk(w.d[j] + par)));
}

// Μείωση αναδίπλωσης (2.4.3.4.10.1) για το ζεύγος σειρών i και m = 17 − i, χωρίς να γραφτεί τίποτα πίσω στο x. Το πρότυπο
// στρέφει το ζεύγος
//   lo = x[17 − i][sb − 1], hi = x[i][sb]   σε   lo·cs − hi·ca,  hi·cs + lo·ca        (sb = 1 … 31, i = 0 … 7).
// Κάθε τιμή μπαίνει σε ΕΝΑ μόνο ζεύγος (οι σειρές 0 … 7 με την από κάτω υποζώνη, οι 10 … 17 με την από πάνω, οι 8 και 9 με
// καμία), άρα η νέα τιμή κάθε θέσης βγαίνει από τις ΠΑΛΙΕΣ τιμές δύο θέσεων:
//   νέο x[i][sb]      = x[i][sb]·cs      + x[17 − i][sb − 1]·ca      (ο γείτονας: ίδια μνήμη, μία θέση πριν)
//   νέο x[17 − i][sb] = x[17 − i][sb]·cs − x[i][sb + 1]·ca           (ο γείτονας: μία θέση μετά)
// Η υποζώνη 0 δεν έχει από κάτω και η 31 δεν έχει από πάνω: εκεί η «μία θέση πριν/μετά» πέφτει σε άλλη σειρά του ίδιου
// πίνακα (πάντα μέσα στα 576: η πιο πίσω ανάγνωση είναι το x[32·10 − 1], η πιο μπροστά το x[32·7 + 32]), και οι πίνακες
// έχουν σε εκείνη τη θέση (aLow στη θέση 0, aHigh στη θέση 31) ca = 0, δηλαδή «χωρίς γείτονα». Γι' αυτό οι συντελεστές
// είναι γραμμένοι ανά υποζώνη, 32 θέσεις ο καθένας: ο πυρήνας διαβάζει απλώς τη θέση του, χωρίς καμία επιλογή πίνακα.
//
// ΟΙ ΣΥΝΤΕΛΕΣΤΕΣ ΕΧΟΥΝ ΜΕΣΑ ΚΑΙ ΤΟ cos b_i ΤΗΣ ΣΤΡΟΦΗΣ ΠΟΥ ΑΚΟΛΟΥΘΕΙ (βλ. HK_ROT): οι πίνακες δίνουν cs·cos b_i και
// ca·cos b_i (και με πλην για τη σειρά m όταν το i είναι περιττό). Οι δύο πολλαπλασιασμοί που ούτως ή άλλως γίνονται εδώ
// κάνουν έτσι και τη μισή δουλειά της στροφής, δωρεάν.
#define HK_ALIAS(i, m)                                                                                                  \
    const V X##i = vmadd(vld(x + 32 * m - 1), vk(t.aLow.ca[i] + oLo), vmul(vld(x + 32 * i), vk(t.aLow.cs[i] + oLo)));    \
    const V X##m = vnmadd(vld(x + 32 * i + 1), vk(t.aHigh.ca[i] + oHi), vmul(vld(x + 32 * m), vk(t.aHigh.cs[i] + oHi)));

// Η στροφή του ζεύγους (x_k, x_m), m = 17 − k, κατά b_k = π(2k + 1)/72 (Βήμα 2β), με t = tan b_k:
//   P  = x_k·cos + x_m·sin = cos·(x_k + t·x_m)        Q = x_m·cos − x_k·sin = cos·(x_m − t·x_k)
// Το cos το έχουν ήδη οι τιμές που βγάζει το HK_ALIAS (X_k = cos·x_k), άρα μένει μία εντολή ανά έξοδο αντί για δύο.
// Το Q' του Βήματος 2δ είναι (−1)^k·Q:
//   k ζυγό:    X_m = +cos·x_m:   P = X_k + t·X_m      Q' = X_m − t·X_k
//   k περιττό: X_m = −cos·x_m:   P = X_k − t·X_m      Q' = X_m + t·X_k  (= x_k·sin − x_m·cos)
// Το t είναι από 0,044 ως 0,77: η στροφή σε αυτή τη μορφή δεν μεγαλώνει το σφάλμα στρογγύλευσης.
#define HK_ROT_EVEN(k, m)                                                                                               \
    const V P##k = vmadd(X##m, vk(t.rotT[k]), X##k);                                                                    \
    const V Q##k = vnmadd(X##k, vk(t.rotT[k]), X##m);
#define HK_ROT_ODD(k, m)                                                                                                \
    const V P##k = vnmadd(X##m, vk(t.rotT[k]), X##k);                                                                   \
    const V Q##k = vmadd(X##k, vk(t.rotT[k]), X##m);

// Ένα διάνυσμα υποζωνών μακριού μπλοκ. x, ov, out δείχνουν στη σειρά 0 της πρώτης υποζώνης του διανύσματος.
// oLo: η θέση στους συντελεστές της μείωσης αναδίπλωσης για τις σειρές 0 … 7 (ζευγάρι με την ΑΠΟ ΚΑΤΩ υποζώνη)· oHi: για τις
// σειρές 10 … 17 (ζευγάρι με την ΑΠΟ ΠΑΝΩ). Κανονικά και οι δύο είναι η πρώτη υποζώνη του διανύσματος. par: 0 στα διανύσματα·
// στον απλό δρόμο η ισοτιμία της υποζώνης.
HK_FN HK_INLINE void long1(const float *HK_RESTRICT x, float *HK_RESTRICT ov, float *HK_RESTRICT out, int oLo, int oHi, const LongWin &w, int par)
{
    const Tables &t = kT;
    // Η ΣΕΙΡΑ ΤΩΝ ΠΡΑΞΕΩΝ είναι διαλεγμένη ώστε να ζουν μαζί όσο γίνεται λιγότερες τιμές (ο επεξεργαστής έχει 16 καταχωρητές
    // και ό,τι δεν χωρά πάει στη στοίβα): οι σειρές μπαίνουν τέσσερις τέσσερις, k, 17 − k, 8 − k, 9 + k, και γίνονται αμέσως
    // τα τέσσερα s και d του διπλώματος των δύο DCT· μετά τρέχει πρώτα το ένα μισό κάθε DCT, γράφεται ό,τι βγήκε, και μετά το άλλο.
    // Βήμα 2: εννέα στροφές (οι οκτώ με μία εντολή ανά έξοδο· η ένατη, των σειρών 8 και 9, ολόκληρη) και το δίπλωμα του Βήματος 3.
    HK_ALIAS(0, 17)
    HK_ROT_EVEN(0, 17)
    const V X8 = vld(x + 32 * 8), X9 = vld(x + 32 * 9);   // οι δύο μεσαίες σειρές: καμία αναδίπλωση
    const V P8 = vmadd(X8, vk(t.rotC8), vmul(X9, vk(t.rotS8)));
    const V Q8 = vnmadd(X8, vk(t.rotS8), vmul(X9, vk(t.rotC8)));
    const V sP0 = vadd(P0, P8), dP0 = vsub(P0, P8), sQ0 = vadd(Q0, Q8), dQ0 = vsub(Q0, Q8);
    HK_ALIAS(1, 16)
    HK_ROT_ODD(1, 16)
    HK_ALIAS(7, 10)
    HK_ROT_ODD(7, 10)
    const V sP1 = vadd(P1, P7), dP1 = vsub(P1, P7), sQ1 = vadd(Q1, Q7), dQ1 = vsub(Q1, Q7);
    HK_ALIAS(2, 15)
    HK_ROT_EVEN(2, 15)
    HK_ALIAS(6, 11)
    HK_ROT_EVEN(6, 11)
    const V sP2 = vadd(P2, P6), dP2 = vsub(P2, P6), sQ2 = vadd(Q2, Q6), dQ2 = vsub(Q2, Q6);
    HK_ALIAS(3, 14)
    HK_ROT_ODD(3, 14)
    HK_ALIAS(5, 12)
    HK_ROT_ODD(5, 12)
    const V sP3 = vadd(P3, P5), dP3 = vsub(P3, P5), sQ3 = vadd(Q3, Q5), dQ3 = vsub(Q3, Q5);
    HK_ALIAS(4, 13)
    HK_ROT_EVEN(4, 13)

    // Βήματα 3, 2(δ) και 4. C_0 = A_0 · C_2l = A_l + B_(9−l) · C_(2l−1) = A_l − B_(9−l) · C_17 = −B_0 (το πλην είναι στο
    // παράθυρο), όπου A ο DCT-II του P και B του Q'. Κάθε C θέλει ένα A και ένα B με δείκτες αντίθετης ισοτιμίας, άρα:
    //   το περιττό μισό του A με το ζυγό του B δίνουν τα C_1, C_2, C_5, C_6 (επικάλυψη) και C_9, C_10, C_13, C_14, C_17 (έξοδος)
    //   το ζυγό μισό του A με το περιττό του B δίνουν τα C_0, C_3, C_4, C_7, C_8 (επικάλυψη) και C_11, C_12, C_15, C_16 (έξοδος)
    // Η νέα επικάλυψη γράφεται πάνω στην παλιά, άρα ΤΕΛΕΥΤΑΙΑ: πρώτα όλη η έξοδος (διαβάζει την παλιά και των 18 σειρών).
    V A1, A3, A5, A7, B0, B2, B4, B6, B8;
    dct9odd(dP0, dP1, dP2, dP3, A1, A3, A5, A7);
    dct9even(sQ0, sQ1, sQ2, sQ3, Q4, B0, B2, B4, B6, B8);
    longOut(vsub(A5, B4), 9, w, par, ov, out);
    longOut(vadd(A5, B4), 10, w, par, ov, out);
    longOut(vsub(A7, B2), 13, w, par, ov, out);
    longOut(vadd(A7, B2), 14, w, par, ov, out);
    longOut(B0, 17, w, par, ov, out);
    const V C1 = vsub(A1, B8), C2 = vadd(A1, B8), C5 = vsub(A3, B6), C6 = vadd(A3, B6);   // περιμένουν τη σειρά τους

    V A0, A2, A4, A6, A8, B1, B3, B5, B7;
    dct9even(sP0, sP1, sP2, sP3, P4, A0, A2, A4, A6, A8);
    dct9odd(dQ0, dQ1, dQ2, dQ3, B1, B3, B5, B7);
    longOut(vsub(A6, B3), 11, w, par, ov, out);
    longOut(vadd(A6, B3), 12, w, par, ov, out);
    longOut(vsub(A8, B1), 15, w, par, ov, out);
    longOut(vadd(A8, B1), 16, w, par, ov, out);

    longOv(A0, 0, w, par, ov);
    longOv(C1, 1, w, par, ov);
    longOv(C2, 2, w, par, ov);
    longOv(vsub(A2, B7), 3, w, par, ov);
    longOv(vadd(A2, B7), 4, w, par, ov);
    longOv(C5, 5, w, par, ov);
    longOv(C6, 6, w, par, ov);
    longOv(vsub(A4, B5), 7, w, par, ov);
    longOv(vadd(A4, B5), 8, w, par, ov);
}
#undef HK_ALIAS
#undef HK_ROT_EVEN
#undef HK_ROT_ODD

// Ο DCT-IV 6 σημείων ενός κοντού παραθύρου (Βήματα 2 και 3 με N = 6, M = 3): x δείχνει στην πρώτη από τις 6 σειρές του.
// Βγάζει g_j = C_j για j = 0 … 4 και g5 = −C_5 (το πλην είναι στον πίνακα του παραθύρου). 16 πολλαπλασιασμοί, 18 προσθέσεις
// (26 εντολές με FMA).
HK_FN HK_INLINE void short6(const float *HK_RESTRICT x, V &g0, V &g1, V &g2, V &g3, V &g4, V &g5)
{
    const Tables &t = kT;
    const V x0 = vld(x), x1 = vld(x + 32), x2 = vld(x + 64), x3 = vld(x + 96), x4 = vld(x + 128), x5 = vld(x + 160);
    // Στροφές: (x0, x5) κατά π/24, (x1, x4) κατά 3π/24, (x2, x3) κατά 5π/24.
    const V c0 = vk(t.srotC[0]), s0 = vk(t.srotS[0]), c1 = vk(t.srotC[1]), s1 = vk(t.srotS[1]), c2 = vk(t.srotC[2]), s2 = vk(t.srotS[2]);
    const V p0 = vmadd(x0, c0, vmul(x5, s0)), q0 = vnmadd(x0, s0, vmul(x5, c0));
    const V p1 = vmadd(x1, c1, vmul(x4, s1)), q1 = vnmadd(x4, c1, vmul(x1, s1));   // Q'_1 = −Q_1
    const V p2 = vmadd(x2, c2, vmul(x3, s2)), q2 = vnmadd(x2, s2, vmul(x3, c2));
    // DCT-II 3 σημείων: A_0 = p0 + p1 + p2 · A_1 = cos 30°·(p0 − p2) · A_2 = (p0 + p2)/2 − p1.
    const V sp = vadd(p0, p2), sq = vadd(q0, q2);
    const V a0 = vadd(sp, p1), a1 = vmul(vsub(p0, p2), vk(t.c30)), a2 = vmsub(sp, vk(t.half), p1);
    const V b0 = vadd(sq, q1), b1 = vmul(vsub(q0, q2), vk(t.c30)), b2 = vmsub(sq, vk(t.half), q1);
    // C_0 = A_0 · C_1 = A_1 − B_2 · C_2 = A_1 + B_2 · C_3 = A_2 − B_1 · C_4 = A_2 + B_1 · C_5 = −B_0.
    g0 = a0;
    g1 = vsub(a1, b2);
    g2 = vadd(a1, b2);
    g3 = vsub(a2, b1);
    g4 = vadd(a2, b1);
    g5 = b0;
}

// Ένα διάνυσμα υποζωνών κοντού μπλοκ: τρία παράθυρα 12 σημείων στις θέσεις 6, 12, 18 των 36 (2.4.3.4.10.3, block_type 2).
// Με h(w, i) = παράθυρο[i]·x_w[i] (i = 0 … 11) και z[6 + 6w + i] += h(w, i):
//   z[0 … 5] = 0 · z[6 … 11] = h(0, 0 … 5) · z[12 … 17] = h(0, 6 … 11) + h(1, 0 … 5)
//   z[18 … 23] = h(1, 6 … 11) + h(2, 0 … 5) · z[24 … 29] = h(2, 6 … 11) · z[30 … 35] = 0
// Τα z[0 … 17] προστίθενται στην παλιά επικάλυψη, τα z[18 … 35] είναι η νέα. Από το Βήμα 1 με N = 6:
//   x[0 … 5] = C_3, C_4, C_5, −C_5, −C_4, −C_3      x[6 … 11] = −C_2, −C_1, −C_0, −C_0, −C_1, −C_2
// Η θέση i του παραθύρου πέφτει πάντα σε σειρά με την ίδια ισοτιμία με το i (οι μετατοπίσεις 6, 12, 18 είναι ζυγές), άρα ο
// ίδιος πίνακας swin[i] (παράθυρο · πρόσημο του x · πρόσημο αντιστροφής συχνότητας) κάνει και για τα τρία παράθυρα.
//
// HK_SHORT_ROW(i, f, s): οι θέσεις i (πρώτο μισό, από το C_f) και 6 + i (δεύτερο μισό, από το C_s) των τριών παραθύρων
// u (πρώτο), v (δεύτερο), w (τρίτο). Γράφει τις σειρές i, 6 + i, 12 + i της εξόδου και τη σειρά i της νέας επικάλυψης (η
// παλιά σειρά i έχει ήδη διαβαστεί στην πρώτη γραμμή).
#define HK_SHORT_ROW(i, f, s)                                                                                           \
    {                                                                                                                   \
        const V k1 = vk(t.swin[i] + par), k2 = vk(t.swin[6 + i] + par);                                                 \
        vst(out + 32 * i, vld(ov + 32 * i));                                                                            \
        vst(out + 32 * (6 + i), vmadd(u##f, k1, vld(ov + 32 * (6 + i))));                                               \
        vst(out + 32 * (12 + i), vmadd(v##f, k1, vmadd(u##s, k2, vld(ov + 32 * (12 + i)))));                            \
        vst(ov + 32 * i, vmadd(w##f, k1, vmul(v##s, k2)));                                                              \
    }
// Μετά από ΟΛΑ τα παραπάνω (που διάβασαν τις παλιές σειρές 6 … 17): οι σειρές 6 + i και 12 + i της νέας επικάλυψης.
#define HK_SHORT_TAIL(i, s)                                                                                             \
    vst(ov + 32 * (6 + i), vmul(w##s, vk(t.swin[6 + i] + par)));                                                        \
    vst(ov + 32 * (12 + i), zero);

HK_FN HK_INLINE void short1(const float *HK_RESTRICT x, float *HK_RESTRICT ov, float *HK_RESTRICT out, int par)
{
    const Tables &t = kT;
    V u0, u1, u2, u3, u4, u5, v0, v1, v2, v3, v4, v5, w0, w1, w2, w3, w4, w5;
    short6(x, u0, u1, u2, u3, u4, u5);
#if defined(MP3_HYBRID_MUTANT) && MP3_HYBRID_MUTANT == 5
    short6(x + 32 * 12, v0, v1, v2, v3, v4, v5);   // σκόπιμο σφάλμα: το δεύτερο και το τρίτο παράθυρο σε λάθος θέση
    short6(x + 32 * 6, w0, w1, w2, w3, w4, w5);
#else
    short6(x + 32 * 6, v0, v1, v2, v3, v4, v5);
    short6(x + 32 * 12, w0, w1, w2, w3, w4, w5);
#endif
    const V zero = vzero();
    HK_SHORT_ROW(0, 3, 2)
    HK_SHORT_ROW(1, 4, 1)
    HK_SHORT_ROW(2, 5, 0)
    HK_SHORT_ROW(3, 5, 0)
    HK_SHORT_ROW(4, 4, 1)
    HK_SHORT_ROW(5, 3, 2)
    HK_SHORT_TAIL(0, 2)
    HK_SHORT_TAIL(1, 1)
    HK_SHORT_TAIL(2, 0)
    HK_SHORT_TAIL(3, 0)
    HK_SHORT_TAIL(4, 1)
    HK_SHORT_TAIL(5, 2)
}
#undef HK_SHORT_ROW
#undef HK_SHORT_TAIL

// Μακρύ μπλοκ στις υποζώνες [0, n). Στους δρόμους με διανύσματα το n είναι πολλαπλάσιο του HK_LANES.
HK_FN void runLong(const float *x, float *ov, float *out, int n, const LongWin *w)
{
    HK_ENTER();
    for (int s = 0; s < n; s += HK_LANES) {
#if HK_LANES == 1
        const int par = s & 1;
#else
        const int par = 0;
#endif
        long1(x + s, ov + s, out + s, s, s, *w, par);
    }
    HK_LEAVE();
}

// Κοντό μπλοκ στις υποζώνες [from, n).
HK_FN void runShort(const float *x, float *ov, float *out, int from, int n)
{
    HK_ENTER();
    for (int s = from; s < n; s += HK_LANES) {
#if HK_LANES == 1
        const int par = s & 1;
#else
        const int par = 0;
#endif
        short1(x + s, ov + s, out + s, par);
    }
    HK_LEAVE();
}

// Οι υποζώνες [done, 32), που δεν είχαν τίποτα να μετασχηματίσουν: η έξοδός τους είναι σκέτη η παλιά επικάλυψη (που έχει
// ήδη μέσα το πρόσημο) και η νέα τους επικάλυψη είναι μηδέν. Πέρα από το live η παλιά επικάλυψη ΕΙΝΑΙ μηδέν (αναλλοίωτο του
// HybridState), άρα εκεί γράφεται μόνο μηδέν στην έξοδο και η κατάσταση δεν αγγίζεται.
HK_FN void runRest(float *HK_RESTRICT ov, float *HK_RESTRICT out, int done, int live)
{
#if HK_LANES == 1
    // Ο απλός δρόμος: ανά σειρά, ένα κομμάτι αντιγραφή και δύο κομμάτια μηδέν.
    const int copy = live > done ? live - done : 0, from = done + copy;
    for (int r = 0; r < 18; ++r) {
        float *const o = out + 32 * r, *const v = ov + 32 * r;
#if defined(MP3_HYBRID_MUTANT) && MP3_HYBRID_MUTANT == 3
        const int clear = copy > 0 ? copy - 1 : 0;   // σκόπιμο σφάλμα: η επικάλυψη της τελευταίας ζωντανής υποζώνης δεν σβήνεται
#else
        const int clear = copy;
#endif
        if (copy) {
            std::memcpy(o + done, v + done, sizeof(float) * size_t(copy));
            std::memset(v + done, 0, sizeof(float) * size_t(clear));
        }
#if defined(MP3_HYBRID_MUTANT) && MP3_HYBRID_MUTANT == 15
        if (r == 11) continue;   // σκόπιμο σφάλμα: μια σειρά της εξόδου μένει άγραφη πέρα από τις ζωντανές υποζώνες
#endif
        if (from < 32) std::memset(o + from, 0, sizeof(float) * size_t(32 - from));
    }
#else
    // Οι δρόμοι με διανύσματα: ανά στήλη διανύσματος, 18 σειρές. (Ανά σειρά ο GCC το έκανε 18 κλήσεις της memset: μετρημένο,
    // 35 ns ανά granule μόνο γι' αυτές.) Το live στρογγυλεύεται προς τα πάνω σε ολόκληρο διάνυσμα: οι θέσεις [live, liveR)
    // είναι μηδέν, η αντιγραφή τους δεν αλλάζει τίποτα.
    const int liveR = (live + HK_LANES - 1) & ~(HK_LANES - 1);
    const V zero = vzero();
    int s = done;
    for (; s < liveR; s += HK_LANES)
        for (int r = 0; r < 18; ++r) {
            vst(out + 32 * r + s, vld(ov + 32 * r + s));
#if defined(MP3_HYBRID_MUTANT) && MP3_HYBRID_MUTANT == 3
            if (s + HK_LANES < liveR)   // σκόπιμο σφάλμα: η επικάλυψη των τελευταίων ζωντανών υποζωνών δεν σβήνεται
#endif
                vst(ov + 32 * r + s, zero);
        }
    for (; s < 32; s += HK_LANES)
        for (int r = 0; r < 18; ++r) {
#if defined(MP3_HYBRID_MUTANT) && MP3_HYBRID_MUTANT == 15
            if (r == 11) continue;   // σκόπιμο σφάλμα: μια σειρά της εξόδου μένει άγραφη πέρα από τις ζωντανές υποζώνες
#endif
            vst(out + 32 * r + s, zero);
        }
#endif
    HK_LEAVE();
}

#if HK_LANES == 1
// Οι δύο μακριές υποζώνες ενός μικτού μπλοκ (2.4.3.4.10.1: μείωση αναδίπλωσης ΜΟΝΟ ανάμεσα στην 0 και την 1).
// Υποζώνη 0: από κάτω τίποτα (η θέση 0 του aLow), από πάνω κανονικά. Υποζώνη 1: από κάτω κανονικά, από πάνω τίποτα: παίρνει
// τους συντελεστές της θέσης 31 του aHigh, που είναι ακριβώς «χωρίς γείτονα από πάνω».
HK_FN void mixedLong(const float *x, float *ov, float *out, const LongWin *w)
{
    long1(x, ov, out, 0, 0, *w, 0);
#if defined(MP3_HYBRID_MUTANT) && MP3_HYBRID_MUTANT == 11
    long1(x + 1, ov + 1, out + 1, 1, 1, *w, 1);   // σκόπιμο σφάλμα: αναδίπλωση και ανάμεσα στην 1 και τη 2
#else
    long1(x + 1, ov + 1, out + 1, 1, 31, *w, 1);
#endif
    HK_LEAVE();
}
#endif

#undef HK_MUTANT_PATH
