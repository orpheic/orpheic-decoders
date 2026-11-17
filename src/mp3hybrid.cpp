// Το υβριδικό στάδιο του αποκωδικοποιητή MP3, γρήγορα. Τι είναι και τι υπόσχεται: mp3hybrid.h. Προδιαγραφή: exact::hybrid.
// Οι αριθμοί ενοτήτων είναι του ISO/IEC 11172-3. Ό,τι ακολουθεί βγαίνει από τους τύπους του προτύπου με τριγωνομετρία·
// κάθε βήμα ελέγχθηκε και αριθμητικά και ολόκληρο απέναντι στο exact::hybrid.
// Οι ταυτότητες δεν είναι καινούριες, μόνο η παραγωγή τους ξαναγράφτηκε εδώ βήμα βήμα για να ελέγχεται: ο μετασχηματισμός είναι ο
// MDCT των J. P. Princen και A. B. Bradley («Analysis/Synthesis Filter Bank Design Based on Time Domain Aliasing Cancellation», IEEE
// Transactions on Acoustics, Speech, and Signal Processing 34(5), 1986), και η σχέση του με τον DCT-IV είναι κοινή γνώση της βιβλιογραφίας.
//
// ═══ Η ΠΑΡΑΓΩΓΗ ═══════════════════════════════════════════════════════════════════════════════════════════════════════
// Ο τύπος (2.4.3.4.10.2):  x_i = Σ_{k<N} X_k·cos(π/(2n)·(2i + 1 + n/2)(2k + 1)),  i = 0 … n − 1,  n = 36 ή 12.
// Γράφω N = n/2 (18 ή 6) και M = N/2 (9 ή 3). Η γωνία γίνεται θ(i, k) = π(2i + 1 + N)(2k + 1)/(4N).
//
// Βήμα 1. ΑΠΟ ΤΙΣ 2N ΤΙΜΕΣ ΟΙ N ΕΙΝΑΙ ΑΝΕΞΑΡΤΗΤΕΣ.
//   Θέτω i = 3M − 1 − j (j = 0 … N − 1). Τότε 2i + 1 + N = 4N − (2j + 1), άρα θ = π(2k + 1) − π(2j + 1)(2k + 1)/(4N) και,
//   αφού cos(π·περιττός − φ) = −cos φ:
//       x_(3M−1−j) = −C_j,   C_j = Σ_{k<N} X_k·cos(π(2j + 1)(2k + 1)/(4N))                        («DCT-IV» N σημείων)
//   Αυτό δίνει τα x_M … x_(3M−1). Τα υπόλοιπα από δύο συμμετρίες του ίδιου τύπου:
//       2(N − 1 − i) + 1 + N = 4N − (2i + 1 + N)  ⇒  θ(N − 1 − i) = π(2k + 1) − θ(i)   ⇒  x_(N−1−i) = −x_i
//       2(3N − 1 − i) + 1 + N = 8N − (2i + 1 + N) ⇒  θ(3N − 1 − i) = 2π(2k + 1) − θ(i) ⇒  x_(3N−1−i) = x_i
//   Συνολικά:   πρώτο μισό:    x_i     =  C_(M+i)      (i = 0 … M − 1)       x_i     = −C_(3M−1−i)   (i = M … N − 1)
//               δεύτερο μισό:  x_(N+m) = −C_(M−1−m)    (m = 0 … M − 1)       x_(N+m) = −C_(m−M)      (m = M … N − 1)
//   Δηλαδή το πρώτο μισό (που βγαίνει τώρα) θέλει μόνο τα C_M … C_(N−1) και το δεύτερο (η επικάλυψη) μόνο τα C_0 … C_(M−1),
//   και κάθε C_j πηγαίνει σε ΔΥΟ θέσεις.
//
// Βήμα 2. Ο DCT-IV N ΣΗΜΕΙΩΝ = M ΣΤΡΟΦΕΣ + ΔΥΟ DCT-II M ΣΗΜΕΙΩΝ.
//   Γράφω φ(j, k) = π(2j + 1)(2k + 1)/(4N).
//   (α) Ζευγάρωμα k ↔ N − 1 − k: 2(N − 1 − k) + 1 = 2N − (2k + 1), άρα φ(j, N − 1 − k) = π(2j + 1)/2 − φ(j, k) και
//       cos φ(j, N − 1 − k) = (−1)^j·sin φ(j, k). Οπότε  C_j = Σ_{k<M} [X_k·cos φ + (−1)^j·X_(N−1−k)·sin φ].
//   (β) j ζυγό, j = 2l: φ(2l, k) = α + β με α = πl(2k + 1)/(2M) και β = b_k = π(2k + 1)/(4N) (δεν εξαρτάται από το l).
//       X·cos(α + β) + X'·sin(α + β) = cos α·(X·cos β + X'·sin β) + sin α·(X'·cos β − X·sin β). Με τη ΣΤΡΟΦΗ
//           P_k = X_k·cos b_k + X_(N−1−k)·sin b_k        Q_k = X_(N−1−k)·cos b_k − X_k·sin b_k
//       γίνεται  C_2l = Σ_{k<M} [P_k·cos α + Q_k·sin α].
//   (γ) j περιττό, γραμμένο j = N − 1 − 2l: 2j + 1 = 2N − (4l + 1), άρα φ(j, k) = π(2k + 1)/2 − φ(2l, k) και
//       cos φ(j, k) = (−1)^k·sin φ(2l, k):  C_(N−1−2l) = Σ_{k<N} (−1)^k·X_k·sin φ(2l, k). Ζευγάρωμα όπως στο (α):
//       sin φ(2l, N − 1 − k) = sin(π(4l + 1)/2 − φ(2l, k)) = cos φ(2l, k) και (−1)^(N−1−k) = −(−1)^k (N ζυγό), άρα
//           C_(N−1−2l) = Σ_{k<M} (−1)^k·[X_k·sin(α + β) − X_(N−1−k)·cos(α + β)] = Σ_{k<M} (−1)^k·[P_k·sin α − Q_k·cos α]
//       (γιατί X·sin(α + β) − X'·cos(α + β) = sin α·(X·cos β + X'·sin β) − cos α·(X'·cos β − X·sin β)).
//   (δ) Τα ημίτονα γίνονται συνημίτονα με l → M − l: πl'(2k + 1)/(2M) με l' = M − l είναι π(2k + 1)/2 − α, άρα
//       cos α(M − l) = (−1)^k·sin α(l) και sin α(M − l) = (−1)^k·cos α(l). Ορίζω Q'_k = (−1)^k·Q_k και τους δύο DCT-II
//           A_l = Σ_{k<M} P_k·cos(πl(2k + 1)/(2M))        B_l = Σ_{k<M} Q'_k·cos(πl(2k + 1)/(2M))        (l = 0 … M − 1)
//       Τότε Σ Q_k·sin α(l) = B_(M−l), Σ (−1)^k·P_k·sin α(l) = A_(M−l), Σ (−1)^k·Q_k·cos α(l) = B_l (με A_M = B_M = 0, αφού
//       cos(π(2k + 1)/2) = 0), και:
//           C_0 = A_0        C_2l = A_l + B_(M−l)        C_(2l−1) = A_l − B_(M−l)   (l = 1 … M − 1)        C_(N−1) = −B_0
//
// Βήμα 3. Ο DCT-II 9 ΣΗΜΕΙΩΝ:  A_l = Σ_{k=0}^{8} p_k·cos(l(2k + 1)·10°). Γράφω cNN = cos NN°.
//   Δίπλωμα k ↔ 8 − k: cos(l(18 − (2k + 1))·10°) = (−1)^l·cos(l(2k + 1)·10°). Με s_k = p_k + p_(8−k), d_k = p_k − p_(8−k)
//   (k = 0 … 3) και s_4 = p_4, τα ζυγά l βλέπουν μόνο τα s και τα περιττά μόνο τα d:
//       A_0 = s0 + s1 + s2 + s3 + s4
//       A_2 = s0·c20 + s1·c60 + s2·c100 + s3·c140 + s4·c180 =  s0·c20 + s1/2 − s2·c80 − s3·c40 − s4
//       A_4 = s0·c40 + s1·c120 + s2·c200 + s3·c280 + s4·c360 =  s0·c40 − s1/2 − s2·c20 + s3·c80 + s4
//       A_6 = s0·c60 + s1·c180 + s2·c300 + s3·c420 + s4·c540 = (s0 + s2 + s3)/2 − s1 − s4
//       A_8 = s0·c80 + s1·c240 + s2·c400 + s3·c560 + s4·c720 =  s0·c80 − s1/2 + s2·c40 − s3·c20 + s4
//       A_1 = d0·c10 + d1·c30 + d2·c50 + d3·c70
//       A_3 = d0·c30 + d1·c90 + d2·c150 + d3·c210            =  c30·(d0 − d2 − d3)
//       A_5 = d0·c50 + d1·c150 + d2·c250 + d3·c350           =  d0·c50 − d1·c30 − d2·c70 + d3·c10
//       A_7 = d0·c70 + d1·c210 + d2·c350 + d3·c490           =  d0·c70 − d1·c30 + d2·c10 − d3·c50
//   Δύο ταυτότητες, και οι δύο από 2·cos 60°·cos θ = cos θ:   c20 = c40 + c80   και   c10 = c50 + c70.
//   Με u = s0 − s2, v = s3 − s2, e = s1/2 − s4 (αντικατάσταση του c20, έλεγχος όρο προς όρο):
//       A_2 =  e + c20·u − c40·v        A_4 = −e + c40·u + c80·v        A_8 = −e + c80·u − c20·v
//   Με u' = d0 + d3, v' = d3 − d2, g = c30·d1 (αντικατάσταση του c10):
//       A_1 =  g + c10·u' − c50·v'      A_5 = −g + c50·u' + c70·v'      A_7 = −g + c70·u' − c10·v'
//   Ο DCT-II 3 ΣΗΜΕΙΩΝ: A_0 = p0 + p1 + p2 · A_1 = c30·(p0 − p2) · A_2 = p0·c60 + p1·c180 + p2·c300 = (p0 + p2)/2 − p1.
//
// Βήμα 4. ΠΑΡΑΘΥΡΟ, ΕΠΙΚΑΛΥΨΗ, ΠΡΟΣΗΜΟ ΣΕ ΜΙΑ ΠΡΑΞΗ ΑΝΑ ΘΕΣΗ. Με z_i = παράθυρο[i]·x_i, η έξοδος είναι z_i + παλιά
//   επικάλυψη (i < N) και η νέα επικάλυψη z_(N+i). Από το Βήμα 1 κάθε θέση είναι ΕΝΑ C_j επί μία σταθερά (± το παράθυρο):
//       έξοδος   σειρά j − 9:   +w[j − 9]·C_j        σειρά 26 − j:  −w[26 − j]·C_j        (j = 9 … 17)
//       επικάλ.  σειρά 8 − j:   −w[26 − j]·C_j       σειρά 9 + j:   −w[27 + j]·C_j        (j = 0 … 8)
//   Το παράθυρο δεν κοστίζει λοιπόν τίποτα παραπάνω από τους 36 πολλαπλασιασμούς του, και οι 18 της εξόδου είναι η ίδια
//   εντολή με την πρόσθεση της επικάλυψης (a·b + c).
//
// ΤΟ ΠΡΟΣΗΜΟ (2.4.3.4.10.5: κάθε δεύτερο δείγμα κάθε δεύτερης υποζώνης αλλάζει πρόσημο). Η αλλαγή εφαρμόζεται στο ΑΘΡΟΙΣΜΑ
//   z + επικάλυψη της σειράς i, και η επικάλυψη της σειράς i προστίθεται πάντα στη σειρά i. Άρα αν η επικάλυψη φυλάγεται ήδη
//   με το πρόσημο της θέσης της, αρκεί κάθε σταθερά του παραθύρου να έχει μέσα το πρόσημο της σειράς ΟΠΟΥ ΠΕΦΤΕΙ:
//   −(z + ov) = (−w)·x + (−ov). Οι σταθερές είναι γι' αυτό μικρά διανύσματα 8 θέσεων (+, ±, +, ±, …): τα διανύσματα
//   αρχίζουν πάντα σε ζυγή υποζώνη, άρα περιττή θέση = περιττή υποζώνη. Καμία πράξη για το πρόσημο.
//
// ΤΟ `subbands`. Με τιμές μόνο στις πρώτες s υποζώνες: στα μακριά μπλοκ η μείωση αναδίπλωσης δίνει τιμές και στην υποζώνη s
//   (το ζεύγος s − 1, s), άρα μετασχηματίζονται οι [0, min(s + 1, 32))· στα κοντά οι [0, s)· στα μικτά οι [0, max(s, 2)) (η
//   μόνη αναδίπλωση είναι ανάμεσα στην 0 και την 1). Με s = 0 καμία. Οι υπόλοιπες υποζώνες βγάζουν σκέτη την παλιά τους
//   επικάλυψη και η νέα τους είναι μηδέν· το HybridState::live θυμάται ως πού υπάρχει επικάλυψη, ώστε το επόμενο granule να
//   μην αγγίξει καν ό,τι είναι μηδέν. Οι δρόμοι με διανύσματα στρογγυλεύουν προς τα πάνω σε ολόκληρο διάνυσμα: οι παραπάνω
//   θέσεις έχουν είσοδο μηδέν (το εγγυάται ο καλών) και ο μετασχηματισμός του μηδενός είναι μηδέν, ακριβώς.
//
// ΜΕΤΡΗΜΑ ΠΡΑΞΕΩΝ ανά υποζώνη (ή ανά διάνυσμα υποζωνών), μακρύ μπλοκ. «π» = πολλαπλασιασμός, «+» = πρόσθεση/αφαίρεση.
//       μείωση αναδίπλωσης   32π 16+      (16 τιμές, καθεμία 2π 1+· οι συντελεστές έχουν μέσα και το cos της στροφής)
//       στροφές              20π 18+      (οκτώ με 2π 2+ χάρη στο παραπάνω, βλ. HK_ROT στον πυρήνα· η ένατη, των σειρών 8 και 9, 4π 2+)
//       δύο DCT-II           32π 64+      (καθένας 16π 32+· από τους 16π οι 2 είναι επί ½)
//       συνδυασμός C_j        0π 16+
//       παράθυρο + επικάλυψη 36π 18+
//       ΣΥΝΟΛΟ              120π 132+     · ο τύπος κατά γράμμα: 32 + 18·36 + 36 = 716π και 16 + 17·36 + 18 = 646+
//   Με FMA οι 80 από τις 132 προσθέσεις είναι η ίδια εντολή με έναν πολλαπλασιασμό: 172 εντολές πράξεων ανά διάνυσμα, και
//   μαζί τους 34 αναγνώσεις του x, 18 της επικάλυψης και 36 γραφές.
//   Κοντό μπλοκ: τρεις DCT-IV 6 σημείων (καθένας 12π 6+ οι στροφές, 4π 8+ οι δύο DCT-II 3 σημείων, 4+ ο συνδυασμός) και
//   36π 24+ για παράθυρο, άθροισμα γειτονικών παραθύρων και επικάλυψη: 84π 78+ (114 εντολές με FMA) αντί 252π 234+.
//
// ΑΚΡΙΒΕΙΑ (η δοκιμή του σταδίου, 135 734 κλήσεις ανά δρόμο): μέγιστη διαφορά από το exact::hybrid 3,1e-7 του μεγίστου του μπλοκ
//   στον απλό δρόμο και στον SSE2, 2,9e-7 με FMA· το όριο είναι 1e-6.
//
// ΤΙ ΑΛΛΟ ΔΟΚΙΜΑΣΤΗΚΕ ΚΑΙ ΜΕΤΡΗΘΗΚΕ (οι παραλλαγές είναι στο alt/, με το alt/compare.sh· ns ανά μακρύ μπλοκ 32 υποζωνών σε
//   i9-14900K, g++ 13 -O2, δρόμος AVX2 / SSE2 / απλός, το μικρότερο από 4 ως 6 τρεξίματα):
//   · Ο πιο φτηνός DCT-IV (προσθέσεις X_k + X_(k−1), DCT-III, και κλιμάκωση 1/(2cos(π(2j + 1)/(4N))) διπλωμένη στο παράθυρο:
//     9 εντολές λιγότερες με FMA, 18 χωρίς). Η κλιμάκωση φτάνει το 11,5 στο j = 17 και μεγαλώνει το σφάλμα στρογγύλευσης όσο
//     και το σήμα: σε float (200 000 τυχαία μπλοκ) μέγιστη διαφορά 4,2e-6 του μεγίστου, έξω από το όριο·
//     ο δρόμος εδώ (μόνο στροφές και αθροίσματα) 3,0e-7 στην ίδια μέτρηση. Ταχύτητα: δεν μετρήθηκε.
//   · Στροφές με cos και sin χωριστά, χωρίς το δίπλωμα στη μείωση αναδίπλωσης (136π): 66,9 / 143 / 584 αντί 63,5 / 136 / 563.
//   · DCT-II με 12π αντί 16π (m = c80·(u + v) − e, c4 = m + (c40 − c80)·u, c8 = m − (c20 + c80)·v, c2 = c4 + c8 + 3e): 63,1 /
//     134 / 551 αντί 63,5 / 136 / 563, αλλά η μέγιστη διαφορά ανεβαίνει από 2,9e-7 σε 4,1e-7 (FMA) και από 3,1e-7 σε 3,9e-7.
//     Κέρδος κάτω από 1% στον κύριο δρόμο με χειρότερη ακρίβεια: δεν κρατήθηκε.
//   · Η μείωση αναδίπλωσης ως χωριστό πέρασμα ΕΠΙ ΤΟΠΟΥ στο x (3π ανά ζεύγος) και μετά ο μετασχηματισμός, με το x να
//     αντιγράφεται φρέσκο πριν από κάθε κλήση και στους δύο: 85,7 / 160 / 576 αντί 75,5 / 149 / 578. Δεν κρατήθηκε.
//   · Η σειρά των πράξεων «πρώτα το ένα μισό κάθε DCT, μετά το άλλο» (όπως είναι τώρα) αντί «ολόκληρος ο ένας, μετά ο
//     άλλος»: 59,5 / 134 / 546 αντί 62,5 / 136 / 560 με g++, 66 αντί 71,5 με clang++ (λιγότερες τιμές στη στοίβα).
// ═════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
#include "mp3hybrid.h"
#include "mp3consts.h"

#include <cmath>
#include <cstring>

// Ανά μεταγλωττιστή αλλάζει μόνο το πώς γράφεται «βάλ' το μέσα», «μην υποθέσεις επικάλυψη δεικτών» και η ερώτηση «τι
// ξέρει αυτός ο επεξεργαστής» (ίδιες συμβάσεις με το flacdecoder.cpp).
#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#define HK_INLINE __forceinline
#define HK_TARGET_AVX2
#else
#define HK_INLINE inline __attribute__((always_inline))
#define HK_TARGET_AVX2 __attribute__((target("avx2,fma")))
#endif
#define HK_RESTRICT __restrict

#if defined(__x86_64__) || defined(_M_X64)
#define HK_X86 1
#if !defined(_MSC_VER) || defined(__clang__)
#include <cpuid.h>
#endif
#include <immintrin.h>
#else
#define HK_X86 0
#endif
#if defined(__aarch64__) && !defined(_MSC_VER)
#define HK_NEON 1
#include <arm_neon.h>
#else
#define HK_NEON 0
#endif

namespace orpheic::mp3 {
namespace {

// ── ΟΙ ΠΙΝΑΚΕΣ ──
// Κάθε σταθερά είναι 8 float (ένα διάνυσμα AVX· οι δρόμοι των 4 διαβάζουν το πρώτο μισό, ο απλός δρόμος μία θέση). Οι πιο
// πολλές έχουν την ίδια τιμή και στις 8 θέσεις· οι σταθερές του παραθύρου έχουν το πρόσημο της αντιστροφής συχνότητας στις
// περιττές θέσεις. Οι συντελεστές της μείωσης αναδίπλωσης είναι γραμμένοι ανά υποζώνη (32 θέσεις), με «χωρίς γείτονα» στη
// θέση της υποζώνης που δεν έχει. (Με χωριστούς πίνακες «πρώτο / μεσαίο / τελευταίο διάνυσμα» και επιλογή ανά διάνυσμα ο
// clang έβγαζε μία επιλογή διεύθυνσης ανά σταθερά: μετρημένο, 97 ns ανά μακρύ μπλοκ αντί 80, και 810 αντί 540 στον απλό δρόμο.)
struct AliasHalf {
    float cs[8][32], ca[8][32];   // ανά υποζώνη
};
struct LongWin {
    float a[9][8], b[9][8];   // έξοδος: C_j (j = 9 … 17) → σειρά j − 9 και σειρά 26 − j
    float c[9][8], d[9][8];   // νέα επικάλυψη: C_j (j = 0 … 8) → σειρά 8 − j και σειρά 9 + j
};
struct Tables {
    alignas(32) float rotT[8][8];        // tan b_k, b_k = π(2k + 1)/72, k = 0 … 7: οι στροφές του μακριού μπλοκ (Βήμα 2)
    float rotC8[8], rotS8[8];            // cos b_8, sin b_8: η στροφή των δύο μεσαίων σειρών, που δεν περνούν από αναδίπλωση
    float srotC[3][8], srotS[3][8];      // cos, sin του b_k = π(2k + 1)/24: οι στροφές του κοντού
    float c10[8], c20[8], c30[8], c40[8], c50[8], c70[8], c80[8], half[8];   // Βήμα 3
    // Πίνακας B.9 → cs, ca, πολλαπλασιασμένα με το cos b_i της στροφής (βλ. HK_ALIAS στον πυρήνα). aLow: οι σειρές 0 … 7,
    // που ζευγαρώνουν με την από κάτω υποζώνη (η 0 δεν έχει: ca = 0 στη θέση 0)· aHigh: οι σειρές 17 … 10, με την από πάνω
    // (η 31 δεν έχει: ca = 0 στη θέση 31), και με πλην για i περιττό.
    AliasHalf aLow, aHigh;
    LongWin win[4];                      // ανά block_type· το [2] δεν χρησιμοποιείται (τα μακριά του μικτού παίρνουν το [0])
    float swin[12][8];                   // το παράθυρο 12 σημείων του κοντού, με τα πρόσημα μέσα
    Tables();
};

void fill(float *lanes, double v)
{
    for (int i = 0; i < 8; ++i) lanes[i] = float(v);
}
// Το πρόσημο της αντιστροφής συχνότητας στη σειρά row, θέση lane του διανύσματος (περιττή θέση = περιττή υποζώνη).
double flip(int row, int lane) { return (row & 1) && (lane & 1) ? -1.0 : 1.0; }

Tables::Tables()
{
    for (int k = 0; k < 8; ++k) fill(rotT[k], consts::tanPi576(8 * (2 * k + 1)));   // tan(π(2k+1)/72)
    fill(rotC8, consts::cosPi576(8 * 17));
    fill(rotS8, consts::sinPi576(8 * 17));
    for (int k = 0; k < 3; ++k) {
        fill(srotC[k], consts::cosPi576(24 * (2 * k + 1)));
        fill(srotS[k], consts::sinPi576(24 * (2 * k + 1)));
    }
    fill(c10, consts::cosPi576(32));   // cos(10°)
    fill(c20, consts::cosPi576(64));   // cos(20°)
    fill(c30, consts::cosPi576(96));   // cos(30°)
    fill(c40, consts::cosPi576(128));   // cos(40°)
    fill(c50, consts::cosPi576(160));   // cos(50°)
    fill(c70, consts::cosPi576(224));   // cos(70°)
    fill(c80, consts::cosPi576(256));   // cos(80°)
    fill(half, 0.5);

    // Πίνακας B.9: cs = 1/√(1 + c²), ca = c/√(1 + c²), εδώ επί cos b_i (και επί −1 στις σειρές 17 − i με i περιττό).
    const double c[8] = {-0.6, -0.535, -0.33, -0.185, -0.095, -0.041, -0.0142, -0.0037};
    for (int i = 0; i < 8; ++i) {
        const double cosb = consts::cosPi576(8 * (2 * i + 1)), sign = (i & 1) ? -1.0 : 1.0;   // cos(π(2i+1)/72)
        const double cs = 1.0 / std::sqrt(1.0 + c[i] * c[i]);
        double ca = c[i] / std::sqrt(1.0 + c[i] * c[i]);
#if defined(MP3_HYBRID_MUTANT) && MP3_HYBRID_MUTANT == 9
        if (i == 7) ca *= 1.01;   // σκόπιμο σφάλμα: ο μικρότερος συντελεστής αναδίπλωσης λάθος κατά 1%
#endif
        for (int sb = 0; sb < 32; ++sb) {
            aLow.cs[i][sb] = float(cs * cosb);
            aLow.ca[i][sb] = float(ca * cosb);
            aHigh.cs[i][sb] = float(sign * cs * cosb);
            aHigh.ca[i][sb] = float(sign * ca * cosb);
        }
        // Χωρίς γείτονα η τιμή μένει όπως είναι (cs = 1, ca = 0), πάλι επί το cos b_i της στροφής.
#if !(defined(MP3_HYBRID_MUTANT) && MP3_HYBRID_MUTANT == 12)   // σκόπιμο σφάλμα 12: η υποζώνη 0 «ζευγαρώνει» με ξένη σειρά
        aLow.cs[i][0] = float(cosb);
        aLow.ca[i][0] = 0.0f;
#endif
#if !(defined(MP3_HYBRID_MUTANT) && MP3_HYBRID_MUTANT == 13)   // σκόπιμο σφάλμα 13: το ίδιο για την υποζώνη 31
        aHigh.cs[i][31] = float(sign * cosb);
        aHigh.ca[i][31] = 0.0f;
#endif
    }

    // Τα παράθυρα 36 σημείων (2.4.3.4.10.3) και η τοποθέτησή τους κατά το Βήμα 4.
    for (int bt = 0; bt < 4; ++bt) {
        double w[36];
        for (int i = 0; i < 36; ++i) {
            const double longSine = consts::sinPi576(8 * (2 * i + 1));   // sin(π/36·(i+½))
            if (bt == 0)
                w[i] = longSine;
            else if (bt == 1)   // start: μακρύ μισό, μετά 1, μετά το δεύτερο μισό του κοντού, μετά 0
                w[i] = i < 18 ? longSine : i < 24 ? 1.0 : i < 30 ? consts::sinPi576(24 * (2 * (i - 18) + 1)) : 0.0;
            else if (bt == 3)   // stop: 0, το πρώτο μισό του κοντού, 1, μακρύ μισό
                w[i] = i < 6 ? 0.0 : i < 12 ? consts::sinPi576(24 * (2 * (i - 6) + 1)) : i < 18 ? 1.0 : longSine;
            else
                w[i] = 0.0;
        }
        LongWin &lw = win[bt];
        for (int lane = 0; lane < 8; ++lane) {
            for (int j = 9; j < 18; ++j) {
                const double g = j == 17 ? -1.0 : 1.0;   // ο πυρήνας δίνει B_0 = −C_17
                lw.a[j - 9][lane] = float(flip(j - 9, lane) * w[j - 9] * g);
                lw.b[j - 9][lane] = float(-flip(26 - j, lane) * w[26 - j] * g);
            }
            for (int j = 0; j < 9; ++j) {
                lw.c[j][lane] = float(-flip(8 - j, lane) * w[26 - j]);
                lw.d[j][lane] = float(-flip(9 + j, lane) * w[27 + j]);
            }
        }
    }
    // Το παράθυρο 12 σημείων. Η θέση i παίρνει (Βήμα 1, N = 6):  i = 0 … 5: +C_3, +C_4, +C_5, −C_5, −C_4, −C_3 ·
    // i = 6 … 11: −C_2, −C_1, −C_0, −C_0, −C_1, −C_2. Ο πυρήνας δίνει −C_5 αντί C_5: το πρόσημο αυτό μπαίνει εδώ.
    for (int i = 0; i < 12; ++i) {
        const double w = consts::sinPi576(24 * (2 * i + 1));   // sin(π/12·(i+½))
        const double sign = (i < 3) ? 1.0 : -1.0;
        const double g = (i == 2 || i == 3) ? -1.0 : 1.0;
        for (int lane = 0; lane < 8; ++lane) swin[i][lane] = float(flip(i, lane) * sign * g * w);
    }

#if defined(MP3_HYBRID_MUTANT)
    // Σκόπιμα σφάλματα σε ΜΙΑ τιμή παραθύρου (μία στις δέκα χιλιάδες): η δοκιμή πρέπει να τα πιάνει.
    const float off = 1.0001f;
    for (int lane = 0; lane < 8; ++lane) {
        if (MP3_HYBRID_MUTANT == 1) win[0].a[4][lane] *= off;    // κανονικό παράθυρο, σειρά 4
        if (MP3_HYBRID_MUTANT == 7) win[1].c[3][lane] *= off;    // παράθυρο start, επικάλυψη σειρά 5
        if (MP3_HYBRID_MUTANT == 8) win[3].a[7][lane] *= off;    // παράθυρο stop, σειρά 7
        if (MP3_HYBRID_MUTANT == 14) swin[9][lane] *= off;       // κοντό παράθυρο, θέση 9
        if (MP3_HYBRID_MUTANT == 19) win[1].d[2][lane] = win[0].d[2][lane];   // start με μία τιμή του κανονικού παραθύρου
    }
#endif
}

const Tables kT;

// ── ΟΙ ΔΡΟΜΟΙ ──
// Ο ίδιος πυρήνας (mp3hybrid_kernel.inc) τέσσερις φορές, καθεμία πάνω σε άλλον τύπο διανύσματος.

// (α) Απλή C++, μία υποζώνη τη φορά: ο δρόμος κάθε άλλου επεξεργαστή, και των μικτών μπλοκ σε όλους.
namespace scalar {
#define HK_LANES 1
#define HK_FN
#define HK_ENTER()
#define HK_LEAVE()
using V = float;
HK_INLINE V vld(const float *p) { return *p; }
HK_INLINE void vst(float *p, V a) { *p = a; }
HK_INLINE V vk(const float *p) { return *p; }
HK_INLINE V vadd(V a, V b) { return a + b; }
HK_INLINE V vsub(V a, V b) { return a - b; }
HK_INLINE V vmul(V a, V b) { return a * b; }
HK_INLINE V vmadd(V a, V b, V c) { return a * b + c; }
HK_INLINE V vmsub(V a, V b, V c) { return a * b - c; }
HK_INLINE V vnmadd(V a, V b, V c) { return c - a * b; }
HK_INLINE V vzero() { return 0.0f; }
#include "mp3hybrid_kernel.inc"
#undef HK_LANES
#undef HK_FN
#undef HK_ENTER
#undef HK_LEAVE
}   // namespace scalar

#if HK_X86
// (β) SSE2, που το έχει κάθε x86-64: 4 υποζώνες τη φορά. Οι σταθερές διαβάζονται στοιχισμένες (οι πίνακες είναι στα 32
// byte και οι μετατοπίσεις πολλαπλάσια των 16), ώστε να μπαίνουν κατευθείαν ως τελεστέος μνήμης.
namespace sse2 {
#define HK_LANES 4
#define HK_FN
#define HK_ENTER()
#define HK_LEAVE()
using V = __m128;
HK_INLINE V vld(const float *p) { return _mm_loadu_ps(p); }
HK_INLINE void vst(float *p, V a) { _mm_storeu_ps(p, a); }
HK_INLINE V vk(const float *p) { return _mm_load_ps(p); }
HK_INLINE V vadd(V a, V b) { return _mm_add_ps(a, b); }
HK_INLINE V vsub(V a, V b) { return _mm_sub_ps(a, b); }
HK_INLINE V vmul(V a, V b) { return _mm_mul_ps(a, b); }
HK_INLINE V vmadd(V a, V b, V c) { return _mm_add_ps(_mm_mul_ps(a, b), c); }
HK_INLINE V vmsub(V a, V b, V c) { return _mm_sub_ps(_mm_mul_ps(a, b), c); }
HK_INLINE V vnmadd(V a, V b, V c) { return _mm_sub_ps(c, _mm_mul_ps(a, b)); }
HK_INLINE V vzero() { return _mm_setzero_ps(); }
#include "mp3hybrid_kernel.inc"
#undef HK_LANES
#undef HK_FN
#undef HK_ENTER
#undef HK_LEAVE
}   // namespace sse2

// (γ) AVX2 + FMA, αν τα έχει ο επεξεργαστής ΚΑΙ τα σώζει το λειτουργικό: 8 υποζώνες τη φορά, a·b ± c σε μία εντολή.
// Κάθε συνάρτηση του δρόμου γράφει target(...)· ο MSVC δεν το χρειάζεται (δέχεται τις εντολές παντού).
namespace avx2 {
#define HK_LANES 8
#define HK_FN HK_TARGET_AVX2
// Ο πυρήνας κρατά πάνω από 16 τιμές ζωντανές, άρα μερικές περνούν από τη στοίβα. Ο clang δεν στοιχίζει τη στοίβα στα 32 byte
// μόνο γι' αυτές, και οι μισές θέσεις πέφτουν τότε πάνω σε σύνορο γραμμής της κρυφής μνήμης (μετρημένο: 67 ns ανά μακρύ μπλοκ
// αντί 63,5· κοντό 48 αντί 41). Μια τοπική μεταβλητή στοιχισμένη στα 32 τον υποχρεώνει· ο GCC το κάνει ήδη μόνος του.
#define HK_ENTER()                                                                                                      \
    alignas(32) volatile float hkAligned[8];                                                                            \
    hkAligned[0] = 0.0f;                                                                                                \
    (void)hkAligned
#define HK_LEAVE() _mm256_zeroupper()
using V = __m256;
HK_FN HK_INLINE V vld(const float *p) { return _mm256_loadu_ps(p); }
HK_FN HK_INLINE void vst(float *p, V a) { _mm256_storeu_ps(p, a); }
HK_FN HK_INLINE V vk(const float *p) { return _mm256_load_ps(p); }
HK_FN HK_INLINE V vadd(V a, V b) { return _mm256_add_ps(a, b); }
HK_FN HK_INLINE V vsub(V a, V b) { return _mm256_sub_ps(a, b); }
HK_FN HK_INLINE V vmul(V a, V b) { return _mm256_mul_ps(a, b); }
HK_FN HK_INLINE V vmadd(V a, V b, V c) { return _mm256_fmadd_ps(a, b, c); }
HK_FN HK_INLINE V vmsub(V a, V b, V c) { return _mm256_fmsub_ps(a, b, c); }
HK_FN HK_INLINE V vnmadd(V a, V b, V c) { return _mm256_fnmadd_ps(a, b, c); }
HK_FN HK_INLINE V vzero() { return _mm256_setzero_ps(); }
#include "mp3hybrid_kernel.inc"
#undef HK_LANES
#undef HK_FN
#undef HK_ENTER
#undef HK_LEAVE
}   // namespace avx2
#endif

#if HK_NEON
// (δ) NEON του aarch64 (το έχει κάθε τέτοιος επεξεργαστής): 4 υποζώνες τη φορά, καθρέφτης του δρόμου SSE2 με a·b ± c σε
// μία εντολή. Γράφτηκε χωρίς τέτοιο μηχάνημα· έτρεξε στις 14/11/2026 σε Raspberry Pi 5 (Cortex-A76, g++ 14.2) μέσα σε ολόκληρο τον
// αποκωδικοποιητή: 32 αρχεία ίδια σε κάθε πακέτο με τον αποκωδικοποιητή του FFmpeg. Ό,τι διαφέρει από τους άλλους δρόμους είναι
// αυτές οι δέκα γραμμές· ο πυρήνας είναι ο ίδιος. Η δοκιμή του ίδιου του σταδίου δεν έχει τρέξει ακόμη σε τέτοιο μηχάνημα.
namespace neon {
#define HK_LANES 4
#define HK_FN
#define HK_ENTER()
#define HK_LEAVE()
using V = float32x4_t;
HK_INLINE V vld(const float *p) { return vld1q_f32(p); }
HK_INLINE void vst(float *p, V a) { vst1q_f32(p, a); }
HK_INLINE V vk(const float *p) { return vld1q_f32(p); }
HK_INLINE V vadd(V a, V b) { return vaddq_f32(a, b); }
HK_INLINE V vsub(V a, V b) { return vsubq_f32(a, b); }
HK_INLINE V vmul(V a, V b) { return vmulq_f32(a, b); }
HK_INLINE V vmadd(V a, V b, V c) { return vfmaq_f32(c, a, b); }               // c + a·b
HK_INLINE V vmsub(V a, V b, V c) { return vnegq_f32(vfmsq_f32(c, a, b)); }    // −(c − a·b)
HK_INLINE V vnmadd(V a, V b, V c) { return vfmsq_f32(c, a, b); }              // c − a·b
HK_INLINE V vzero() { return vdupq_n_f32(0.0f); }
#include "mp3hybrid_kernel.inc"
#undef HK_LANES
#undef HK_FN
#undef HK_ENTER
#undef HK_LEAVE
}   // namespace neon
#endif

// ── ΠΟΙΟΣ ΔΡΟΜΟΣ ──
struct Kernel {
    int lanes;
    void (*runLong)(const float *, float *, float *, int, const LongWin *);
    void (*runShort)(const float *, float *, float *, int, int);
    void (*runRest)(float *, float *, int, int);
    const char *text;
};
const Kernel kScalar = {1, scalar::runLong, scalar::runShort, scalar::runRest, "hybrid: plain c++, one subband at a time"};
#if HK_X86
const Kernel kSse2 = {4, sse2::runLong, sse2::runShort, sse2::runRest, "hybrid: sse2, 4 subbands at a time"};
const Kernel kAvx2 = {8, avx2::runLong, avx2::runShort, avx2::runRest, "hybrid: avx2 + fma, 8 subbands at a time"};

// AVX2 και FMA στον επεξεργαστή, και το λειτουργικό να σώζει τους καταχωρητές των 256 μπιτ (αλλιώς οι εντολές πέφτουν).
bool cpuHasAvx2Fma()
{
#if defined(_MSC_VER) && !defined(__clang__)
    int r[4];
    __cpuid(r, 0);
    if (r[0] < 7) return false;
    __cpuid(r, 1);
    const int need = (1 << 12) | (1 << 27) | (1 << 28);   // FMA, OSXSAVE, AVX
    if ((r[2] & need) != need) return false;
    if ((_xgetbv(0) & 6) != 6) return false;
    __cpuidex(r, 7, 0);
    return (r[1] & (1 << 5)) != 0;
#else
    unsigned a = 0, b = 0, c = 0, d = 0;
    const unsigned need = (1u << 12) | (1u << 27) | (1u << 28);   // FMA, OSXSAVE, AVX
    if (!__get_cpuid(1, &a, &b, &c, &d) || (c & need) != need) return false;
    unsigned lo = 0, hi = 0;
    __asm__ volatile("xgetbv" : "=a"(lo), "=d"(hi) : "c"(0));   // επιτρέπεται, αφού το OSXSAVE είναι αναμμένο
    if ((lo & 6) != 6) return false;
    return __get_cpuid_count(7, 0, &a, &b, &c, &d) && (b & (1u << 5));
#endif
}
#elif HK_NEON
const Kernel kNeon = {4, neon::runLong, neon::runShort, neon::runRest, "hybrid: neon, 4 subbands at a time"};
#endif

const Kernel *pick(bool basicOnly, bool scalarOnly)
{
    if (scalarOnly) return &kScalar;
#if HK_X86
    return (!basicOnly && cpuHasAvx2Fma()) ? &kAvx2 : &kSse2;
#elif HK_NEON
    (void)basicOnly;
    return &kNeon;
#else
    (void)basicOnly;
    return &kScalar;
#endif
}

bool gBasicOnly = false, gScalarOnly = false;
const Kernel *gKernel = pick(false, false);   // αλλάζει μόνο με τους διακόπτες δοκιμής

}  // namespace

void hybridUseBasicOnly(bool on)
{
    gBasicOnly = on;
    gKernel = pick(gBasicOnly, gScalarOnly);
}

void hybridUseScalar(bool on)
{
    gScalarOnly = on;
    gKernel = pick(gBasicOnly, gScalarOnly);
}

const char *hybridKernelText() { return gKernel->text; }

void hybridReset(HybridState &st)
{
#if !(defined(MP3_HYBRID_MUTANT) && MP3_HYBRID_MUTANT == 10)   // σκόπιμο σφάλμα 10: η παλιά επικάλυψη μένει
    std::memset(st.overlap, 0, sizeof st.overlap);
#endif
    st.live = 0;
}

void hybrid(HybridState &st, float *x, int blockType, bool mixed, int subbands, float *out)
{
    const Kernel &k = *gKernel;
    const int nz = subbands < 0 ? 0 : subbands > 32 ? 32 : subbands;
    blockType &= 3;
#if defined(MP3_HYBRID_MUTANT) && MP3_HYBRID_MUTANT == 6
    mixed = false;   // σκόπιμο σφάλμα: το μικτό μπλοκ σαν σκέτο κοντό
#endif
    int live;   // πόσες πρώτες υποζώνες έχουν κάτι να μετασχηματίσουν (βλ. «ΤΟ `subbands`» στην κορυφή)
    if (blockType != 2) {
#if defined(MP3_HYBRID_MUTANT) && MP3_HYBRID_MUTANT == 4
        live = nz;   // σκόπιμο σφάλμα: ξεχνά την υποζώνη που γεμίζει η μείωση αναδίπλωσης
#else
        live = nz == 0 ? 0 : nz < 32 ? nz + 1 : 32;
#endif
        const int done = (live + k.lanes - 1) & ~(k.lanes - 1);
        k.runLong(x, st.overlap, out, done, &kT.win[blockType]);
        k.runRest(st.overlap, out, done, st.live);
    } else if (!mixed) {
        live = nz;
        const int done = (live + k.lanes - 1) & ~(k.lanes - 1);
        k.runShort(x, st.overlap, out, 0, done);
        k.runRest(st.overlap, out, done, st.live);
    } else {
        // Μικτό μπλοκ: σπάνιο, πάει ολόκληρο από τον απλό δρόμο, μία υποζώνη τη φορά.
        live = nz == 0 ? 0 : nz < 2 ? 2 : nz;
        if (live) {
            scalar::mixedLong(x, st.overlap, out, &kT.win[0]);
            scalar::runShort(x, st.overlap, out, 2, live);
        }
        scalar::runRest(st.overlap, out, live, st.live);
    }
    st.live = live;
#if defined(MP3_HYBRID_MUTANT) && MP3_HYBRID_MUTANT == 2
    for (int r = 1; r < 18; r += 2) out[32 * r + 5] = -out[32 * r + 5];   // σκόπιμο σφάλμα: η υποζώνη 5 χωρίς αλλαγή προσήμου
#endif
}

}  // namespace orpheic::mp3
