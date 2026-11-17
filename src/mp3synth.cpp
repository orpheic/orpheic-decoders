// Η τράπεζα σύνθεσης του MP3 σε απλή ακρίβεια, γρήγορα. Τι υπόσχεται: mp3synth.h. Ο μάρτυράς της: exact::synth (mp3exact.h).
//
// Όλα εδώ βγαίνουν από τους τύπους του σχήματος A.2 του ISO/IEC 11172-3, όπως είναι γραμμένοι στο exact::synth, με τριγωνομετρία
// σχολείου. Κάθε βήμα είναι γραμμένο ώστε να ελέγχεται με το χέρι, και όλα μαζί ελέγχονται αριθμητικά, δείγμα προς δείγμα,
// απέναντι στο exact::synth.
// Τα αποτελέσματα δεν είναι καινούρια, μόνο η παραγωγή τους ξαναγράφτηκε εδώ βήμα βήμα για να ελέγχεται. Ότι ο πίνακας της σύνθεσης
// είναι ένας DCT 32 σημείων το έδειξε ο K. Konstantinides («Fast Subband Filtering in MPEG Audio Coding», IEEE Signal Processing
// Letters 1(2), 1994). Για τον DCT-IV με μιγαδικό FFT μισού μήκους και μία περιστροφή πριν και μία μετά, η συνηθισμένη αναφορά είναι
// P. Duhamel, Y. Mahieux και J.-P. Petit («A Fast Algorithm for the Implementation of Filter Banks Based on Time Domain Aliasing
// Cancellation», ICASSP 1991). Το βήμα ζυγά/μονά του FFT είναι των J. W. Cooley και J. W. Tukey («An Algorithm for the Machine
// Calculation of Complex Fourier Series», Mathematics of Computation 19, 1965).
//
// ── ΤΟ ΔΙΑΓΡΑΜΜΑ ──
// Για κάθε χρονική θέση (32 δείγματα υποζωνών S[0 … 31]):
//   (1) τα 1024 V γλιστρούν 64 θέσεις · (2) V[i] = Σ_k cos((16+i)(2k+1)π/64)·S[k], i = 0 … 63 ·
//   (3) out[j] = Σ_{i=0…7} ( V[128i + j]·D[64i + j] + V[128i + 96 + j]·D[64i + 32 + j] ), j = 0 … 31.
// Γράφω V_t για τις 64 τιμές που υπολογίστηκαν πριν από t θέσεις (t = 0 η τωρινή). Τότε V[128i + j] = V_{2i}[j] και
// V[128i + 96 + j] = V_{2i+1}[32 + j], άρα το (3) είναι
//   out[j] = Σ_{t=0…15} D[32t + j]·V_t[j + (t μονό ? 32 : 0)].                                                        (Π)
// Κάθε δείγμα εξόδου είναι 16 όροι, ένας από καθεμία από τις 16 τελευταίες θέσεις: από τις ζυγές t το πρώτο μισό του V, από
// τις μονές το δεύτερο.
//
// ── V ΑΠΟ C: ΜΟΝΟ 32 ΚΑΙΝΟΥΡΙΟΙ ΑΡΙΘΜΟΙ ΑΝΑ ΘΕΣΗ ──
// Έστω C[n] = Σ_k S[k]·cos((2k+1)·n·π/64) για κάθε ακέραιο n. Τότε V[i] = C[16 + i], και επειδή το 2k+1 είναι περιττό:
//   C[−n] = C[n]              (το συνημίτονο είναι άρτιο)
//   C[64 − n] = −C[n]         (cos((2k+1)π − x) = −cos x)
//   C[64 + n] = −C[n]         (τα δύο προηγούμενα μαζί)
//   C[32] = 0                 (cos((2k+1)π/2) = 0)
// άρα όλα τα V γράφονται με τα C[0 … 31], που είναι ο DCT-II 32 σημείων του S:
//   V[i]  =  C[16 + i]     i = 0 … 15
//   V[16] =  0
//   V[i]  = −C[48 − i]     i = 17 … 48        (C[16 + i] = C[64 − (48 − i)])
//   V[i]  = −C[i − 48]     i = 48 … 63        (C[16 + i] = C[64 + (i − 48)])
//
// ── ΤΟ ΠΑΡΑΘΥΡΟ ──
// Βάζω τα παραπάνω στο (Π), για j = 1 … 16:
//   out[j]      =  Σ_{t ζυγό} D[32t + j]·C_t[16 + j]       −  Σ_{t μονό} D[32t + j]·C_t[16 − j]
//   out[32 − j] = −Σ_{t ζυγό} D[32t + 32 − j]·C_t[16 + j]  −  Σ_{t μονό} D[32t + 32 − j]·C_t[16 − j]
// (στο j = 16 οι δύο γραμμές είναι το ίδιο δείγμα και συμφωνούν, γιατί C[32] = 0), και χωριστά
//   out[0]      =  Σ_{t ζυγό} D[32t]·C_t[16]  −  Σ_{t μονό} D[32t]·C_t[16].
// Δηλαδή τα out[j] και out[32 − j] θέλουν ΤΙΣ ΙΔΙΕΣ τιμές από κάθε θέση: από τις ζυγές τα C[17 … 32], από τις μονές τα
// C[15 … 0]. Γι' αυτό η σειρά μιας θέσης στην ιστορία είναι E[m] = C[17 + m] (με E[15] = C[32] = 0) και O[m] = C[15 − m],
// m = j − 1 = 0 … 15, και το C[16] πάει χωριστά (SynthState::mid).
//
// Το παράθυρο είναι συμμετρικό: D[512 − i] = −D[i] για κάθε i που δεν είναι πολλαπλάσιο του 64. Ο πίνακας B.3 το έχει ακριβώς,
// ακέραιο προς ακέραιο: το ελέγχει ρητά η δοκιμή (και αν έλειπε, η σύγκριση με το exact::synth θα αποτύχαινε). Με G[t][m] = D[32t + m + 1]:
//   D[32t + 32 − j] = D[512 − (32(15−t) + j)] = −G[15−t][m]
// και η δεύτερη γραμμή γίνεται   out[32 − j] = Σ_{t ζυγό} G[15−t][m]·E_t[m] + Σ_{t μονό} G[15−t][m]·O_t[m].
// Άρα φτάνει μισός πίνακας παραθύρου (16 × 16 τιμές), και οι όροι t και 15 − t μοιράζονται και τις δύο σειρές του παραθύρου και
// τις δύο σειρές της ιστορίας (mp3synth_kernel.inc, ΜΕΡΟΣ 2).
// Οι συντελεστές είναι ακέραιοι ως 75038 διά 65536: χωρούν ΑΚΡΙΒΩΣ σε float, άρα το παράθυρο δεν προσθέτει δικό του σφάλμα.
//
// ── Ο DCT ──
// (α) DCT-II N σημείων, C[n] = Σ_{k<N} x[k]·cos((2k+1)nπ/(2N)). Ο όρος του N−1−k έχει cos(nπ − (2k+1)nπ/(2N)) =
//     (−1)^n·cos((2k+1)nπ/(2N)), άρα C[n] = Σ_{k<N/2} (x[k] + (−1)^n·x[N−1−k])·cos((2k+1)nπ/(2N)):
//       C[2m]   = Σ_{k<N/2} (x[k] + x[N−1−k])·cos((2k+1)mπ/N)              DCT-II N/2 σημείων των αθροισμάτων
//       C[2m+1] = Σ_{k<N/2} (x[k] − x[N−1−k])·cos((2k+1)(2m+1)π/(2N))      DCT-IV N/2 σημείων των διαφορών
//     Για N = 2: C[0] = x0 + x1, C[1] = (x0 − x1)·cos(π/4).
// (β) DCT-IV M σημείων, X[m] = Σ_{k<M} d[k]·cos θ(k,m), με θ(k,m) = (2k+1)(2m+1)π/(4M). Για p = 0 … M/2−1 ορίζω
//       z[p] = Σ_{n<M/2} (d[2n] − j·d[M−1−2n])·e^(jθ(2n,2p)).
//     Επειδή θ(M−1−2n, 2p) = (4p+1)π/2 − θ(2n,2p), είναι cos θ(M−1−2n,2p) = sin θ(2n,2p) και sin θ(M−1−2n,2p) = cos θ(2n,2p). Άρα
//       Re z[p] = Σ_n d[2n]·cos θ(2n,2p) + d[M−1−2n]·sin θ(2n,2p) = Σ_k d[k]·cos θ(k,2p) = X[2p]
//     (τα k = 2n είναι όλα τα ζυγά k, τα k = M−1−2n όλα τα μονά). Και επειδή θ(k, M−1−2p) = (2k+1)π/2 − θ(k,2p), είναι
//     cos θ(k, M−1−2p) = (−1)^k·sin θ(k,2p), άρα
//       X[M−1−2p] = Σ_n d[2n]·sin θ(2n,2p) − d[M−1−2n]·sin θ(M−1−2n,2p) = Σ_n d[2n]·sin θ(2n,2p) − d[M−1−2n]·cos θ(2n,2p) = Im z[p].
//     Τέλος θ(2n,2p) = (4n+1)(4p+1)π/(4M) = (4n+1)π/(4M) + 2π·n·p/(M/2) + pπ/M: ο εκθέτης χωρίζει σε τρία κομμάτια,
//       u[n] = (d[2n] − j·d[M−1−2n])·e^(jπ(4n+1)/(4M))      μία περιστροφή ανά n
//       U[p] = Σ_n u[n]·e^(2πj·n·p/(M/2))                     DFT M/2 σημείων
//       z[p] = U[p]·e^(jπp/M)                                 μία περιστροφή ανά p (καμία για p = 0)
// (γ) DFT L σημείων, U[p] = Σ_n u[n]·w^(np), w = e^(2πj/L). Χωρίζοντας ζυγά και μονά n: U[p] = E[p] + w^p·O[p], όπου E, O οι DFT
//     L/2 σημείων των u[2n'] και u[2n'+1], και επειδή w^(L/2) = −1: U[p + L/2] = E[p] − w^p·O[p].
// Για N = 32 αυτά δίνουν 236 πράξεις όπου υπάρχει FMA (289 χωρίς) για τις 32 τιμές, όλες αθροίσματα, διαφορές και περιστροφές.
//
// ── Η ΔΙΑΤΑΞΗ ΤΗΣ ΔΟΥΛΕΙΑΣ ──
// Πρώτα ο DCT ΟΛΩΝ των 18 θέσεων, μετά το παράθυρο όλων. Οι 18 σειρές γράφονται στο μισό της ιστορίας που ΔΕΝ έγραψε η
// προηγούμενη κλήση· το παράθυρο διαβάζει τις 15 παλιές σειρές από το άλλο μισό, και καμία σειρά δεν αντιγράφεται (mp3synth.h).
//   · Ο DCT δουλεύει «κατά μήκος του χρόνου»: ένα διάνυσμα κρατά την ίδια ποσότητα για 8 (AVX2) ή 4 (SSE2, NEON) διαδοχικές
//     θέσεις. Έτσι ο βαθμωτός αλγόριθμος τρέχει αυτούσιος, χωρίς καμία ανακατάταξη μέσα σε διάνυσμα, και μόνο η είσοδος και η
//     έξοδός του θέλουν μια αναστροφή 4×4 (gather, scatter). Οι 18 θέσεις είναι 8 + 8 και μία τετράδα για τις 14 … 17 (AVX2), ή
//     πέντε τετράδες 0, 4, 8, 12, 14 (SSE2, NEON)· η τελευταία τετράδα ξαναϋπολογίζει τις θέσεις 14 και 15.
//   · Το παράθυρο δουλεύει «κατά πλάτος της εξόδου»: 8 (ή 4) γειτονικά δείγματα εξόδου μαζί, δύο χρονικές θέσεις σε κάθε γύρο,
//     και όλες οι αναγνώσεις είναι συνεχόμενες και στοιχισμένες.
//   · Επειδή όλη η είσοδος έχει διαβαστεί πριν γραφτεί η πρώτη έξοδος, in και out μπορούν να είναι ο ίδιος πίνακας.
//
// ── ΤΙ ΜΕΤΡΗΘΗΚΕ ──
// i9-14900K, πυρήνας P στα 5,6 GHz περίπου, g++ 13.3 -O2, ένα νήμα, η πιο γρήγορη από πολλές κοντές μετρήσεις. Το μηχάνημα ήταν
// φορτωμένο: οι αριθμοί παίζουν ±3%. «ns» = νανοδευτερόλεπτα ανά κλήση (ένα granule ενός καναλιού), εκτός αν λέει αλλιώς.
//   Ολόκληρο:      AVX2+FMA 203 · SSE2 500 · σκέτη C++ 1530 · κατά γράμμα (exact::synth, double) 8500
//   Μέρη, AVX2:    ο DCT των 18 θέσεων 86 · το παράθυρο των 18 θέσεων μαζί με το κατέβασμα της ιστορίας 121
// Τι δοκιμάστηκε και ΔΕΝ κρατήθηκε:
//   · DCT με την αναδρομή «διαφορά επί 1/(2cos), και ξανά DCT-II» αντί για περιστροφές: ίδιος χρόνος (35,0 αντί για 33,5 ns
//     η ομάδα των 8 θέσεων), αλλά σφάλμα 1,6e-6 ως 2,7e-6 όταν η ενέργεια είναι στις υποζώνες 15 και 16 (εκεί ο συντελεστής της
//     φτάνει το 10,2 και μετά αφαιρούνται μεγάλοι αριθμοί). Πάνω από την απαίτηση· οι περιστροφές δίνουν 2,9e-7 στην ίδια είσοδο.
//   · DCT ως πίνακας 32×32 μέσα σε κάθε θέση: 260 ns για τις 18 θέσεις αντί για 86.
//   · Παράθυρο κατά γράμμα πάνω σε ολόκληρο το V (64 τιμές ανά θέση, ολόκληρος ο πίνακας D): 181 ns αντί για 121, και θέλει
//     διπλάσιο scatter.
//   · Παράθυρο μία θέση τη φορά, με την ιστορία να κατεβαίνει σε ξεχωριστό πέρασμα στο τέλος: ολόκληρο 242 αντί για 225 ns.
//   · Τρίτη ομάδα των 8 (θέσεις 10 … 17) αντί για την τετράδα: ολόκληρο 220 αντί για 211 ns.
//   · Αναστροφή 4×4 με unpcklps/unpckhps αντί για shufps: ο DCT των 18 θέσεων 100 αντί για 94 ns. Σε αυτούς τους πυρήνες
//     το unpcklps τρέχει μία φορά ανά κύκλο και το shufps δύο.
//   · Οι πράξεις του DCT ανακατεμένες ώστε περιστροφές και αθροίσματα να εναλλάσσονται, και κάθε τετράδα του scatter μόλις είναι
//     έτοιμη: 204 κύκλοι η ομάδα αντί για 187, γιατί περισσότερες τιμές καταλήγουν στη στοίβα.
//   · Οι αφαιρέσεις μέσα από τη μονάδα του FMA (a − b·1), για να μοιραστεί το φορτίο ανάμεσα στις μονάδες: καμία διαφορά (±2%).
//   · Το scatter μιας ομάδας και το gather της επόμενης «φιλοξενούμενα» μέσα στον βρόχο του παραθύρου, με τον DCT να περνά από
//     τη μνήμη (περνά τη δοκιμή): 231 ως 242 ns όταν το παραδοτέο της ίδιας ώρας έδινε 216 ως 229.
//   · 14/11/2026, ΑΛΛΑΞΕ: η ιστορία ήταν ένα κομμάτι 33 σειρών (15 παλιές, 18 νέες) και οι 15 τελευταίες κατέβαιναν στην αρχή
//     μέσα στο πέρασμα του παραθύρου. Τώρα είναι δύο μισά των 18 που εναλλάσσονται, χωρίς αντιγραφή σειρών. Η έξοδος είναι η ίδια
//     ΜΠΙΤ ΠΡΟΣ ΜΠΙΤ (20.000 granules ανά δρόμο, μαζί με in = out, μηδενισμό και αντιγραφή της κατάστασης). Χρόνος ολόκληρου, νέο
//     προς παλιό, δεμένο σε έναν πυρήνα, πέντε τρεξίματα: i9-14900K AVX2+FMA ίσο (από −1% ως +1%), SSE2 ίσο μέσα στον θόρυβο
//     (από −2% ως +5%), σκέτη C++ +3%· Raspberry Pi 5 NEON 1321 αντί για 1495 ns (−12%), σκέτη C++ −13%. Οι αριθμοί «Ολόκληρο»
//     και «Μέρη» πιο πάνω είναι της παλιάς διάταξης.
//   · Τα 33 C[16] αντιγραμμένα σε τοπικό πίνακα ΑΚΡΙΒΩΣ ΠΡΙΝ από το παράθυρο (αντί πριν από τον DCT): ολόκληρο 206 αντί για 198 ns.
//     Οι φαρδιές αναγνώσεις του παραθύρου έπεφταν πάνω σε εγγραφές που δεν είχαν τελειώσει.
//   · ΜΕΣΑ ΣΕ ΟΛΟΚΛΗΡΟ ΤΟΝ ΑΠΟΚΩΔΙΚΟΠΟΙΗΤΗ μετρά ο ΟΓΚΟΣ του κώδικα μηχανής, που το στάδιο μόνο του δεν τον δείχνει (36 αρχεία,
//     δεμένο σε έναν πυρήνα, ο καλύτερος από πέντε ή έξι γύρους εναλλάξ· perf stat για τις αστοχίες της μνήμης εντολών L1):
//       παλιά διάταξη (δρόμος AVX2 15 kB)                                   299 ms · 447 εκατ. αστοχίες
//       δύο μισά, εννέα αντίγραφα του παραθύρου, ένα ανά s (22 kB)          306 ms
//       το ίδιο, με συνεχόμενο βρόχο στις μισές κλήσεις (22 kB σε δύο συναρτήσεις)  302 ms · 579 εκατ.
//       δύο αντίγραφα, ένα ανά μισό του granule, η επιλογή στις μισές αναγνώσεις (18 kB)  301,2 ms · 494 εκατ.
//       ΕΝΑ αντίγραφο, η επιλογή μισού σε κάθε ανάγνωση (17 kB)             300,3 ms · 445 εκατ.   ← αυτό κρατήθηκε
//     Το στάδιο μόνο του με το ένα αντίγραφο: AVX2 +2% ως +3%, SSE2 −6%, σκέτη C++ −6% απέναντι στην παλιά διάταξη.
//     Οι χρονομετρήσεις μόνες τους παραπλανούν εδώ: ο παλιός κώδικας με μια άχρηστη συνάρτηση δίπλα του έδωσε 302,5, σαν να ήταν
//     «θέμα διάταξης του εκτελέσιμου». Οι μετρητές του επεξεργαστή έδειξαν ότι εκείνο το χτίσιμο είχε τις ίδιες αστοχίες και
//     τους ίδιους κύκλους με το παλιό: ήταν θόρυβος. Το αληθινό ήταν οι +30% αστοχίες της μνήμης εντολών.
//   · Πίνακας 17 δεικτών ανά ζεύγος θέσεων αντί για επιλογή ανά ανάγνωση: το στάδιο +17% ως +33% (AVX2). Δεν κρατήθηκε.
//   · Το ίδιο το synth με την έξοδο σε 16 διαφορετικές αποστάσεις από την κατάσταση: 205 ns παντού. Ο χρόνος δεν εξαρτάται
//     από το πού βρίσκονται οι πίνακες του καλούντος.
#include "mp3synth.h"
#include "mp3consts.h"

#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

// Άλλοι μεταγλωττιστές: ο ίδιος κώδικας παντού· αλλάζει μόνο ο τρόπος που ζητάμε «βάλ' το μέσα» και «βγάλε εντολές AVX2 σε
// αυτή τη συνάρτηση» (ο MSVC δεν θέλει το δεύτερο: δέχεται κάθε intrinsic παντού).
// Ο MSVC ΔΕΝ ΕΧΕΙ ΔΟΚΙΜΑΣΤΕΙ (δεν υπήρχε στο μηχάνημα όπου γράφτηκε). Δύο πράγματα να κοιταχτούν την πρώτη φορά εκεί: ότι χτίζει,
// και τι βγάζει για τις εντολές των 128 bit ΜΕΣΑ στις συναρτήσεις του δρόμου AVX2 (τα _mm_loadu_ps / _mm_storeu_ps των αναστροφών):
// αν τις κωδικοποιεί χωρίς VEX δίπλα σε εντολές των 256 bit, οι παλιότεροι επεξεργαστές Intel πληρώνουν ακριβά κάθε εναλλαγή. Τότε
// η θεραπεία είναι ο δρόμος AVX2 σε δικό του αρχείο χτισμένο με /arch:AVX2 (είναι ένα #include του mp3synth_kernel.inc).
#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#define SYN_INLINE static __forceinline
#define SYN_TARGET_AVX2
#else
#define SYN_INLINE static inline __attribute__((always_inline))
#define SYN_TARGET_AVX2 __attribute__((target("avx2,fma")))
#endif
#if defined(__x86_64__) || defined(_M_X64)
#define SYN_X86 1
#if !defined(_MSC_VER) || defined(__clang__)
#include <cpuid.h>
#endif
#include <immintrin.h>
#else
#define SYN_X86 0
#endif
#if defined(__aarch64__) && !defined(_MSC_VER)
#define SYN_NEON 1
#include <arm_neon.h>
#else
#define SYN_NEON 0
#endif

namespace orpheic::mp3 {
namespace {

namespace tab {
#include "mp3tables.inc"
}

// Οι σταθερές. Κάθε σταθερά των περιστροφών είναι γραμμένη 8 φορές στη σειρά, ώστε ο πυρήνας να τη διαβάζει ως έτοιμο
// διάνυσμα κατευθείαν μέσα στην εντολή του πολλαπλασιασμού (ο βαθμωτός δρόμος διαβάζει την πρώτη, ο SSE2 τις 4 πρώτες).
struct SynthConsts {
    alignas(32) float pre16[16][8];    // cos, sin του π(4n+1)/64, n = 0 … 7      η περιστροφή πριν από τον DFT, M = 16
    alignas(32) float pre8[8][8];      // … του π(4n+1)/32, n = 0 … 3             M = 8
    alignas(32) float pre4[4][8];      // … του π(4n+1)/16, n = 0, 1              M = 4
    alignas(32) float pre2[2][8];      // … του π/8                               M = 2
    alignas(32) float post16[16][8];   // cos, sin του πp/16, p = 0 … 7           η περιστροφή μετά τον DFT, M = 16
    alignas(32) float post8[8][8];     // … του πp/8, p = 0 … 3                   M = 8
    alignas(32) float half[8];         // √½ = cos(π/4)
    alignas(32) float g[16][16];       // G[t][m] = D[32t + m + 1]                το μισό παράθυρο
    alignas(32) float g0[16][8];       // g0[u] = ±D[32(15−u)], «+» για 15−u ζυγό: οι συντελεστές του out[0] στη σειρά της μνήμης

    static void put(float (*dst)[8], long m)   // η γωνία π·m/576 (mp3consts.h)
    {
        for (int i = 0; i < 8; ++i) {
            dst[0][i] = float(consts::cosPi576(m));
            dst[1][i] = float(consts::sinPi576(m));
        }
    }
    SynthConsts()
    {
        for (int n = 0; n < 8; ++n) put(pre16 + 2 * n, 9 * (4 * n + 1));    // π(4n+1)/64
        for (int n = 0; n < 4; ++n) put(pre8 + 2 * n, 18 * (4 * n + 1));    // π(4n+1)/32
        for (int n = 0; n < 2; ++n) put(pre4 + 2 * n, 36 * (4 * n + 1));    // π(4n+1)/16
        put(pre2, 72);                                                       // π/8
        for (int p = 0; p < 8; ++p) put(post16 + 2 * p, 36 * p);           // πp/16
        for (int p = 0; p < 4; ++p) put(post8 + 2 * p, 72 * p);             // πp/8
        for (int i = 0; i < 8; ++i) half[i] = float(std::sqrt(0.5));
        for (int t = 0; t < 16; ++t)
            for (int m = 0; m < 16; ++m) g[t][m] = float(tab::kWindow[32 * t + m + 1]) / 65536.0f;   // ακριβές: |ακέραιος| < 2^24
        for (int u = 0; u < 16; ++u) {
            const int t = 15 - u;
            for (int i = 0; i < 8; ++i) g0[u][i] = float((t & 1) ? -tab::kWindow[32 * t] : tab::kWindow[32 * t]) / 65536.0f;
        }
    }
};

const SynthConsts &consts()
{
    static const SynthConsts k;
    return k;
}

// Τα C[16] (ένας αριθμός ανά θέση) μένουν στη σειρά του χρόνου, γιατί το παράθυρο τα διαβάζει κατά μήκος του χρόνου: τα 15 τελευταία
// του προηγούμενου granule μπαίνουν μπροστά ΠΡΙΝ από τον DCT, ο DCT γράφει τα 18 νέα δίπλα τους (το `mc`), και στο τέλος τα 18 νέα
// φυλάγονται για την επόμενη κλήση. Η αντιγραφή γίνεται νωρίς επίτηδες: αν γινόταν ακριβώς πριν από το παράθυρο, οι φαρδιές
// αναγνώσεις του θα έπεφταν πάνω σε εγγραφές που δεν έχουν ακόμη τελειώσει (μετρημένο: +4% σε ολόκληρο το στάδιο).
#define SYN_MID_BEGIN                                                                                                    \
    float mid[40];                                                                                                      \
    std::memcpy(mid, st.mid + 3, 15 * sizeof(float));                                                                   \
    float *const mc = mid + 15;
#define SYN_MID_END std::memcpy(st.mid, mc, 18 * sizeof(float));

// ── Ο ΔΡΟΜΟΣ ΤΗΣ ΣΚΕΤΗΣ C++ ── ένα «διάνυσμα» είναι ένας float: ο DCT τρέχει μία θέση τη φορά, το παράθυρο ένα ζεύγος δειγμάτων.
#define SYN_NAME(x) x##Plain
#define SYN_ATTR
#define SYN_LANES 1
#define SYN_V float
#define SYN_ADD(a, b) ((a) + (b))
#define SYN_SUB(a, b) ((a) - (b))
#define SYN_MUL(a, b) ((a) * (b))
#define SYN_MADD(a, b, c) ((a) * (b) + (c))
#define SYN_NMADD(a, b, c) ((c) - (a) * (b))
#define SYN_ZERO() 0.0f
#define SYN_TW(p) ((p)[0])
#define SYN_LOAD(p) (*(p))
#define SYN_LOADU(p) (*(p))
#define SYN_STORE(p, v) (*(p) = (v))
#define SYN_STOREU(p, v) (*(p) = (v))
#define SYN_REV(v) (v)
#define SYN_GATHER4(in, o, x0, x1, x2, x3)                                                                              \
    {                                                                                                                   \
        x0 = (in)[(o)];                                                                                                 \
        x1 = (in)[(o) + 1];                                                                                             \
        x2 = (in)[(o) + 2];                                                                                             \
        x3 = (in)[(o) + 3];                                                                                             \
    }
#define SYN_SCATTER4(row, o, q0, q1, q2, q3)                                                                            \
    {                                                                                                                   \
        (row)[(o)] = (q0);                                                                                              \
        (row)[(o) + 1] = (q1);                                                                                          \
        (row)[(o) + 2] = (q2);                                                                                          \
        (row)[(o) + 3] = (q3);                                                                                          \
    }
#define SYN_WANT_DCT 1
#define SYN_WANT_WINDOW 1
#include "mp3synth_kernel.inc"

void synthPlain(SynthState &st, const float *in, float *out)
{
    const SynthConsts &k = consts();
    float *cur = st.rows[st.last ^ 1];
    SYN_MID_BEGIN
    for (int s = 0; s < 18; ++s) dctPlain(k, in + 32 * s, cur + 32 * s, mc + s);
    windowPlain(k, st.rows[st.last], cur, mid, out);
    SYN_MID_END
    st.last ^= 1;
}

#if SYN_X86
// ── SSE2 ── ό,τι έχει κάθε x86-64. Τέσσερις χρονικές θέσεις ανά διάνυσμα στον DCT, τέσσερα δείγματα εξόδου στο παράθυρο.
// Η μνήμη μας διαβάζεται με στοιχισμένες εντολές: στους επεξεργαστές του 2008 η ανάγνωση χωρίς στοίχιση είναι δύο μικροεντολές
// αντί για μία, και μόνο η στοιχισμένη μπαίνει μέσα στην εντολή του πολλαπλασιασμού.
SYN_INLINE __m128 revSse(__m128 v) { return _mm_shuffle_ps(v, v, 0x1B); }
// Αναστροφή 4×4: από r0 … r3 (μία θέση το καθένα, 4 υποζώνες) σε x0 … x3 (μία υποζώνη το καθένα, 4 θέσεις), και ανάποδα.
// Πρώτα ζεύγη γραμμών (t0 = r0[0] r0[1] r1[0] r1[1], t1 = r0[2] r0[3] r1[2] r1[3], το ίδιο για r2, r3), μετά από κάθε ζεύγος τα
// ζυγά ή τα μονά στοιχεία. Και τα οκτώ βήματα είναι η ΙΔΙΑ εντολή (shufps), επίτηδες: η συνηθισμένη αναστροφή με unpcklps και
// unpckhps τρέχει στους πυρήνες P της Intel μία εντολή ανά κύκλο, ενώ το shufps δύο (μετρημένο).
#define SYN_T4(r0, r1, r2, r3, x0, x1, x2, x3)                                                                          \
    {                                                                                                                   \
        const SYN_V t0 = SYN_SHUF(r0, r1, 0x44), t1 = SYN_SHUF(r0, r1, 0xEE), t2 = SYN_SHUF(r2, r3, 0x44), t3 = SYN_SHUF(r2, r3, 0xEE); \
        x0 = SYN_SHUF(t0, t2, 0x88);                                                                                    \
        x1 = SYN_SHUF(t0, t2, 0xDD);                                                                                    \
        x2 = SYN_SHUF(t1, t3, 0x88);                                                                                    \
        x3 = SYN_SHUF(t1, t3, 0xDD);                                                                                    \
    }
#define SYN_GATHER4_128(in, o, x0, x1, x2, x3)                                                                          \
    {                                                                                                                   \
        const __m128 g0 = _mm_loadu_ps((in) + (o)), g1 = _mm_loadu_ps((in) + 32 + (o));                                 \
        const __m128 g2 = _mm_loadu_ps((in) + 64 + (o)), g3 = _mm_loadu_ps((in) + 96 + (o));                            \
        SYN_T4(g0, g1, g2, g3, x0, x1, x2, x3)                                                                          \
    }
#define SYN_SCATTER4_128(row, o, q0, q1, q2, q3)                                                                        \
    {                                                                                                                   \
        __m128 y0, y1, y2, y3;                                                                                          \
        SYN_T4(q0, q1, q2, q3, y0, y1, y2, y3)                                                                          \
        _mm_store_ps((row) + (o), y0);                                                                                  \
        _mm_store_ps((row) + 32 + (o), y1);                                                                             \
        _mm_store_ps((row) + 64 + (o), y2);                                                                             \
        _mm_store_ps((row) + 96 + (o), y3);                                                                             \
    }

#define SYN_NAME(x) x##Sse2
#define SYN_ATTR
#define SYN_LANES 4
#define SYN_V __m128
#define SYN_ADD(a, b) _mm_add_ps(a, b)
#define SYN_SUB(a, b) _mm_sub_ps(a, b)
#define SYN_MUL(a, b) _mm_mul_ps(a, b)
#define SYN_MADD(a, b, c) _mm_add_ps(_mm_mul_ps(a, b), c)
#define SYN_NMADD(a, b, c) _mm_sub_ps(c, _mm_mul_ps(a, b))
#define SYN_ZERO() _mm_setzero_ps()
#define SYN_TW(p) _mm_load_ps(p)
#define SYN_LOAD(p) _mm_load_ps(p)
#define SYN_LOADU(p) _mm_loadu_ps(p)
#define SYN_STORE(p, v) _mm_store_ps(p, v)
#define SYN_STOREU(p, v) _mm_storeu_ps(p, v)
#define SYN_REV(v) revSse(v)
#define SYN_SHUF(a, b, i) _mm_shuffle_ps(a, b, i)
#define SYN_GATHER4 SYN_GATHER4_128
#define SYN_SCATTER4 SYN_SCATTER4_128
#define SYN_WANT_DCT 1
#define SYN_WANT_WINDOW 1
#include "mp3synth_kernel.inc"

void synthSse2(SynthState &st, const float *in, float *out)
{
    const SynthConsts &k = consts();
    // 18 θέσεις σε πέντε τετράδες: 0, 4, 8, 12 και 14 (η τελευταία ξαναγράφει τις 14 και 15 με τις ίδιες τιμές).
    float *cur = st.rows[st.last ^ 1];
    SYN_MID_BEGIN
    for (int s = 0; s < 16; s += 4) dctSse2(k, in + 32 * s, cur + 32 * s, mc + s);
    dctSse2(k, in + 32 * 14, cur + 32 * 14, mc + 14);
    windowSse2(k, st.rows[st.last], cur, mid, out);
    SYN_MID_END
    st.last ^= 1;
}

// ── AVX2 + FMA ── οκτώ χρονικές θέσεις ανά διάνυσμα στον DCT, οκτώ δείγματα εξόδου στο παράθυρο, και a·b + c σε μία εντολή.
// Εδώ οι αναγνώσεις είναι «χωρίς στοίχιση»: με την κωδικοποίηση VEX κοστίζουν το ίδιο και μπαίνουν μέσα στις εντολές.
//
// Η σειρά των κομματιών του DCT ανάμεσα στις τρεις ομάδες είναι μετρημένη, και εξαρτάται από τον μεταγλωττιστή (το αποτέλεσμα
// είναι ακριβώς το ίδιο σε κάθε περίπτωση):
//   · g++: οι πράξεις της επόμενης ομάδας ΠΡΙΝ από το σκόρπισμα της προηγούμενης, όλα μέσα σε μία συνάρτηση. Περίπου 3%
//     λιγότερος χρόνος από το «η μία ομάδα μετά την άλλη»: οι ανακατατάξεις του σκορπίσματος και οι πολλαπλασιασμοί
//     των πράξεων τρέχουν σε διαφορετικές μονάδες του επεξεργαστή και έτσι βρίσκονται μαζί μπροστά του.
//   · clang++ και ό,τι άλλο: η μία ομάδα μετά την άλλη, καθεμία σε δική της κλήση. Ο clang 18, όταν τα δει όλα σε μία συνάρτηση,
//     κρατά πολύ περισσότερες τιμές στη στοίβα (534 εγγραφές αντί για 210 του g++) και βγαίνει πιο αργός.
#ifndef SYN_AVX2_INTERLEAVED
#if defined(__GNUC__) && !defined(__clang__)
#define SYN_AVX2_INTERLEAVED 1
#else
#define SYN_AVX2_INTERLEAVED 0
#endif
#endif
SYN_INLINE SYN_TARGET_AVX2 __m256 revAvx(__m256 v) { return _mm256_permutevar8x32_ps(v, _mm256_setr_epi32(7, 6, 5, 4, 3, 2, 1, 0)); }
// Η αναστροφή γίνεται μέσα σε κάθε μισό των 128 bit χωριστά (οι εντολές unpack/shuffle των 256 bit έτσι δουλεύουν). Γι' αυτό
// το κάτω μισό κάθε διανύσματος εισόδου φορτώνεται από τη θέση i και το πάνω από τη θέση i + 4: μετά την αναστροφή 4×4 κάθε
// διάνυσμα έχει μία υποζώνη για τις θέσεις 0 … 3 κάτω και 4 … 7 πάνω, δηλαδή 8 διαδοχικές θέσεις, χωρίς καμία εντολή που να
// περνά από το ένα μισό στο άλλο. Το ίδιο ανάποδα στην έξοδο.
#define SYN_GATHER4_256(in, o, x0, x1, x2, x3)                                                                          \
    {                                                                                                                   \
        const __m256 g0 = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_loadu_ps((in) + (o))), _mm_loadu_ps((in) + 128 + (o)), 1);        \
        const __m256 g1 = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_loadu_ps((in) + 32 + (o))), _mm_loadu_ps((in) + 160 + (o)), 1);   \
        const __m256 g2 = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_loadu_ps((in) + 64 + (o))), _mm_loadu_ps((in) + 192 + (o)), 1);   \
        const __m256 g3 = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_loadu_ps((in) + 96 + (o))), _mm_loadu_ps((in) + 224 + (o)), 1);   \
        SYN_T4(g0, g1, g2, g3, x0, x1, x2, x3)                                                                          \
    }
#define SYN_SCATTER4_256(row, o, q0, q1, q2, q3)                                                                        \
    {                                                                                                                   \
        __m256 y0, y1, y2, y3;                                                                                          \
        SYN_T4(q0, q1, q2, q3, y0, y1, y2, y3)                                                                          \
        _mm_storeu_ps((row) + (o), _mm256_castps256_ps128(y0));                                                         \
        _mm_storeu_ps((row) + 32 + (o), _mm256_castps256_ps128(y1));                                                    \
        _mm_storeu_ps((row) + 64 + (o), _mm256_castps256_ps128(y2));                                                    \
        _mm_storeu_ps((row) + 96 + (o), _mm256_castps256_ps128(y3));                                                    \
        _mm_storeu_ps((row) + 128 + (o), _mm256_extractf128_ps(y0, 1));                                                 \
        _mm_storeu_ps((row) + 160 + (o), _mm256_extractf128_ps(y1, 1));                                                 \
        _mm_storeu_ps((row) + 192 + (o), _mm256_extractf128_ps(y2, 1));                                                 \
        _mm_storeu_ps((row) + 224 + (o), _mm256_extractf128_ps(y3, 1));                                                 \
    }

#define SYN_NAME(x) x##Avx2
#define SYN_ATTR SYN_TARGET_AVX2
#define SYN_LANES 8
#define SYN_V __m256
#define SYN_ADD(a, b) _mm256_add_ps(a, b)
#define SYN_SUB(a, b) _mm256_sub_ps(a, b)
#define SYN_MUL(a, b) _mm256_mul_ps(a, b)
#define SYN_MADD(a, b, c) _mm256_fmadd_ps(a, b, c)
#define SYN_NMADD(a, b, c) _mm256_fnmadd_ps(a, b, c)
#define SYN_ZERO() _mm256_setzero_ps()
#define SYN_TW(p) _mm256_load_ps(p)
#define SYN_LOAD(p) _mm256_loadu_ps(p)
#define SYN_LOADU(p) _mm256_loadu_ps(p)
#define SYN_STORE(p, v) _mm256_storeu_ps(p, v)
#define SYN_STOREU(p, v) _mm256_storeu_ps(p, v)
#define SYN_REV(v) revAvx(v)
#define SYN_SHUF(a, b, i) _mm256_shuffle_ps(a, b, i)
#define SYN_GATHER4 SYN_GATHER4_256
#define SYN_SCATTER4 SYN_SCATTER4_256
// Στο παράθυρο κάθε σειρά της ιστορίας πολλαπλασιάζεται με δύο σειρές του παραθύρου. Ο g++ τη βάζει τότε ως τελεστή μνήμης και
// στις δύο εντολές, δηλαδή τη διαβάζει δύο φορές, και οι αναγνώσεις γίνονται το στενό σημείο (τρεις ανά κύκλο σε αυτούς τους
// πυρήνες). Η άδεια εντολή assembly παρακάτω δεν βγάζει κώδικα· λέει μόνο «αυτές οι τέσσερις τιμές είναι σε καταχωρητή», άρα
// διαβάζονται μία φορά. Μετρημένο (i9-14900K): το παράθυρο από 137,6 σε 121,5 ns, με τα ίδια ακριβώς bit στην έξοδο.
// Μόνο για τον g++: ο clang δεν διπλασιάζει αναγνώσεις, ο MSVC δεν έχει τέτοια εντολή, και στον δρόμο SSE2 η ίδια αλλαγή μετρήθηκε
// λίγο ΠΙΟ ΑΡΓΗ (358 → 369 ns), γι' αυτό εκεί δεν μπαίνει.
#if defined(__GNUC__) && !defined(__clang__)
#define SYN_KEEP(x0, y0, x1, y1) __asm__ volatile("" : "+x"(x0), "+x"(y0), "+x"(x1), "+x"(y1));
#endif
#define SYN_WANT_DCT (!SYN_AVX2_INTERLEAVED)
#define SYN_WANT_WINDOW 1
#include "mp3synth_kernel.inc"

// Ο ίδιος DCT σε μισό πλάτος (4 θέσεις) αλλά με FMA: για τις θέσεις που περισσεύουν μετά τις δύο οκτάδες (18 = 8 + 8 + 2).
#define SYN_NAME(x) x##Avx2Half
#define SYN_ATTR SYN_TARGET_AVX2
#define SYN_LANES 4
#define SYN_V __m128
#define SYN_ADD(a, b) _mm_add_ps(a, b)
#define SYN_SUB(a, b) _mm_sub_ps(a, b)
#define SYN_MUL(a, b) _mm_mul_ps(a, b)
#define SYN_MADD(a, b, c) _mm_fmadd_ps(a, b, c)
#define SYN_NMADD(a, b, c) _mm_fnmadd_ps(a, b, c)
#define SYN_ZERO() _mm_setzero_ps()
#define SYN_TW(p) _mm_load_ps(p)
#define SYN_LOAD(p) _mm_loadu_ps(p)
#define SYN_LOADU(p) _mm_loadu_ps(p)
#define SYN_STORE(p, v) _mm_storeu_ps(p, v)
#define SYN_STOREU(p, v) _mm_storeu_ps(p, v)
#define SYN_SHUF(a, b, i) _mm_shuffle_ps(a, b, i)
#define SYN_GATHER4 SYN_GATHER4_128
#define SYN_SCATTER4 SYN_SCATTER4_128
#define SYN_WANT_DCT (!SYN_AVX2_INTERLEAVED)
#define SYN_WANT_WINDOW 0
#include "mp3synth_kernel.inc"

// 18 θέσεις = 8 + 8 και μία τετράδα σε μισό πλάτος για τις 14 … 17 (οι 14 και 15 ξαναγράφονται με τις ίδιες τιμές).
SYN_TARGET_AVX2 void synthAvx2(SynthState &st, const float *in, float *out)
{
    const SynthConsts &k = consts();
    float *cur = st.rows[st.last ^ 1];
    SYN_MID_BEGIN
#if SYN_AVX2_INTERLEAVED
    __m256 xa[32], xb[32], ca[32], cb[32];
    __m128 xh[32], ch[32];
    gatherAvx2(in, xa);                                        // θέσεις 0 … 7
    coreAvx2(k, xa, ca);
    gatherAvx2(in + 32 * 8, xb);                               // 8 … 15
    coreAvx2(k, xb, cb);
    scatterAvx2(ca, cur, mc);
    gatherAvx2Half(in + 32 * 14, xh);                          // 14 … 17
    coreAvx2Half(k, xh, ch);
    scatterAvx2(cb, cur + 32 * 8, mc + 8);
    scatterAvx2Half(ch, cur + 32 * 14, mc + 14);
#else
    dctAvx2(k, in, cur, mc);                    // θέσεις 0 … 7
    dctAvx2(k, in + 32 * 8, cur + 32 * 8, mc + 8);           // 8 … 15
    dctAvx2Half(k, in + 32 * 14, cur + 32 * 14, mc + 14);      // 14 … 17
#endif
    windowAvx2(k, st.rows[st.last], cur, mid, out);
    SYN_MID_END
    st.last ^= 1;
}

// Έχει αυτός ο επεξεργαστής AVX2 και FMA, ΚΑΙ σώζει το λειτουργικό τους καταχωρητές των 256 bit όταν αλλάζει νήμα;
bool cpuHasAvx2Fma()
{
#if defined(_MSC_VER) && !defined(__clang__)
    int r[4];
    __cpuid(r, 0);
    if (r[0] < 7) return false;
    __cpuid(r, 1);
    const bool fma = (r[2] & (1 << 12)) != 0, osxsave = (r[2] & (1 << 27)) != 0, avx = (r[2] & (1 << 28)) != 0;
    if (!fma || !osxsave || !avx) return false;
    if ((_xgetbv(0) & 6) != 6) return false;
    __cpuidex(r, 7, 0);
    return (r[1] & (1 << 5)) != 0;
#else
    unsigned a = 0, b = 0, c = 0, d = 0;
    if (!__get_cpuid(1, &a, &b, &c, &d)) return false;
    const bool fma = (c & (1u << 12)) != 0, osxsave = (c & (1u << 27)) != 0, avx = (c & (1u << 28)) != 0;
    if (!fma || !osxsave || !avx) return false;
    unsigned lo = 0, hi = 0;
    __asm__ volatile(".byte 0x0f, 0x01, 0xd0" : "=a"(lo), "=d"(hi) : "c"(0u));   // xgetbv: ο καταχωρητής XCR0
    if ((lo & 6) != 6) return false;                                             // μπιτ 1 και 2: το λειτουργικό σώζει XMM και YMM
    if (!__get_cpuid_count(7, 0, &a, &b, &c, &d)) return false;
    return (b & (1u << 5)) != 0;
#endif
}
#undef SYN_T4
#undef SYN_GATHER4_128
#undef SYN_SCATTER4_128
#undef SYN_GATHER4_256
#undef SYN_SCATTER4_256
#endif   // SYN_X86

#if SYN_NEON
// ── NEON (ARM64) ── καθρέφτης του δρόμου SSE2, με a·b + c σε μία εντολή (κάθε ARM64 την έχει).
// Γράφτηκε σε μηχάνημα x86-64. Τι έχει γίνει και τι όχι:
//   · μεταγλωττίζεται χωρίς προειδοποίηση με clang++ --target=aarch64-linux-gnu·
//   · η λογική του περνά ολόκληρη τη δοκιμή, με τις 14 εντολές NEON που χρησιμοποιεί γραμμένες σε σκέτη C++ όπως τις ορίζει η ARM·
//     τα μαθηματικά είναι έτσι κι αλλιώς το ίδιο κείμενο (mp3synth_kernel.inc) με τους
//     άλλους δρόμους, και δικοί του είναι μόνο οι παρακάτω ορισμοί·
//   · σε αληθινό επεξεργαστή ARM64 έτρεξε στις 14/11/2026 (Raspberry Pi 5, Cortex-A76, g++ 14.2): δίνει τα ίδια δείγματα μπιτ προς
//     μπιτ με την προηγούμενη διάταξη της ιστορίας (20.000 granules), και μέσα σε ολόκληρο τον αποκωδικοποιητή 32 αρχεία βγήκαν
//     ίδια σε κάθε πακέτο με τον αποκωδικοποιητή του FFmpeg. Η δοκιμή του ίδιου του σταδίου απέναντι στον ακριβή δρόμο δεν έχει
//     τρέξει ακόμη εκεί.
// Ο MSVC για ARM64 δεν παίρνει αυτόν τον δρόμο (όπως και στο flacdecoder.cpp): εκεί τρέχει η σκέτη C++.
SYN_INLINE float32x4_t revNeon(float32x4_t v)
{
    const float32x4_t h = vrev64q_f32(v);   // 1 0 3 2
    return vextq_f32(h, h, 2);              // 3 2 1 0
}
// Αναστροφή 4×4 όπως στο SSE2: zip1/zip2 = unpacklo/unpackhi, και μετά τα κάτω (ή τα πάνω) μισά των δύο ενώνονται.
#define SYN_T4(r0, r1, r2, r3, x0, x1, x2, x3)                                                                          \
    {                                                                                                                   \
        const float32x4_t t0 = vzip1q_f32(r0, r1), t1 = vzip2q_f32(r0, r1), t2 = vzip1q_f32(r2, r3), t3 = vzip2q_f32(r2, r3); \
        x0 = vcombine_f32(vget_low_f32(t0), vget_low_f32(t2));                                                          \
        x1 = vcombine_f32(vget_high_f32(t0), vget_high_f32(t2));                                                        \
        x2 = vcombine_f32(vget_low_f32(t1), vget_low_f32(t3));                                                          \
        x3 = vcombine_f32(vget_high_f32(t1), vget_high_f32(t3));                                                        \
    }
#define SYN_NAME(x) x##Neon
#define SYN_ATTR
#define SYN_LANES 4
#define SYN_V float32x4_t
#define SYN_ADD(a, b) vaddq_f32(a, b)
#define SYN_SUB(a, b) vsubq_f32(a, b)
#define SYN_MUL(a, b) vmulq_f32(a, b)
#define SYN_MADD(a, b, c) vfmaq_f32(c, a, b)
#define SYN_NMADD(a, b, c) vfmsq_f32(c, a, b)
#define SYN_ZERO() vdupq_n_f32(0.0f)
#define SYN_TW(p) vld1q_f32(p)
#define SYN_LOAD(p) vld1q_f32(p)
#define SYN_LOADU(p) vld1q_f32(p)
#define SYN_STORE(p, v) vst1q_f32(p, v)
#define SYN_STOREU(p, v) vst1q_f32(p, v)
#define SYN_REV(v) revNeon(v)
#define SYN_GATHER4(in, o, x0, x1, x2, x3)                                                                              \
    {                                                                                                                   \
        const float32x4_t g0 = vld1q_f32((in) + (o)), g1 = vld1q_f32((in) + 32 + (o));                                  \
        const float32x4_t g2 = vld1q_f32((in) + 64 + (o)), g3 = vld1q_f32((in) + 96 + (o));                             \
        SYN_T4(g0, g1, g2, g3, x0, x1, x2, x3)                                                                          \
    }
#define SYN_SCATTER4(row, o, q0, q1, q2, q3)                                                                            \
    {                                                                                                                   \
        float32x4_t y0, y1, y2, y3;                                                                                     \
        SYN_T4(q0, q1, q2, q3, y0, y1, y2, y3)                                                                          \
        vst1q_f32((row) + (o), y0);                                                                                     \
        vst1q_f32((row) + 32 + (o), y1);                                                                                \
        vst1q_f32((row) + 64 + (o), y2);                                                                                \
        vst1q_f32((row) + 96 + (o), y3);                                                                                \
    }
#define SYN_WANT_DCT 1
#define SYN_WANT_WINDOW 1
#include "mp3synth_kernel.inc"
#undef SYN_T4

void synthNeon(SynthState &st, const float *in, float *out)
{
    const SynthConsts &k = consts();
    float *cur = st.rows[st.last ^ 1];
    SYN_MID_BEGIN
    for (int s = 0; s < 16; s += 4) dctNeon(k, in + 32 * s, cur + 32 * s, mc + s);
    dctNeon(k, in + 32 * 14, cur + 32 * 14, mc + 14);
    windowNeon(k, st.rows[st.last], cur, mid, out);
    SYN_MID_END
    st.last ^= 1;
}
#endif   // SYN_NEON

// ── ΠΟΙΟΣ ΔΡΟΜΟΣ ──
using SynthFn = void (*)(SynthState &, const float *, float *);
bool gBasicOnly = false, gScalar = false;

struct Choice {
    SynthFn fn;
    const char *text;
};

Choice pick()
{
    if (gScalar) return {synthPlain, "plain c++"};
#if SYN_X86
    if (!gBasicOnly && cpuHasAvx2Fma()) return {synthAvx2, "avx2+fma"};
    return {synthSse2, "sse2"};
#elif SYN_NEON
    return {synthNeon, "neon"};
#else
    return {synthPlain, "plain c++"};
#endif
}

// Η πρώτη κλήση διαλέγει τον δρόμο και τον θυμάται. Ο δείκτης έχει από την αρχή έγκυρη τιμή (χωρίς κώδικα εκκίνησης), άρα το
// synth δουλεύει και αν το καλέσει κάποιος πριν τρέξουν οι αρχικοποιήσεις αυτού του αρχείου· και δύο νήματα που φτάνουν
// μαζί γράφουν την ίδια τιμή.
void synthFirst(SynthState &st, const float *in, float *out);
std::atomic<SynthFn> gFn{synthFirst};

void synthFirst(SynthState &st, const float *in, float *out)
{
    const SynthFn fn = pick().fn;
    gFn.store(fn, std::memory_order_relaxed);
    fn(st, in, out);
}

}  // namespace

void synthReset(SynthState &st)
{
    std::memset(st.rows, 0, sizeof st.rows);
    std::memset(st.mid, 0, sizeof st.mid);
    st.last = 0;
}

void synth(SynthState &st, const float *in, float *out) { gFn.load(std::memory_order_relaxed)(st, in, out); }

const char *synthKernelText() { return pick().text; }

void synthUseBasicOnly(bool on)
{
    gBasicOnly = on;
    gFn.store(pick().fn, std::memory_order_relaxed);
}

void synthUseScalar(bool on)
{
    gScalar = on;
    gFn.store(pick().fn, std::memory_order_relaxed);
}

}  // namespace orpheic::mp3
