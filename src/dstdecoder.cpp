// Ο ΑΠΟΣΥΜΠΙΕΣΤΗΣ DST, ο γρήγορος. Δες dstdecoder.h για το τι βγάζει, με τι συγκρίνεται και από τι
// γράφτηκε (το πρότυπο και ο προηγούμενος δικός μας αποσυμπιεστής).
//
// «Ο ΠΑΛΙΟΣ» σε όλα τα σχόλια αυτού του αρχείου είναι ο προηγούμενος δικός μας αποσυμπιεστής, που σήμερα μένει μόνο ως
// μέτρο σύγκρισης.
//
// ΑΠΟ ΤΙ ΓΡΑΦΤΗΚΕ. Μόνο από τον παλιό και από ό,τι ορίζει το πρότυπο όπως είναι ήδη υλοποιημένο εκεί
// (ISO/IEC 14496-3 Υποενότητα 10, 10.5 έως 10.7). Η ανάγνωση της επικεφαλίδας (τμήματα, αντιστοιχίσεις, φίλτρα, πίνακες
// πιθανοτήτων, μαζί με ΟΛΟΥΣ τους ελέγχους για χαλασμένα πλαίσια) είναι αντίγραφο εκείνου, με τις αλλαγές του διαιτητή
// (ενός απλού αποσυμπιεστή αναφοράς γραμμένου για τον έλεγχο) σημειωμένες εκεί που γίνονται. Αλλάζει η ανασύνθεση (10.7.2.3 και 10.7.2.4), με
// γενικές τεχνικές: πίνακες αναζήτησης, ιστορικό σε μπιτ, παράθυρο μπιτ στον αριθμητικό αποκωδικοποιητή, μέτρημα
// μηδενικών αντί για βρόχο. Για τον γρήγορο αποσυμπιεστή δεν διαβάστηκε κανένας άλλος αποσυμπιεστής DST.
//
// ΓΙΑΤΙ Η ΕΞΟΔΟΣ ΕΙΝΑΙ Η ΙΔΙΑ, ΜΠΙΤ ΠΡΟΣ ΜΠΙΤ (η απόδειξη είναι η σύγκριση της εξόδου, εδώ το επιχείρημα):
//
//  α. Η ΠΡΟΒΛΕΨΗ, ΟΚΤΩ ΜΠΙΤ ΤΗ ΦΟΡΑ. Ο παλιός υπολογίζει για κάθε μπιτ n: Z = Σ_{i<N} c_i · Y[n-1-i], με Y = ±1 και
//     Y[m] = (-1)^m για m < 0 (|Z| <= 128 * 256, χωρίς υπερχείλιση). Εδώ το ιστορικό R έχει στο μπιτ k το b[n-1-k] (1 όταν
//     Y = +1), και στην αρχή κάθε byte (n πολλαπλάσιο του 8) υπολογίζονται ΚΑΙ ΟΙ ΟΚΤΩ προβλέψεις του byte μαζί. Για το μπιτ
//     n+p (p = 0..7), χωρίζοντας το άθροισμα σε «παλιά» μπιτ (του R, y_k = ±1) και «νέα» (X_0..X_{p-1} του ίδιου byte):
//         Z_p = Σ_{k=0}^{N-1-p} c_{k+p} · y_k   +   Σ_{j<p} c_{p-1-j} · Y(X_j).
//     Το πρώτο άθροισμα είναι το άθροισμα ενός πίνακα φτιαγμένου από τους συντελεστές ΜΕΤΑΤΟΠΙΣΜΕΝΟΥΣ κατά p. Οι πίνακες
//     έχουν δείκτη ένα ΤΕΤΡΑΜΠΙΤΟ του ιστορικού (32 γραμμές των 4 συντελεστών), και κάθε θέση τους κρατά και τις οκτώ
//     μετατοπίσεις μαζί, σε οκτώ «λωρίδες» των 16 μπιτ: μία πρόσθεση 128 μπιτ ανά γραμμή δίνει και τις οκτώ προβλέψεις.
//     Το δεύτερο άθροισμα μπαίνει μπιτ προς μπιτ: κάθε λωρίδα ξεκινά σαν όλα τα νέα μπιτ να είναι 0 (Y = -1), και όταν βγει
//     ένα μπιτ X_j = 1 προστίθεται 2 c_{p-1-j} στις λωρίδες p > j (μία πρόσθεση: οι τιμές 2c σε μία λέξη, μία ανά λωρίδα).
//     Ακέραια αθροίσματα και τα δύο, με άλλη σειρά πρόσθεσης: ο ΙΔΙΟΣ αριθμός.
//     ΓΙΑΤΙ ΟΙ ΛΩΡΙΔΕΣ ΔΕΝ ΠΕΦΤΟΥΝ Η ΜΙΑ ΣΤΗΝ ΑΛΛΗ. Κάθε λωρίδα κρατά U = Z + S, με S = Σ|c_i|: κάθε όρος μπαίνει ως
//     c · (±1) + |c|, που είναι 0 ή 2|c|, ποτέ αρνητικός. Άρα κάθε λωρίδα, και κάθε ενδιάμεσο άθροισμά της, είναι στο
//     [0, 2S]. Με S <= 32767 χωρά σε 16 μπιτ χωρίς κρατούμενο προς τη διπλανή. Οι συντελεστές είναι στο [-256, 255], άρα
//     S <= 32768, και το 32768 το φτάνει ΜΟΝΟ το φίλτρο με 128 συντελεστές όλους -256: πλαίσιο με τέτοιο φίλτρο πηγαίνει από
//     τον αργό δρόμο (reconstructSlow: το άθροισμα όπως το ορίζει το πρότυπο, ένας όρος τη φορά). Υγιής δίσκος δεν το έχει.
//     ΟΙ ΠΡΟΣΘΕΣΕΙΣ ΤΩΝ 2c ΕΙΝΑΙ ΠΡΟΣΗΜΑΣΜΕΝΕΣ μέσα σε λέξη 64 μπιτ (τέσσερις λωρίδες): η λέξη, ως ακέραιος, είναι
//     Σ λωρίδα_t · 2^{16t}, και η πρόσθεση του Σ 2c_t · 2^{16t} δίνει Σ (λωρίδα_t + 2c_t) · 2^{16t}· αφού κάθε νέα λωρίδα
//     είναι πάλι στο [0, 65535], αυτό ΕΙΝΑΙ η σωστή λέξη (το modulo 2^64 δεν αλλάζει τίποτα, η τιμή είναι μέσα).
//     Όπου αλλάζει φίλτρο (όριο τμήματος) οι οκτώ προβλέψεις ξαναϋπολογίζονται ολόκληρες από το R με το νέο φίλτρο.
//  β. Αριθμητικός αποκωδικοποιητής. Ο παλιός: A στο [2048, 4095], P στο [1, 128] (όλοι οι δρόμοι το εξασφαλίζουν: 7 μπιτ + 1,
//     ή περιορισμός στο [1, 128], ή 128), K = (A>>8)|((A>>7)&1) στο [8, 15], Q = K·P στο [8, 1920], A - Q >= 1024. Αν στην
//     αρχή 0 <= C < A, τότε 0 <= C < A ΓΙΑ ΠΑΝΤΑ: αν C >= A-Q, νέο C = C-(A-Q) < Q = νέο A· αλλιώς C < A-Q = νέο A· και κάθε
//     διπλασιασμός δίνει 2C+b <= 2A-1 < 2A. Άρα το C χωρά σε 12 μπιτ και η λέξη V = C·2^51 + (τα επόμενα μπιτ της εισόδου)
//     δίνει: C >= D ακριβώς όταν V >= D·2^51, και η αφαίρεση και η ολίσθηση γίνονται πάνω στο V χωρίς να χαθεί μπιτ. Η
//     ολίσθηση είναι πόσες φορές θα έτρεχε ο βρόχος του παλιού, δηλαδή 11 - (θέση του υψηλότερου μπιτ του A), και το πολύ 8
//     (A >= 8)· γι' αυτό το παράθυρο ξαναγεμίζει όταν μείνουν λιγότερα από 8 μπιτ. Πέρα από το τέλος της εισόδου μπαίνουν
//     μηδενικά, όπως επιστρέφει ο αναγνώστης του παλιού.
//     ΑΝ ΤΑ ΠΡΩΤΑ 13 ΜΠΙΤ ΔΙΝΟΥΝ C >= 4095 (= το αρχικό A) η ανισότητα δεν ισχύει: το πλαίσιο ΑΠΟΡΡΙΠΤΕΤΑΙ (όπως στον
//     διαιτητή). Το πρότυπο θέλει το πρώτο μπιτ 0 και η τιμή ενός κώδικα είναι πάντα μέσα στο διάστημα, άρα κανένας
//     κωδικοποιητής δεν το γράφει· ο παλιός το κάνει ένα πλαίσιο θορύβου (και, εκτός αν είναι ακριβώς 4095 με μόνο
//     μηδενικά μετά, ολισθαίνει έναν αρνητικό int στη διαδρομή).
//  γ. Έξοδος. Ο παλιός μηδενίζει και μετά βάζει 1 << (7 - (n&7)) για κάθε X = 1. Εδώ τα 8 μπιτ του byte μαζεύονται με τη
//     σειρά (το πρώτο ψηλότερα) και γράφονται μαζί. Ο μηδενισμός στην αρχή μένει, για τους δρόμους που επιστρέφουν false.
//  δ. Τμήματα. Ο παλιός ρωτά το segOf σε κάθε byte. Το segOf συγκρίνει μόνο με τις αρχές των τμημάτων, άρα είναι σταθερό ανάμεσα
//     σε δύο διαδοχικές αρχές: εδώ ρωτιέται μία φορά στην αρχή κάθε τέτοιου διαστήματος (με την ΙΔΙΑ συνάρτηση), και το
//     διάστημα κόβεται και στο byte 16, ώστε ο έλεγχος «n < N» του Half_Probability να μένει μόνο στα πρώτα 128 μπιτ.
#include "dstdecoder.h"

#include <cstddef>
#include <cstring>
#include <utility>

// Οι οκτώ λωρίδες σε έναν καταχωρητή 128 μπιτ όπου υπάρχει (SSE2: κάθε x86 των 64 μπιτ· NEON: κάθε ARM των 64 μπιτ), αλλιώς
// σε δύο λέξεις των 64 (κάθε μηχάνημα των 32 μπιτ, όπως το armeabi-v7a, και κάθε big-endian: η διάταξη στη μνήμη, λωρίδα 0
// στα χαμηλά 16 μπιτ της πρώτης λέξης, ταιριάζει με τους καταχωρητές μόνο σε little-endian). Δύο διακόπτες ΜΟΝΟ για τις
// δοκιμές, ώστε να χτίζεται και να τρέχει εδώ και ο κώδικας που αλλιώς θα έβλεπε μόνο ένα άλλο μηχάνημα·
// ένα κανονικό χτίσιμο δεν ορίζει κανέναν:
//   DST_NO_SIMD   οι δύο λέξεις αντί για τον καταχωρητή 128 μπιτ (ό,τι θα χτίσει ένα ARM των 32 μπιτ).
//   DST_PORTABLE  το ίδιο, και επιπλέον ο βρόχος αντί για το μέτρημα μηδενικών (ό,τι θα χτίσει ένας άγνωστος μεταγλωττιστής).
#if defined(DST_PORTABLE) && !defined(DST_NO_SIMD)
#define DST_NO_SIMD 1
#endif
#if !defined(DST_NO_SIMD) && defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__)
#if __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
#define DST_NO_SIMD 1
#endif
#endif
// Το __SSE2__ ορίζεται και σε x86 των 32 μπιτ, όπου το _mm_cvtsi128_si64 δεν υπάρχει: μόνο στα 64.
#if !defined(DST_NO_SIMD) && ((defined(__SSE2__) && defined(__x86_64__)) || defined(_M_X64) || defined(_M_AMD64))
#define DST_SSE2 1
#include <emmintrin.h>
#elif !defined(DST_NO_SIMD) && (defined(__aarch64__) || defined(_M_ARM64))
#define DST_NEON 1
#include <arm_neon.h>
#endif
#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#endif

#if defined(__GNUC__) || defined(__clang__)
#define DST_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define DST_INLINE __forceinline
#else
#define DST_INLINE inline
#endif

namespace {
namespace dstd {   // δεύτερο όνομα: σε χτίσιμο «όλα σε ένα αρχείο» να μη συγκρουστεί με τα ίδια ονόματα του παλιού και του διαιτητή

constexpr int kMaxCh = 6;
constexpr int kMaxSeg = 8;                   // MAXNRSEGS (Ptable) = 8
constexpr int kMaxPredOrder = 128;           // 7 μπιτ + 1
constexpr int kMaxPtableLen = 64;            // 6 μπιτ + 1
constexpr int kMaxElem = kMaxCh * kMaxSeg;   // το πολύ ένα νέο στοιχείο ανά τμήμα: 48 φίλτρα, 48 πίνακες
constexpr int kRows = kMaxPredOrder / 4;     // 32 γραμμές των 4 συντελεστών
constexpr int kTabStride = kRows * 16;       // θέσεις (των 16 bytes) ανά φίλτρο: 8 KB
constexpr int kHeadBytes = kMaxPredOrder / 8;   // n < N μόνο στα πρώτα 16 bytes (128 μπιτ)
constexpr int kMaxAbsSum = 32767;            // Σ|c| ως εδώ χωρά σε λωρίδα 16 μπιτ (α πάνω)

using Entry = DstLanes;                      // δύο λέξεις: λωρίδες 0..3 και 4..7, η χαμηλότερη λωρίδα στα χαμηλά μπιτ

// ── ΑΝΤΙΓΡΑΦΟ ΤΟΥ ΠΑΛΙΟΥ: πίνακες του προτύπου, αναγνώστης μπιτ, τμήματα, αντιστοιχίσεις ──
const int kCcpc8[3][3] = {{-8, 0, 0}, {-16, 8, 0}, {-9, -5, 6}};
const int kPcpc[3][3] = {{-1, 0, 0}, {-2, 1, 0}, {-3, 3, -1}};
const int kCcpo[3] = {1, 2, 3};
const int kPcpo[3] = {1, 2, 3};

struct BitReader {
  const uint8_t* p;
  int64_t nbits;   // 64 μπιτ: ο long του παλιού είναι 32 μπιτ στον MSVC
  int64_t pos = 0;
  BitReader(const uint8_t* d, int bytes) : p(d), nbits(int64_t(bytes) * 8) {}
  bool overrun() const { return pos > nbits; }
  int bit() {
    if (pos >= nbits) { ++pos; return 0; }
    int b = (p[pos >> 3] >> (7 - (pos & 7))) & 1;
    ++pos;
    return b;
  }
  uint32_t bits(int k) {
    uint32_t v = 0;
    while (k-- > 0) v = (v << 1) | (uint32_t)bit();
    return v;
  }
};

int nbitsFor(int v) {
  int n = 0;
  while ((1 << n) <= v) ++n;
  return n;
}

struct Seg {
  int nr = 1;
  int start[kMaxSeg + 1] = {0};
  int elem[kMaxSeg] = {0};
};

void readChannelSeg(BitReader& br, Seg& s, int frameLen, int minSeg, int& res, bool& resRead) {
  s.nr = 1;
  s.start[0] = 0;
  int eoc = br.bit();
  while (eoc == 0 && s.nr < kMaxSeg) {
    if (!resRead) { res = (int)br.bits(13); resRead = true; }
    int range = (frameLen - s.start[s.nr - 1] - minSeg) / (res > 0 ? res : 1);
    int sl = (int)br.bits(nbitsFor(range));
    int segLen = res * sl;
    s.start[s.nr] = s.start[s.nr - 1] + segLen;
    ++s.nr;
    eoc = br.bit();
  }
  s.start[s.nr] = frameLen;
}

void readSegAlloc(BitReader& br, Seg seg[], int C, int frameLen, int minSeg) {
  int res = 1;
  bool resRead = false;
  int same = br.bit();
  if (same == 0) {
    for (int ch = 1; ch <= C; ++ch) readChannelSeg(br, seg[ch], frameLen, minSeg, res, resRead);
  } else {
    readChannelSeg(br, seg[1], frameLen, minSeg, res, resRead);
    for (int ch = 2; ch <= C; ++ch) seg[ch] = seg[1];
  }
}

void readChannelMap(BitReader& br, Seg& s, int& nrElem) {
  for (int sg = 0; sg < s.nr; ++sg) {
    int e = (int)br.bits(nbitsFor(nrElem));
    s.elem[sg] = e;
    if (e == nrElem) ++nrElem;
  }
}

int readMaps(BitReader& br, Seg seg[], int C) {
  int nrElem = 0;
  int same = br.bit();
  if (same == 0) {
    for (int ch = 1; ch <= C; ++ch) readChannelMap(br, seg[ch], nrElem);
  } else {
    readChannelMap(br, seg[1], nrElem);
    for (int ch = 2; ch <= C; ++ch)
      for (int sg = 0; sg < seg[ch].nr; ++sg) seg[ch].elem[sg] = seg[1].elem[sg];
  }
  return nrElem;
}

int segOf(const Seg& s, int n) {
  int byte = n >> 3;
  for (int sg = 0; sg < s.nr; ++sg)
    if (byte >= s.start[sg] && byte < s.start[sg + 1]) return sg;
  return s.nr - 1;
}
// ── τέλος του αντιγράφου ──

// ── ΟΙ ΛΩΡΙΔΕΣ ─────────────────────────────────────────────────────────────────────────────────────────────────────
#if defined(DST_SSE2)
using Lanes = __m128i;
DST_INLINE Lanes lanesZero() { return _mm_setzero_si128(); }
DST_INLINE Lanes lanesAdd(Lanes a, const Entry* e) {
  return _mm_add_epi16(a, _mm_load_si128(reinterpret_cast<const __m128i*>(e)));
}
DST_INLINE uint64_t lanesLo(Lanes a) { return uint64_t(_mm_cvtsi128_si64(a)); }
DST_INLINE uint64_t lanesHi(Lanes a) { return uint64_t(_mm_cvtsi128_si64(_mm_unpackhi_epi64(a, a))); }
#elif defined(DST_NEON)
using Lanes = uint16x8_t;
DST_INLINE Lanes lanesZero() { return vdupq_n_u16(0); }
DST_INLINE Lanes lanesAdd(Lanes a, const Entry* e) { return vaddq_u16(a, vreinterpretq_u16_u64(vld1q_u64(e->w))); }
DST_INLINE uint64_t lanesLo(Lanes a) { return vgetq_lane_u64(vreinterpretq_u64_u16(a), 0); }
DST_INLINE uint64_t lanesHi(Lanes a) { return vgetq_lane_u64(vreinterpretq_u64_u16(a), 1); }
#else
// Δύο λέξεις, τέσσερις λωρίδες η καθεμία. Η πρόσθεση ολόκληρης λέξης είναι πρόσθεση ανά λωρίδα επειδή καμία λωρίδα δεν
// περνά τα 16 μπιτ (α πάνω).
struct Lanes {
  uint64_t lo, hi;
};
DST_INLINE Lanes lanesZero() { return Lanes{0, 0}; }
DST_INLINE Lanes lanesAdd(Lanes a, const Entry* e) { return Lanes{a.lo + e->w[0], a.hi + e->w[1]}; }
DST_INLINE uint64_t lanesLo(Lanes a) { return a.lo; }
DST_INLINE uint64_t lanesHi(Lanes a) { return a.hi; }
#endif

// Πόσες φορές θα έτρεχε ο βρόχος κανονικοποίησης του παλιού (`while (A < 2048) A <<= 1`), για 8 <= A <= 4095: 11 μείον
// τη θέση του υψηλότερου μπιτ. Μία εντολή όπου υπάρχει (bsr/clz), αλλιώς ο ίδιος βρόχος.
DST_INLINE unsigned normShift(unsigned A) {
#if !defined(DST_PORTABLE) && (defined(__GNUC__) || defined(__clang__))
  // 31 ^ clz είναι η θέση του υψηλότερου μπιτ (για 32 μπιτ, 31 - x == 31 ^ x όταν 0 <= x <= 31): γραμμένο έτσι ο
  // μεταγλωττιστής βγάζει σκέτο bsr σε x86, μία πράξη λιγότερη πάνω στην αλυσίδα από A σε A από το «clz - 20».
  return 11u - (31u ^ unsigned(__builtin_clz(A)));
#elif !defined(DST_PORTABLE) && defined(_MSC_VER)
  unsigned long top;
  _BitScanReverse(&top, A);
  return 11u - unsigned(top);
#else
  unsigned k = 0;
  while ((A << k) < 2048u) ++k;
  return k;
#endif
}

// 8 bytes, το πρώτο ψηλότερα (όπως διαβάζει ο αναγνώστης μπιτ: MSB πρώτο). Χωρίς υπόθεση για τη σειρά bytes του
// μηχανήματος· οι μεταγλωττιστές το κάνουν ένα φόρτωμα και μία αντιστροφή.
DST_INLINE uint64_t loadBE64(const uint8_t* p) {
  return (uint64_t(p[0]) << 56) | (uint64_t(p[1]) << 48) | (uint64_t(p[2]) << 40) | (uint64_t(p[3]) << 32) |
         (uint64_t(p[4]) << 24) | (uint64_t(p[5]) << 16) | (uint64_t(p[6]) << 8) | uint64_t(p[7]);
}

// Ο αριθμητικός αποκωδικοποιητής, για 0 <= C < A (β πάνω· αλλιώς το πλαίσιο έχει ήδη απορριφθεί).
// V: μπιτ 63 πάντα 0, μπιτ 62..51 = C, από κάτω `avail` (<= 51) μπιτ της εισόδου που ακολουθούν, από κάτω μηδενικά.
// Το μπιτ 63 μένει άδειο ώστε το «V >= D·2^51» να είναι το πρόσημο μιας αφαίρεσης: αριθμητική, όχι σύγκριση που ο
// μεταγλωττιστής θα έκανε άλμα (το gcc το έκανε). Το E είναι σχεδόν τυχαίο, και ένα άλμα που μαντεύεται λάθος τη μισή
// φορά κοστίζει περισσότερο από όλη την υπόλοιπη αποκωδικοποίηση του μπιτ.
struct FastArith {
  uint64_t V;
  int avail;
  unsigned A;
  const uint8_t* in;
  int64_t inBytes;
  int64_t next;   // το επόμενο byte της εισόδου που δεν έχει μπει στο V

  // Γεμίζει με ολόκληρα bytes, όσο χωράνε κάτω από το bit 51. Καλείται με 0 <= avail < 8, άρα χωράνε 5 ή 6.
  // Κάθε byte διαβάζεται μόνο αν είναι μέσα στην είσοδο· πέρα από αυτήν μηδέν, όπως ο αναγνώστης του παλιού.
  DST_INLINE void refill() {
    if (next + 8 <= inBytes) {
      const int nb = (51 - avail) >> 3;
      const uint64_t w = loadBE64(in + next);
      V |= (w >> (64 - 8 * nb)) << (51 - avail - 8 * nb);
      avail += 8 * nb;
      next += nb;
    } else {
      while (avail <= 43) {
        const uint64_t b = next < inBytes ? uint64_t(in[next]) : uint64_t(0);
        ++next;
        V |= b << (43 - avail);
        avail += 8;
      }
    }
  }

  // Επιστρέφει 1 όταν E = 0 (C < A-Q) και 0 όταν E = 1: το αντίστροφο του E, γιατί αυτό βγαίνει κατευθείαν από το πρόσημο.
  DST_INLINE unsigned decode(unsigned P) {
    const unsigned K = (A >> 8) | ((A >> 7) & 1u);
    const unsigned Q = K * P;
    const unsigned D = A - Q;
    const uint64_t t = uint64_t(D) << 51;
    const uint64_t diff = V - t;                 // V, t < 2^63: το μπιτ 63 της διαφοράς είναι 1 ακριβώς όταν V < t
    const uint64_t lt = diff >> 63;              // 1 = E 0 (C < A-Q), 0 = E 1
    V = diff + (t & (uint64_t(0) - lt));         // E 1: C -= A-Q· E 0: όπως ήταν
    A = Q + ((D - Q) & (0u - unsigned(lt)));     // E 1: A = Q· E 0: A = A-Q
    const unsigned k = normShift(A);
    A <<= k;
    V <<= k;
    avail -= int(k);
    if (avail < 8) refill();
    return unsigned(lt);
  }
};

// Ό,τι χρειάζεται η κρίσιμη διαδρομή από ένα φίλτρο, εκτός από τον πίνακά του.
struct Filt {
  uint64_t p4[5];   // p4[k]: οι τιμές 2c_k, 2c_{k+1}, 2c_{k+2}, 2c_{k+3} σε τέσσερις λωρίδες, προσημασμένες
  int S;            // Σ|c_i|
  int N;            // τάξη
  int rows;         // γραμμές που έχουν συντελεστές: (N + 3) / 4
  int steps;        // πόσα από τα 8 βήματα ενός byte έχουν αναζητήσεις «υπολοίπου» (δες lookups)
};

// Ένα κανάλι.
struct Chan {
  Lanes acc;            // τα αθροίσματα για το ΕΠΟΜΕΝΟ byte, όσα έχουν μαζευτεί ως τώρα
  uint64_t rlo, rhi;    // το ιστορικό στην ΑΡΧΗ του τρέχοντος byte: μπιτ k του (rhi:rlo) = b[n-1-k], 1 όταν Y = +1
  uint64_t wlo, whi;    // οι προβλέψεις του τρέχοντος byte, U = Z + S, λωρίδες 0..3 και 4..7 (δες step)
  const Entry* T;       // ο πίνακας του φίλτρου
  const uint8_t* pt;    // ο πίνακας πιθανοτήτων (1..128), συμπληρωμένος ως τις 64 θέσεις με την τελευταία του τιμή
  uint64_t p4[5];       // του φίλτρου
  int S;                // του φίλτρου
  int halfUntil;        // Half_Probability: P = 128 για n < halfUntil (η τάξη του φίλτρου, ή 0)
  int steps;            // του φίλτρου
  unsigned r;           // τα μπιτ του τρέχοντος byte ως τώρα, το πρώτο ψηλότερα
};

// Η θέση του πίνακα για το τετράμπιτο K του ιστορικού (K = 0 το πιο πρόσφατο), στη γραμμή Row. Η μετατόπιση μέσα στη γραμμή
// βγαίνει κατευθείαν σε bytes (τετράμπιτο επί 16): μία ολίσθηση και μία μάσκα.
template <int K, int Row>
DST_INLINE const Entry* cell(const Chan& c) {
  const uint64_t w = K < 16 ? c.rlo : c.rhi;
  size_t off;
  if constexpr ((K & 15) == 0) off = size_t((w << 4) & 0xF0u);
  else off = size_t((w >> (4 * (K & 15) - 4)) & 0xF0u);
  return reinterpret_cast<const Entry*>(reinterpret_cast<const char*>(c.T) + size_t(Row) * 256 + off);
}

// ΤΑ ΑΘΡΟΙΣΜΑΤΑ ΓΙΑ ΤΟ ΕΠΟΜΕΝΟ BYTE, μοιρασμένα στα 8 βήματα του τρέχοντος ώστε ο επεξεργαστής να τα κάνει παράλληλα με
// τον αριθμητικό αποκωδικοποιητή. Το ιστορικό του επόμενου byte είναι (R << 8) | (τα 8 νέα μπιτ): το τετράμπιτο K του R
// γίνεται το τετράμπιτο K+2, άρα η γραμμή K+2 με δείκτη το τετράμπιτο K του R (γνωστό από την αρχή του byte)· οι γραμμές 1
// και 0 με δείκτη τα νέα μπιτ, στο τελευταίο βήμα. Το βήμα J παίρνει τις γραμμές 4J+2 .. 4J+5 (το 7: τις 30 και 31), και
// μόνο αν το φίλτρο φτάνει ως εκεί (c.steps): ένα φίλτρο μικρής τάξης δεν πληρώνει για γραμμές που είναι μηδέν.
// ΔΥΟ ΑΛΛΕΣ ΔΙΑΤΑΞΕΙΣ ΔΟΚΙΜΑΣΤΗΚΑΝ ΚΑΙ ΑΠΟΡΡΙΦΘΗΚΑΝ ΜΕ ΜΕΤΡΗΣΗ (12/11), για να μην ξαναδοκιμαστούν: όλες οι αναζητήσεις
// μαζεμένες στο τελευταίο μπιτ (ο αθροιστής σε καταχωρητή αντί για τη μνήμη) βγήκαν 7 έως 16% πιο αργές σε σύγχρονο
// επεξεργαστή· ανά δύο μπιτ, 1 έως 5% πιο αργές εκεί και 4 έως 5% πιο αργές σε Core 2 με χτίσιμο MSVC (όπου περίμενα κέρδος,
// γιατί εκείνος ο μεταγλωττιστής γράφει τον αθροιστή με τη μη στοιχισμένη εντολή). Η δουλειά απλωμένη ομοιόμορφα κερδίζει.
template <int J>
DST_INLINE void lookups(Chan& c) {
  if constexpr (J == 0) c.acc = lanesZero();
  if (J < c.steps) {
    Lanes a = c.acc;
    a = lanesAdd(a, cell<4 * J, 4 * J + 2>(c));
    a = lanesAdd(a, cell<4 * J + 1, 4 * J + 3>(c));
    if constexpr (J < 7) {
      a = lanesAdd(a, cell<4 * J + 2, 4 * J + 4>(c));
      a = lanesAdd(a, cell<4 * J + 3, 4 * J + 5>(c));
    }
    c.acc = a;
  }
  if constexpr (J == 7) {
    // Τα νέα μπιτ: X0..X3 (γραμμή 1) και X4..X6 με το X7 ως 0 (γραμμή 0· το X7 δεν έχει βγει ακόμη, το μερίδιό του μπαίνει
    // στο step). Το c.r έχει εδώ τα 7 πρώτα μπιτ.
    Lanes a = c.acc;
    a = lanesAdd(a, c.T + 16 + ((c.r >> 3) & 15u));
    a = lanesAdd(a, c.T + ((c.r << 1) & 15u));
    c.acc = a;
  }
}

// Ένα μπιτ ενός καναλιού, το J-οστό του byte: πρόβλεψη -> πιθανότητα -> αριθμητικός αποκωδικοποιητής -> μπιτ DSD.
// ΠΩΣ ΚΑΤΑΝΑΛΩΝΟΝΤΑΙ ΟΙ ΛΩΡΙΔΕΣ. Τα μπιτ 0..3 διαβάζουν τη χαμηλή λωρίδα του wlo, που μετά κατεβαίνει κατά 16· τα 4..7 το
// ίδιο με το whi. Όταν βγει το μπιτ J με τιμή 1, οι λωρίδες p > J παίρνουν 2 c_{p-1-J}: στη λέξη που κατεβαίνει αυτό είναι
// πάντα το p4[0] (2c_0 στη νέα χαμηλή λωρίδα, 2c_1 στην επόμενη, ...)· στο whi, όσο ακόμη καταναλώνεται το wlo, είναι οι
// λωρίδες 4..7, δηλαδή 2c_{3-J} .. 2c_{6-J} = p4[3-J]. Ό,τι περισσεύει πάνω από τις λωρίδες που θα διαβαστούν είναι
// σκουπίδι που δεν διαβάζεται ποτέ, και κρατούμενο πηγαίνει μόνο προς τα πάνω.
template <int J, bool Head, class Arith>
DST_INLINE void step(Chan& c, Arith& ad, int n) {
  const uint64_t w = J < 4 ? c.wlo : c.whi;
  const int Z = int(unsigned(w) & 0xFFFFu) - c.S;
  const unsigned neg = unsigned(Z) >> 31;          // 1 όταν Z < 0, δηλαδή F = 0
  // Ο παλιός: pone[min(|Z| >> 3, μήκος - 1)]. Εδώ min(.., 63) σε πίνακα που από το μήκος και πέρα επαναλαμβάνει
  // την τελευταία τιμή: ίδια τιμή, και ο περιορισμός είναι σταθερά αντί για ακόμη μία μεταβλητή ανά κανάλι.
  unsigned idx = unsigned(Z < 0 ? -Z : Z) >> 3;
  if (idx > unsigned(kMaxPtableLen - 1)) idx = unsigned(kMaxPtableLen - 1);
  unsigned P = c.pt[idx];
  if constexpr (Head) {
    if (n < c.halfUntil) P = 128;
  } else {
    (void)n;
  }
  lookups<J>(c);
  // X = E xor F, με E = 1 - lt και F = 1 - neg: τα δύο «1 -» αλληλοαναιρούνται.
  const unsigned X = ad.decode(P) ^ neg;
  const uint64_t m = uint64_t(0) - uint64_t(X);
  if constexpr (J < 3) {
    c.wlo = (w >> 16) + (c.p4[0] & m);
    c.whi += c.p4[3 - J] & m;
  } else if constexpr (J == 3) {
    c.whi += c.p4[0] & m;
  } else if constexpr (J < 7) {
    c.whi = (w >> 16) + (c.p4[0] & m);
  }
  if constexpr (J == 0) {
    c.r = X;
  } else if constexpr (J < 7) {
    c.r = (c.r << 1) | X;
  } else {
    // Το byte τελείωσε. Οι προβλέψεις του επόμενου: τα αθροίσματα που μαζεύτηκαν (με το X7 ως 0) και το μερίδιο του X7,
    // που είναι 2c_p σε κάθε λωρίδα p. Και το ιστορικό προχωρά κατά 8.
    const uint64_t byte = uint64_t((c.r << 1) | X);
    c.wlo = lanesLo(c.acc) + (c.p4[0] & m);
    c.whi = lanesHi(c.acc) + (c.p4[4] & m);
    c.rhi = (c.rhi << 8) | (c.rlo >> 56);
    c.rlo = (c.rlo << 8) | byte;
  }
}

// Οι οκτώ προβλέψεις από την αρχή, από όλο το ιστορικό: στην αρχή κάθε διαστήματος (το φίλτρο μπορεί να άλλαξε).
void startStretch(Chan& c, int rows) {
  Lanes a = lanesZero();
  for (int g = 0; g < rows; ++g) {
    const uint64_t w = g < 16 ? c.rlo : c.rhi;
    a = lanesAdd(a, c.T + g * 16 + int((w >> (4 * (g & 15))) & 15u));
  }
  c.wlo = lanesLo(a);
  c.whi = lanesHi(a);
}

// Οι πίνακες ενός φίλτρου και οι σταθερές του. Γραμμή g, θέση x (ένα τετράμπιτο του ιστορικού), λωρίδα p:
//   Σ_{t<4} ( c_{4g+t+p} · (μπιτ t του x ? +1 : -1) + |c_{4g+t+p}| ),   με c_i = 0 για i >= N,
// και στη γραμμή 0 επιπλέον Σ_{i<p} (|c_i| - c_i): το μερίδιο των νέων μπιτ του byte σαν να είναι όλα 0 (α στην κορυφή).
// Κάθε τιμή είναι άθροισμα μη αρνητικών όρων, το πολύ 4 * 512 + 7 * 512.
void buildFilter(Entry* tab, const int* cf, int N, Filt& q) {
  int c[kMaxPredOrder + 16];
  int S = 0;
  for (int i = 0; i < kMaxPredOrder + 16; ++i) {
    c[i] = i < N ? cf[i] : 0;
    S += c[i] < 0 ? -c[i] : c[i];
  }
  const int rows = (N + 3) >> 2;
  const int steps = rows <= 2 ? 0 : (rows + 1) >> 2;   // ceil((rows - 2) / 4), το πολύ 8
  q.S = S;
  q.N = N;
  q.rows = rows;
  q.steps = steps;
  for (int k = 0; k < 5; ++k) {
    uint64_t v = 0;
    for (int t = 0; t < 4; ++t) v += uint64_t(int64_t(2 * c[k + t])) << (16 * t);   // προσημασμένο, modulo 2^64
    q.p4[k] = v;
  }
  if (S > kMaxAbsSum) return;   // αργός δρόμος: οι πίνακες δεν θα χρησιμοποιηθούν
  int k0[8];
  k0[0] = 0;
  for (int p = 1; p < 8; ++p) k0[p] = k0[p - 1] + (c[p - 1] < 0 ? -2 * c[p - 1] : 0);   // |c| - c: 0 ή 2|c|
  for (int g = 0; g < rows; ++g) {
    int e[16][8];
    for (int p = 0; p < 8; ++p) {
      int v = g == 0 ? k0[p] : 0;
      for (int t = 0; t < 4; ++t) {
        const int x = c[4 * g + t + p];
        v += x < 0 ? -2 * x : 0;
      }
      e[0][p] = v;
    }
    for (int t = 0; t < 4; ++t) {
      const int half = 1 << t;
      for (int x = 0; x < half; ++x)
        for (int p = 0; p < 8; ++p) e[half + x][p] = e[x][p] + 2 * c[4 * g + t + p];
    }
    Entry* row = tab + g * 16;
    for (int x = 0; x < 16; ++x) {
      uint64_t lo = 0, hi = 0;
      for (int p = 0; p < 4; ++p) {
        lo |= uint64_t(uint32_t(e[x][p])) << (16 * p);
        hi |= uint64_t(uint32_t(e[x][p + 4])) << (16 * p);
      }
      row[x].w[0] = lo;
      row[x].w[1] = hi;
    }
  }
  // Γραμμές που διαβάζονται χωρίς να έχουν συντελεστές: η 1 (πάντα), και όσες πιάνει το τελευταίο βήμα με αναζητήσεις.
  const int pad = steps == 8 ? kRows : (steps > 0 ? 4 * steps + 2 : 2);
  for (int g = rows; g < pad; ++g) std::memset(static_cast<void*>(tab + g * 16), 0, 16 * sizeof(Entry));
}

// Τα 8 μπιτ ενός byte για όλα τα κανάλια, με τη σειρά μπιτ 0: κανάλι 0, 1, ..· μπιτ 1: κανάλι 0, 1, .. (η σειρά μετρά: ένας
// αριθμητικός αποκωδικοποιητής για όλα).
template <bool Head, class Arith, std::size_t... J>
DST_INLINE void byteOfTwo(Chan& a, Chan& b, Arith& ad, int n, std::index_sequence<J...>) {
  ((step<int(J), Head>(a, ad, n + int(J)), step<int(J), Head>(b, ad, n + int(J))), ...);
}
template <int J, bool Head, class Arith>
DST_INLINE void stepAll(Chan* s, int C, Arith& ad, int n) {
  for (int ch = 0; ch < C; ++ch) step<J, Head>(s[ch], ad, n);
}
template <bool Head, class Arith, std::size_t... J>
DST_INLINE void byteOfAll(Chan* s, int C, Arith& ad, int n, std::index_sequence<J...>) {
  (stepAll<int(J), Head>(s, C, ad, n + int(J)), ...);
}

// Bytes [b0, b1) όπου κανένα κανάλι δεν αλλάζει φίλτρο ή πίνακα πιθανοτήτων. CC == 2: τα δύο κανάλια ξεδιπλωμένα, η
// κατάστασή τους σε τοπικές μεταβλητές (σχεδόν όλοι οι δίσκοι έχουν στερεοφωνική περιοχή)· CC == 0: όσα λέει το Crt, με
// την κατάσταση στη μνήμη. Head: τα πρώτα 16 bytes, όπου ισχύει ο έλεγχος «n < N» του Half_Probability.
template <int CC, bool Head, class Arith>
void runBytes(Arith& arith, Chan* chans, int Crt, uint8_t* out, int b0, int b1) {
  // Αντίγραφο σε τοπική μεταβλητή: η έξοδος γράφεται μέσα από uint8_t*, που κατά τη γλώσσα μπορεί να δείχνει οπουδήποτε·
  // ό,τι ζει σε τοπική μεταβλητή που δεν δόθηκε σε κανέναν, ο μεταγλωττιστής ξέρει ότι δεν το αγγίζει.
  Arith ad = arith;
  using Seq = std::make_index_sequence<8>;
  if constexpr (CC == 2) {
    (void)Crt;
    Chan a = chans[0], b = chans[1];
    for (int by = b0; by < b1; ++by) {
      byteOfTwo<Head>(a, b, ad, by * 8, Seq{});
      out[size_t(by) * 2] = uint8_t(a.rlo);
      out[size_t(by) * 2 + 1] = uint8_t(b.rlo);
    }
    chans[0] = a;
    chans[1] = b;
  } else {
    const int C = Crt;
    for (int by = b0; by < b1; ++by) {
      byteOfAll<Head>(chans, C, ad, by * 8, Seq{});
      uint8_t* o = out + size_t(by) * size_t(C);
      for (int ch = 0; ch < C; ++ch) o[ch] = uint8_t(chans[ch].rlo);
    }
  }
  arith = ad;
}

struct Frame {
  int C;
  int frameLen;
  int nrFilters;
  const Seg* fseg;
  const Seg* pseg;
  const int* predOrder;
  const int* coef;        // kMaxPredOrder ανά φίλτρο
  const Filt* filt;       // ανά φίλτρο
  const Entry* tab;       // kTabStride ανά φίλτρο
  const uint8_t* p8;      // kMaxPtableLen ανά πίνακα
  const int* halfProb;    // [1..C]
  const int* brk;         // όρια διαστημάτων σε bytes, αύξοντα, από 0 ως frameLen
  int nbrk;
};

template <class Arith>
void reconstruct(Arith& ad, const Frame& f, uint8_t* out) {
  const int C = f.C;
  Chan s[kMaxCh];
  for (int ch = 0; ch < C; ++ch) {
    s[ch].rlo = s[ch].rhi = 0xAAAAAAAAAAAAAAAAull;   // Y[m] = (-1)^m για m < 0
    s[ch].acc = lanesZero();
    s[ch].r = 0;
  }
  for (int r = 0; r + 1 < f.nbrk; ++r) {
    const int b0 = f.brk[r], b1 = f.brk[r + 1];
    for (int ch = 0; ch < C; ++ch) {
      // Η ΙΔΙΑ ερώτηση με τον παλιό (segOf στο πρώτο μπιτ του byte), μία φορά για όλο το διάστημα (δ στην κορυφή).
      const int fs = f.fseg[ch + 1].elem[segOf(f.fseg[ch + 1], b0 * 8)];
      const int ps = f.pseg[ch + 1].elem[segOf(f.pseg[ch + 1], b0 * 8)];
      const Filt& q = f.filt[fs];
      Chan& c = s[ch];
      c.T = f.tab + size_t(fs) * kTabStride;
      c.pt = f.p8 + size_t(ps) * kMaxPtableLen;
      for (int k = 0; k < 5; ++k) c.p4[k] = q.p4[k];
      c.S = q.S;
      c.halfUntil = f.halfProb[ch + 1] == 1 ? q.N : 0;
      c.steps = q.steps;
      startStretch(c, q.rows);   // ολόκληρο: το φίλτρο μπορεί να άλλαξε
    }
    if (b0 < kHeadBytes) runBytes<0, true>(ad, s, C, out, b0, b1);
    else if (C == 2) runBytes<2, false>(ad, s, C, out, b0, b1);
    else runBytes<0, false>(ad, s, C, out, b0, b1);
  }
}

// Ο ΑΡΓΟΣ ΔΡΟΜΟΣ: πλαίσιο με φίλτρο που το Σ|c| του δεν χωρά σε λωρίδα (128 συντελεστές όλοι -256· υγιής δίσκος δεν το
// έχει). Η πρόβλεψη όπως την ορίζει το πρότυπο, ένας όρος τη φορά πάνω στο ίδιο ιστορικό των 128 μπιτ, και το φίλτρο και
// ο πίνακας ρωτιούνται σε κάθε byte όπως στον παλιό. Τίποτα εδώ δεν χρειάζεται να είναι γρήγορο.
template <class Arith>
void reconstructSlow(Arith& ad, const Frame& f, uint8_t* out) {
  const int C = f.C;
  uint64_t lo[kMaxCh], hi[kMaxCh];
  for (int ch = 0; ch < C; ++ch) lo[ch] = hi[ch] = 0xAAAAAAAAAAAAAAAAull;
  for (int b = 0; b < f.frameLen; ++b) {
    const int* cf[kMaxCh];
    const uint8_t* pt[kMaxCh];
    int N[kMaxCh];
    for (int ch = 0; ch < C; ++ch) {
      const int fs = f.fseg[ch + 1].elem[segOf(f.fseg[ch + 1], b * 8)];
      const int ps = f.pseg[ch + 1].elem[segOf(f.pseg[ch + 1], b * 8)];
      cf[ch] = f.coef + size_t(fs) * kMaxPredOrder;
      pt[ch] = f.p8 + size_t(ps) * kMaxPtableLen;
      N[ch] = f.predOrder[fs];
    }
    for (int bit = 0; bit < 8; ++bit) {
      const int n = b * 8 + bit;
      for (int ch = 0; ch < C; ++ch) {
        int Z = 0;
        for (int i = 0; i < N[ch]; ++i) {
          const uint64_t w = i < 64 ? lo[ch] : hi[ch];
          Z += ((w >> (i & 63)) & 1u) ? cf[ch][i] : -cf[ch][i];
        }
        const unsigned neg = Z < 0 ? 1u : 0u;
        unsigned idx = unsigned(Z < 0 ? -Z : Z) >> 3;
        if (idx > unsigned(kMaxPtableLen - 1)) idx = unsigned(kMaxPtableLen - 1);
        unsigned P = pt[ch][idx];
        if (f.halfProb[ch + 1] == 1 && n < N[ch]) P = 128;
        const unsigned X = ad.decode(P) ^ neg;
        hi[ch] = (hi[ch] << 1) | (lo[ch] >> 63);
        lo[ch] = (lo[ch] << 1) | uint64_t(X);
      }
    }
    for (int ch = 0; ch < C; ++ch) out[size_t(b) * size_t(C) + size_t(ch)] = uint8_t(lo[ch]);
  }
}

// ΜΟΝΟ ΓΙΑ ΤΗ ΜΕΤΡΗΣΗ: τι είδους πλαίσιο ήταν αυτό. Καλείται μία φορά ανά συμπιεσμένο πλαίσιο, πριν από την ανασύνθεση,
// και μόνο όταν το εργαλείο έδωσε χώρο για μετρητές· η κρίσιμη διαδρομή δεν έχει ούτε μία εντολή από εδώ.
void noteFrame(DstStats& st, const Frame& f, bool slow) {
  ++st.coded;
  if (slow) ++st.slow;
  bool multi = false, anyHalf = false, narrow = true;
  for (int ch = 1; ch <= f.C; ++ch) {
    if (f.fseg[ch].nr > 1 || f.pseg[ch].nr > 1) multi = true;
    if (f.halfProb[ch] == 1) anyHalf = true;
  }
  if (multi) ++st.multiSeg;
  if (anyHalf) ++st.halfProb;
  if (f.nrFilters > st.maxFilters) st.maxFilters = f.nrFilters;
  // Τα όρια είναι αύξοντα και το 16 είναι πάντα ένα από αυτά όταν το πλαίσιο είναι μακρύτερο: όριο μέσα στο (0, 16) = κάποιο
  // τμήμα αρχίζει μέσα στα πρώτα 16 bytes.
  if (f.nbrk > 1 && f.brk[1] < kHeadBytes && f.brk[1] < f.frameLen) ++st.headCut;
  st.stretches += f.nbrk > 0 ? f.nbrk - 1 : 0;
  for (int r = 0; r + 1 < f.nbrk; ++r) {
    const int b0 = f.brk[r];
    const long long len = f.brk[r + 1] - b0;
    int used[kMaxCh];
    int nUsed = 0, nLong = 0;
    for (int ch = 0; ch < f.C; ++ch) {
      const int fs = f.fseg[ch + 1].elem[segOf(f.fseg[ch + 1], b0 * 8)];
      const int N = f.predOrder[fs];
      st.orderBytes[(N - 1) >> 4] += len;
      if (N > 64) narrow = false;
      bool seen = false;
      for (int i = 0; i < nUsed; ++i) seen = seen || used[i] == fs;
      if (!seen) {
        used[nUsed++] = fs;
        if (N > 64) ++nLong;
      }
    }
    st.stretchBytes += len;
    st.distinctBytes[nUsed] += len;
    st.longBytes[nLong] += len;
  }
  if (narrow) ++st.narrow;
}

// Το σώμα του DstDecoder::decodeFrame, με τον πρόχειρο χώρο του ως ορίσματα. Μέσα στο dstd και όχι με «using namespace dstd»
// μέσα στη μέθοδο: εκείνο, σε χτίσιμο «όλα σε ένα αρχείο» μαζί με τον παλιό ή τον διαιτητή, θα έκανε διφορούμενα τα ίδια
// ονόματα (BitReader, Seg, kMaxCh, ...) των αρχείων. Εδώ η αναζήτηση βρίσκει πρώτα τα ονόματα του dstd.
bool decode(int channels_, DstRefusal& why, DstStats* st, std::vector<int>& coef_, std::vector<int>& pone_,
            std::vector<uint8_t>& p8_, std::vector<DstLanes>& tab_, const uint8_t* in, int inBytes, uint8_t* dsdOut,
            int bytesPerChannel) {
  why = DstRefusal::Shape;
  const int C = channels_;
  if (C < 1 || C > kMaxCh) return false;
  const int frameLen = bytesPerChannel;
  // Αρνητικό μήκος θα μηδένιζε «αρνητικά» bytes, και από 2^28 bytes ανά κανάλι το «byte * 8» δεν χωρά σε int. Κανείς δεν
  // καλεί έτσι (ένας αναγνώστης SACD δίνει πάντα 4704)· άρνηση, ώστε η ασφάλεια να μη στηρίζεται στους καλούντες.
  if (frameLen < 0 || frameLen > kDstMaxBytesPerChannel) return false;
  if (inBytes < 0) inBytes = 0;
  std::memset(dsdOut, 0, size_t(frameLen) * size_t(C));

  BitReader br(in, inBytes);

  const int processingMode = br.bit();          // Table 10.4
  if (processingMode == 0) {                     // plain DSD
    // 1 + 1 + 6 μπιτ, και μετά frameLen * C bytes DSD, με τη σειρά: τα bytes 1.. της εισόδου ως έχουν. Αν είναι λιγότερα,
    // το πλαίσιο κόπηκε: άρνηση, όχι το κομμάτι που έφτασε και μετά μηδενικά (το 0x00 στο DSD δεν είναι σιωπή αλλά άλμα
    // συνεχούς τάσης σε πλήρη κλίμακα).
    const int64_t want = int64_t(frameLen) * C;
    if (inBytes < 1 || int64_t(inBytes) - 1 < want) { why = DstRefusal::ShortPlain; return false; }
    if (want > 0) std::memcpy(dsdOut, in + 1, size_t(want));
    why = DstRefusal::None;
    if (st) ++st->plain;
    return true;
  }

  // ---------- DST coded: ΑΝΤΙΓΡΑΦΟ του παλιού, με σταθερούς πίνακες αντί για vectors ----------
  // Κάθε `return false` ως το τέλος του αντιγράφου είναι έλεγχος που έχει και ο παλιός: ο λόγος είναι «επικεφαλίδα».
  why = DstRefusal::Header;
  Seg fseg[kMaxCh + 1], pseg[kMaxCh + 1];

  const int sameSeg = br.bit();
  if (sameSeg == 0) {
    readSegAlloc(br, fseg, C, frameLen, 128);
    readSegAlloc(br, pseg, C, frameLen, 4);
  } else {
    readSegAlloc(br, fseg, C, frameLen, 128);
    for (int ch = 1; ch <= C; ++ch) pseg[ch] = fseg[ch];
  }

  int nrFilters, nrPtables;
  const int sameMap = br.bit();
  if (sameMap == 0) {
    nrFilters = readMaps(br, fseg, C);
    nrPtables = readMaps(br, pseg, C);
  } else {
    nrFilters = readMaps(br, fseg, C);
    for (int ch = 1; ch <= C; ++ch)
      for (int sg = 0; sg < fseg[ch].nr; ++sg) pseg[ch].elem[sg] = fseg[ch].elem[sg];
    nrPtables = nrFilters;
  }

  // Ο έλεγχος του παλιού: κάθε τμήμα να δείχνει στοιχείο που δηλώθηκε.
  for (int ch = 1; ch <= C; ++ch) {
    for (int sg = 0; sg < fseg[ch].nr; ++sg)
      if (fseg[ch].elem[sg] < 0 || fseg[ch].elem[sg] >= nrFilters) return false;
    for (int sg = 0; sg < pseg[ch].nr; ++sg)
      if (pseg[ch].elem[sg] < 0 || pseg[ch].elem[sg] >= nrPtables) return false;
  }
  // Οι σταθεροί πίνακες εδώ χωρούν kMaxElem στοιχεία. Το readChannelMap αυξάνει το πλήθος το πολύ μία φορά ανά τμήμα, άρα
  // το πλήθος δεν περνά ποτέ τα C * 8 <= 48 και αυτός ο έλεγχος δεν πιάνει ποτέ τίποτα (ίδια απάντηση με τον παλιό). Μένει
  // για να μην εξαρτάται η ασφάλεια της μνήμης από εκείνο το επιχείρημα.
  if (nrFilters < 1 || nrFilters > kMaxElem || nrPtables < 1 || nrPtables > kMaxElem) return false;

  int halfProb[kMaxCh + 1] = {0};
  for (int ch = 1; ch <= C; ++ch) halfProb[ch] = br.bit();

  // Filter_Coef_Sets (Table 10.13).
  if ((int)coef_.size() < kMaxElem * kMaxPredOrder) coef_.assign(size_t(kMaxElem) * kMaxPredOrder, 0);
  int predOrder[kMaxElem];
  for (int f = 0; f < nrFilters; ++f) {
    int* cf = coef_.data() + size_t(f) * kMaxPredOrder;
    int po = (int)br.bits(7) + 1;
    predOrder[f] = po;
    for (int i = 0; i < po; ++i) cf[i] = 0;
    int coded = br.bit();
    if (coded == 0) {
      for (int i = 0; i < po; ++i) {
        int v = (int)br.bits(9);
        cf[i] = (v & 0x100) ? v - 512 : v;
      }
    } else {
      int method = (int)br.bits(2);
      if (method > 2) return false;              // δύο μπιτ, τρεις μέθοδοι (όπως ο παλιός)
      int ccpo = kCcpo[method];
      if (po < ccpo) return false;               // ccpo θέσεις σε πίνακα po (όπως ο παλιός)
      for (int i = 0; i < ccpo; ++i) {
        int v = (int)br.bits(9);
        cf[i] = (v & 0x100) ? v - 512 : v;
      }
      int ccm = (int)br.bits(3);
      for (int i = ccpo; i < po; ++i) {
        int run = 0;
        while (br.bit() == 0) { if (++run > 1024 || br.overrun()) return false; }
        int lsbs = (int)br.bits(ccm);
        int delta = (run << ccm) + lsbs;
        if (delta != 0 && br.bit() == 1) delta = -delta;
        // ΟΧΙ ΑΝΤΙΓΡΑΦΟ ΤΟΥ ΠΑΛΙΟΥ. Ο παλιός αθροίζει το d8 σε int, πάνω σε συντελεστές που
        // ακόμη δεν έχουν περιοριστεί· με τη μέθοδο 1 είναι c = δ + 2c' - c'', που μεγαλώνει με το τετράγωνο του i, και ένα
        // ΦΤΙΑΓΜΕΝΟ πλαίσιο (τάξη 128, παράμετρος Rice 7, μεγάλα υπόλοιπα) περνά το όριο του int: απροσδιόριστη συμπεριφορά,
        // και κάθε μεταγλωττιστής βγάζει άλλο αποτέλεσμα. Εδώ η ΙΔΙΑ πράξη σε 64 μπιτ, ακριβής, όπως στον διαιτητή
        // (το αποτέλεσμα χωρά πάντα σε int). Όπου ο παλιός δεν υπερχειλίζει,
        // ο ίδιος αριθμός· υγιής δίσκος δεν φτάνει ποτέ κοντά (οι συντελεστές του είναι 9 μπιτ).
        int64_t d8 = 0;
        for (int t = 0; t < ccpo; ++t) d8 += int64_t(kCcpc8[method][t]) * cf[i - t - 1];
        int64_t c = delta;
        if (d8 >= 0) c -= (d8 + 4) / 8;
        else c += (-d8 + 3) / 8;
        cf[i] = int(c > INT32_MAX ? int64_t(INT32_MAX) : (c < INT32_MIN ? int64_t(INT32_MIN) : c));
      }
    }
  }
  // Ο ίδιος περιορισμός στο [-256, 255] (όπως ο παλιός).
  for (int f = 0; f < nrFilters; ++f) {
    int* cf = coef_.data() + size_t(f) * kMaxPredOrder;
    for (int i = 0; i < predOrder[f]; ++i) {
      const int c = cf[i];
      cf[i] = (c > 255) ? 255 : (c < -256 ? -256 : c);
    }
  }

  // Probability_Tables (Table 10.14).
  if ((int)pone_.size() < kMaxElem * kMaxPtableLen) pone_.assign(size_t(kMaxElem) * kMaxPtableLen, 128);
  int ptLen[kMaxElem];
  for (int t = 0; t < nrPtables; ++t) {
    int* pt = pone_.data() + size_t(t) * kMaxPtableLen;
    int len = (int)br.bits(6) + 1;
    ptLen[t] = len;
    for (int e = 0; e < len; ++e) pt[e] = 128;
    if (len == 1) { pt[0] = 128; continue; }
    int coded = br.bit();
    if (coded == 0) {
      for (int e = 0; e < len; ++e) pt[e] = (int)br.bits(7) + 1;
    } else {
      int method = (int)br.bits(2);
      if (method > 2) return false;              // όπως ο παλιός
      int pcpo = kPcpo[method];
      if (len < pcpo) return false;              // όπως ο παλιός
      for (int e = 0; e < pcpo; ++e) pt[e] = (int)br.bits(7) + 1;
      int pcm = (int)br.bits(3);
      for (int e = pcpo; e < len; ++e) {
        int run = 0;
        while (br.bit() == 0) { if (++run > 1024 || br.overrun()) return false; }
        int lsbs = (int)br.bits(pcm);
        int delta = (run << pcm) + lsbs;
        if (delta != 0 && br.bit() == 1) delta = -delta;
        pt[e] = delta;
        for (int tap = 0; tap < pcpo; ++tap) pt[e] -= kPcpc[method][tap] * pt[e - tap - 1];
        if (pt[e] < 1) pt[e] = 1;                // [1, 128], όπως ο παλιός
        else if (pt[e] > 128) pt[e] = 128;
      }
    }
  }

  // ΑΥΣΤΗΡΟΤΕΡΟ ΑΠΟ ΤΟΝ ΠΑΛΙΟ (όπως ο διαιτητής): ολόκληρη η επικεφαλίδα πρέπει να ήταν μέσα στο πλαίσιο. Πέρα από το τέλος
  // ο αναγνώστης δίνει μηδενικά, άρα ένα πλαίσιο κομμένο μέσα στην επικεφαλίδα του θα έδινε φίλτρα και πίνακες από μηδενικά,
  // δηλαδή ένα πλαίσιο θορύβου. Τίποτα δεν έχει γραφτεί ακόμη στην έξοδο.
  if (br.overrun()) { why = DstRefusal::CutHeader; return false; }

  // P(DST_X_Bit) από το Coef[0][0] (10.7.2.4), όπως ο παλιός. Ο παλιός ρωτά «nrFilters > 0 και το
  // coef[0] όχι άδειο»: εδώ ισχύουν πάντα (πάνω: nrFilters >= 1, και κάθε φίλτρο έχει τάξη >= 1).
  int pDstX = 1;
  {
    const int c00 = coef_[0] & 0x1FF;
    int v = 0;
    for (int b = 0; b <= 6; ++b) v = (v << 1) | ((c00 >> b) & 1);
    pDstX = v + 1;
  }
  // ---------- τέλος του αντιγράφου ----------

  // Τα πρώτα 13 μπιτ του αριθμητικού κώδικα, όπως ο ArithDec του παλιού.
  uint32_t C0 = 0;
  for (int i = 0; i <= 12; ++i) C0 = (C0 << 1) + uint32_t(br.bit());
  // ΑΥΣΤΗΡΟΤΕΡΟ ΑΠΟ ΤΟΝ ΠΑΛΙΟ (όπως ο διαιτητής): C >= A από την αρχή (β στην κορυφή του αρχείου). Τίποτα δεν έχει γραφτεί
  // ακόμη στην έξοδο.
  if (C0 >= 4095u) { why = DstRefusal::CodeStart; return false; }

  // Οι πιθανότητες σε bytes (όλες στο [1, 128]), κάθε πίνακας ως τις 64 θέσεις με την τελευταία του τιμή (δες step), και οι
  // πίνακες και οι σταθερές κάθε φίλτρου. Ο χώρος μεγαλώνει μόνο.
  if ((int)p8_.size() < kMaxElem * kMaxPtableLen) p8_.assign(size_t(kMaxElem) * kMaxPtableLen, 128);
  for (int t = 0; t < nrPtables; ++t) {
    const int* pt = pone_.data() + size_t(t) * kMaxPtableLen;
    uint8_t* q = p8_.data() + size_t(t) * kMaxPtableLen;
    const int len = ptLen[t];
    for (int e = 0; e < kMaxPtableLen; ++e) q[e] = uint8_t(pt[e < len ? e : len - 1]);
  }
  if (tab_.size() < size_t(nrFilters) * kTabStride) tab_.resize(size_t(nrFilters) * kTabStride);
  Filt filt[kMaxElem];
  bool slow = false;
  for (int f = 0; f < nrFilters; ++f) {
    buildFilter(tab_.data() + size_t(f) * kTabStride, coef_.data() + size_t(f) * kMaxPredOrder, predOrder[f], filt[f]);
    if (filt[f].S > kMaxAbsSum) slow = true;
  }

  // Τα όρια των διαστημάτων (δ στην κορυφή): 0, 16, frameLen και κάθε αρχή τμήματος ανάμεσα, ταξινομημένα, χωρίς διπλά.
  int brk[3 + 2 * kMaxCh * (kMaxSeg + 1)];
  int nbrk = 0;
  if (frameLen > 0) {
    brk[nbrk++] = 0;
    brk[nbrk++] = frameLen;
    if (kHeadBytes < frameLen) brk[nbrk++] = kHeadBytes;
    for (int ch = 1; ch <= C; ++ch) {
      for (int sg = 0; sg <= fseg[ch].nr; ++sg)
        if (fseg[ch].start[sg] > 0 && fseg[ch].start[sg] < frameLen) brk[nbrk++] = fseg[ch].start[sg];
      for (int sg = 0; sg <= pseg[ch].nr; ++sg)
        if (pseg[ch].start[sg] > 0 && pseg[ch].start[sg] < frameLen) brk[nbrk++] = pseg[ch].start[sg];
    }
    for (int i = 1; i < nbrk; ++i)   // ταξινόμηση με εισαγωγή: το πολύ 111 τιμές, μία φορά ανά πλαίσιο
      for (int j = i; j > 0 && brk[j - 1] > brk[j]; --j) { const int t = brk[j]; brk[j] = brk[j - 1]; brk[j - 1] = t; }
    int u = 1;
    for (int i = 1; i < nbrk; ++i)
      if (brk[i] != brk[u - 1]) brk[u++] = brk[i];
    nbrk = u;
  }

  Frame fr;
  fr.C = C;
  fr.frameLen = frameLen;
  fr.nrFilters = nrFilters;
  fr.fseg = fseg;
  fr.pseg = pseg;
  fr.predOrder = predOrder;
  fr.coef = coef_.data();
  fr.filt = filt;
  fr.tab = tab_.data();
  fr.p8 = p8_.data();
  fr.halfProb = halfProb;
  fr.brk = brk;
  fr.nbrk = nbrk;
  if (st) noteFrame(*st, fr, slow);

  // ---------- Reconstruction ----------
  FastArith fa;
  fa.A = 4095u;
  fa.V = uint64_t(C0) << 51;
  fa.avail = 0;
  while (br.pos & 7) {                           // ως το επόμενο ολόκληρο byte, ένα μπιτ τη φορά
    fa.V |= uint64_t(br.bit()) << (50 - fa.avail);
    ++fa.avail;
  }
  fa.in = in;
  fa.inBytes = inBytes;
  fa.next = br.pos >> 3;
  fa.refill();
  fa.decode(unsigned(pDstX));                    // E[0] = DST_X_Bit (αγνοείται)
  if (slow) reconstructSlow(fa, fr, dsdOut);
  else reconstruct(fa, fr, dsdOut);
  why = DstRefusal::None;
  return true;
}

}  // namespace dstd
}  // namespace

bool DstDecoder::decodeFrame(const uint8_t* in, int inBytes, uint8_t* dsdOut, int bytesPerChannel) {
  const bool ok = dstd::decode(channels_, why_, stats_, coef_, pone_, p8_, tab_, in, inBytes, dsdOut, bytesPerChannel);
  if (stats_) ++stats_->refused[int(why_)];
  return ok;
}
