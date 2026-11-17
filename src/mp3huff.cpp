// Ο γρήγορος πυρήνας των κωδίκων Huffman του MP3. Τι υπόσχεται: mp3huff.h. Η προδιαγραφή του είναι το exact::huffman του
// mp3exact.h (ISO/IEC 11172-3, 2.4.1.7 και 2.4.3.4.6)· εδώ εξηγείται ΠΩΣ βγαίνει το ίδιο αποτέλεσμα γρήγορα και ΓΙΑΤΙ είναι
// το ίδιο. Οι πίνακες αναζήτησης χτίζονται μία φορά από τους πίνακες του προτύπου (mp3tables.inc)· το σχέδιό τους είναι δικό
// μας, βγαλμένο από δύο ιδιότητες των ίδιων των κωδίκων: κανένας κώδικας δεν είναι αρχή άλλου, και κάθε σειρά μπιτ είναι
// κάποιος κώδικας (και τα δύο ελέγχονται στο χτίσιμο των πινάκων).
//
// ΔΥΟ ΔΡΟΜΟΙ, ΜΙΑ ΑΛΗΘΕΙΑ.
//   · Ο ΠΡΟΣΕΚΤΙΚΟΣ δρόμος (finish) είναι ολόκληρος αποκωδικοποιητής: κάνει ό,τι και το exact::huffman, βήμα προς βήμα, με
//     τους ίδιους ελέγχους στην ίδια σειρά. Η μόνη διαφορά του είναι ότι βρίσκει τον κώδικα με πίνακα αντί για δέντρο και
//     παίρνει τα μπιτ από ένα παράθυρο 64 μπιτ (peek) όπου ό,τι βρίσκεται από το endBits και πέρα είναι ήδη μηδέν. Όλες οι
//     ιδιαιτερότητες ζουν ΜΟΝΟ εκεί: ζεύγος ή τετράδα που περνά το τέλος του granule, οι δύο πολιτικές της τελευταίας
//     τετράδας, ανάγνωση πέρα από το τέλος των δεδομένων, πίνακας που «δεν χρησιμοποιείται», οι σημειώσεις, το μηδένισμα.
//   · Ο ΓΡΗΓΟΡΟΣ δρόμος (mp3huff_kernel.inc) είναι επιταχυντής: αποκωδικοποιεί όσο ΚΑΝΕΝΑ από τα παραπάνω δεν μπορεί να
//     συμβεί, και παραδίδει στον προσεκτικό την κατάσταση (γραμμή, θέση) εκεί που σταμάτησε. Ένα βήμα του γίνεται δεκτό
//     μόνο αν τελειώνει ως το όριο min(gend, 8·⌈endBits/8⌉ − 63). Άρα: (α) δεν πέρασε το τέλος του granule, οπότε το exact
//     θα το κρατούσε· (β) όλα τα μπιτ που κατανάλωσε είναι πριν από το endBits, δηλαδή αληθινά, οπότε δεν χρειάζεται μάσκα·
//     (γ) το επόμενο γέμισμα του παραθύρου διαβάζει μέσα στα 8 bytes περιθωρίου (βλ. HUFF_REFILL). Βήμα που δεν γίνεται
//     δεκτό ξεγίνεται, και ο προσεκτικός το ξαναδιαβάζει από την αρχή του: ο γρήγορος δρόμος δεν παίρνει ποτέ απόφαση για
//     άκρη. Για τα τελευταία 63 μπιτ των δεδομένων ο ίδιος γρήγορος δρόμος τρέχει δεύτερη φορά πάνω σε αντίγραφο του
//     τέλους με μηδενικά μετά το endBits (βλ. huffman()).
//
// ΟΙ ΠΙΝΑΚΕΣ. Οι 32 αριθμοί πίνακα μοιράζονται 15 σύνολα κωδίκων (οι 16 … 23 και οι 24 … 31 διαφέρουν μόνο στα linbits).
//   G («γενικός», για κάθε σύνολο): ρίζα 2^n θέσεων με δείκτη τα πρώτα n μπιτ. Κώδικας ως n μπιτ πιάνει όλες τις θέσεις που
//      αρχίζουν με αυτόν (φύλλο: x, y, μήκος κώδικα, και μήκος κώδικα + μπιτ προσήμου). Οι μακρύτεροι κώδικες με κοινά πρώτα
//      n μπιτ πάνε σε υποπίνακα 2^s θέσεων, s = το μακρύτερο υπόλοιπο κάτω από αυτή την αρχή· η θέση της ρίζας γράφει πού
//      είναι και το s. Δεύτερο επίπεδο υποπίνακα δεν υπάρχει: ένας κώδικας θέλει το πολύ δύο αναζητήσεις.
//   D («έτοιμος», για τα σύνολα χωρίς linbits όπου αξίζει): 2^kDBits θέσεις με δείκτη τα επόμενα kDBits μπιτ της ροής.
//      Επειδή τα μπιτ προσήμου ακολουθούν αμέσως τον κώδικα, ένα ζεύγος που χωρά ολόκληρο (κώδικας + πρόσημα) στον δείκτη
//      είναι γνωστό ΜΑΖΙ με τα πρόσημά του από τον δείκτη και μόνο. Αν στα μπιτ που περισσεύουν χωρά ολόκληρο και δεύτερο
//      ζεύγος, η θέση τα έχει και τα δύο: τέσσερις τιμές με πρόσημο, πόσα μπιτ είναι όλα μαζί, και πόσες γραμμές γράφονται
//      (2 ή 4). Θέση 0 = το πρώτο ζεύγος δεν χωρά στον δείκτη: πήγαινε από τον G. Έτσι στα μικρά σύνολα (ζεύγη των 2 … 6
//      μπιτ) μία αναζήτηση δίνει συνήθως δύο ζεύγη. Ο χρόνος ορίζεται από την αλυσίδα «μπιτ → αναζήτηση → μήκος → ολίσθηση»,
//      που τρέχει μία φορά ανά αναζήτηση, όχι ανά ζεύγος.
//      Ποια σύνολα παίρνουν D το αποφασίζει ο ίδιος ο κώδικας του συνόλου (βλ. Tables()): εκείνα όπου, με τις πιθανότητες
//      που υπονοούν τα μήκη των κωδίκων (2^−μήκος), τουλάχιστον 88% των ζευγών χωράνε στον δείκτη, ώστε η διακλάδωση «έτοιμο
//      ή όχι» να προβλέπεται. Με kDBits = 10 αυτά είναι οι πίνακες 1 … 13· ο 15 (83%) και όσοι έχουν linbits πάνε από τον G.
//      Τα μεγέθη (δείκτης 10 μπιτ στον D, ρίζα 10 μπιτ στον G των μεγάλων συνόλων) είναι μετρημένα: με 11 μπιτ οι πίνακες
//      δεν χωράνε πια μαζί στην πιο κοντινή κρυφή μνήμη και ο πυρήνας αργεί· με 9 αυξάνουν οι δεύτερες αναζητήσεις.
//   Τετράδες: ρίζα 64 θέσεων (ο πίνακας A έχει κώδικες ως 6 μπιτ, ο B πάντα 4): μήκος με τα πρόσημα, μήκος κώδικα, ποιες
//      από τις v, w, x, y είναι 1. Οι τέσσερις τιμές με πρόσημο βγαίνουν από δεύτερο, κοινό πίνακα 256 θέσεων με δείκτη
//      (ποιες είναι 1, τα επόμενα 4 μπιτ)· αυτή η δεύτερη αναζήτηση είναι έξω από την αλυσίδα που ορίζει τον χρόνο.
//   Σε κάθε θέση το μήκος που καταναλώνεται είναι ολόκληρο το πρώτο byte: βγαίνει με μία πράξη, και είναι το μόνο πεδίο
//   πάνω στην αλυσίδα που ορίζει τον χρόνο. Τα υπόλοιπα πεδία είναι βαλμένα ώστε να βγαίνουν με όσο γίνεται λιγότερες
//   ολισθήσεις, γιατί στους x86-64 οι ολισθήσεις τρέχουν μόνο σε δύο από τις μονάδες εκτέλεσης και εκεί στενεύει ο βρόχος.
//   Μεγέθη: huffmanTableBytes() και huffmanSelfCheck().
//
// ΤΑ LINBITS. Στους πίνακες 16 … 31 μια τιμή 15 ακολουθείται από linbits μπιτ που προστίθενται. Η θέση του G γράφει ποιες από
// τις δύο τιμές είναι 15 και, σε δικό του byte, κατά πόσο ολισθαίνει το 2·linbits για να δώσει τα επιπλέον μπιτ του ζεύγους
// (0, linbits ή 2·linbits). Έτσι το συνολικό μήκος βγαίνει από τη θέση και μόνο (το επόμενο ζεύγος δεν περιμένει τις τιμές
// αυτού) και ολόκληρο το ζεύγος διαβάζεται από το ίδιο παράθυρο χωρίς διακλάδωση: οι τιμές 15 είναι συχνές εκεί (στο
// σύνολο 24 σχεδόν ένα ζεύγος στα τέσσερα) και απρόβλεπτες.
#include "mp3huff.h"

#include <cstdio>
#include <cstring>

// Άλλα OS και άλλοι μεταγλωττιστές: ο ίδιος κώδικας παντού· αλλάζει μόνο το όνομα λίγων εντολών (αντιστροφή bytes, «τι
// ξέρει αυτός ο επεξεργαστής», ιδιότητες συναρτήσεων). Ό,τι δεν είναι ούτε MSVC ούτε GCC/Clang παίρνει τη σκέτη C++, και
// το ίδιο κάνει το MP3HUFF_PLAIN_CPP (για να δοκιμάζεται εδώ αυτός ο δρόμος).
#if defined(MP3HUFF_PLAIN_CPP)
#define HUFF_MSVC 0
#define HUFF_GNU 0
#elif defined(_MSC_VER) && !defined(__clang__)
#define HUFF_MSVC 1
#define HUFF_GNU 0
#elif defined(__GNUC__) || defined(__clang__)
#define HUFF_MSVC 0
#define HUFF_GNU 1
#else
#define HUFF_MSVC 0
#define HUFF_GNU 0
#endif

#if HUFF_MSVC
#include <intrin.h>
#include <stdlib.h>
#define HUFF_LIKELY(x) (x)
#define HUFF_UNLIKELY(x) (x)
#define HUFF_INLINE __forceinline
#define HUFF_NOINLINE __declspec(noinline)
#define HUFF_TARGET_BMI2
#elif HUFF_GNU
#define HUFF_LIKELY(x) __builtin_expect(!!(x), 1)
#define HUFF_UNLIKELY(x) __builtin_expect(!!(x), 0)
#define HUFF_INLINE inline __attribute__((always_inline))
#define HUFF_NOINLINE __attribute__((noinline))
#define HUFF_TARGET_BMI2 __attribute__((target("bmi2")))
#else
#define HUFF_LIKELY(x) (x)
#define HUFF_UNLIKELY(x) (x)
#define HUFF_INLINE inline
#define HUFF_NOINLINE
#define HUFF_TARGET_BMI2
#endif

// Ο επεξεργαστής γράφει το μικρό άκρο πρώτο; (όλοι οι στόχοι του MSVC ναι· στους GCC/Clang το λέει ο μεταγλωττιστής)
#if HUFF_MSVC || (HUFF_GNU && defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
#define HUFF_LITTLE 1
#else
#define HUFF_LITTLE 0
#endif

#if (defined(__x86_64__) || defined(_M_X64)) && (HUFF_MSVC || HUFF_GNU)
#define HUFF_X86 1
#if !HUFF_MSVC
#include <cpuid.h>
#endif
#include <immintrin.h>
#else
#define HUFF_X86 0
#endif
#if !HUFF_X86 && defined(__aarch64__) && HUFF_GNU && HUFF_LITTLE
#include <arm_neon.h>
#define HUFF_NEON 1
#else
#define HUFF_NEON 0
#endif

// Οι σκόπιμες βλάβες της δοκιμής μετάλλαξης: κάθε αριθμός χαλάει ΕΝΑ πράγμα, για να φανεί ότι η δοκιμή
// το πιάνει. Στο κανονικό χτίσιμο το MP3HUFF_MUTATION είναι 0 και κάθε HUFF_MUT(n) είναι η σταθερά «ψευδές».
#ifndef MP3HUFF_MUTATION
#define MP3HUFF_MUTATION 0
#endif
#define HUFF_MUT(n) (MP3HUFF_MUTATION == (n))

// Μέτρημα χρήσης κάθε θέσης κάθε πίνακα, μόνο για τη δοκιμή κάλυψης (που χτίζει ξεχωριστό εκτελέσιμο με αυτό).
#ifdef MP3HUFF_COVERAGE
#define HUFF_COV(counter) (void)(++(counter))
#else
#define HUFF_COV(counter) (void)0
#endif

namespace orpheic::mp3 {
namespace {

namespace tab = exact::tab;

// Τα 8 bytes από το p, το πρώτο στο πιο σημαντικό byte.
HUFF_INLINE uint64_t be64(const uint8_t *p)
{
#if HUFF_MSVC
    uint64_t v;
    std::memcpy(&v, p, 8);
    return _byteswap_uint64(v);
#elif HUFF_GNU && HUFF_LITTLE
    uint64_t v;
    std::memcpy(&v, p, 8);
    return __builtin_bswap64(v);
#else
    uint64_t v = 0;
    for (int i = 0; i < 8; ++i) v = (v << 8) | p[i];
    return v;
#endif
}

// ── ΟΙ ΠΙΝΑΚΕΣ ─────────────────────────────────────────────────────────────────────────────────────────────────────────

#ifndef MP3HUFF_DBITS
#define MP3HUFF_DBITS 10
#endif
#ifndef MP3HUFF_GBITS
#define MP3HUFF_GBITS 10
#endif
#ifndef MP3HUFF_GSMALL
#define MP3HUFF_GSMALL 8
#endif
#ifndef MP3HUFF_DIRECT_PERCENT
#define MP3HUFF_DIRECT_PERCENT 88
#endif
constexpr unsigned kDBits = MP3HUFF_DBITS;        // δείκτης του D
constexpr unsigned kGBits = MP3HUFF_GBITS;        // ρίζα του G στα σύνολα με μακριούς κώδικες
constexpr unsigned kGSmallBits = MP3HUFF_GSMALL;  // ρίζα του G (το πολύ) στα σύνολα με κώδικες ως 11 μπιτ: εκεί ο G είναι
                                                  // μόνο η εφεδρεία του D και ο προσεκτικός δρόμος
constexpr uint32_t kGCap = (5u << kGBits) + 4096u;   // θέσεις του G για όλα τα σύνολα μαζί (αν δεν φτάνουν, ok = ψευδές)
constexpr int kMaxSets = 15, kMaxDirect = 13;
// Τα όρια που κάνουν σωστό το «τόσα βήματα ανά γέμισμα» του πυρήνα (εκεί εξηγείται το καθένα).
static_assert(kDBits >= 4 && 4u * kDBits <= 56u && 21u + 3u * kDBits <= 56u, "D: τρία βήματα ανά γέμισμα");
static_assert(kGBits >= 1 && kGBits <= 14, "G: δύο ζεύγη των 21 μπιτ ανά γέμισμα αφήνουν 14");
static_assert(kGSmallBits >= 1 && kGSmallBits <= 11, "G των μικρών συνόλων");

// Θέση του G (32 μπιτ).
//   φύλλο:      byte 0   μήκος κώδικα + ένα μπιτ προσήμου για κάθε μη μηδενική τιμή (χωρίς linbits)
//               byte 1   μπιτ 8-12 μήκος κώδικα · μπιτ 13 «το x είναι 15» · μπιτ 14 «το y είναι 15» (μόνο στα σύνολα με linbits)
//               byte 2   x
//               byte 3   μπιτ 24-27 y · μπιτ 28-30 (μόνο στα σύνολα με linbits) κατά πόσο ολισθαίνει δεξιά το 2·linbits για
//                        να δώσει τα μπιτ linbits του ζεύγους: 7 → καμία τιμή 15 (2·13 >> 7 = 0) · 1 → μία · 0 → δύο
//               Στα σύνολα χωρίς linbits τα bytes 1, 2, 3 είναι σκέτα το μήκος κώδικα, το x και το y.
//   παραπομπή:  μπιτ 15 · byte 0 το s του υποπίνακα · μπιτ 16-31 πού αρχίζει ο υποπίνακας, μετρώντας από τη ρίζα του συνόλου.
//               Η αρχή είναι κάτω από 2^14, ώστε τα μπιτ 28-31 μιας παραπομπής να δίνουν μετρητή ολίσθησης κάτω από 4: ο
//               πυρήνας των linbits τον χρησιμοποιεί πριν μάθει ότι η θέση είναι παραπομπή (και πετά το αποτέλεσμα).
// Το φύλλο δεν είναι ποτέ 0 (μήκος ≥ 1), οπότε το 0 σημαίνει «άγραφη θέση» όσο χτίζεται ο πίνακας.
constexpr uint32_t kLink = 0x8000u;

constexpr uint32_t leafEntry(int x, int y, int len, bool linSet)
{
    const int k15 = linSet ? (x == 15) + (y == 15) : 0;
    return uint32_t(len + (x != 0) + (y != 0)) | uint32_t(len) << 8 | uint32_t(linSet && x == 15) << 13 | uint32_t(linSet && y == 15) << 14 |
           uint32_t(x) << 16 | uint32_t(y) << 24 | (!linSet ? 0u : k15 == 0 ? 7u : k15 == 1 ? 1u : 0u) << 28;
}

// Οι μάσκες προσήμου ενός ζεύγους χωρίς linbits, για τον γρήγορο δρόμο: δείκτης 4·d + t, όπου d = πόσες από τις δύο τιμές
// δεν είναι 0 (μήκος με πρόσημα − μήκος κώδικα) και t = τα δύο μπιτ αμέσως μετά τον κώδικα (το πρώτο στο μπιτ 1). Τα κάτω
// 16 μπιτ είναι η μάσκα του x, τα πάνω του y· όλο άσοι = αρνητικό. Με δύο μη μηδενικές τιμές το πρώτο μπιτ είναι το πρόσημο
// του x και το δεύτερο του y. Με μία, το πρώτο μπιτ είναι το πρόσημό της, όποια κι αν είναι· μπαίνει και στις δύο μάσκες,
// γιατί η αλλαγή προσήμου του 0 δίνει 0. Με καμία οι μάσκες δεν παίζουν ρόλο.
constexpr uint32_t kSignMasks[12] = {0, 0, 0, 0, 0, 0, 0xFFFFFFFFu, 0xFFFFFFFFu, 0, 0xFFFF0000u, 0x0000FFFFu, 0xFFFFFFFFu};

enum Kind : uint8_t {
    kZero,     // πίνακας 0: όλη η περιοχή μηδέν, χωρίς κανένα μπιτ
    kBad,      // πίνακας 4 ή 14 (ή αριθμός έξω από 0 … 31): δεν χρησιμοποιείται
    kDirect,   // χωρίς linbits, με πίνακα D (και εφεδρεία τον G)
    kPlain,    // χωρίς linbits, μόνο με τον G
    kLin,      // με linbits (16 … 31): G
};

struct Set {
    uint32_t gBase = 0;   // πού αρχίζει η ρίζα του στο g[]
    uint32_t gSize = 0;   // ρίζα + υποπίνακες
    uint8_t gBits = 0;    // n της ρίζας
    int8_t d = -1;        // ποιος D, ή −1
    uint8_t firstSel = 0; // ο μικρότερος αριθμός πίνακα που το χρησιμοποιεί (για τις αναφορές)
    uint8_t directPercent = 0;   // πόσα στα εκατό των ζευγών χωράνε ολόκληρα στον δείκτη του D (με πιθανότητες 2^−μήκος)
};

struct Tables {
    uint32_t g[kGCap];
    // Θέση του D (32 μπιτ): byte 0 τα μπιτ όλων των ζευγών της θέσης · μπιτ 11-15, 16-20, 21-25, 26-30 οι τέσσερις τιμές
    // (5 μπιτ η καθεμία, συμπλήρωμα ως προς δύο: −15 … 15) · μπιτ 31 «η θέση έχει δύο ζεύγη». Θέση χωρίς ζεύγος: 0.
    uint32_t d[kMaxDirect][1u << kDBits];
    // Τετράδες. quad[πίνακας][πρώτα 6 μπιτ]: byte 0 μήκος κώδικα + ένα μπιτ προσήμου για κάθε 1 · μπιτ 8-11 μήκος κώδικα ·
    // μπιτ 12-15 ποιες τιμές είναι 1 (v το πιο σημαντικό). quadVal[(ποιες είναι 1) · 16 + τα 4 μπιτ μετά τον κώδικα]: οι
    // τέσσερις τιμές με πρόσημο, όπως γράφονται στο is.
    uint16_t quad[2][64];
    alignas(8) int16_t quadVal[256][4];
    Set set[kMaxSets];
    uint8_t kind[33], setOf[33], lin[33];   // η θέση 32 είναι ο «άγνωστος αριθμός πίνακα»
    int sets = 0, directs = 0;
    uint32_t gUsed = 0;
    bool ok = true;   // ψευδές μόνο αν οι πίνακες του προτύπου δεν είναι αυτό που περιμένει το σχέδιο (το ελέγχει η δοκιμή)

    // Ο G ενός συνόλου, στο out[0 … ). Επιστρέφει πόσες θέσεις έπιασε (0 = δεν χωρά).
    uint32_t buildG(uint32_t *out, uint32_t cap, const tab::HuffCode *codes, int count, unsigned n, bool linSet)
    {
        const uint32_t rootSize = 1u << n;
        if (rootSize > cap) return 0;
        for (uint32_t i = 0; i < rootSize; ++i) out[i] = 0;
        // Πρώτα πόσο μακρύ είναι το μακρύτερο υπόλοιπο κάτω από κάθε αρχή n μπιτ (προσωρινά μέσα στη ρίζα).
        for (int c = 0; c < count; ++c) {
            const tab::HuffCode &h = codes[c];
            if (h.len <= n) continue;
            const uint32_t pre = h.code >> (h.len - n);
            if (h.len - n > out[pre]) out[pre] = h.len - n;
        }
        uint32_t used = rootSize;
        for (uint32_t pre = 0; pre < rootSize; ++pre) {
            const uint32_t s = out[pre];
            if (!s) continue;
            if (used + (1u << s) > cap || used >= 0x4000u) return 0;
            out[pre] = kLink | used << 16 | s;
            for (uint32_t i = 0; i < (1u << s); ++i) out[used + i] = 0;
            used += 1u << s;
        }
        // Κάθε κώδικας γράφεται σε όλες τις θέσεις που αρχίζουν με αυτόν. Επειδή κανένας κώδικας δεν είναι αρχή άλλου, δύο
        // κώδικες δεν διεκδικούν ποτέ την ίδια θέση· επειδή κάθε σειρά μπιτ είναι κάποιος κώδικας, δεν μένει θέση άγραφη.
        // Και τα δύο ελέγχονται εδώ (ok), ώστε λάθος στους πίνακες του προτύπου να φανεί αντί να γίνει λάθος ήχος.
        for (int c = 0; c < count; ++c) {
            const tab::HuffCode &h = codes[c];
            if (h.len == 0) continue;
            const uint32_t e = leafEntry(h.x, h.y, h.len, linSet);
            uint32_t at, span;
            if (h.len <= n) {
                span = 1u << (n - h.len);
                at = h.code << (n - h.len);
            } else {
                const uint32_t link = out[h.code >> (h.len - n)];
                const uint32_t s = link & 0xFFu, rest = h.len - n;
                span = 1u << (s - rest);
                at = (link >> 16) + ((h.code & ((1u << rest) - 1u)) << (s - rest));
            }
            for (uint32_t i = 0; i < span; ++i) {
                if (out[at + i] != 0) ok = false;
                out[at + i] = e;
            }
        }
        for (uint32_t i = 0; i < used; ++i)
            if (out[i] == 0) ok = false;
        return used;
    }

    // Ο D ενός συνόλου χωρίς linbits. Για κάθε δείκτη (τα επόμενα kDBits μπιτ της ροής) γίνεται εδώ, μία φορά, ό,τι θα έκανε
    // το exact::huffman για ως δύο ζεύγη: ο κώδικας που αρχίζει εκεί, ένα μπιτ προσήμου για το x αν δεν είναι 0, ένα για το y.
    // Ένα ζεύγος μπαίνει μόνο αν ΟΛΑ τα μπιτ του είναι μέσα στον δείκτη: τότε το αποτέλεσμα εξαρτάται μόνο από μπιτ που
    // καταναλώνονται, και ό,τι ακολουθεί στον δείκτη δεν το αλλάζει.
    void buildD(uint32_t *out, const tab::HuffCode *codes, int count)
    {
        for (uint32_t i = 0; i < (1u << kDBits); ++i) {
            int v[4] = {0, 0, 0, 0};
            unsigned used = 0, pairs = 0;
            for (unsigned k = 0; k < 2; ++k) {
                const unsigned room = kDBits - used;
                const tab::HuffCode *hit = nullptr;
                for (int c = 0; c < count && !hit; ++c)
                    if (codes[c].len != 0 && codes[c].len <= room &&
                        ((i >> (room - codes[c].len)) & ((1u << codes[c].len) - 1u)) == codes[c].code)
                        hit = &codes[c];
                if (!hit) break;   // ο κώδικας που αρχίζει εδώ είναι μακρύτερος από όσα μπιτ έμειναν
                const unsigned t = hit->len + (hit->x != 0) + (hit->y != 0);
                // ο κώδικας χωρά, τα πρόσημά του όχι: το ζεύγος δεν μπαίνει
                if (t > room && !(HUFF_MUT(17) && k == 1)) break;   // ΒΛΑΒΗ 17: το δεύτερο ζεύγος μπαίνει και χωρίς τα πρόσημά του
                int at = int(room - hit->len);   // πόσα μπιτ του δείκτη μένουν μετά τον κώδικα
                int x = hit->x, y = hit->y;
                if (x) {
                    --at;
                    if (at >= 0 && ((i >> at) & 1u)) x = -x;
                }
                if (y) {
                    --at;
                    if (at >= 0 && ((i >> at) & 1u)) y = -y;
                }
                v[2 * k] = x;
                v[2 * k + 1] = y;
                used += t;
                ++pairs;
            }
            out[i] = !pairs ? 0u
                            : used | (uint32_t(v[0]) & 31u) << 11 | (uint32_t(v[1]) & 31u) << 16 | (uint32_t(v[2]) & 31u) << 21 |
                                  (uint32_t(v[3]) & 31u) << 26 | (pairs - 1u) << 31;
        }
        if (HUFF_MUT(2)) out[0x2A5 & ((1u << kDBits) - 1u)] ^= 1u << 11;   // ΒΛΑΒΗ 2: μία θέση ενός D δίνει άλλο x
    }

    Tables()
    {
        const tab::HuffCode *seen[kMaxSets] = {};
        bool linSet[kMaxSets] = {};
        int maxLen[kMaxSets] = {}, maxVal[kMaxSets] = {};
        for (int i = 0; i < 33; ++i) {
            kind[i] = kBad;
            setOf[i] = 0;
            lin[i] = 0;
        }
        kind[0] = kZero;
        // Ποια σύνολα υπάρχουν, και αν κάποιος αριθμός πίνακα τα χρησιμοποιεί με linbits.
        for (int sel = 1; sel < 32; ++sel) {
            const tab::HuffTable &ht = tab::kHuffTable[sel];
            if (!ht.codes) continue;
            int s = 0;
            while (s < sets && seen[s] != ht.codes) ++s;
            if (s == sets) {
                if (sets == kMaxSets) {
                    ok = false;
                    continue;
                }
                seen[sets++] = ht.codes;
                set[s].firstSel = uint8_t(sel);
                for (int c = 0; c < ht.count; ++c) {
                    if (ht.codes[c].len > maxLen[s]) maxLen[s] = ht.codes[c].len;
                    if (ht.codes[c].x > maxVal[s]) maxVal[s] = ht.codes[c].x;
                    if (ht.codes[c].y > maxVal[s]) maxVal[s] = ht.codes[c].y;
                }
            }
            if (ht.linbits) linSet[s] = true;
            setOf[sel] = uint8_t(s);
            lin[sel] = uint8_t(ht.linbits);
            // Ό,τι υποθέτει το παράθυρο: κώδικας ως 19 μπιτ, linbits ως 13, τιμές ως 15· και το μακρύτερο ζεύγος αυτού του
            // αριθμού πίνακα (κώδικας, linbits για κάθε 15, πρόσημα) αφήνει από τα 56 μπιτ αρκετά για την επόμενη αναζήτηση.
            if (ht.linbits < 0 || ht.linbits > 13 || maxVal[s] > 15 || maxLen[s] > 19) ok = false;
            for (int c = 0; c < ht.count; ++c) {
                const tab::HuffCode &h = ht.codes[c];
                const unsigned bits = h.len + unsigned(h.x != 0) + unsigned(h.y != 0) +
                                      unsigned(ht.linbits) * (unsigned(h.x == 15) + unsigned(h.y == 15));
                if (ht.linbits && bits + kGBits > 56u) ok = false;
            }
        }
        for (int s = 0; s < sets; ++s) {
            const tab::HuffTable &ht = tab::kHuffTable[set[s].firstSel];
            // Πόσα από τα ζεύγη χωράνε ολόκληρα (κώδικας + πρόσημα) στον δείκτη του D, αν κάθε κώδικας έρχεται με πιθανότητα
            // 2^−μήκος (αυτή για την οποία σχεδιάστηκε). Σε μονάδες 2^−19, γιατί ο μακρύτερος κώδικας είναι 19 μπιτ.
            uint32_t fit = 0;
            for (int c = 0; c < ht.count; ++c) {
                const tab::HuffCode &h = ht.codes[c];
                if (h.len && h.len <= 19 && h.len + unsigned(h.x != 0) + unsigned(h.y != 0) <= kDBits) fit += 1u << (19 - h.len);
            }
            set[s].directPercent = uint8_t(uint64_t(fit) * 100u >> 19);
            const bool direct = !linSet[s] && set[s].directPercent >= MP3HUFF_DIRECT_PERCENT && directs < kMaxDirect;
            const unsigned n = !direct || maxLen[s] > 11 ? kGBits : unsigned(maxLen[s]) < kGSmallBits ? unsigned(maxLen[s]) : kGSmallBits;
            set[s].gBase = gUsed;
            set[s].gBits = uint8_t(n);
            // ΒΛΑΒΗ 21: τα σύνολα χωρίς linbits σημειώνουν κι αυτά ποιες τιμές είναι 15 (και ο πυρήνας τις μετρά)
            set[s].gSize = buildG(g + gUsed, kGCap - gUsed, ht.codes, ht.count, n, linSet[s] || (HUFF_MUT(21) && !direct));
            if (!set[s].gSize) ok = false;
            gUsed += set[s].gSize;
            if (direct) {
                set[s].d = int8_t(directs);
                buildD(d[directs++], ht.codes, ht.count);
            }
        }
        for (int sel = 1; sel < 32; ++sel) {
            if (!tab::kHuffTable[sel].codes) continue;
            const Set &st = set[setOf[sel]];
            kind[sel] = st.d >= 0 ? kDirect : (lin[sel] || HUFF_MUT(21)) ? kLin : kPlain;
        }
        if (HUFF_MUT(1)) {
            // ΒΛΑΒΗ 1: μία θέση ενός υποπίνακα (ο μακρύτερος κώδικας του πίνακα 13, 19 μπιτ όλο μηδέν) δίνει άλλο y
            const Set &st = set[setOf[13]];
            g[st.gBase + (g[st.gBase] >> 16)] ^= 1u << 24;
        }
        // Τετράδες.
        for (int t = 0; t < 2; ++t) {
            const tab::HuffCode *codes = t ? tab::kHuffQuadB : tab::kHuffQuadA;
            for (uint32_t i = 0; i < 64; ++i) {
                quad[t][i] = 0;
                for (int c = 0; c < 16; ++c) {
                    const tab::HuffCode &h = codes[c];
                    if (h.len == 0 || h.len > 6 || (i >> (6 - h.len)) != h.code) continue;
                    if (quad[t][i]) ok = false;
                    const unsigned ones = (h.x >> 3 & 1u) + (h.x >> 2 & 1u) + (h.x >> 1 & 1u) + (h.x & 1u);
                    quad[t][i] = uint16_t((h.len + ones) | unsigned(h.len) << 8 | unsigned(h.x) << 12);
                }
                if (!quad[t][i]) ok = false;
            }
        }
        // Τα πρόσημα της τετράδας: ένα μπιτ για κάθε τιμή που είναι 1, με τη σειρά v, w, x, y (2.4.1.7).
        for (unsigned m = 0; m < 16; ++m)
            for (unsigned s = 0; s < 16; ++s) {
                unsigned next = 3;   // το επόμενο αδιάβαστο από τα 4 μπιτ, το πιο σημαντικό πρώτο
                for (unsigned k = 0; k < 4; ++k) {
                    const unsigned who = HUFF_MUT(11) ? k ^ 1u : k;   // ΒΛΑΒΗ 11: τα πρόσημα των v, w (και των x, y) ανάποδα
                    int val = 0;
                    if ((m >> (3 - who)) & 1u) {
                        val = (s >> next) & 1u ? -1 : 1;
                        --next;
                    }
                    quadVal[m << 4 | s][who] = int16_t(val);
                }
            }
    }
};

const Tables &tables()
{
    static const Tables t;
    return t;
}

#ifdef MP3HUFF_COVERAGE
struct Coverage {
    uint32_t d[kMaxDirect][1u << kDBits];
    uint32_t gFast[kGCap], gCareful[kGCap];
    uint32_t quadFast[2][64], quadCareful[2][64];
    uint32_t quadValFast[256], quadValCareful[256];
};
Coverage gCov;
#endif

// ── ΤΙ ΠΡΕΠΕΙ ΝΑ ΓΙΝΕΙ ΚΑΙ ΠΟΥ ΒΡΙΣΚΟΜΑΣΤΕ ──────────────────────────────────────────────────────────────────────────────

// Το exact::HuffmanIn, λυμένο μία φορά: το exact ρωτά για κάθε ζεύγος «l < r1 ? 0 : l < r2 ? 1 : 2»· εδώ η ίδια ερώτηση
// έχει απαντηθεί ως τρία τέλη περιοχών.
struct Plan {
    int end[3];       // η γραμμή όπου τελειώνει κάθε περιοχή: ζυγή, αύξουσα, ≤ 576· end[2] = οι γραμμές των ζευγών
    uint8_t sel[3];   // ο πίνακας κάθε περιοχής, 0 … 31· 32 = αριθμός έξω από αυτά
    uint8_t count1;   // 0 = πίνακας A των τετράδων, 1 = B
};

struct State {
    int l;           // πόσες γραμμές έχουν γραφτεί
    int region;      // σε ποια περιοχή είμαστε (0 … 2)· 3 = τα ζεύγη τελείωσαν όλα, σειρά έχουν οι τετράδες
    size_t pos;      // το επόμενο αδιάβαστο μπιτ
    unsigned nLin;   // τιμές που διάβασαν linbits ως εδώ
};

// ── Ο ΓΡΗΓΟΡΟΣ ΔΡΟΜΟΣ ──────────────────────────────────────────────────────────────────────────────────────────────────

// Οι τέσσερις τιμές μιας θέσης του D ως τέσσερις int16, με μία αποθήκευση 8 bytes. Η θέση είναι δύο λέξεις 16 μπιτ: στην
// κάτω η πρώτη τιμή είναι ήδη στα 5 πάνω μπιτ· η πάνω λέξη έχει τις άλλες τρεις στα μπιτ 0-4, 5-9, 10-14. Τέσσερις λωρίδες
// (κάτω, πάνω, πάνω, πάνω), πολλαπλασιασμός με 1, 2^11, 2^6, 2 ώστε η τιμή κάθε λωρίδας να φτάσει στα 5 πάνω μπιτ της (ό,τι
// ξεχειλίζει χάνεται, μαζί και το μπιτ 31), και αριθμητική ολίσθηση δεξιά κατά 11 που απλώνει το πρόσημο. Το SSE2 το έχει κάθε x86-64, το NEON
// κάθε ARM64· αλλού, τέσσερις αποθηκεύσεις.
HUFF_INLINE void expand4(int16_t *o, uint32_t m)
{
#if HUFF_X86
    __m128i v = _mm_shufflelo_epi16(_mm_cvtsi32_si128(int(m)), 0x54);
    v = _mm_srai_epi16(_mm_mullo_epi16(v, _mm_set_epi16(0, 0, 0, 0, 2, 64, 2048, 1)), 11);
    _mm_storel_epi64(reinterpret_cast<__m128i *>(o), v);
#elif HUFF_NEON
    static const int16_t up[4] = {0, 11, 6, 1};
    const int16x4_t v = vset_lane_s16(int16_t(uint16_t(m)), vdup_n_s16(int16_t(uint16_t(m >> 16))), 0);
    vst1_s16(o, vshr_n_s16(vshl_s16(v, vld1_s16(up)), 11));
#else
    for (unsigned k = 0; k < 4; ++k) o[k] = int16_t(int(((m >> (11u + 5u * k)) & 31u) ^ 16u) - 16);   // 5 μπιτ με πρόσημο
#endif
}

// Ένα ζεύγος χωρίς linbits στο o[0], o[1]: οι τιμές x, y της θέσης `e` του G με τα πρόσημά τους. `w`: το παράθυρο
// ολισθημένο ώστε να αρχίζει αμέσως μετά τον κώδικα· `d`: πόσες από τις δύο τιμές δεν είναι 0. Οι μάσκες από το kSignMasks
// (εκεί εξηγείται γιατί αρκούν τα δύο πρώτα μπιτ μετά τον κώδικα)· τιμή με μάσκα «όλο άσοι» αλλάζει πρόσημο: (v ^ m) − m.
// Στον x86-64 τα δύο μαζί ως δύο λωρίδες 16 μπιτ: έτσι το ζεύγος θέλει τρεις ολισθήσεις αντί για επτά.
// ΒΛΑΒΗ 3: το πρόσημο του y ανάποδα.
HUFF_INLINE void storePair(int16_t *o, uint32_t e, uint64_t w, unsigned d)
{
    const uint32_t m = kSignMasks[4u * d + unsigned(w >> 62)] ^ (HUFF_MUT(3) ? 0xFFFF0000u : 0u);
#if HUFF_X86
    const __m128i mask = _mm_cvtsi32_si128(int(m));
    const __m128i xy = _mm_unpacklo_epi8(_mm_cvtsi32_si128(int(e >> 16)), _mm_setzero_si128());
    const int both = _mm_cvtsi128_si32(_mm_sub_epi16(_mm_xor_si128(xy, mask), mask));
    std::memcpy(o, &both, 4);
#else
    const uint32_t x = (e >> 16) & 0xFFu, y = e >> 24, mx = m & 0xFFFFu, my = m >> 16;
    o[0] = int16_t(uint16_t((x ^ mx) - mx));
    o[1] = int16_t(uint16_t((y ^ my) - my));
#endif
}

// Ο πυρήνας σε δύο μορφές στον x86-64: με τις ολισθήσεις που έχει κάθε επεξεργαστής (μετρητής στον CL: τρεις μικροεντολές
// στους Intel), και με SHLX/SHRX (μία μικροεντολή, όπου υπάρχει BMI2: από το 2013 και μετά). Διαλέγεται μία φορά.
#define HUFF_KERNEL_NAME fastPlain
#define HUFF_KERNEL_ATTR
#define HUFF_SHL(v, n) ((v) << (n))
#define HUFF_SHR(v, n) ((v) >> (n))
#include "mp3huff_kernel.inc"
#undef HUFF_KERNEL_NAME
#undef HUFF_KERNEL_ATTR
#undef HUFF_SHL
#undef HUFF_SHR

#if HUFF_X86
#define HUFF_KERNEL_NAME fastBmi2
#define HUFF_KERNEL_ATTR HUFF_TARGET_BMI2
#if HUFF_MSVC
#define HUFF_SHL(v, n) _shlx_u64((v), unsigned(n))
#define HUFF_SHR(v, n) _shrx_u64((v), unsigned(n))
#else   // μέσα σε συνάρτηση με target("bmi2") ο μεταγλωττιστής βγάζει μόνος του SHLX/SHRX
#define HUFF_SHL(v, n) ((v) << (n))
#define HUFF_SHR(v, n) ((v) >> (n))
#endif
#include "mp3huff_kernel.inc"
#undef HUFF_KERNEL_NAME
#undef HUFF_KERNEL_ATTR
#undef HUFF_SHL
#undef HUFF_SHR

bool cpuHasBmi2()
{
#if HUFF_MSVC
    int r[4];
    __cpuid(r, 0);
    if (r[0] < 7) return false;
    __cpuidex(r, 7, 0);
    return (r[1] & (1 << 8)) != 0;
#else
    unsigned a = 0, b = 0, c = 0, d = 0;
    return __get_cpuid_count(7, 0, &a, &b, &c, &d) && (b & (1u << 8));
#endif
}
#endif

using FastFn = void (*)(const Tables &, const uint8_t *, int64_t, const Plan &, int16_t *, State &);
#if HUFF_X86
FastFn gFast = cpuHasBmi2() ? fastBmi2 : fastPlain;   // αλλάζει μόνο με το huffmanUseBasicOnly (δοκιμή του βασικού δρόμου)
#else
FastFn gFast = fastPlain;
#endif
bool gCarefulOnly = false;

// ── Ο ΠΡΟΣΕΚΤΙΚΟΣ ΔΡΟΜΟΣ ───────────────────────────────────────────────────────────────────────────────────────────────

// Τα μπιτ της ροής από το `pos`, το πρώτο στο πιο σημαντικό μπιτ, με μηδέν σε κάθε θέση από το endBits και πέρα: ό,τι θα
// έδιναν 64 διαδοχικά exact::Bits::get1(). Τα πρώτα 57 τουλάχιστον είναι σωστά (το διάβασμα αρχίζει σε όριο byte, άρα ως 7
// μπιτ χάνονται από κάτω και μπαίνουν μηδενικά)· κανένα ζεύγος δεν θέλει πάνω από 19 + 2·13 + 2 = 47.
// Μνήμη: με pos < endBits το byte pos/8 είναι byte δεδομένων, άρα τα 8 bytes από εκεί τελειώνουν το πολύ στο 7ο byte του
// περιθωρίου. Με pos ≥ endBits δεν διαβάζεται τίποτα.
HUFF_INLINE uint64_t peek(const uint8_t *data, size_t endBits, size_t pos)
{
    if (pos >= endBits) return 0;
    const uint64_t v = be64(data + (pos >> 3)) << (pos & 7);
    const size_t have = endBits - pos;
    if (HUFF_MUT(10)) return v;   // ΒΛΑΒΗ 10: τα μπιτ μετά το endBits δεν μηδενίζονται (το περιθώριο επηρεάζει το αποτέλεσμα)
    return have >= 64 ? v : v & (~uint64_t(0) << (64 - have));
}

// Συνεχίζει από την κατάσταση `st` ως το τέλος και κάνει όλο το κλείσιμο. Με st = {0, 0, pos, 0} είναι ολόκληρο το
// exact::huffman· τα «(Ε)» στα σχόλια δείχνουν το αντίστοιχο σημείο εκεί.
HUFF_NOINLINE int finish(const Tables &T, const uint8_t *data, size_t endBits, size_t gend, const Plan &plan, int16_t *is,
                         const State &st, size_t &posOut, exact::HuffmanNotes *notes)
{
    int l = st.l;
    size_t pos = st.pos;
    unsigned nLin = st.nLin;
    bool overrun = false;
    for (int region = st.region; region < 3 && !overrun; ++region) {
        const int end = plan.end[region];
        if (l >= end) continue;
        const unsigned sel = plan.sel[region];
        if (T.kind[sel] == kZero) {
            // (Ε) ο πίνακας 0 δεν έχει κώδικες: μηδενικά χωρίς κανένα μπιτ, και χωρίς έλεγχο τέλους (δεν διαβάστηκε τίποτα)
            if (HUFF_MUT(19) && pos > gend) {   // ΒΛΑΒΗ 19: η περιοχή με πίνακα 0 ελέγχει το τέλος του granule
                overrun = true;
                break;
            }
            std::memset(is + l, 0, size_t(end - l) * sizeof(int16_t));
            l = end;
            continue;
        }
        if (T.kind[sel] == kBad) {
            // (Ε) bad: is = όλο μηδέν, −1. Η θέση μένει εκεί που ήταν (αυτό το ζεύγος δεν διάβασε τίποτα), και όσα linbits
            // μετρήθηκαν πριν μένουν μετρημένα.
            std::memset(is, 0, 576 * sizeof(int16_t));
            if (notes) notes->linbitsValues += nLin;
            posOut = HUFF_MUT(16) ? pos + 1 : pos;   // ΒΛΑΒΗ 16: λάθος θέση όταν ο πίνακας δεν χρησιμοποιείται
            return -1;
        }
        const Set &set = T.set[T.setOf[sel]];
        const uint32_t *const G = T.g + set.gBase;
        const unsigned gb = set.gBits, lin = T.lin[sel];
        while (l < end) {
            uint64_t w = peek(data, endBits, pos);
            const uint32_t i0 = uint32_t(w >> (64 - gb));
            uint32_t e = G[i0];
            HUFF_COV(gCov.gCareful[set.gBase + i0]);
            if (e & kLink) {
                const uint32_t i1 = (e >> 16) + uint32_t((w << gb) >> (64 - (e & 0xFFu)));
                e = G[i1];
                HUFF_COV(gCov.gCareful[set.gBase + i1]);
            }
            unsigned used = (e >> 8) & 31u;
            int x = int((e >> 16) & 15u), y = int((e >> 24) & 15u);
            w <<= used;
            // (Ε) x: linbits αν είναι 15 και ο πίνακας έχει· μετά ένα μπιτ προσήμου αν δεν είναι 0. Το ίδιο για το y.
            if (x == 15 && (lin || HUFF_MUT(5))) {   // ΒΛΑΒΗ 5: τα linbits «εφαρμόζονται» και όταν είναι 0 (μετριούνται)
                if (lin) {
                    x += int(w >> (64 - lin));
                    w <<= lin;
                    used += lin;
                }
                ++nLin;
            }
            if (x) {
                if (w >> 63) x = -x;
                w <<= 1;
                ++used;
            }
            if (y == 15 && (lin || HUFF_MUT(5))) {
                if (lin) {
                    y += int(w >> (64 - lin));
                    w <<= lin;
                    used += lin;
                }
                ++nLin;
            }
            if (y) {
                if (bool(w >> 63) != bool(HUFF_MUT(4))) y = -y;   // ΒΛΑΒΗ 4: το πρόσημο του y ανάποδα (προσεκτικός δρόμος)
                ++used;
            }
            pos += used;
            if (pos > gend) {
                // (Ε) το ζεύγος πέρασε το τέλος των μπιτ του granule: δεν υπάρχει, ούτε όσα θα ακολουθούσαν. Η θέση μένει
                // μετά από αυτό, και τα linbits του έχουν ήδη μετρηθεί (το exact τα μετρά πριν από τον έλεγχο).
                overrun = true;
                if (HUFF_MUT(6)) {   // ΒΛΑΒΗ 6: το ζεύγος που πέρασε το τέλος κρατιέται
                    is[l] = int16_t(x);
                    is[l + 1] = int16_t(y);
                    l += 2;
                }
                if (HUFF_MUT(9)) --pos;   // ΒΛΑΒΗ 9: η τελική θέση ένα μπιτ πριν
                break;
            }
            is[l] = int16_t(x);
            is[l + 1] = int16_t(y);
            l += 2;
        }
    }
    if (!overrun) {
        // (Ε) οι τετράδες: ώσπου να τελειώσουν τα μπιτ του granule ή οι 576 γραμμές, ό,τι έρθει πρώτο.
        const uint16_t *const Q = T.quad[plan.count1];
        while (l + 4 <= 576 && (pos < gend || (HUFF_MUT(12) && pos == gend))) {   // ΒΛΑΒΗ 12: οι τετράδες συνεχίζουν και με θέση = gend
            const uint64_t w = peek(data, endBits, pos);
            const unsigned q = Q[w >> 58];
            const unsigned vi = ((q >> 8) & 0xF0u) | unsigned((w << ((q >> 8) & 15u)) >> 60);
            HUFF_COV(gCov.quadCareful[plan.count1][w >> 58]);
            HUFF_COV(gCov.quadValCareful[vi]);
            const int16_t *const v = T.quadVal[vi];
            pos += q & 0xFFu;
            if (pos > gend) {
                if (notes) {
                    notes->quadOverrun = true;
                    notes->quadOverrunLine = HUFF_MUT(15) ? l + 4 : l;   // ΒΛΑΒΗ 15: λάθος γραμμή στη σημείωση
                }
                break;
            }
            std::memcpy(is + l, v, 4 * sizeof(int16_t));
            l += 4;
        }
        if (notes && pos < gend) notes->stuffingBits += gend - pos - (HUFF_MUT(14) ? 1 : 0);   // ΒΛΑΒΗ 14: ένα μπιτ γεμίσματος λιγότερο
    } else if (notes) {
        notes->bigOverrun = true;
    }
    std::memset(is + l, 0, size_t(576 - l) * sizeof(int16_t));
    if (notes) notes->linbitsValues += nLin;
    posOut = pos;
    return l;
}

}  // namespace

int huffman(const uint8_t *data, size_t endBits, size_t &pos, size_t gend, const exact::HuffmanIn &in, int16_t *is,
            exact::HuffmanNotes *notes)
{
    const Tables &T = tables();
    // Τα τέλη των περιοχών. Το exact προχωρά ανά δύο γραμμές από το 0, άρα ρωτά μόνο για ζυγές γραμμές: «l < r» για ζυγό l
    // είναι το ίδιο με «l < το r στρογγυλεμένο προς τα πάνω σε ζυγό». Περιοχή 0: l < r1. Περιοχή 1: r1 ≤ l < r2 (κενή αν
    // r2 ≤ r1). Περιοχή 2: τα υπόλοιπα ως το big. Το big πάνω από 576 θα έγραφε έξω από το is· κόβεται.
    const auto even = [](int v) { return v <= 0 ? 0 : v >= 576 ? 576 : (v + 1) & ~1; };
    Plan plan;
    const int big = even(in.big);
    const int r1 = HUFF_MUT(13) ? (in.r1 <= 0 ? 0 : in.r1 >= 576 ? 576 : in.r1 & ~1) : even(in.r1);   // ΒΛΑΒΗ 13: μονό r1 προς τα κάτω
    const int r2 = even(in.r2);
    plan.end[0] = r1 < big ? r1 : big;
    plan.end[1] = r2 < plan.end[0] ? plan.end[0] : r2 < big ? r2 : big;
    plan.end[2] = big;
    for (int i = 0; i < 3; ++i) plan.sel[i] = uint8_t(unsigned(in.tableSelect[i]) < 32u ? in.tableSelect[i] : 32);
    plan.count1 = uint8_t(in.count1Table & 1);

    State st{0, 0, pos, 0};
    if (!gCarefulOnly) {
        // 1η φάση: πάνω στα ίδια τα δεδομένα. Ένα βήμα επιτρέπεται να τελειώνει ως το min(gend, 8·⌈endBits/8⌉ − 63) (βλ. την
        // κορυφή του αρχείου). Όριο αρνητικό ή μικρότερο από το pos: η φάση δεν ξεκινά (π.χ. pos ήδη μετά το gend).
        const uint64_t nBytes = uint64_t((endBits + 7) >> 3);
        // ΒΛΑΒΗ 23: ο γρήγορος δρόμος φτάνει ένα byte πιο κοντά στο τέλος των δεδομένων (το γέμισμα διαβάζει έξω από το περιθώριο)
        const int64_t byData = nBytes * 8u <= uint64_t(INT64_MAX) ? int64_t(nBytes * 8u) - (HUFF_MUT(23) ? 55 : 63) : -1;
        const int64_t byGranule = gend <= uint64_t(INT64_MAX) ? int64_t(gend) : -1;
        const int64_t limit = byGranule < byData ? byGranule : byData;
        if (limit >= 0 && pos <= uint64_t(limit)) gFast(T, data, limit, plan, is, st);
        // 2η φάση: αν αυτό που σταμάτησε την 1η ήταν το τέλος των ΔΕΔΟΜΕΝΩΝ και όχι του granule (το granule τελειώνει στα
        // τελευταία 63 μπιτ: συμβαίνει στο τελευταίο granule κάθε καρέ όταν η αποθήκη μπιτ είναι άδεια), τα λίγα bytes που
        // μένουν αντιγράφονται σε δικό μας χώρο 32 bytes, με μηδέν σε κάθε μπιτ από το endBits και πέρα και μηδενικά bytes
        // μετά. Εκεί η ροή ΕΙΝΑΙ αυτό που ορίζει το exact::Bits (μηδενικά μετά το τέλος) και το περιθώριο είναι δικό μας,
        // οπότε ο ίδιος γρήγορος δρόμος συνεχίζει ως το gend. Οι θέσεις μετριούνται από το πρώτο byte του αντιγράφου.
        // Γιατί μένουν το πολύ 16 bytes: η 1η φάση σταμάτησε είτε πριν ξεκινήσει (pos > όριο: κάτω από 8 bytes ως το
        // τέλος), είτε σε βήμα που πέρασε το όριο (ένα βήμα είναι ως 47 μπιτ: κάτω από 14 bytes). Αν σταμάτησε για άλλο
        // λόγο (πίνακας που δεν χρησιμοποιείται, 576 γραμμές γεμάτες) μπορεί να μένουν περισσότερα: τότε δεν υπάρχει
        // τίποτα να κερδηθεί και η φάση παραλείπεται.
        const uint64_t firstByte = uint64_t(st.pos >> 3);
        if (byGranule > byData && st.pos < gend && firstByte <= nBytes && nBytes - firstByte <= 16u) {
            uint8_t tail[32] = {};
            const size_t n = size_t(nBytes - firstByte);
            if (n) std::memcpy(tail, data + firstByte, n);
            if (n && (endBits & 7u) && !HUFF_MUT(18))   // ΒΛΑΒΗ 18: το αντίγραφο κρατά τα μπιτ του τελευταίου byte μετά το endBits
                tail[n - 1] = uint8_t(tail[n - 1] & (0xFFu << (8u - unsigned(endBits & 7u))));
            // Όριο: το gend, και οπωσδήποτε ως εκεί που το γέμισμα μένει μέσα στα 32 bytes (8·(32 − 8) − 63 = 129 μπιτ·
            // με gend ≤ endBits το gend είναι ≤ 8·16 = 128, άρα δεν κόβεται ποτέ από αυτό).
            const uint64_t gendHere = uint64_t(gend) - firstByte * 8u;
            const int64_t limitHere = gendHere < 129u ? int64_t(gendHere) : 129;
            State here{st.l, 0, size_t(st.pos & 7u), st.nLin};
            gFast(T, tail, limitHere, plan, is, here);
            st.l = here.l;
            st.region = here.region;
            st.pos = size_t(firstByte * 8u) + here.pos;
            st.nLin = here.nLin;
        }
    }
    return finish(T, data, endBits, gend, plan, is, st, pos, notes);
}

const char *huffmanKernelText()
{
    if (gCarefulOnly) return "huffman: careful path only (test switch)";
#if HUFF_X86
    return gFast == fastBmi2 ? "huffman: two-pair tables, 64-bit window · shlx/shrx" : "huffman: two-pair tables, 64-bit window · plain shifts";
#else
    return "huffman: two-pair tables, 64-bit window · plain shifts";
#endif
}

void huffmanUseBasicOnly(bool on)
{
#if HUFF_X86
    gFast = (!on && cpuHasBmi2()) ? fastBmi2 : fastPlain;
#else
    (void)on;
#endif
}

void huffmanUseCarefulOnly(bool on) { gCarefulOnly = on; }

size_t huffmanTableBytes()
{
    const Tables &T = tables();
    return size_t(T.gUsed) * sizeof(uint32_t) + size_t(T.directs) * sizeof T.d[0] + sizeof T.quad + sizeof T.quadVal;
}

bool huffmanSelfCheck(char *text, size_t textSize)
{
    const Tables &T = tables();
    if (text && textSize) {
        size_t at = 0;
        text[0] = 0;
        for (int s = 0; s < T.sets && at + 1 < textSize; ++s) {
            const Set &st = T.set[s];
            char dText[64] = "no D";
            if (st.d >= 0) std::snprintf(dText, sizeof dText, "D 2^%u = %u bytes", kDBits, unsigned(sizeof T.d[0]));
            const int n = std::snprintf(text + at, textSize - at,
                                        "  code set of table %2d: G root 2^%d + %u in subtables = %u bytes · %s (%d%% of pairs fit in %u bits)\n",
                                        st.firstSel, st.gBits, st.gSize - (1u << st.gBits), st.gSize * 4u, dText, st.directPercent, kDBits);
            if (n < 0 || size_t(n) >= textSize - at) break;
            at += size_t(n);
        }
    }
    return T.ok;
}

#ifdef MP3HUFF_COVERAGE
void huffmanCoverageReset() { std::memset(&gCov, 0, sizeof gCov); }

size_t huffmanCoverage(HuffmanCoverageRow *rows, size_t max)
{
    const Tables &T = tables();
    size_t n = 0;
    const auto add = [&](const char *what, int sel, const uint32_t *uses, uint32_t count) {
        if (n == max) return;
        HuffmanCoverageRow &r = rows[n++];
        if (sel >= 0)
            std::snprintf(r.name, sizeof r.name, "table %2d: %s", sel, what);
        else
            std::snprintf(r.name, sizeof r.name, "%s", what);
        r.entries = count;
        r.used = 0;
        r.minUses = count ? UINT32_MAX : 0;
        r.uses = 0;
        r.firstUnused = count;
        for (uint32_t i = 0; i < count; ++i) {
            r.uses += uses[i];
            r.used += uses[i] != 0;
            if (!uses[i] && r.firstUnused == count) r.firstUnused = i;
            if (uses[i] < r.minUses) r.minUses = uses[i];
        }
    };
    for (int s = 0; s < T.sets; ++s) {
        const Set &st = T.set[s];
        if (st.d >= 0) add("D, fast path", st.firstSel, gCov.d[st.d], 1u << kDBits);
        add("G, fast path", st.firstSel, gCov.gFast + st.gBase, st.gSize);
        add("G, careful path", st.firstSel, gCov.gCareful + st.gBase, st.gSize);
    }
    add("quadruples A, fast path", -1, gCov.quadFast[0], 64);
    add("quadruples B, fast path", -1, gCov.quadFast[1], 64);
    add("quadruple values, fast path", -1, gCov.quadValFast, 256);
    add("quadruples A, careful path", -1, gCov.quadCareful[0], 64);
    add("quadruples B, careful path", -1, gCov.quadCareful[1], 64);
    add("quadruple values, careful path", -1, gCov.quadValCareful, 256);
    return n;
}
#endif

}  // namespace orpheic::mp3
