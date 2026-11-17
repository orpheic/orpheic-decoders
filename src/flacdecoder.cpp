// Ο δικός μας αποκωδικοποιητής FLAC από το πρότυπο (RFC 9639). Τι είναι και τι υπόσχεται: flacdecoder.h.
// Οι αριθμοί ενοτήτων στα σχόλια είναι του RFC 9639.
#include "flacdecoder.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <type_traits>

// Άλλα OS: ο ίδιος κώδικας παντού· ανά μεταγλωττιστή αλλάζει μόνο το όνομα λίγων εντολών του επεξεργαστή (αντιστροφή
// bytes, μέτρημα μηδενικών μπροστά, ερώτηση «τι ξέρει αυτός ο επεξεργαστής»). Όλοι οι στόχοι μας (x86-64, ARM64) γράφουν το
// μικρό άκρο πρώτο.
#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#include <stdlib.h>
#define FLAC_LIKELY(x) (x)
#define FLAC_UNLIKELY(x) (x)
#define FLAC_INLINE __forceinline
#define FLAC_NOINLINE __declspec(noinline)
#define FLAC_TARGET_BMI2
#else
#define FLAC_LIKELY(x) __builtin_expect(!!(x), 1)
#define FLAC_UNLIKELY(x) __builtin_expect(!!(x), 0)
#define FLAC_INLINE inline __attribute__((always_inline))
#define FLAC_NOINLINE __attribute__((noinline))
#define FLAC_TARGET_BMI2 __attribute__((target("bmi2,lzcnt")))
#endif
#if defined(__x86_64__) || defined(_M_X64)
#define FLAC_X86 1
#if !defined(_MSC_VER) || defined(__clang__)
#include <cpuid.h>
#endif
#include <immintrin.h>
#else
#define FLAC_X86 0
#if defined(__aarch64__) && !defined(_MSC_VER)
#include <arm_acle.h>
#include <arm_neon.h>
#endif
#endif
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
#error "flacdecoder.cpp: written for little-endian processors"
#endif

namespace orpheic::flac {
namespace {

#if defined(_MSC_VER) && !defined(__clang__)
inline uint64_t bswap64(uint64_t v) { return _byteswap_uint64(v); }
inline int clz64(uint64_t v)   // v != 0
{
    unsigned long i;
    _BitScanReverse64(&i, v);
    return 63 - int(i);
}
#else
inline uint64_t bswap64(uint64_t v) { return __builtin_bswap64(v); }
inline int clz64(uint64_t v) { return __builtin_clzll(v); }   // v != 0
#endif

inline double nowNs() { return std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

inline uint64_t be64(const uint8_t *p)
{
    uint64_t v;
    std::memcpy(&v, p, 8);
    return bswap64(v);
}

// ── CRC (9.1.8, 9.3): και τα δύο αρχίζουν από 0, το πιο σημαντικό μπιτ πρώτο ──
// Το CRC-16 με πολλούς πίνακες (έναν για κάθε byte από όσα διαβάζονται μαζί) είναι δημοσιευμένη μέθοδος, το «slicing-by-8»:
// M. E. Kounavis και F. L. Berry, «A Systematic Approach to Building High Performance Software-Based CRC Generators», IEEE ISCC 2005,
// και «Novel Table Lookup-Based Algorithms for High-Performance CRC Generation», IEEE Transactions on Computers 57(11), 2008.
// Τα 16 bytes τη φορά εδώ είναι η ίδια ιδέα με διπλάσιους πίνακες.
struct CrcTables {
    uint8_t c8[256];
    uint16_t c16[16][256];   // c16[k][v] = το CRC του byte v ακολουθούμενου από k μηδενικά bytes
    CrcTables()
    {
        for (int i = 0; i < 256; ++i) {
            uint8_t a = uint8_t(i);
            uint16_t b = uint16_t(i << 8);
            for (int bit = 0; bit < 8; ++bit) {
                a = uint8_t((a << 1) ^ ((a & 0x80) ? 0x07 : 0));       // x^8 + x^2 + x + 1
                b = uint16_t((b << 1) ^ ((b & 0x8000) ? 0x8005 : 0));  // x^16 + x^15 + x^2 + 1
            }
            c8[i] = a;
            c16[0][i] = b;
        }
        for (int k = 1; k < 16; ++k)
            for (int i = 0; i < 256; ++i) {
                const uint16_t c = c16[k - 1][i];
                c16[k][i] = uint16_t((c << 8) ^ c16[0][c >> 8]);
            }
    }
};
const CrcTables kCrc;

inline uint8_t crc8(const uint8_t *p, size_t n)
{
    uint8_t c = 0;
    for (size_t i = 0; i < n; ++i) c = kCrc.c8[c ^ p[i]];
    return c;
}

uint16_t crc16Table(const uint8_t *p, size_t n, uint16_t c = 0)
{
    const auto &t = kCrc.c16;
    // 16 bytes μαζί: το CRC είναι γραμμικό, κάθε byte συνεισφέρει ό,τι θα έδινε μόνο του με τα μηδενικά πίσω του. Τα bytes
    // διαβάζονται ως δύο λέξεις των 64 μπιτ (2 διαβάσματα αντί για 16): σε παλιό επεξεργαστή με μία θύρα ανάγνωσης, τα
    // διαβάσματα είναι αυτά που ορίζουν τον χρόνο. Το πρώτο byte της ροής είναι το κάτω byte της λέξης (μικρό άκρο πρώτο).
    while (n >= 16) {
        uint64_t a, b;
        std::memcpy(&a, p, 8);
        std::memcpy(&b, p + 8, 8);
        a ^= uint64_t(uint16_t((c >> 8) | (c << 8)));
        c = uint16_t(t[15][a & 0xFF] ^ t[14][(a >> 8) & 0xFF] ^ t[13][(a >> 16) & 0xFF] ^ t[12][(a >> 24) & 0xFF] ^ t[11][(a >> 32) & 0xFF]
                     ^ t[10][(a >> 40) & 0xFF] ^ t[9][(a >> 48) & 0xFF] ^ t[8][a >> 56] ^ t[7][b & 0xFF] ^ t[6][(b >> 8) & 0xFF]
                     ^ t[5][(b >> 16) & 0xFF] ^ t[4][(b >> 24) & 0xFF] ^ t[3][(b >> 32) & 0xFF] ^ t[2][(b >> 40) & 0xFF] ^ t[1][(b >> 48) & 0xFF]
                     ^ t[0][b >> 56]);
        p += 16;
        n -= 16;
    }
    while (n >= 8) {
        c = uint16_t(t[7][(c >> 8) ^ p[0]] ^ t[6][(c & 0xFF) ^ p[1]] ^ t[5][p[2]] ^ t[4][p[3]] ^ t[3][p[4]] ^ t[2][p[5]]
                     ^ t[1][p[6]] ^ t[0][p[7]]);
        p += 8;
        n -= 8;
    }
    while (n--) c = uint16_t((c << 8) ^ t[0][(c >> 8) ^ *p++]);
    return c;
}

// Το ίδιο CRC-16 με πολλαπλασιασμό πολυωνύμων χωρίς κρατούμενο (PCLMULQDQ στον x86, PMULL στον ARM64), 64 bytes τη φορά.
// Η ιδέα: τα bytes του καρέ είναι ένα πολυώνυμο M(x) και το CRC είναι M(x)·x^16 mod P. Ένα κομμάτι S των 128 μπιτ που το
// ακολουθούν 128·m μπιτ ακόμη αξίζει S·x^(128m) = S_πάνω·x^(128m+64) + S_κάτω·x^(128m), και κάθε x^k μπορεί να αντικατασταθεί
// από το x^k mod P (16 μπιτ). Έτσι τέσσερα κομμάτια «διπλώνονται» πάνω στα επόμενα 64 bytes με δύο πολλαπλασιασμούς το
// καθένα, και στο τέλος ό,τι έμεινε (16 bytes που είναι ισοδύναμα mod P με όλη την αρχή, και η ουρά) περνά από τον πίνακα.
// Οι σταθερές x^k mod P υπολογίστηκαν με το χέρι του προγράμματος (βλ. kFold) και η μέθοδος ελέγχεται στην εκκίνηση
// απέναντι στον πίνακα· αν δεν συμφωνούσε, μένει ο πίνακας.
// Η μέθοδος είναι δημοσιευμένη: V. Gopal, E. Ozturk, J. Guilford, G. Wolrich, W. Feghali, M. Dixon και D. Karakoyunlu, «Fast CRC
// Computation for Generic Polynomials Using PCLMULQDQ Instruction», Intel, Δεκέμβριος 2009.
constexpr uint64_t kFold[8] = {0x0106, 0x1666, 0x8011, 0x9323, 0x926F, 0x4AAA, 0x8107, 0x1446};   // x^128, 192, 256, 320, 384, 448, 512, 576 mod P

#if FLAC_X86
#if defined(_MSC_VER) && !defined(__clang__)
#define FLAC_TARGET_CLMUL
#else
#define FLAC_TARGET_CLMUL __attribute__((target("pclmul,ssse3")))
#endif
#define FLAC_FOLD(s, k) _mm_xor_si128(_mm_clmulepi64_si128((s), (k), 0x11), _mm_clmulepi64_si128((s), (k), 0x00))
FLAC_NOINLINE FLAC_TARGET_CLMUL uint16_t crc16Fold(const uint8_t *p, size_t n)
{
    if (n < 64) return crc16Table(p, n);
    const __m128i rev = _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);   // το πρώτο byte γίνεται το πιο σημαντικό
    const __m128i k1 = _mm_set_epi64x(int64_t(kFold[1]), int64_t(kFold[0]));
    const __m128i k2 = _mm_set_epi64x(int64_t(kFold[3]), int64_t(kFold[2]));
    const __m128i k3 = _mm_set_epi64x(int64_t(kFold[5]), int64_t(kFold[4]));
    const __m128i k4 = _mm_set_epi64x(int64_t(kFold[7]), int64_t(kFold[6]));
#define FLAC_LOAD(q) _mm_shuffle_epi8(_mm_loadu_si128(reinterpret_cast<const __m128i *>(q)), rev)
    __m128i s0 = FLAC_LOAD(p), s1 = FLAC_LOAD(p + 16), s2 = FLAC_LOAD(p + 32), s3 = FLAC_LOAD(p + 48);
    p += 64;
    n -= 64;
    while (n >= 64) {
        s0 = _mm_xor_si128(FLAC_FOLD(s0, k4), FLAC_LOAD(p));
        s1 = _mm_xor_si128(FLAC_FOLD(s1, k4), FLAC_LOAD(p + 16));
        s2 = _mm_xor_si128(FLAC_FOLD(s2, k4), FLAC_LOAD(p + 32));
        s3 = _mm_xor_si128(FLAC_FOLD(s3, k4), FLAC_LOAD(p + 48));
        p += 64;
        n -= 64;
    }
    __m128i s = _mm_xor_si128(_mm_xor_si128(s3, FLAC_FOLD(s2, k1)), _mm_xor_si128(FLAC_FOLD(s1, k2), FLAC_FOLD(s0, k3)));
    while (n >= 16) {
        s = _mm_xor_si128(FLAC_FOLD(s, k1), FLAC_LOAD(p));
        p += 16;
        n -= 16;
    }
    uint8_t head[16];
    _mm_storeu_si128(reinterpret_cast<__m128i *>(head), _mm_shuffle_epi8(s, rev));
#undef FLAC_LOAD
    return crc16Table(p, n, crc16Table(head, 16));
}
#undef FLAC_FOLD
bool cpuHasClmul()
{
#if defined(_MSC_VER) && !defined(__clang__)
    int r[4];
    __cpuid(r, 1);
    return (r[2] & (1 << 1)) && (r[2] & (1 << 9));
#else
    unsigned a = 0, b = 0, c = 0, d = 0;
    return __get_cpuid(1, &a, &b, &c, &d) && (c & (1u << 1)) && (c & (1u << 9));
#endif
}
#define FLAC_HAVE_FOLD 1
#elif defined(__aarch64__) && defined(__ARM_FEATURE_AES) && !defined(_MSC_VER)
inline uint8x16_t revBytes(uint8x16_t v)
{
    v = vrev64q_u8(v);
    return vextq_u8(v, v, 8);
}
inline uint8x16_t foldNeon(uint8x16_t s, uint64_t kHigh, uint64_t kLow)
{
    const uint64x2_t q = vreinterpretq_u64_u8(s);
    const poly128_t a = vmull_p64(poly64_t(vgetq_lane_u64(q, 1)), poly64_t(kHigh));
    const poly128_t b = vmull_p64(poly64_t(vgetq_lane_u64(q, 0)), poly64_t(kLow));
    return veorq_u8(vreinterpretq_u8_p128(a), vreinterpretq_u8_p128(b));
}
FLAC_NOINLINE uint16_t crc16Fold(const uint8_t *p, size_t n)
{
    if (n < 64) return crc16Table(p, n);
    uint8x16_t s0 = revBytes(vld1q_u8(p)), s1 = revBytes(vld1q_u8(p + 16)), s2 = revBytes(vld1q_u8(p + 32)), s3 = revBytes(vld1q_u8(p + 48));
    p += 64;
    n -= 64;
    while (n >= 64) {
        s0 = veorq_u8(foldNeon(s0, kFold[7], kFold[6]), revBytes(vld1q_u8(p)));
        s1 = veorq_u8(foldNeon(s1, kFold[7], kFold[6]), revBytes(vld1q_u8(p + 16)));
        s2 = veorq_u8(foldNeon(s2, kFold[7], kFold[6]), revBytes(vld1q_u8(p + 32)));
        s3 = veorq_u8(foldNeon(s3, kFold[7], kFold[6]), revBytes(vld1q_u8(p + 48)));
        p += 64;
        n -= 64;
    }
    uint8x16_t s = veorq_u8(veorq_u8(s3, foldNeon(s2, kFold[1], kFold[0])), veorq_u8(foldNeon(s1, kFold[3], kFold[2]), foldNeon(s0, kFold[5], kFold[4])));
    while (n >= 16) {
        s = veorq_u8(foldNeon(s, kFold[1], kFold[0]), revBytes(vld1q_u8(p)));
        p += 16;
        n -= 16;
    }
    uint8_t head[16];
    vst1q_u8(head, revBytes(s));
    return crc16Table(p, n, crc16Table(head, 16));
}
inline bool cpuHasClmul() { return true; }   // κάθε ARM64 που χτίζεται με αυτή τη δυνατότητα την έχει
#define FLAC_HAVE_FOLD 1
#else
#define FLAC_HAVE_FOLD 0
#endif

using Crc16Fn = uint16_t (*)(const uint8_t *, size_t);
uint16_t crc16Plain(const uint8_t *p, size_t n) { return crc16Table(p, n); }

// Διαλέγει τον γρήγορο δρόμο μόνο αν ο επεξεργαστής τον έχει ΚΑΙ δίνει τα ίδια με τον πίνακα σε μια δοκιμή με κάθε μήκος
// από 0 ως 300 bytes και ένα μεγάλο.
Crc16Fn pickCrc16()
{
#if FLAC_HAVE_FOLD
    if (!cpuHasClmul()) return crc16Plain;
    uint8_t data[5000];
    uint64_t x = 0x9E3779B97F4A7C15ULL;
    for (uint8_t &b : data) {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        b = uint8_t(x >> 32);
    }
    for (size_t n = 0; n <= 300; ++n)
        if (crc16Fold(data + (n & 7), n) != crc16Table(data + (n & 7), n)) return crc16Plain;
    if (crc16Fold(data, sizeof data) != crc16Table(data, sizeof data)) return crc16Plain;
    return crc16Fold;
#else
    return crc16Plain;
#endif
}
Crc16Fn kCrc16 = pickCrc16();
inline uint16_t crc16(const uint8_t *p, size_t n) { return kCrc16(p, n); }

// ── Ο ΑΝΑΓΝΩΣΤΗΣ ΜΠΙΤ ──
// `buf`: τα επόμενα μπιτ της ροής, το πρώτο στο πιο σημαντικό μπιτ. Τα πάνω `cnt` είναι «δικά μας» (έχουν μετρηθεί στο
// `pos`)· ό,τι βρίσκεται πιο κάτω είναι τα αμέσως επόμενα μπιτ της ροής, που θα ξαναμπούν ίδια στο επόμενο γέμισμα.
// Το γέμισμα διαβάζει πάντα 8 bytes από τη θέση `pos` και προχωρά όσα ολόκληρα bytes χωράνε: μετά από αυτό cnt ≥ 56.
// Πέρα από το τέλος δεν διαβάζει: μπαίνουν μηδενικά, το `pos` συνεχίζει να μετρά, και το overrun() λέει ότι ξεπεράστηκε.
// Το γέμισμα χωρίς διακλάδωση (τρεις πράξεις: «ή» με τα επόμενα 8 bytes, η θέση προχωρά (63 − cnt)/8 bytes, cnt |= 56) είναι
// δημοσιευμένη τεχνική: Fabian Giesen, «Reading bits in far too many ways (part 2)», ιστολόγιο «The ryg blog», Φεβρουάριος 2018.
struct BitReader {
    const uint8_t *base;
    size_t n;        // bytes δεδομένων· διαβάζονται ως και base[n + 7]
    size_t pos = 0;  // bytes που έχουν μπει στο buf
    uint64_t buf = 0;
    int cnt = 0;

    BitReader(const uint8_t *p, size_t size) : base(p), n(size) {}

    inline void refill()
    {
        if (FLAC_LIKELY(pos <= n)) buf |= be64(base + pos) >> cnt;
        pos += size_t(63 - cnt) >> 3;
        cnt |= 56;
    }
    inline uint64_t read(int bits)   // 0..56
    {
        if (!bits) return 0;
        refill();
        const uint64_t v = buf >> (64 - bits);
        buf <<= bits;
        cnt -= bits;
        return v;
    }
    inline int64_t readSigned(int bits)   // 0..56, συμπλήρωμα ως προς δύο
    {
        if (!bits) return 0;
        refill();
        const int64_t v = int64_t(buf) >> (64 - bits);
        buf <<= bits;
        cnt -= bits;
        return v;
    }
    bool overrun() const { return pos * 8 - size_t(cnt) > n * 8; }
    size_t bitsUsed() const { return pos * 8 - size_t(cnt); }

    // Μηδενικά ως το πρώτο 1 (που καταναλώνεται). false = τα δεδομένα τελείωσαν.
    bool unary(uint64_t &zeros)
    {
        zeros = 0;
        for (;;) {
            refill();
            if (buf >> (64 - cnt)) {   // υπάρχει 1 μέσα στα δικά μας μπιτ
                const int z = clz64(buf);
                zeros += uint64_t(z);
                buf <<= z + 1;
                cnt -= z + 1;
                return true;
            }
            zeros += uint64_t(cnt);
            buf <<= cnt;   // cnt ≤ 63
            cnt = 0;
            if (overrun()) return false;
        }
    }
};

inline int32_t unfold(uint32_t u) { return int32_t((u >> 1) ^ (0u - (u & 1))); }   // 9.2.7.2: ζιγκ-ζαγκ

// Ένα σύμβολο Rice από τον προσεκτικό δρόμο (ελέγχει το τέλος): για μακριές σειρές μηδενικών και για την άκρη των δεδομένων.
FLAC_NOINLINE bool riceSlow(BitReader &br, unsigned k, uint32_t &u)
{
    uint64_t q;
    if (!br.unary(q)) return false;
    u = uint32_t(q << k) | uint32_t(br.read(int(k)));
    return true;
}

// Ο πυρήνας Rice (flacdecoder_rice.inc) σε τρεις μορφές. Σε x86-64 υπάρχουν δύο για τα 32 μπιτ: η βασική
// (ολίσθηση με μετρητή στον CL, BSR) και, όπου ο επεξεργαστής τις έχει (από το 2013 και μετά), με SHLX/SHRX/LZCNT, που
// είναι μία μικροεντολή η καθεμία· διαλέγεται μία φορά στην εκκίνηση. Στον ARM64 η βασική είναι ήδη μία εντολή η καθεμία.
// Το ξεδίπλωμα του προσήμου (9.2.7.2, ζιγκ-ζαγκ) για ένα ολόκληρο τμήμα, επιτόπου, 4 τιμές τη φορά όπου γίνεται.
FLAC_NOINLINE void unfoldRun(int32_t *p, uint32_t n)
{
    uint32_t i = 0;
#if FLAC_X86
    const __m128i one = _mm_set1_epi32(1), zero = _mm_setzero_si128();
    for (; i + 4 <= n; i += 4) {
        const __m128i x = _mm_loadu_si128(reinterpret_cast<const __m128i *>(p + i));
        const __m128i sign = _mm_sub_epi32(zero, _mm_and_si128(x, one));
        _mm_storeu_si128(reinterpret_cast<__m128i *>(p + i), _mm_xor_si128(_mm_srli_epi32(x, 1), sign));
    }
#elif defined(__aarch64__) && !defined(_MSC_VER)
    for (; i + 4 <= n; i += 4) {
        const uint32x4_t x = vreinterpretq_u32_s32(vld1q_s32(p + i));
        const uint32x4_t sign = vreinterpretq_u32_s32(vnegq_s32(vreinterpretq_s32_u32(vandq_u32(x, vdupq_n_u32(1)))));
        vst1q_s32(p + i, vreinterpretq_s32_u32(veorq_u32(vshrq_n_u32(x, 1), sign)));
    }
#endif
    for (; i < n; ++i) p[i] = unfold(uint32_t(p[i]));
}
inline void unfoldRun(int64_t *p, uint32_t n)
{
    for (uint32_t i = 0; i < n; ++i) p[i] = unfold(uint32_t(uint64_t(p[i])));
}

// «Πόσα μηδενικά έχει μπροστά». Στον x86-64 χωρίς LZCNT η εντολή (BSR) δεν ορίζεται στο μηδέν: το μηδέν ελέγχεται χωριστά,
// έξω από την αλυσίδα. Στον ARM64 και με LZCNT το μηδέν δίνει 64, που ο πυρήνας το βλέπει ως «δεν χωράει».
#if FLAC_X86 && defined(_MSC_VER) && !defined(__clang__)
#define FLAC_LEAD_PLAIN(v, z, onZero)                              \
    {                                                              \
        unsigned long top_;                                        \
        if (FLAC_UNLIKELY(!_BitScanReverse64(&top_, (v)))) { onZero; } \
        z = 63u - unsigned(top_);                                  \
    }
#elif FLAC_X86
#define FLAC_LEAD_PLAIN(v, z, onZero)                              \
    {                                                              \
        if (FLAC_UNLIKELY((v) == 0)) { onZero; }                   \
        z = unsigned(__builtin_clzll(v));                          \
    }
#elif defined(__aarch64__) && !defined(_MSC_VER)
#define FLAC_LEAD_PLAIN(v, z, onZero) { z = unsigned(__clzll(v)); }
#else
#define FLAC_LEAD_PLAIN(v, z, onZero)                              \
    {                                                              \
        if (FLAC_UNLIKELY((v) == 0)) { onZero; }                   \
        z = unsigned(clz64(v));                                    \
    }
#endif

#define FLAC_RICE_ATTR
#define FLAC_RICE_LEAD(v, z, onZero) FLAC_LEAD_PLAIN(v, z, onZero)
#define FLAC_RICE_SHL(v, n) ((v) << (n))
#define FLAC_RICE_SHR(v, n) ((v) >> (n))
#define FLAC_RICE_UNFOLD(p, count) unfoldRun(p, count)
#if !FLAC_X86   // στον x86-64 τη θέση της παίρνουν οι δύο μορφές πιο κάτω
#define FLAC_RICE_NAME riceRunPlain
#define FLAC_RICE_T int32_t
#include "flacdecoder_rice.inc"
#undef FLAC_RICE_NAME
#undef FLAC_RICE_T
#endif

#define FLAC_RICE_NAME riceRunWide
#define FLAC_RICE_T int64_t
#include "flacdecoder_rice.inc"
#undef FLAC_RICE_NAME
#undef FLAC_RICE_ATTR
#undef FLAC_RICE_T
#undef FLAC_RICE_LEAD
#undef FLAC_RICE_SHL
#undef FLAC_RICE_SHR
#undef FLAC_RICE_UNFOLD

#if FLAC_X86
// Η βασική μορφή στον x86-64, μία συνάρτηση για κάθε k (0..30): η ολίσθηση με μεταβλητό μετρητή κοστίζει εκεί δύο-τρεις
// μικροεντολές και δένει τον καταχωρητή CL, ενώ με σταθερό μετρητή είναι μία.
#define FLAC_RICE_NAME riceRunPlainK
#define FLAC_RICE_ATTR
#define FLAC_RICE_T int32_t
#if defined(_MSC_VER) && !defined(__clang__)
#define FLAC_RICE_TOP(v, top, onZero)                                  \
    {                                                                  \
        unsigned long top_;                                            \
        if (FLAC_UNLIKELY(!_BitScanReverse64(&top_, (v)))) { onZero; } \
        top = unsigned(top_);                                          \
    }
#else
#define FLAC_RICE_TOP(v, top, onZero)                                  \
    {                                                                  \
        if (FLAC_UNLIKELY((v) == 0)) { onZero; }                       \
        top = 63u - unsigned(__builtin_clzll(v));                      \
    }
#endif
#define FLAC_RICE_SHL(v, n) ((v) << (n))
#define FLAC_RICE_SHR(v, n) ((v) >> (n))
#define FLAC_RICE_UNFOLD(p, count) unfoldRun(p, count)
#define FLAC_RICE_KCONST
#include "flacdecoder_rice.inc"
#undef FLAC_RICE_KCONST
#undef FLAC_RICE_NAME
#undef FLAC_RICE_ATTR
#undef FLAC_RICE_T
#undef FLAC_RICE_TOP
#undef FLAC_RICE_SHL
#undef FLAC_RICE_SHR
#undef FLAC_RICE_UNFOLD
using RiceKFn = bool (*)(BitReader &, int32_t *, uint32_t, uint64_t *);
#define FLAC_K(n) riceRunPlainK<n>
const RiceKFn kRicePlainK[31] = {FLAC_K(0),  FLAC_K(1),  FLAC_K(2),  FLAC_K(3),  FLAC_K(4),  FLAC_K(5),  FLAC_K(6),  FLAC_K(7),
                                 FLAC_K(8),  FLAC_K(9),  FLAC_K(10), FLAC_K(11), FLAC_K(12), FLAC_K(13), FLAC_K(14), FLAC_K(15),
                                 FLAC_K(16), FLAC_K(17), FLAC_K(18), FLAC_K(19), FLAC_K(20), FLAC_K(21), FLAC_K(22), FLAC_K(23),
                                 FLAC_K(24), FLAC_K(25), FLAC_K(26), FLAC_K(27), FLAC_K(28), FLAC_K(29), FLAC_K(30)};
#undef FLAC_K
bool riceRunBasic(BitReader &br, int32_t *out, uint32_t count, unsigned k, uint64_t *slow) { return kRicePlainK[k](br, out, count, slow); }

#define FLAC_RICE_NAME riceRunBmi2
#define FLAC_RICE_ATTR FLAC_TARGET_BMI2
#define FLAC_RICE_T int32_t
#define FLAC_RICE_LEAD(v, z, onZero) { z = unsigned(_lzcnt_u64(v)); }
#define FLAC_RICE_UNFOLD(p, count) unfoldRun(p, count)
#if defined(_MSC_VER) && !defined(__clang__)
#define FLAC_RICE_SHL(v, n) _shlx_u64((v), unsigned(n))
#define FLAC_RICE_SHR(v, n) _shrx_u64((v), unsigned(n))
#else   // μέσα σε συνάρτηση με target("bmi2,lzcnt") ο μεταγλωττιστής βγάζει μόνος του SHLX/SHRX
#define FLAC_RICE_SHL(v, n) ((v) << (n))
#define FLAC_RICE_SHR(v, n) ((v) >> (n))
#endif
#include "flacdecoder_rice.inc"
#undef FLAC_RICE_NAME
#undef FLAC_RICE_ATTR
#undef FLAC_RICE_T
#undef FLAC_RICE_LEAD
#undef FLAC_RICE_SHL
#undef FLAC_RICE_SHR
#undef FLAC_RICE_UNFOLD

bool cpuHasBmi2()
{
#if defined(_MSC_VER) && !defined(__clang__)
    int r[4];
    __cpuid(r, 0);
    if (r[0] < 7) return false;
    __cpuidex(r, 7, 0);
    const bool bmi2 = (r[1] & (1 << 8)) != 0;
    __cpuid(r, int(0x80000001));
    return bmi2 && (r[2] & (1 << 5)) != 0;
#else
    unsigned a = 0, b = 0, c = 0, d = 0;
    if (!__get_cpuid_count(7, 0, &a, &b, &c, &d) || !(b & (1u << 8))) return false;
    return __get_cpuid(0x80000001u, &a, &b, &c, &d) && (c & (1u << 5));
#endif
}
using RiceFn = bool (*)(BitReader &, int32_t *, uint32_t, unsigned, uint64_t *);
RiceFn kRiceRun = cpuHasBmi2() ? riceRunBmi2 : riceRunBasic;   // αλλάζει μόνο με το useBasicOnly (δοκιμή του βασικού δρόμου)
#else
constexpr auto kRiceRun = riceRunPlain;
#endif

inline bool riceRun(BitReader &br, int32_t *out, uint32_t count, int k, uint64_t *slow) { return kRiceRun(br, out, count, unsigned(k), slow); }
inline bool riceRun(BitReader &br, int64_t *out, uint32_t count, int k, uint64_t *slow) { return riceRunWide(br, out, count, unsigned(k), slow); }

// Το κωδικοποιημένο υπόλοιπο (9.2.7): τα δείγματα order .. n-1 του υποκαρέ.
template <class T>
Status residual(BitReader &br, T *out, uint32_t n, uint32_t order, Stats *st)
{
    const int method = int(br.read(2));
    if (method > 1) return Status::BadData;
    const int paramBits = method ? 5 : 4;
    const int escape = (1 << paramBits) - 1;
    const int po = int(br.read(4));
    const uint32_t psize = n >> po;
    if (po && (n & ((1u << po) - 1))) return Status::BadData;   // το μπλοκ δεν μοιράζεται σε ίσα τμήματα
    if (psize < order) return Status::BadData;                  // το πρώτο τμήμα δεν χωράει τα δείγματα προθέρμανσης
    if (st) ++st->partitionOrder[po];
    T *p = out + order;
    uint32_t count = psize - order;
    for (uint32_t part = 0; part < (1u << po); ++part, count = psize) {
        const int k = int(br.read(paramBits));
        if (k == escape) {   // 9.2.7.1: αυτούσια, με σταθερό πλάτος
            const int bits = int(br.read(5));
            for (uint32_t i = 0; i < count; ++i) p[i] = T(br.readSigned(bits));
            if (st) st->escapedSamples += count;
        } else {
            if (!riceRun(br, p, count, k, st ? &st->slowRice : nullptr)) return Status::Truncated;
            if (st) st->riceSamples[k] += count;
        }
        if (br.overrun()) return Status::Truncated;
        p += count;
    }
    return Status::Ok;
}

// Σταθεροί προβλέπτες (9.2.5), επιτόπου: στο s[] είναι τα δείγματα προθέρμανσης και μετά το υπόλοιπο. Modulo 2^πλάτος.
// Τα τελευταία δείγματα κρατιούνται σε μεταβλητές (καταχωρητές): αν το καθένα ξαναδιαβαζόταν από τη μνήμη αμέσως μετά το
// γράψιμό του, η αλυσίδα «δείγμα → επόμενο δείγμα» θα περνούσε κάθε φορά από τη μνήμη.
template <class T>
FLAC_NOINLINE void fixedPredict(T *s, uint32_t n, uint32_t order)
{
    using U = std::make_unsigned_t<T>;
    if (order == 0 || n <= order) return;
    U a = U(s[order - 1]);                      // s[i-1]
    U b = order > 1 ? U(s[order - 2]) : 0;      // s[i-2]
    U c = order > 2 ? U(s[order - 3]) : 0;
    U d = order > 3 ? U(s[order - 4]) : 0;
    switch (order) {
    case 1:
        for (uint32_t i = 1; i < n; ++i) { a += U(s[i]); s[i] = T(a); }
        break;
    case 2:
        for (uint32_t i = 2; i < n; ++i) { const U x = U(s[i]) + 2 * a - b; s[i] = T(x); b = a; a = x; }
        break;
    case 3:
        for (uint32_t i = 3; i < n; ++i) { const U x = U(s[i]) + 3 * a - 3 * b + c; s[i] = T(x); c = b; b = a; a = x; }
        break;
    default:
        for (uint32_t i = 4; i < n; ++i) { const U x = U(s[i]) + 4 * a - 6 * b + 4 * c - d; s[i] = T(x); d = c; c = b; b = a; a = x; }
        break;
    }
}

// Γραμμικός προβλέπτης (9.2.6), επιτόπου. A = ο τύπος του αθροίσματος (χωρίς πρόσημο, για να μην υπάρχει απροσδιόριστη
// υπερχείλιση σε χαλασμένα δεδομένα): uint32_t όταν το πρότυπο εγγυάται ότι χωράει, αλλιώς uint64_t. Τα δύο πιο πρόσφατα
// δείγματα κρατιούνται σε καταχωρητές (αυτά ορίζουν πόσο γρήγορα βγαίνει το επόμενο)· τα παλιότερα διαβάζονται από τη
// μνήμη, όπου έχουν προλάβει να γραφτούν. Η τάξη είναι σταθερά του προτύπου (template): ο μεταγλωττιστής ξετυλίγει το άθροισμα.
template <int ORDER, class T, class A>
FLAC_NOINLINE void lpcOrder(T *s, uint32_t n, const int32_t *coef, int shift)
{
    using SA = std::make_signed_t<A>;
    using U = std::make_unsigned_t<T>;
    A c[ORDER];
    for (int j = 0; j < ORDER; ++j) c[j] = A(SA(coef[j]));
    A p1 = A(SA(s[ORDER - 1]));
    A p2 = ORDER > 1 ? A(SA(s[ORDER - 2])) : A(0);
    for (uint32_t i = ORDER; i < n; ++i) {
        A sum = 0;
        for (int j = 2; j < ORDER; ++j) sum += c[j] * A(SA(s[i - 1 - j]));
        if (ORDER > 1) sum += c[1] * p2;
        sum += c[0] * p1;
        const U x = U(s[i]) + U(SA(sum) >> shift);
        s[i] = T(x);
        p2 = p1;
        p1 = A(SA(T(x)));
    }
}

template <class T, class A>
FLAC_NOINLINE void lpcAny(T *s, uint32_t n, uint32_t order, const int32_t *coef, int shift)
{
    using SA = std::make_signed_t<A>;
    using U = std::make_unsigned_t<T>;
    for (uint32_t i = order; i < n; ++i) {
        A sum = 0;
        for (uint32_t j = 0; j < order; ++j) sum += A(SA(coef[j])) * A(SA(s[i - 1 - j]));
        s[i] = T(U(s[i]) + U(SA(sum) >> shift));
    }
}

template <class T, class A>
void lpcPredict(T *s, uint32_t n, uint32_t order, const int32_t *coef, int shift)
{
    switch (order) {
    case 1: lpcOrder<1, T, A>(s, n, coef, shift); break;
    case 2: lpcOrder<2, T, A>(s, n, coef, shift); break;
    case 3: lpcOrder<3, T, A>(s, n, coef, shift); break;
    case 4: lpcOrder<4, T, A>(s, n, coef, shift); break;
    case 5: lpcOrder<5, T, A>(s, n, coef, shift); break;
    case 6: lpcOrder<6, T, A>(s, n, coef, shift); break;
    case 7: lpcOrder<7, T, A>(s, n, coef, shift); break;
    case 8: lpcOrder<8, T, A>(s, n, coef, shift); break;
    case 9: lpcOrder<9, T, A>(s, n, coef, shift); break;
    case 10: lpcOrder<10, T, A>(s, n, coef, shift); break;
    case 11: lpcOrder<11, T, A>(s, n, coef, shift); break;
    case 12: lpcOrder<12, T, A>(s, n, coef, shift); break;
    default: lpcAny<T, A>(s, n, order, coef, shift); break;
    }
}

inline int ceilLog2(uint32_t v)   // v ≥ 1
{
    int r = 0;
    while ((1u << r) < v) ++r;
    return r;
}

// Τι βρέθηκε σε ένα υποκαρέ. Η πρόβλεψη γίνεται ΜΕΤΑ, όταν έχουν διαβαστεί όλα τα κανάλια του καρέ και έχει περάσει το
// CRC-16: έτσι δύο κανάλια προβλέπονται μαζί (δύο ανεξάρτητες αλυσίδες στον ίδιο βρόχο), και ένα χαλασμένο καρέ δεν κοστίζει.
struct Sub {
    enum Kind : uint8_t { Plain, Fixed, Lpc };   // Plain = σταθερό ή αυτούσιο: τα δείγματα είναι ήδη στη θέση τους
    Kind kind = Plain;
    uint8_t order = 0;
    uint8_t shift = 0;
    uint8_t wasted = 0;
    uint8_t bits = 0;       // το βάθος του υποκαρέ: του καρέ, +1 στο κανάλι διαφοράς, μείον τα αχρησιμοποίητα μπιτ
    bool wideSum = false;   // ο γραμμικός προβλέπτης θέλει άθροισμα 64 μπιτ (παράρτημα A.3)
    int32_t coef[32];
};

// Ένα υποκαρέ (9.2), χωρίς την πρόβλεψη: στο out[] μένουν τα δείγματα προθέρμανσης και μετά το υπόλοιπο.
// `bits` = το βάθος του καρέ, +1 αν είναι κανάλι διαφοράς.
template <class T>
Status parseSub(BitReader &br, T *out, uint32_t n, int bits, Sub &s, Stats *st)
{
    const uint32_t h = uint32_t(br.read(8));
    if (h & 0x80) return Status::BadData;
    const uint32_t type = (h >> 1) & 0x3F;
    s.kind = Sub::Plain;
    s.order = 0;
    s.wasted = 0;
    if (h & 1) {   // 9.2.2
        uint64_t z;
        if (!br.unary(z)) return Status::Truncated;
        if (z >= uint64_t(bits - 1)) return Status::BadData;   // δεν θα έμενε ούτε ένα μπιτ
        s.wasted = uint8_t(z + 1);
        bits -= int(s.wasted);
        if (st) ++st->wasted;
    }
    s.bits = uint8_t(bits);
    if (st) ++st->subframes;
    if (type == 0) {   // 9.2.3
        const T v = T(br.readSigned(bits));
        for (uint32_t i = 0; i < n; ++i) out[i] = v;
        if (st) ++st->constant;
    } else if (type == 1) {   // 9.2.4
        for (uint32_t i = 0; i < n; ++i) {
            out[i] = T(br.readSigned(bits));
            if ((i & 255) == 255 && br.overrun()) return Status::Truncated;
        }
        if (st) ++st->verbatim;
    } else if (type >= 8 && type <= 12) {   // 9.2.5
        const uint32_t order = type - 8;
        if (order > n) return Status::BadData;
        for (uint32_t i = 0; i < order; ++i) out[i] = T(br.readSigned(bits));
        const Status r = residual(br, out, n, order, st);
        if (r != Status::Ok) return r;
        s.kind = Sub::Fixed;
        s.order = uint8_t(order);
        if (st) ++st->fixed[order];
    } else if (type >= 32) {   // 9.2.6
        const uint32_t order = type - 31;
        if (order > n) return Status::BadData;
        for (uint32_t i = 0; i < order; ++i) out[i] = T(br.readSigned(bits));
        const int precision = int(br.read(4)) + 1;
        if (precision == 16) return Status::BadData;
        const int shift = int(br.readSigned(5));
        if (shift < 0) return Status::BadData;
        for (uint32_t i = 0; i < order; ++i) s.coef[i] = int32_t(br.readSigned(precision));
        const Status r = residual(br, out, n, order, st);
        if (r != Status::Ok) return r;
        s.kind = Sub::Lpc;
        s.order = uint8_t(order);
        s.shift = uint8_t(shift);
        s.wideSum = sizeof(T) > 4 || bits + precision + ceilLog2(order) > 32;   // παράρτημα A.3
        if (st) {
            ++st->lpc[order];
            ++st->precision[precision];
            ++st->shift[shift];
            if (s.wideSum) ++st->lpcWide;
        }
    } else {
        return Status::BadData;
    }
    return br.overrun() ? Status::Truncated : Status::Ok;
}

// Ο σταθερός προβλέπτης τάξης K είναι K τρέχοντα αθροίσματα το ένα πάνω στο άλλο: το υπόλοιπο είναι η K-οστή διαφορά των
// δειγμάτων (π.χ. τάξη 2: x[n] - 2x[n-1] + x[n-2] = (x[n] - x[n-1]) - (x[n-1] - x[n-2])), άρα αθροίζοντας K φορές
// ξαναβγαίνουν τα δείγματα. Ένα τρέχον άθροισμα γίνεται 4 δείγματα τη φορά (δύο ολισθήσεις θέσεων και δύο προσθέσεις), και
// μεταξύ των τετράδων περνά μόνο ένα «κρατούμενο» ανά επίπεδο. Modulo 2^32, όπως και ο απλός δρόμος: ίδια μπιτ.
// Αρχικές τιμές: οι διαφορές των δειγμάτων προθέρμανσης στο τελευταίο τους σημείο.
// Το τρέχον άθροισμα με ολισθήσεις που διπλασιάζονται (μία θέση, μετά δύο) είναι το γνωστό παράλληλο άθροισμα προθεμάτων:
// W. D. Hillis και G. L. Steele, «Data Parallel Algorithms», Communications of the ACM 29(12), 1986.
#if FLAC_X86
#define FLAC_SUMS_VEC __m128i
#define FLAC_SUMS_LOAD(p) _mm_loadu_si128(reinterpret_cast<const __m128i *>(p))
#define FLAC_SUMS_STORE(p, v) _mm_storeu_si128(reinterpret_cast<__m128i *>(p), v)
#define FLAC_SUMS_SET(x) _mm_set1_epi32(int32_t(x))
#define FLAC_SUMS_RUN(v) v = _mm_add_epi32(v, _mm_slli_si128(v, 4)); v = _mm_add_epi32(v, _mm_slli_si128(v, 8))
#define FLAC_SUMS_ADD(a, b) _mm_add_epi32(a, b)
#define FLAC_SUMS_LAST(v) _mm_shuffle_epi32(v, 0xFF)
#define FLAC_SUMS_FIRST(v) uint32_t(_mm_cvtsi128_si32(v))
#define FLAC_HAVE_SUMS 1
#elif defined(__aarch64__) && !defined(_MSC_VER)
#define FLAC_SUMS_VEC int32x4_t
#define FLAC_SUMS_LOAD(p) vld1q_s32(p)
#define FLAC_SUMS_STORE(p, v) vst1q_s32(p, v)
#define FLAC_SUMS_SET(x) vdupq_n_s32(int32_t(x))
#define FLAC_SUMS_RUN(v) v = vaddq_s32(v, vextq_s32(vdupq_n_s32(0), v, 3)); v = vaddq_s32(v, vextq_s32(vdupq_n_s32(0), v, 2))
#define FLAC_SUMS_ADD(a, b) vaddq_s32(a, b)
#define FLAC_SUMS_LAST(v) vdupq_laneq_s32(v, 3)
#define FLAC_SUMS_FIRST(v) uint32_t(vgetq_lane_s32(v, 0))
#define FLAC_HAVE_SUMS 1
#else
#define FLAC_HAVE_SUMS 0
#endif

#if FLAC_HAVE_SUMS
template <int K>
FLAC_NOINLINE void fixedSums(int32_t *s, uint32_t n)
{
    if (n <= uint32_t(K)) return;
    uint32_t w[4] = {};
    for (int j = 0; j < K; ++j) w[j] = uint32_t(s[K - 1 - j]);
    const uint32_t d[4] = {w[0], w[0] - w[1], w[0] - 2 * w[1] + w[2], w[0] - 3 * w[1] + 3 * w[2] - w[3]};
    FLAC_SUMS_VEC carry[K];
    for (int j = 0; j < K; ++j) carry[j] = FLAC_SUMS_SET(d[j]);
    uint32_t i = uint32_t(K);
    for (; i + 4 <= n; i += 4) {
        FLAC_SUMS_VEC v = FLAC_SUMS_LOAD(s + i);
        for (int level = K - 1; level >= 0; --level) {
            FLAC_SUMS_RUN(v);
            v = FLAC_SUMS_ADD(v, carry[level]);
            carry[level] = FLAC_SUMS_LAST(v);
        }
        FLAC_SUMS_STORE(s + i, v);
    }
    uint32_t c[4] = {};
    for (int j = 0; j < K; ++j) c[j] = FLAC_SUMS_FIRST(carry[j]);
    for (; i < n; ++i) {
        uint32_t v = uint32_t(s[i]);
        for (int level = K - 1; level >= 0; --level) {
            v += c[level];
            c[level] = v;
        }
        s[i] = int32_t(v);
    }
}
#endif

inline void fixedFast(int32_t *s, uint32_t n, uint32_t order)
{
#if FLAC_HAVE_SUMS
    switch (order) {
    case 1: fixedSums<1>(s, n); break;
    case 2: fixedSums<2>(s, n); break;
    case 3: fixedSums<3>(s, n); break;
    case 4: fixedSums<4>(s, n); break;
    default: break;
    }
#else
    fixedPredict(s, n, order);
#endif
}
inline void fixedFast(int64_t *s, uint32_t n, uint32_t order) { fixedPredict(s, n, order); }

// Η πρόβλεψη ενός υποκαρέ με τον απλό δρόμο (ένα δείγμα τη φορά), επιτόπου.
template <class T>
void predictScalar(T *s, uint32_t n, const Sub &sub)
{
    if (sub.kind == Sub::Fixed) {
        fixedFast(s, n, sub.order);
    } else if (sub.kind == Sub::Lpc) {
        if (sub.wideSum) lpcPredict<T, uint64_t>(s, n, sub.order, sub.coef, sub.shift);
        else lpcPredict<T, uint32_t>(s, n, sub.order, sub.coef, sub.shift);
    }
}

// Ένα δείγμα του γραμμικού προβλέπτη με ακριβές άθροισμα 64 μπιτ: για τις άκρες των γρήγορων δρόμων.
inline int32_t lpcOne(const int32_t *at, int32_t r, const int32_t *coef, uint32_t order, int shift)
{
    int64_t sum = 0;
    for (uint32_t j = 0; j < order; ++j) sum += int64_t(coef[j]) * int64_t(at[-1 - int(j)]);
    return int32_t(uint32_t(r) + uint32_t(uint64_t(sum >> shift)));
}

#if FLAC_X86
// ── Ο ΓΡΑΜΜΙΚΟΣ ΠΡΟΒΛΕΠΤΗΣ «ΣΕ ΔΑΧΤΥΛΙΔΙ» (x86-64) ──
// Ο απλός δρόμος κάνει `τάξη` πολλαπλασιασμούς ανά δείγμα, και ο επεξεργαστής έχει ΕΝΑΝ πολλαπλασιαστή ακεραίων: 8 ως 12
// κύκλοι ανά δείγμα, όσο κι αν περισσεύει ο υπόλοιπος. Εδώ το άθροισμα γυρίζει ανάποδα: αντί κάθε δείγμα να μαζεύει τα N
// προηγούμενα, κάθε δείγμα που βγαίνει ΜΟΙΡΑΖΕΙ τη συνεισφορά του (δείγμα × συντελεστής) στα N επόμενα, με μία ή δύο
// εντολές που πολλαπλασιάζουν 4 θέσεις μαζί. Οι N «υπό κατασκευή» προβλέψεις ζουν σε N θέσεις καταχωρητών που δεν
// μετακινούνται ποτέ: η θέση h ανήκει στα δείγματα h, h+N, h+2N… (δαχτυλίδι), και οι συντελεστές είναι γραμμένοι εκ των
// προτέρων σε N περιστροφές, μία για κάθε βήμα. Το υπόλοιπο του δείγματος προστίθεται στη θέση του ολισθημένο αριστερά,
// ώστε το δείγμα να είναι σκέτη ολίσθηση δεξιά: floor((άθροισμα + υπόλοιπο·2^s) / 2^s) = υπόλοιπο + floor(άθροισμα / 2^s).
// Η αλυσίδα «δείγμα → επόμενο δείγμα» είναι πιο μακριά από του απλού δρόμου (8 κύκλοι αντί 6), αλλά κάθε βήμα είναι λίγες
// εντολές που δεν περιμένουν τον πολλαπλασιαστή: δύο κανάλια τρέχουν πλεγμένα στον ίδιο βρόχο σχεδόν στον χρόνο του ενός.
//
// ΑΚΡΙΒΕΙΑ, μορφή 16 μπιτ (ήχος ως 16 μπιτ, SSE2). Το PMADDWD πολλαπλασιάζει τα κάτω 16 μπιτ του δείγματος: ακριβές όσο
// κάθε δείγμα χωράει σε 16 μπιτ. Οι θέσεις είναι 32 μπιτ και αθροίζουν modulo 2^32· το τελικό ποσό ενός δείγματος x είναι
// στο [x·2^s, (x+1)·2^s), άρα χωράει σε 32 μπιτ όσο το x χωράει σε 17 (s ≤ 15), όπως κι αν τυλίχτηκαν τα ενδιάμεσα. Το
// κανάλι διαφοράς μπορεί να θέλει 17 μπιτ: κάθε δείγμα που βγαίνει ελέγχεται, και αν έστω ένα δεν χωράει σε 16 το υποκαρέ
// ξαναγίνεται από την αρχή με τη μορφή 64 μπιτ ή με τον απλό δρόμο (το υπόλοιπο έχει μείνει άθικτο: η έξοδος γράφεται αλλού).
// ΑΚΡΙΒΕΙΑ, μορφή 64 μπιτ (κάθε βάθος, SSE4.1). Το PMULDQ δίνει ακριβές γινόμενο 64 μπιτ και οι θέσεις είναι 64 μπιτ: είναι
// η αριθμητική του προτύπου αυτούσια, χωρίς καμία προϋπόθεση. Από το ολισθημένο ποσό κρατιούνται τα κάτω 32 μπιτ.
template <int L>
FLAC_INLINE __m128i laneAll(__m128i v) { return _mm_shuffle_epi32(v, L * 0x55); }
template <int L>
FLAC_INLINE __m128i laneOff() { return _mm_set_epi32(L == 3 ? 0 : -1, L == 2 ? 0 : -1, L == 1 ? 0 : -1, L == 0 ? 0 : -1); }

template <int N>
struct Ring16 {
    __m128i acc[N / 4];
    __m128i count;
    const int32_t *res;
    int32_t *out;
    const Sub *sub;
    alignas(16) int16_t c[N][N / 4][8];   // c[βήμα][καταχωρητής]: ζεύγη (συντελεστής, 0) για το PMADDWD
};

// Τα πρώτα N δείγματα με τον ακριβή απλό δρόμο, οι N περιστροφές των συντελεστών, και οι N προβλέψεις «υπό κατασκευή» για
// τα δείγματα N .. 2N-1 (ό,τι τους χρωστούν τα δείγματα 0 .. N-1).
template <int N>
void ring16Start(Ring16<N> &r, const int32_t *res, int32_t *out, const Sub &s)
{
    r.res = res;
    r.out = out;
    r.sub = &s;
    r.count = _mm_cvtsi32_si128(s.shift);
    const uint32_t order = s.order;
    for (uint32_t i = 0; i < order; ++i) out[i] = res[i];
    for (uint32_t i = order; i < uint32_t(N); ++i) out[i] = lpcOne(out + i, res[i], s.coef, order, s.shift);
    for (int h = 0; h < N; ++h)
        for (int m = 0; m < N; ++m) {   // το βήμα h δίνει στη θέση m τον συντελεστή του δείγματος που απέχει (m - h) από αυτό
            const uint32_t k = uint32_t((m - h - 1 + N) % N);
            r.c[h][m / 4][(m % 4) * 2] = k < order ? int16_t(s.coef[k]) : int16_t(0);
            r.c[h][m / 4][(m % 4) * 2 + 1] = 0;
        }
    alignas(16) int32_t a[N];
    for (uint32_t k = 0; k < uint32_t(N); ++k) {
        uint32_t sum = 0;
        for (uint32_t j = k; j < order; ++j) sum += uint32_t(s.coef[j]) * uint32_t(out[uint32_t(N) + k - 1 - j]);
        a[k] = int32_t(sum);
    }
    for (int v = 0; v < N / 4; ++v) r.acc[v] = _mm_load_si128(reinterpret_cast<const __m128i *>(a + 4 * v));
}

template <int N, int H>
FLAC_INLINE void ring16Step(__m128i *acc, __m128i count, __m128i &seen, const int16_t (*c)[8], int32_t *dst)
{
    const __m128i x = laneAll<H & 3>(_mm_sra_epi32(acc[H >> 2], count));   // το δείγμα, σε όλες τις θέσεις
    *dst = _mm_cvtsi128_si32(x);
    seen = _mm_or_si128(seen, _mm_add_epi32(x, _mm_set1_epi32(0x8000)));   // μένει κάτω από 2^16 μόνο αν όλα χωράνε σε 16 μπιτ
    acc[H >> 2] = _mm_and_si128(acc[H >> 2], laneOff<H & 3>());            // η θέση αδειάζει: τώρα ανήκει στο δείγμα N θέσεις μετά
    for (int w = 0; w < N / 4; ++w)
        acc[w] = _mm_add_epi32(acc[w], _mm_madd_epi16(x, _mm_load_si128(reinterpret_cast<const __m128i *>(c[w]))));
}

template <int N, int H>
FLAC_INLINE void ring16Steps2(const Ring16<N> &a, const Ring16<N> &b, __m128i *aa, __m128i *ba, __m128i ac, __m128i bc, __m128i &as, __m128i &bs,
                              uint32_t i)
{
    if constexpr ((H & 3) == 0) {   // τα 4 επόμενα υπόλοιπα μπαίνουν στις θέσεις τους, ολισθημένα
        aa[H >> 2] = _mm_add_epi32(aa[H >> 2], _mm_sll_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i *>(a.res + i + H)), ac));
        ba[H >> 2] = _mm_add_epi32(ba[H >> 2], _mm_sll_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i *>(b.res + i + H)), bc));
    }
    ring16Step<N, H>(aa, ac, as, a.c[H], a.out + i + H);
    ring16Step<N, H>(ba, bc, bs, b.c[H], b.out + i + H);
    if constexpr (H + 1 < N) ring16Steps2<N, H + 1>(a, b, aa, ba, ac, bc, as, bs, i);
}

template <int N, int H>
FLAC_INLINE void ring16Steps1(const Ring16<N> &a, __m128i *aa, __m128i ac, __m128i &as, uint32_t i)
{
    if constexpr ((H & 3) == 0)
        aa[H >> 2] = _mm_add_epi32(aa[H >> 2], _mm_sll_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i *>(a.res + i + H)), ac));
    ring16Step<N, H>(aa, ac, as, a.c[H], a.out + i + H);
    if constexpr (H + 1 < N) ring16Steps1<N, H + 1>(a, aa, ac, as, i);
}

// Η ουρά (λιγότερα από N δείγματα) με τον ακριβή απλό δρόμο, και ο έλεγχος «όλα χωράνε σε 16 μπιτ» για αρχή, μέση και ουρά.
template <int N>
bool ring16Finish(const Ring16<N> &r, __m128i seen, uint32_t i, uint32_t n)
{
    uint32_t bits = uint32_t(_mm_cvtsi128_si32(seen));
    for (uint32_t k = 0; k < uint32_t(N); ++k) bits |= uint32_t(r.out[k]) + 0x8000u;
    for (; i < n; ++i) {
        r.out[i] = lpcOne(r.out + i, r.res[i], r.sub->coef, r.sub->order, r.sub->shift);
        bits |= uint32_t(r.out[i]) + 0x8000u;
    }
    return (bits >> 16) == 0;
}

template <int N>
FLAC_NOINLINE void lpcRing16x2(const int32_t *resA, int32_t *outA, const Sub &sa, const int32_t *resB, int32_t *outB, const Sub &sb, uint32_t n,
                               bool &okA, bool &okB)
{
    Ring16<N> a, b;
    ring16Start(a, resA, outA, sa);
    ring16Start(b, resB, outB, sb);
    __m128i aa[N / 4], ba[N / 4];
    for (int v = 0; v < N / 4; ++v) {
        aa[v] = a.acc[v];
        ba[v] = b.acc[v];
    }
    __m128i as = _mm_setzero_si128(), bs = _mm_setzero_si128();
    uint32_t i = uint32_t(N);
    for (; i + uint32_t(N) <= n; i += uint32_t(N)) ring16Steps2<N, 0>(a, b, aa, ba, a.count, b.count, as, bs, i);
    okA = ring16Finish(a, as, i, n);
    okB = ring16Finish(b, bs, i, n);
}

template <int N>
FLAC_NOINLINE bool lpcRing16x1(const int32_t *res, int32_t *out, const Sub &s, uint32_t n)
{
    Ring16<N> a;
    ring16Start(a, res, out, s);
    __m128i aa[N / 4];
    for (int v = 0; v < N / 4; ++v) aa[v] = a.acc[v];
    __m128i as = _mm_setzero_si128();
    uint32_t i = uint32_t(N);
    for (; i + uint32_t(N) <= n; i += uint32_t(N)) ring16Steps1<N, 0>(a, aa, a.count, as, i);
    return ring16Finish(a, as, i, n);
}

// Η μορφή 64 μπιτ: ίδια ιδέα, δύο θέσεις ανά καταχωρητή.
#if defined(_MSC_VER) && !defined(__clang__)
#define FLAC_TARGET_SSE41
#else
#define FLAC_TARGET_SSE41 __attribute__((target("sse4.1")))
#endif

template <int N>
struct Ring64 {
    __m128i acc[N / 2];
    __m128i count;
    const int32_t *res;
    int32_t *out;
    const Sub *sub;
    alignas(16) int64_t c[N][N / 2][2];
};

template <int N>
void ring64Start(Ring64<N> &r, const int32_t *res, int32_t *out, const Sub &s)
{
    r.res = res;
    r.out = out;
    r.sub = &s;
    r.count = _mm_cvtsi32_si128(s.shift);
    const uint32_t order = s.order;
    for (uint32_t i = 0; i < order; ++i) out[i] = res[i];
    for (uint32_t i = order; i < uint32_t(N); ++i) out[i] = lpcOne(out + i, res[i], s.coef, order, s.shift);
    for (int h = 0; h < N; ++h)
        for (int m = 0; m < N; ++m) {
            const uint32_t k = uint32_t((m - h - 1 + N) % N);
            r.c[h][m / 2][m % 2] = k < order ? int64_t(s.coef[k]) : 0;
        }
    alignas(16) int64_t a[N];
    for (uint32_t k = 0; k < uint32_t(N); ++k) {
        int64_t sum = 0;
        for (uint32_t j = k; j < order; ++j) sum += int64_t(s.coef[j]) * int64_t(out[uint32_t(N) + k - 1 - j]);
        a[k] = sum;
    }
    for (int v = 0; v < N / 2; ++v) r.acc[v] = _mm_load_si128(reinterpret_cast<const __m128i *>(a + 2 * v));
}

template <int N, int H>
FLAC_INLINE FLAC_TARGET_SSE41 void ring64Step(__m128i *acc, __m128i count, const int64_t (*c)[2], int32_t *dst)
{
    const __m128i x = _mm_shuffle_epi32(_mm_srl_epi64(acc[H >> 1], count), (H & 1) ? 0xAA : 0x00);   // τα κάτω 32 μπιτ, παντού
    *dst = _mm_cvtsi128_si32(x);
    acc[H >> 1] = _mm_and_si128(acc[H >> 1], (H & 1) ? _mm_set_epi64x(0, -1) : _mm_set_epi64x(-1, 0));
    for (int w = 0; w < N / 2; ++w)
        acc[w] = _mm_add_epi64(acc[w], _mm_mul_epi32(x, _mm_load_si128(reinterpret_cast<const __m128i *>(c[w]))));
}

template <int N, int H>
FLAC_INLINE FLAC_TARGET_SSE41 void ring64Steps2(const Ring64<N> &a, const Ring64<N> &b, __m128i *aa, __m128i *ba, __m128i ac, __m128i bc, uint32_t i)
{
    if constexpr ((H & 1) == 0) {
        aa[H >> 1] = _mm_add_epi64(
            aa[H >> 1], _mm_sll_epi64(_mm_cvtepi32_epi64(_mm_loadl_epi64(reinterpret_cast<const __m128i *>(a.res + i + H))), ac));
        ba[H >> 1] = _mm_add_epi64(
            ba[H >> 1], _mm_sll_epi64(_mm_cvtepi32_epi64(_mm_loadl_epi64(reinterpret_cast<const __m128i *>(b.res + i + H))), bc));
    }
    ring64Step<N, H>(aa, ac, a.c[H], a.out + i + H);
    ring64Step<N, H>(ba, bc, b.c[H], b.out + i + H);
    if constexpr (H + 1 < N) ring64Steps2<N, H + 1>(a, b, aa, ba, ac, bc, i);
}

template <int N, int H>
FLAC_INLINE FLAC_TARGET_SSE41 void ring64Steps1(const Ring64<N> &a, __m128i *aa, __m128i ac, uint32_t i)
{
    if constexpr ((H & 1) == 0)
        aa[H >> 1] = _mm_add_epi64(
            aa[H >> 1], _mm_sll_epi64(_mm_cvtepi32_epi64(_mm_loadl_epi64(reinterpret_cast<const __m128i *>(a.res + i + H))), ac));
    ring64Step<N, H>(aa, ac, a.c[H], a.out + i + H);
    if constexpr (H + 1 < N) ring64Steps1<N, H + 1>(a, aa, ac, i);
}

template <int N>
void ring64Finish(const Ring64<N> &r, uint32_t i, uint32_t n)
{
    for (; i < n; ++i) r.out[i] = lpcOne(r.out + i, r.res[i], r.sub->coef, r.sub->order, r.sub->shift);
}

template <int N>
FLAC_NOINLINE FLAC_TARGET_SSE41 void lpcRing64x2(const int32_t *resA, int32_t *outA, const Sub &sa, const int32_t *resB, int32_t *outB,
                                                 const Sub &sb, uint32_t n)
{
    Ring64<N> a, b;
    ring64Start(a, resA, outA, sa);
    ring64Start(b, resB, outB, sb);
    __m128i aa[N / 2], ba[N / 2];
    for (int v = 0; v < N / 2; ++v) {
        aa[v] = a.acc[v];
        ba[v] = b.acc[v];
    }
    uint32_t i = uint32_t(N);
    for (; i + uint32_t(N) <= n; i += uint32_t(N)) ring64Steps2<N, 0>(a, b, aa, ba, a.count, b.count, i);
    ring64Finish(a, i, n);
    ring64Finish(b, i, n);
}

template <int N>
FLAC_NOINLINE FLAC_TARGET_SSE41 void lpcRing64x1(const int32_t *res, int32_t *out, const Sub &s, uint32_t n)
{
    Ring64<N> a;
    ring64Start(a, res, out, s);
    __m128i aa[N / 2];
    for (int v = 0; v < N / 2; ++v) aa[v] = a.acc[v];
    uint32_t i = uint32_t(N);
    for (; i + uint32_t(N) <= n; i += uint32_t(N)) ring64Steps1<N, 0>(a, aa, a.count, i);
    ring64Finish(a, i, n);
}

// Η φαρδιά μορφή ΧΩΡΙΣ SSE4.1 (επεξεργαστές πριν από το 2008, και οι AMD ως το 2011). Το SSE2 έχει μόνο πολλαπλασιασμό
// χωρίς πρόσημο (PMULUDQ), οπότε οι δύο παράγοντες μετατοπίζονται για να γίνουν μη αρνητικοί: x' = x + 2^31, c' = c + 2^15.
// Τότε x·c = x'·c' - 2^15·x' - 2^31·c' + 2^46. Το x'·c' το δίνει το PMULUDQ ακριβώς· το 2^15·x' αφαιρείται από κάθε θέση σε
// κάθε βήμα· και το υπόλοιπο (-2^31·c' + 2^46) είναι σταθερά: κάθε θέση, από τη στιγμή που αδειάζει ως τη στιγμή που
// διαβάζεται, δέχεται ακριβώς μία φορά κάθε συντελεστή, άρα το άθροισμα αυτών των σταθερών είναι το ίδιο για όλες τις
// θέσεις (G) και μπαίνει μια κι έξω όταν η θέση αδειάζει. Όλα modulo 2^64: ίδια μπιτ με τη μορφή του SSE4.1.
template <int N>
struct RingB64 {
    __m128i acc[N / 2];
    __m128i count;
    __m128i fresh;   // το G και στις δύο θέσεις
    const int32_t *res;
    int32_t *out;
    const Sub *sub;
    alignas(16) uint64_t c[N][N / 2][2];   // c' = c + 2^15 στα κάτω 32 μπιτ κάθε θέσης
};

template <int N>
void ringB64Start(RingB64<N> &r, const int32_t *res, int32_t *out, const Sub &s)
{
    r.res = res;
    r.out = out;
    r.sub = &s;
    r.count = _mm_cvtsi32_si128(s.shift);
    const uint32_t order = s.order;
    for (uint32_t i = 0; i < order; ++i) out[i] = res[i];
    for (uint32_t i = order; i < uint32_t(N); ++i) out[i] = lpcOne(out + i, res[i], s.coef, order, s.shift);
    uint64_t fix[N];   // η σταθερά που αναλογεί στον συντελεστή k
    uint64_t g = 0;
    for (uint32_t k = 0; k < uint32_t(N); ++k) {
        const uint64_t cb = uint64_t(int64_t(k < order ? s.coef[k] : 0) + 32768);
        fix[k] = (uint64_t(1) << 46) - (cb << 31);
        g += fix[k];
    }
    r.fresh = _mm_set1_epi64x(int64_t(g));
    for (int h = 0; h < N; ++h)
        for (int m = 0; m < N; ++m) {
            const uint32_t k = uint32_t((m - h - 1 + N) % N);
            r.c[h][m / 2][m % 2] = uint64_t(int64_t(k < order ? s.coef[k] : 0) + 32768);
        }
    alignas(16) uint64_t a[N];
    for (uint32_t k = 0; k < uint32_t(N); ++k) {
        uint64_t sum = 0;   // ό,τι χρωστούν στο δείγμα N+k τα δείγματα 0..N-1, και οι σταθερές των συντελεστών που θα έρθουν (0..k-1)
        for (uint32_t j = k; j < order; ++j) sum += uint64_t(int64_t(s.coef[j]) * int64_t(out[uint32_t(N) + k - 1 - j]));
        for (uint32_t j = 0; j < k; ++j) sum += fix[j];
        a[k] = sum;
    }
    for (int v = 0; v < N / 2; ++v) r.acc[v] = _mm_load_si128(reinterpret_cast<const __m128i *>(a + 2 * v));
}

template <int N, int H>
FLAC_INLINE void ringB64Step(__m128i *acc, __m128i count, __m128i fresh, const uint64_t (*c)[2], int32_t *dst)
{
    const __m128i x = _mm_shuffle_epi32(_mm_srl_epi64(acc[H >> 1], count), (H & 1) ? 0xAA : 0x00);   // τα κάτω 32 μπιτ, παντού
    *dst = _mm_cvtsi128_si32(x);
    const __m128i xb = _mm_xor_si128(x, _mm_set1_epi32(int32_t(0x80000000u)));                        // x' = x + 2^31
    const __m128i xs = _mm_mul_epu32(xb, _mm_set1_epi64x(32768));                                     // 2^15 · x'
    const __m128i keep = (H & 1) ? _mm_set_epi64x(0, -1) : _mm_set_epi64x(-1, 0);
    acc[H >> 1] = _mm_or_si128(_mm_and_si128(acc[H >> 1], keep), _mm_andnot_si128(keep, fresh));     // η θέση αδειάζει: ξεκινά από το G
    for (int w = 0; w < N / 2; ++w)
        acc[w] = _mm_sub_epi64(_mm_add_epi64(acc[w], _mm_mul_epu32(xb, _mm_load_si128(reinterpret_cast<const __m128i *>(c[w])))), xs);
}

FLAC_INLINE __m128i signExtend2(const int32_t *p)   // δύο δείγματα 32 μπιτ σε δύο θέσεις 64 μπιτ, με το πρόσημό τους
{
    const __m128i v = _mm_loadl_epi64(reinterpret_cast<const __m128i *>(p));
    return _mm_unpacklo_epi32(v, _mm_srai_epi32(v, 31));
}

template <int N, int H>
FLAC_INLINE void ringB64Steps2(const RingB64<N> &a, const RingB64<N> &b, __m128i *aa, __m128i *ba, uint32_t i)
{
    if constexpr ((H & 1) == 0) {
        aa[H >> 1] = _mm_add_epi64(aa[H >> 1], _mm_sll_epi64(signExtend2(a.res + i + H), a.count));
        ba[H >> 1] = _mm_add_epi64(ba[H >> 1], _mm_sll_epi64(signExtend2(b.res + i + H), b.count));
    }
    ringB64Step<N, H>(aa, a.count, a.fresh, a.c[H], a.out + i + H);
    ringB64Step<N, H>(ba, b.count, b.fresh, b.c[H], b.out + i + H);
    if constexpr (H + 1 < N) ringB64Steps2<N, H + 1>(a, b, aa, ba, i);
}

template <int N, int H>
FLAC_INLINE void ringB64Steps1(const RingB64<N> &a, __m128i *aa, uint32_t i)
{
    if constexpr ((H & 1) == 0) aa[H >> 1] = _mm_add_epi64(aa[H >> 1], _mm_sll_epi64(signExtend2(a.res + i + H), a.count));
    ringB64Step<N, H>(aa, a.count, a.fresh, a.c[H], a.out + i + H);
    if constexpr (H + 1 < N) ringB64Steps1<N, H + 1>(a, aa, i);
}

template <int N>
void ringB64Finish(const RingB64<N> &r, uint32_t i, uint32_t n)
{
    for (; i < n; ++i) r.out[i] = lpcOne(r.out + i, r.res[i], r.sub->coef, r.sub->order, r.sub->shift);
}

template <int N>
FLAC_NOINLINE void lpcRingB64x2(const int32_t *resA, int32_t *outA, const Sub &sa, const int32_t *resB, int32_t *outB, const Sub &sb, uint32_t n)
{
    RingB64<N> a, b;
    ringB64Start(a, resA, outA, sa);
    ringB64Start(b, resB, outB, sb);
    __m128i aa[N / 2], ba[N / 2];
    for (int v = 0; v < N / 2; ++v) {
        aa[v] = a.acc[v];
        ba[v] = b.acc[v];
    }
    uint32_t i = uint32_t(N);
    for (; i + uint32_t(N) <= n; i += uint32_t(N)) ringB64Steps2<N, 0>(a, b, aa, ba, i);
    ringB64Finish(a, i, n);
    ringB64Finish(b, i, n);
}

template <int N>
FLAC_NOINLINE void lpcRingB64x1(const int32_t *res, int32_t *out, const Sub &s, uint32_t n)
{
    RingB64<N> a;
    ringB64Start(a, res, out, s);
    __m128i aa[N / 2];
    for (int v = 0; v < N / 2; ++v) aa[v] = a.acc[v];
    uint32_t i = uint32_t(N);
    for (; i + uint32_t(N) <= n; i += uint32_t(N)) ringB64Steps1<N, 0>(a, aa, i);
    ringB64Finish(a, i, n);
}

bool cpuHasSse41()
{
#if defined(_MSC_VER) && !defined(__clang__)
    int r[4];
    __cpuid(r, 1);
    return (r[2] & (1 << 19)) != 0;
#else
    unsigned a = 0, b = 0, c = 0, d = 0;
    return __get_cpuid(1, &a, &b, &c, &d) && (c & (1u << 19));
#endif
}
bool kHasSse41 = cpuHasSse41();

inline int ringSize(uint32_t order) { return order <= 8 ? 8 : order <= 12 ? 12 : order <= 16 ? 16 : order <= 24 ? 24 : 32; }

void ringNarrowPair(int size, const int32_t *resA, int32_t *outA, const Sub &sa, const int32_t *resB, int32_t *outB, const Sub &sb, uint32_t n, bool &okA,
                bool &okB)
{
    if (size == 8) lpcRing16x2<8>(resA, outA, sa, resB, outB, sb, n, okA, okB);
    else if (size == 12) lpcRing16x2<12>(resA, outA, sa, resB, outB, sb, n, okA, okB);
    else if (size == 16) lpcRing16x2<16>(resA, outA, sa, resB, outB, sb, n, okA, okB);
    else if (size == 24) lpcRing16x2<24>(resA, outA, sa, resB, outB, sb, n, okA, okB);
    else lpcRing16x2<32>(resA, outA, sa, resB, outB, sb, n, okA, okB);
}
bool ringNarrowOne(int size, const int32_t *res, int32_t *out, const Sub &s, uint32_t n)
{
    return size == 8 ? lpcRing16x1<8>(res, out, s, n)
           : size == 12 ? lpcRing16x1<12>(res, out, s, n)
           : size == 16 ? lpcRing16x1<16>(res, out, s, n) : size == 24 ? lpcRing16x1<24>(res, out, s, n) : lpcRing16x1<32>(res, out, s, n);
}
void ringWidePair(int size, const int32_t *resA, int32_t *outA, const Sub &sa, const int32_t *resB, int32_t *outB, const Sub &sb, uint32_t n)
{
    if (!kHasSse41) {
        if (size == 8) lpcRingB64x2<8>(resA, outA, sa, resB, outB, sb, n);
        else if (size == 12) lpcRingB64x2<12>(resA, outA, sa, resB, outB, sb, n);
        else if (size == 16) lpcRingB64x2<16>(resA, outA, sa, resB, outB, sb, n);
        else if (size == 24) lpcRingB64x2<24>(resA, outA, sa, resB, outB, sb, n);
        else lpcRingB64x2<32>(resA, outA, sa, resB, outB, sb, n);
        return;
    }
    if (size == 8) lpcRing64x2<8>(resA, outA, sa, resB, outB, sb, n);
    else if (size == 12) lpcRing64x2<12>(resA, outA, sa, resB, outB, sb, n);
    else if (size == 16) lpcRing64x2<16>(resA, outA, sa, resB, outB, sb, n);
    else if (size == 24) lpcRing64x2<24>(resA, outA, sa, resB, outB, sb, n);
    else lpcRing64x2<32>(resA, outA, sa, resB, outB, sb, n);
}
void ringWideOne(int size, const int32_t *res, int32_t *out, const Sub &s, uint32_t n)
{
    if (!kHasSse41) {
        if (size == 8) lpcRingB64x1<8>(res, out, s, n);
        else if (size == 12) lpcRingB64x1<12>(res, out, s, n);
        else if (size == 16) lpcRingB64x1<16>(res, out, s, n);
        else if (size == 24) lpcRingB64x1<24>(res, out, s, n);
        else lpcRingB64x1<32>(res, out, s, n);
        return;
    }
    if (size == 8) lpcRing64x1<8>(res, out, s, n);
    else if (size == 12) lpcRing64x1<12>(res, out, s, n);
    else if (size == 16) lpcRing64x1<16>(res, out, s, n);
    else if (size == 24) lpcRing64x1<24>(res, out, s, n);
    else lpcRing64x1<32>(res, out, s, n);
}
inline bool ringNarrowFits(const Sub &, int frameBits) { return frameBits <= 16; }   // και ελέγχεται κάθε δείγμα όσο τρέχει
#define FLAC_RING_WIDE true
constexpr uint32_t kRingAloneFrom = 9;   // ένα κανάλι μόνο του: ο απλός δρόμος είναι εξίσου γρήγορος ως την τάξη 8
#define FLAC_HAVE_RING 1

#elif defined(__aarch64__) && !defined(_MSC_VER)
// ── Ο ΓΡΑΜΜΙΚΟΣ ΠΡΟΒΛΕΠΤΗΣ «ΣΕ ΔΑΧΤΥΛΙΔΙ» (ARM64) ──
// Η ίδια ιδέα με του x86-64 (βλ. πιο πάνω), με τις εντολές του NEON: «πολλαπλασίασε 4 θέσεις με μία θέση άλλου καταχωρητή
// και πρόσθεσε» (MLA) είναι ΜΙΑ εντολή, άρα το δείγμα δεν χρειάζεται καν να απλωθεί πρώτα.
// ΑΚΡΙΒΕΙΑ, στενή μορφή (θέσεις 32 μπιτ). Εδώ ο πολλαπλασιασμός είναι πλήρης 32×32 modulo 2^32, όχι 16×16: δεν υπάρχει
// προϋπόθεση για τα δείγματα. Το τελικό ποσό ενός δείγματος x είναι στο [x·2^s, (x+1)·2^s), άρα χωράει σε 32 μπιτ όταν
// (βάθος υποκαρέ + ολίσθηση) ≤ 32, και τότε το αποτέλεσμα είναι ακριβώς του προτύπου όπως κι αν τυλίχτηκαν τα ενδιάμεσα.
// ΑΚΡΙΒΕΙΑ, φαρδιά μορφή (θέσεις 64 μπιτ, SMLAL): η αριθμητική του προτύπου αυτούσια, χωρίς προϋπόθεση.
template <int N>
struct Ring32 {
    int32x4_t acc[N / 4];
    int32x4_t down, up;   // -ολίσθηση (για το δείγμα) και +ολίσθηση (για το υπόλοιπο)
    const int32_t *res;
    int32_t *out;
    const Sub *sub;
    alignas(16) int32_t c[N][N / 4][4];
};

template <int N>
void ring32Start(Ring32<N> &r, const int32_t *res, int32_t *out, const Sub &s)
{
    r.res = res;
    r.out = out;
    r.sub = &s;
    r.down = vdupq_n_s32(-int32_t(s.shift));
    r.up = vdupq_n_s32(int32_t(s.shift));
    const uint32_t order = s.order;
    for (uint32_t i = 0; i < order; ++i) out[i] = res[i];
    for (uint32_t i = order; i < uint32_t(N); ++i) out[i] = lpcOne(out + i, res[i], s.coef, order, s.shift);
    for (int h = 0; h < N; ++h)
        for (int m = 0; m < N; ++m) {
            const uint32_t k = uint32_t((m - h - 1 + N) % N);
            r.c[h][m / 4][m % 4] = k < order ? s.coef[k] : 0;
        }
    alignas(16) int32_t a[N];
    for (uint32_t k = 0; k < uint32_t(N); ++k) {
        uint32_t sum = 0;
        for (uint32_t j = k; j < order; ++j) sum += uint32_t(s.coef[j]) * uint32_t(out[uint32_t(N) + k - 1 - j]);
        a[k] = int32_t(sum);
    }
    for (int v = 0; v < N / 4; ++v) r.acc[v] = vld1q_s32(a + 4 * v);
}

template <int N, int H>
FLAC_INLINE void ring32Step(int32x4_t *acc, int32x4_t down, const int32_t (*c)[4], int32_t *dst)
{
    const int32x4_t t = vshlq_s32(acc[H >> 2], down);          // στη θέση H & 3: το δείγμα
    vst1q_lane_s32(dst, t, H & 3);
    acc[H >> 2] = vsetq_lane_s32(0, acc[H >> 2], H & 3);       // η θέση αδειάζει: τώρα ανήκει στο δείγμα N θέσεις μετά
    for (int w = 0; w < N / 4; ++w) acc[w] = vmlaq_laneq_s32(acc[w], vld1q_s32(c[w]), t, H & 3);
}

template <int N, int H>
FLAC_INLINE void ring32Steps2(const Ring32<N> &a, const Ring32<N> &b, int32x4_t *aa, int32x4_t *ba, uint32_t i)
{
    if constexpr ((H & 3) == 0) {
        aa[H >> 2] = vaddq_s32(aa[H >> 2], vshlq_s32(vld1q_s32(a.res + i + H), a.up));
        ba[H >> 2] = vaddq_s32(ba[H >> 2], vshlq_s32(vld1q_s32(b.res + i + H), b.up));
    }
    ring32Step<N, H>(aa, a.down, a.c[H], a.out + i + H);
    ring32Step<N, H>(ba, b.down, b.c[H], b.out + i + H);
    if constexpr (H + 1 < N) ring32Steps2<N, H + 1>(a, b, aa, ba, i);
}

template <int N, int H>
FLAC_INLINE void ring32Steps1(const Ring32<N> &a, int32x4_t *aa, uint32_t i)
{
    if constexpr ((H & 3) == 0) aa[H >> 2] = vaddq_s32(aa[H >> 2], vshlq_s32(vld1q_s32(a.res + i + H), a.up));
    ring32Step<N, H>(aa, a.down, a.c[H], a.out + i + H);
    if constexpr (H + 1 < N) ring32Steps1<N, H + 1>(a, aa, i);
}

template <int N>
void ring32Finish(const Ring32<N> &r, uint32_t i, uint32_t n)
{
    for (; i < n; ++i) r.out[i] = lpcOne(r.out + i, r.res[i], r.sub->coef, r.sub->order, r.sub->shift);
}

template <int N>
FLAC_NOINLINE void lpcRing32x2(const int32_t *resA, int32_t *outA, const Sub &sa, const int32_t *resB, int32_t *outB, const Sub &sb, uint32_t n)
{
    Ring32<N> a, b;
    ring32Start(a, resA, outA, sa);
    ring32Start(b, resB, outB, sb);
    int32x4_t aa[N / 4], ba[N / 4];
    for (int v = 0; v < N / 4; ++v) {
        aa[v] = a.acc[v];
        ba[v] = b.acc[v];
    }
    uint32_t i = uint32_t(N);
    for (; i + uint32_t(N) <= n; i += uint32_t(N)) ring32Steps2<N, 0>(a, b, aa, ba, i);
    ring32Finish(a, i, n);
    ring32Finish(b, i, n);
}

template <int N>
FLAC_NOINLINE void lpcRing32x1(const int32_t *res, int32_t *out, const Sub &s, uint32_t n)
{
    Ring32<N> a;
    ring32Start(a, res, out, s);
    int32x4_t aa[N / 4];
    for (int v = 0; v < N / 4; ++v) aa[v] = a.acc[v];
    uint32_t i = uint32_t(N);
    for (; i + uint32_t(N) <= n; i += uint32_t(N)) ring32Steps1<N, 0>(a, aa, i);
    ring32Finish(a, i, n);
}

template <int N>
struct Ring64 {
    int64x2_t acc[N / 2];
    int64x2_t down, up;
    const int32_t *res;
    int32_t *out;
    const Sub *sub;
    alignas(16) int32_t c[N][N / 2][2];
};

template <int N>
void ring64Start(Ring64<N> &r, const int32_t *res, int32_t *out, const Sub &s)
{
    r.res = res;
    r.out = out;
    r.sub = &s;
    r.down = vdupq_n_s64(-int64_t(s.shift));
    r.up = vdupq_n_s64(int64_t(s.shift));
    const uint32_t order = s.order;
    for (uint32_t i = 0; i < order; ++i) out[i] = res[i];
    for (uint32_t i = order; i < uint32_t(N); ++i) out[i] = lpcOne(out + i, res[i], s.coef, order, s.shift);
    for (int h = 0; h < N; ++h)
        for (int m = 0; m < N; ++m) {
            const uint32_t k = uint32_t((m - h - 1 + N) % N);
            r.c[h][m / 2][m % 2] = k < order ? s.coef[k] : 0;
        }
    alignas(16) int64_t a[N];
    for (uint32_t k = 0; k < uint32_t(N); ++k) {
        int64_t sum = 0;
        for (uint32_t j = k; j < order; ++j) sum += int64_t(s.coef[j]) * int64_t(out[uint32_t(N) + k - 1 - j]);
        a[k] = sum;
    }
    for (int v = 0; v < N / 2; ++v) r.acc[v] = vld1q_s64(a + 2 * v);
}

template <int N, int H>
FLAC_INLINE void ring64Step(int64x2_t *acc, int64x2_t down, const int32_t (*c)[2], int32_t *dst)
{
    // λογική ολίσθηση δεξιά: τα κάτω 32 μπιτ είναι ίδια με της αριθμητικής, και μόνο αυτά κρατιούνται
    const int32x4_t t = vreinterpretq_s32_u64(vshlq_u64(vreinterpretq_u64_s64(acc[H >> 1]), down));
    vst1q_lane_s32(dst, t, 2 * (H & 1));
    acc[H >> 1] = vsetq_lane_s64(0, acc[H >> 1], H & 1);
    for (int w = 0; w < N / 2; ++w) acc[w] = vmlal_laneq_s32(acc[w], vld1_s32(c[w]), t, 2 * (H & 1));
}

template <int N, int H>
FLAC_INLINE void ring64Steps2(const Ring64<N> &a, const Ring64<N> &b, int64x2_t *aa, int64x2_t *ba, uint32_t i)
{
    if constexpr ((H & 1) == 0) {
        aa[H >> 1] = vaddq_s64(aa[H >> 1], vshlq_s64(vmovl_s32(vld1_s32(a.res + i + H)), a.up));
        ba[H >> 1] = vaddq_s64(ba[H >> 1], vshlq_s64(vmovl_s32(vld1_s32(b.res + i + H)), b.up));
    }
    ring64Step<N, H>(aa, a.down, a.c[H], a.out + i + H);
    ring64Step<N, H>(ba, b.down, b.c[H], b.out + i + H);
    if constexpr (H + 1 < N) ring64Steps2<N, H + 1>(a, b, aa, ba, i);
}

template <int N, int H>
FLAC_INLINE void ring64Steps1(const Ring64<N> &a, int64x2_t *aa, uint32_t i)
{
    if constexpr ((H & 1) == 0) aa[H >> 1] = vaddq_s64(aa[H >> 1], vshlq_s64(vmovl_s32(vld1_s32(a.res + i + H)), a.up));
    ring64Step<N, H>(aa, a.down, a.c[H], a.out + i + H);
    if constexpr (H + 1 < N) ring64Steps1<N, H + 1>(a, aa, i);
}

template <int N>
void ring64Finish(const Ring64<N> &r, uint32_t i, uint32_t n)
{
    for (; i < n; ++i) r.out[i] = lpcOne(r.out + i, r.res[i], r.sub->coef, r.sub->order, r.sub->shift);
}

template <int N>
FLAC_NOINLINE void lpcRing64x2(const int32_t *resA, int32_t *outA, const Sub &sa, const int32_t *resB, int32_t *outB, const Sub &sb, uint32_t n)
{
    Ring64<N> a, b;
    ring64Start(a, resA, outA, sa);
    ring64Start(b, resB, outB, sb);
    int64x2_t aa[N / 2], ba[N / 2];
    for (int v = 0; v < N / 2; ++v) {
        aa[v] = a.acc[v];
        ba[v] = b.acc[v];
    }
    uint32_t i = uint32_t(N);
    for (; i + uint32_t(N) <= n; i += uint32_t(N)) ring64Steps2<N, 0>(a, b, aa, ba, i);
    ring64Finish(a, i, n);
    ring64Finish(b, i, n);
}

template <int N>
FLAC_NOINLINE void lpcRing64x1(const int32_t *res, int32_t *out, const Sub &s, uint32_t n)
{
    Ring64<N> a;
    ring64Start(a, res, out, s);
    int64x2_t aa[N / 2];
    for (int v = 0; v < N / 2; ++v) aa[v] = a.acc[v];
    uint32_t i = uint32_t(N);
    for (; i + uint32_t(N) <= n; i += uint32_t(N)) ring64Steps1<N, 0>(a, aa, i);
    ring64Finish(a, i, n);
}

inline int ringSize(uint32_t order) { return order <= 8 ? 8 : order <= 12 ? 12 : order <= 16 ? 16 : order <= 24 ? 24 : 32; }

void ringNarrowPair(int size, const int32_t *resA, int32_t *outA, const Sub &sa, const int32_t *resB, int32_t *outB, const Sub &sb, uint32_t n,
                    bool &okA, bool &okB)
{
    if (size == 8) lpcRing32x2<8>(resA, outA, sa, resB, outB, sb, n);
    else if (size == 12) lpcRing32x2<12>(resA, outA, sa, resB, outB, sb, n);
    else if (size == 16) lpcRing32x2<16>(resA, outA, sa, resB, outB, sb, n);
    else if (size == 24) lpcRing32x2<24>(resA, outA, sa, resB, outB, sb, n);
    else lpcRing32x2<32>(resA, outA, sa, resB, outB, sb, n);
    okA = okB = true;
}
bool ringNarrowOne(int size, const int32_t *res, int32_t *out, const Sub &s, uint32_t n)
{
    if (size == 8) lpcRing32x1<8>(res, out, s, n);
    else if (size == 12) lpcRing32x1<12>(res, out, s, n);
    else if (size == 16) lpcRing32x1<16>(res, out, s, n);
    else if (size == 24) lpcRing32x1<24>(res, out, s, n);
    else lpcRing32x1<32>(res, out, s, n);
    return true;
}
void ringWidePair(int size, const int32_t *resA, int32_t *outA, const Sub &sa, const int32_t *resB, int32_t *outB, const Sub &sb, uint32_t n)
{
    if (size == 8) lpcRing64x2<8>(resA, outA, sa, resB, outB, sb, n);
    else if (size == 12) lpcRing64x2<12>(resA, outA, sa, resB, outB, sb, n);
    else if (size == 16) lpcRing64x2<16>(resA, outA, sa, resB, outB, sb, n);
    else if (size == 24) lpcRing64x2<24>(resA, outA, sa, resB, outB, sb, n);
    else lpcRing64x2<32>(resA, outA, sa, resB, outB, sb, n);
}
void ringWideOne(int size, const int32_t *res, int32_t *out, const Sub &s, uint32_t n)
{
    if (size == 8) lpcRing64x1<8>(res, out, s, n);
    else if (size == 12) lpcRing64x1<12>(res, out, s, n);
    else if (size == 16) lpcRing64x1<16>(res, out, s, n);
    else if (size == 24) lpcRing64x1<24>(res, out, s, n);
    else lpcRing64x1<32>(res, out, s, n);
}
inline bool ringNarrowFits(const Sub &s, int) { return int(s.bits) + int(s.shift) <= 32; }
#define FLAC_RING_WIDE true
constexpr uint32_t kRingAloneFrom = 2;
#define FLAC_HAVE_RING 1

#else
#define FLAC_HAVE_RING 0
#endif

// 4.2: αποσυσχέτιση στέρεο, επιτόπου, modulo 2^32. Με μέσο m και διαφορά d: δεξί = m - (d >> 1), αριστερό = δεξί + d (το
// ίδιο με «διπλασίασε το μέσο, πρόσθεσε το κάτω μπιτ της διαφοράς, πρόσθεσε ή αφαίρεσε τη διαφορά και μισό», χωρίς το
// ενδιάμεσο που ξεχειλίζει). Γραμμένο ρητά 4 δείγματα τη φορά: δεν το αφήνουμε στην καλή διάθεση του μεταγλωττιστή.
FLAC_NOINLINE void stereoUndo(int32_t *pa, int32_t *pb, uint32_t n, int assignment)
{
    uint32_t i = 0;
#if FLAC_X86
    for (; i + 4 <= n; i += 4) {
        const __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i *>(pa + i));
        const __m128i b = _mm_loadu_si128(reinterpret_cast<const __m128i *>(pb + i));
        if (assignment == 8) {
            _mm_storeu_si128(reinterpret_cast<__m128i *>(pb + i), _mm_sub_epi32(a, b));
        } else if (assignment == 9) {
            _mm_storeu_si128(reinterpret_cast<__m128i *>(pa + i), _mm_add_epi32(a, b));
        } else {
            const __m128i right = _mm_sub_epi32(a, _mm_srai_epi32(b, 1));
            _mm_storeu_si128(reinterpret_cast<__m128i *>(pb + i), right);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(pa + i), _mm_add_epi32(right, b));
        }
    }
#elif defined(__aarch64__) && !defined(_MSC_VER)
    for (; i + 4 <= n; i += 4) {
        const int32x4_t a = vld1q_s32(pa + i), b = vld1q_s32(pb + i);
        if (assignment == 8) {
            vst1q_s32(pb + i, vsubq_s32(a, b));
        } else if (assignment == 9) {
            vst1q_s32(pa + i, vaddq_s32(a, b));
        } else {
            const int32x4_t right = vsubq_s32(a, vshrq_n_s32(b, 1));
            vst1q_s32(pb + i, right);
            vst1q_s32(pa + i, vaddq_s32(right, b));
        }
    }
#endif
    uint32_t *a = reinterpret_cast<uint32_t *>(pa), *b = reinterpret_cast<uint32_t *>(pb);
    for (; i < n; ++i) {
        if (assignment == 8) {
            b[i] = a[i] - b[i];
        } else if (assignment == 9) {
            a[i] += b[i];
        } else {
            const uint32_t d = b[i];
            const uint32_t right = a[i] - uint32_t(int32_t(d) >> 1);
            b[i] = right;
            a[i] = right + d;
        }
    }
}

// Δύο κανάλια πλεγμένα στην έξοδο, 16 ή 32 μπιτ, με ολίσθηση αριστερά: η μορφή που ζητά ένας παίκτης για στέρεο.
// Ίδια μπιτ με τον απλό βρόχο (κάθε δείγμα ολισθαίνει και κόβεται στο πλάτος της εξόδου, χωρίς κορεσμό).
FLAC_NOINLINE void interleave2(const int32_t *l, const int32_t *r, uint32_t n, int shift, bool wide, void *dst)
{
    uint32_t i = 0;
    if (wide) {
        int32_t *o = static_cast<int32_t *>(dst);
#if FLAC_X86
        const __m128i count = _mm_cvtsi32_si128(shift);
        for (; i + 4 <= n; i += 4) {
            const __m128i a = _mm_sll_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i *>(l + i)), count);
            const __m128i b = _mm_sll_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i *>(r + i)), count);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(o + 2 * size_t(i)), _mm_unpacklo_epi32(a, b));
            _mm_storeu_si128(reinterpret_cast<__m128i *>(o + 2 * size_t(i) + 4), _mm_unpackhi_epi32(a, b));
        }
#elif defined(__aarch64__) && !defined(_MSC_VER)
        const int32x4_t count = vdupq_n_s32(shift);
        for (; i + 4 <= n; i += 4) {
            int32x4x2_t v;
            v.val[0] = vshlq_s32(vld1q_s32(l + i), count);
            v.val[1] = vshlq_s32(vld1q_s32(r + i), count);
            vst2q_s32(o + 2 * size_t(i), v);
        }
#endif
        for (; i < n; ++i) {
            o[2 * size_t(i)] = int32_t(uint32_t(l[i]) << shift);
            o[2 * size_t(i) + 1] = int32_t(uint32_t(r[i]) << shift);
        }
    } else {
        int16_t *o = static_cast<int16_t *>(dst);
#if FLAC_X86
        // (x << (ολίσθηση + 16)) >> 16 με πρόσημο: τα κάτω 16 μπιτ του ολισθημένου δείγματος, απλωμένα σε 32, ώστε το
        // «πακετάρισμα με κορεσμό» που ακολουθεί να μην έχει τίποτα να κορέσει.
        const __m128i count = _mm_cvtsi32_si128(shift + 16);
        for (; i + 8 <= n; i += 8) {
            const __m128i a0 = _mm_srai_epi32(_mm_sll_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i *>(l + i)), count), 16);
            const __m128i a1 = _mm_srai_epi32(_mm_sll_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i *>(l + i + 4)), count), 16);
            const __m128i b0 = _mm_srai_epi32(_mm_sll_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i *>(r + i)), count), 16);
            const __m128i b1 = _mm_srai_epi32(_mm_sll_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i *>(r + i + 4)), count), 16);
            const __m128i a = _mm_packs_epi32(a0, a1), b = _mm_packs_epi32(b0, b1);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(o + 2 * size_t(i)), _mm_unpacklo_epi16(a, b));
            _mm_storeu_si128(reinterpret_cast<__m128i *>(o + 2 * size_t(i) + 8), _mm_unpackhi_epi16(a, b));
        }
#elif defined(__aarch64__) && !defined(_MSC_VER)
        const int32x4_t count = vdupq_n_s32(shift);
        for (; i + 8 <= n; i += 8) {
            int16x8x2_t v;
            v.val[0] = vcombine_s16(vmovn_s32(vshlq_s32(vld1q_s32(l + i), count)), vmovn_s32(vshlq_s32(vld1q_s32(l + i + 4), count)));
            v.val[1] = vcombine_s16(vmovn_s32(vshlq_s32(vld1q_s32(r + i), count)), vmovn_s32(vshlq_s32(vld1q_s32(r + i + 4), count)));
            vst2q_s16(o + 2 * size_t(i), v);
        }
#endif
        for (; i < n; ++i) {
            o[2 * size_t(i)] = int16_t(uint16_t(uint32_t(l[i]) << shift));
            o[2 * size_t(i) + 1] = int16_t(uint16_t(uint32_t(r[i]) << shift));
        }
    }
}

const uint32_t kRates[12] = {0, 88200, 176400, 192000, 8000, 16000, 22050, 24000, 32000, 44100, 48000, 96000};
const int kDepths[8] = {0, 8, 12, -1, 16, 20, 24, 32};

}  // namespace

void useBasicOnly(bool on)
{
#if FLAC_X86
    kRiceRun = (!on && cpuHasBmi2()) ? riceRunBmi2 : riceRunBasic;
    kHasSse41 = !on && cpuHasSse41();
#endif
    kCrc16 = on ? crc16Plain : pickCrc16();
}

const char *kernelText()
{
#if FLAC_X86
    const bool bmi = kRiceRun == riceRunBmi2;
#else
    const bool bmi = false;
#endif
    const bool fold = kCrc16 != crc16Plain;
#if FLAC_X86
    const bool wide = kHasSse41;
    static char text[160];
    std::snprintf(text, sizeof text, "rice: %s · crc-16: %s · prediction: ring, sse2%s", bmi ? "shlx/lzcnt" : "plain", fold ? "carry-less multiply" : "table",
                  wide ? " + sse4.1 for more than 16 bits" : " only");
    return text;
#elif FLAC_HAVE_RING
    (void)bmi;
    return fold ? "rice: plain · crc-16: carry-less multiply · prediction: ring, neon" : "rice: plain · crc-16: table · prediction: ring, neon";
#else
    (void)bmi;
    return fold ? "rice: plain · crc-16: carry-less multiply · prediction: one sample at a time" : "rice: plain · crc-16: table · prediction: one sample at a time";
#endif
}

// ── ΤΑ ΜΕΤΑΔΕΔΟΜΕΝΑ (8) ──
namespace {
// Τα 34 bytes του STREAMINFO (8.2), όπως είναι γραμμένα: τίποτα δεν κρίνεται εδώ.
void fillStreamInfo(const uint8_t *s, StreamInfo &info)
{
    info.minBlock = (uint32_t(s[0]) << 8) | s[1];
    info.maxBlock = (uint32_t(s[2]) << 8) | s[3];
    info.minFrame = (uint32_t(s[4]) << 16) | (uint32_t(s[5]) << 8) | s[6];
    info.maxFrame = (uint32_t(s[7]) << 16) | (uint32_t(s[8]) << 8) | s[9];
    const uint64_t v = be64(s + 10);
    info.sampleRate = uint32_t(v >> 44);
    info.channels = int((v >> 41) & 7) + 1;
    info.bits = int((v >> 36) & 31) + 1;
    info.totalSamples = v & 0xFFFFFFFFFULL;
    std::memcpy(info.md5, s + 18, 16);
    info.present = true;
}
}  // namespace

bool parseStreamInfo(const uint8_t *body, size_t size, StreamInfo &info)
{
    info = StreamInfo();
    if (!body || size < 34) return false;
    fillStreamInfo(body, info);
    // 8.2: ρυθμός 0 σημαίνει «όχι ήχος»· μπλοκ κάτω από 16 δείγματα δεν επιτρέπεται, ούτε ελάχιστο πάνω από το μέγιστο.
    if (info.sampleRate == 0 || info.minBlock < 16 || info.maxBlock < info.minBlock) { info = StreamInfo(); return false; }
    return true;
}

bool readHead(const uint8_t *d, size_t size, StreamInfo &info, size_t &audioStart)
{
    info = StreamInfo();
    size_t at = 0;
    // Ετικέτα ID3v2 μπροστά (δεν είναι του προτύπου, αλλά υπάρχουν τέτοια αρχεία): 10 bytes κεφαλίδα, μέγεθος σε 4 × 7 μπιτ.
    while (size - at >= 10 && d[at] == 'I' && d[at + 1] == 'D' && d[at + 2] == '3' && d[at + 3] != 0xFF && d[at + 4] != 0xFF
           && !((d[at + 6] | d[at + 7] | d[at + 8] | d[at + 9]) & 0x80)) {
        size_t len = (size_t(d[at + 6]) << 21) | (size_t(d[at + 7]) << 14) | (size_t(d[at + 8]) << 7) | size_t(d[at + 9]);
        len += (d[at + 5] & 0x10) ? 20 : 10;
        if (len > size - at) return false;
        at += len;
    }
    if (size - at >= 4 && std::memcmp(d + at, "fLaC", 4) == 0) {
        at += 4;
        bool broken = false;
        for (;;) {
            if (size - at < 4) { broken = true; break; }
            const bool last = (d[at] & 0x80) != 0;
            const int type = d[at] & 0x7F;
            const size_t len = (size_t(d[at + 1]) << 16) | (size_t(d[at + 2]) << 8) | size_t(d[at + 3]);
            if (type == 127 || len > size - at - 4) { broken = true; break; }
            at += 4;
            if (type == 0 && !info.present && len >= 34) fillStreamInfo(d + at, info);   // 8.2
            at += len;
            if (last) break;
        }
        // Μεταδεδομένα που δεν βγάζουν νόημα: ο ήχος είναι εκεί που αρχίζει το πρώτο καρέ που διαβάζεται.
        audioStart = broken ? findFrame(d, size, at, info) : at;
        return audioStart < size || !broken;
    }
    // Χωρίς «fLaC»: δεκτό μόνο αν βρεθεί καρέ που περιγράφει μόνο του τον ήχο του.
    const size_t f = findFrame(d, size, at, info);
    if (f >= size) return false;
    audioStart = f;
    return true;
}

// ── Η ΚΕΦΑΛΙΔΑ ΤΟΥ ΚΑΡΕ (9.1) ──
Status parseFrameHead(const uint8_t *p, size_t n, const StreamInfo &si, FrameHead &h)
{
    if (n < 6) return Status::Truncated;
    if (p[0] != 0xFF || (p[1] & 0xFE) != 0xF8) return Status::NotAFrame;
    const int bsCode = p[2] >> 4, srCode = p[2] & 15, chCode = p[3] >> 4, depthCode = (p[3] >> 1) & 7;
    if ((p[3] & 1) || bsCode == 0 || srCode == 15 || chCode > 10 || depthCode == 3) return Status::NotAFrame;
    h.variable = (p[1] & 1) != 0;

    // 9.1.5: ο αριθμός, γραμμένος όπως το UTF-8 αλλά ως 36 μπιτ
    size_t i = 4;
    const uint8_t b = p[i++];
    int extra = 0;
    uint64_t v = b;
    if (b >= 0x80) {
        if (b < 0xC0 || b == 0xFF) return Status::NotAFrame;
        extra = 1;
        while (extra < 6 && (b & (0x40 >> extra))) ++extra;
        v = b & (0x3F >> extra);
    }
    const size_t need = i + size_t(extra) + (bsCode == 6 ? 1 : bsCode == 7 ? 2 : 0) + (srCode == 12 ? 1 : srCode >= 13 ? 2 : 0) + 1;
    if (need > n) return Status::Truncated;
    for (int k = 0; k < extra; ++k) {
        const uint8_t c = p[i++];
        if ((c & 0xC0) != 0x80) return Status::NotAFrame;
        v = (v << 6) | (c & 0x3F);
    }
    h.number = v;

    if (bsCode == 1) h.blockSize = 192;                                  // 9.1.1
    else if (bsCode <= 5) h.blockSize = 144u << bsCode;
    else if (bsCode == 6) h.blockSize = uint32_t(p[i++]) + 1;            // 9.1.6
    else if (bsCode == 7) { h.blockSize = ((uint32_t(p[i]) << 8) | p[i + 1]) + 1; i += 2; }
    else h.blockSize = 1u << bsCode;
    if (h.blockSize > 65535) return Status::NotAFrame;                   // το 65536 απαγορεύεται

    if (srCode == 0) h.sampleRate = si.sampleRate;                       // 9.1.2
    else if (srCode < 12) h.sampleRate = kRates[srCode];
    else if (srCode == 12) h.sampleRate = uint32_t(p[i++]) * 1000;       // 9.1.7
    else {
        h.sampleRate = (uint32_t(p[i]) << 8) | p[i + 1];
        if (srCode == 14) h.sampleRate *= 10;
        i += 2;
    }

    if (chCode < 8) { h.channels = chCode + 1; h.assignment = 0; }       // 9.1.3
    else { h.channels = 2; h.assignment = chCode; }

    if (crc8(p, i) != p[i]) return Status::NotAFrame;                    // 9.1.8
    h.headBytes = int(i) + 1;

    h.bits = depthCode ? kDepths[depthCode] : si.bits;                   // 9.1.4
    if (h.bits < 4 || h.bits > 32) return Status::Unsupported;           // «όπως το STREAMINFO», και δεν υπάρχει
    return Status::Ok;
}

size_t findFrame(const uint8_t *p, size_t n, size_t from, const StreamInfo &si)
{
    size_t i = from;
    while (i + 6 <= n) {
        const void *hit = std::memchr(p + i, 0xFF, n - i - 5);
        if (!hit) break;
        i = size_t(static_cast<const uint8_t *>(hit) - p);
        FrameHead h;
        if ((p[i + 1] & 0xFE) == 0xF8 && parseFrameHead(p + i, n - i, si, h) == Status::Ok) return i;
        ++i;
    }
    return n;
}

Decoder::Decoder(const StreamInfo &si) : si_(si) {}

// Η πρόβλεψη όλων των καναλιών ενός καρέ (όχι του δρόμου 64 μπιτ για ήχο 32 μπιτ με διαφορά). Αφήνει στο cur_[c] πού
// βρίσκονται τα δείγματα του καναλιού c.
void Decoder::predict(const void *subs, const FrameHead &h)
{
    const Sub *sub = static_cast<const Sub *>(subs);
    const uint32_t n = h.blockSize;
    bool done[8] = {};
    for (int c = 0; c < h.channels; ++c) cur_[c] = plane_[c].data();
#if FLAC_HAVE_RING
    // Ο δρόμος «σε δαχτυλίδι» για τους γραμμικούς προβλέπτες: τα κανάλια ανά δύο, πλεγμένα. Στενή μορφή όπου χωράει, αλλιώς
    // φαρδιά (αν ο επεξεργαστής την έχει)· ό,τι δεν παίρνει κανένα από τα δύο μένει στον απλό δρόμο.
    if (n >= 128) {
        int list[8], m = 0;
        bool fits[8] = {};
        for (int c = 0; c < h.channels; ++c)
            if (sub[c].kind == Sub::Lpc) {
                fits[c] = ringNarrowFits(sub[c], h.bits);
                if (fits[c] || FLAC_RING_WIDE) list[m++] = c;
            }
        auto fit = [&](int c) {
            if (alt_[c].size() < size_t(n) + 16) alt_[c].resize(size_t(n) + 16);
        };
        auto take = [&](int c) {
            cur_[c] = alt_[c].data();
            done[c] = true;
        };
        auto alone = [&](int c, bool narrowOk) {
            if (sub[c].order < kRingAloneFrom) return;   // μόνο του, σε μικρή τάξη, ο απλός δρόμος είναι εξίσου γρήγορος
            const int size = ringSize(sub[c].order);
            fit(c);
            if (narrowOk) {
                if (ringNarrowOne(size, plane_[c].data(), alt_[c].data(), sub[c], n)) {
                    take(c);
                    if (stats_) ++stats_->ringOne;
                    return;
                }
                if (stats_) ++stats_->ringRedone;
            }
            if (FLAC_RING_WIDE) {
                ringWideOne(size, plane_[c].data(), alt_[c].data(), sub[c], n);
                take(c);
                if (stats_) ++stats_->ringOne;
            }
        };
        for (int k = 0; k + 1 < m; k += 2) {
            const int ca = list[k], cb = list[k + 1];
            const int size = ringSize(std::max(sub[ca].order, sub[cb].order));
            fit(ca);
            fit(cb);
            if (fits[ca] && fits[cb]) {
                bool okA = false, okB = false;
                ringNarrowPair(size, plane_[ca].data(), alt_[ca].data(), sub[ca], plane_[cb].data(), alt_[cb].data(), sub[cb], n, okA, okB);
                for (int which = 0; which < 2; ++which) {
                    const int c = which ? cb : ca;
                    if (which ? okB : okA) {
                        take(c);
                    } else {   // ένα δείγμα δεν χώρεσε σε 16 μπιτ: το κανάλι ξαναγίνεται
                        if (stats_) ++stats_->ringRedone;
                        alone(c, false);
                    }
                }
            } else {
                ringWidePair(size, plane_[ca].data(), alt_[ca].data(), sub[ca], plane_[cb].data(), alt_[cb].data(), sub[cb], n);
                take(ca);
                take(cb);
            }
            if (stats_) ++stats_->ringPairs;
        }
        if (m & 1) alone(list[m - 1], fits[list[m - 1]]);
    }
#endif
    for (int c = 0; c < h.channels; ++c)
        if (!done[c]) predictScalar(plane_[c].data(), n, sub[c]);
}

// ── ΕΝΑ ΚΑΡΕ (9) ──
Status Decoder::decode(const uint8_t *p, size_t n, size_t &used)
{
    used = 0;
    span_ = n < 16 ? n : 16;   // η κεφαλίδα θέλει ως 16 bytes
    filled_ = 0;
    const double tFrame = stats_ ? nowNs() : 0;
    FrameHead h;
    const Status hs = parseFrameHead(p, n, si_, h);
    if (hs != Status::Ok) return hs;
    const uint32_t bs = h.blockSize;
    for (int c = 0; c < h.channels; ++c)
        if (plane_[c].size() < size_t(bs) + 16) plane_[c].resize(size_t(bs) + 16);

    BitReader br(p + h.headBytes, n - size_t(h.headBytes));
    const bool wide = h.bits == 32 && h.assignment != 0;   // το κανάλι διαφοράς θέλει 33 μπιτ
    if (wide)
        for (int c = 0; c < 2; ++c)
            if (wide_[c].size() < size_t(bs) + 16) wide_[c].resize(size_t(bs) + 16);
    Sub sub[8];
    for (int c = 0; c < h.channels; ++c) {
        const bool side = (h.assignment == 8 && c == 1) || (h.assignment == 9 && c == 0) || (h.assignment == 10 && c == 1);
        filled_ += bs;
        const Status r = wide ? parseSub<int64_t>(br, wide_[c].data(), bs, h.bits + (side ? 1 : 0), sub[c], stats_)
                              : parseSub<int32_t>(br, plane_[c].data(), bs, h.bits + (side ? 1 : 0), sub[c], stats_);
        if (r != Status::Ok) {
            const size_t at = size_t(h.headBytes) + (br.bitsUsed() + 7) / 8;
            span_ = at < n ? at : n;
            return r;
        }
    }

    // 9.3: μηδενικά ως το όριο του byte, και το CRC-16 όλου του καρέ
    const size_t body = (br.bitsUsed() + 7) / 8;
    const size_t end = size_t(h.headBytes) + body;
    span_ = end + 2 < n ? end + 2 : n;
    if (br.overrun() || end + 2 > n) return Status::Truncated;
    const double tCrc = stats_ ? nowNs() : 0;
    if (checkCrc_) {
        const uint16_t want = uint16_t((uint16_t(p[end]) << 8) | p[end + 1]);
        if (crc16(p, end) != want) return Status::BadCrc;
    }
    const double tPredict = stats_ ? nowNs() : 0;

    // Η πρόβλεψη, και μετά τα αχρησιμοποίητα μπιτ (9.2.2): πριν από την αποσυσχέτιση.
    if (wide) {
        for (int c = 0; c < 2; ++c) {
            predictScalar(wide_[c].data(), bs, sub[c]);
            if (sub[c].wasted)
                for (uint32_t i = 0; i < bs; ++i) wide_[c][i] = int64_t(uint64_t(wide_[c][i]) << sub[c].wasted);
            cur_[c] = plane_[c].data();
        }
    } else {
        predict(sub, h);
        for (int c = 0; c < h.channels; ++c)
            if (sub[c].wasted) {
                int32_t *s = cur_[c];
                for (uint32_t i = 0; i < bs; ++i) s[i] = int32_t(uint32_t(s[i]) << sub[c].wasted);
            }
    }
    const double tStereo = stats_ ? nowNs() : 0;

    // 4.2: αποσυσχέτιση στέρεο. Με μέσο m και διαφορά d: δεξί = m - (d >> 1), αριστερό = δεξί + d (το ίδιο με «διπλασίασε το
    // μέσο, πρόσθεσε το κάτω μπιτ της διαφοράς, πρόσθεσε ή αφαίρεσε τη διαφορά και μισό», χωρίς το ενδιάμεσο που ξεχειλίζει).
    if (wide) {
        const int64_t *a = wide_[0].data(), *b = wide_[1].data();
        int32_t *l = cur_[0], *r = cur_[1];
        for (uint32_t i = 0; i < bs; ++i) {
            int64_t left, right;
            if (h.assignment == 8) { left = a[i]; right = a[i] - b[i]; }
            else if (h.assignment == 9) { right = b[i]; left = a[i] + b[i]; }
            else { right = a[i] - (b[i] >> 1); left = right + b[i]; }
            l[i] = int32_t(uint32_t(uint64_t(left)));
            r[i] = int32_t(uint32_t(uint64_t(right)));
        }
    } else if (h.assignment) {
        stereoUndo(cur_[0], cur_[1], bs, h.assignment);
    }
    if (stats_) {
        ++stats_->frames;
        ++stats_->assignment[h.assignment];
        const double tEnd = nowNs();
        stats_->nsFrame += tEnd - tFrame;
        stats_->nsResidual += tCrc - tFrame;
        stats_->nsCrc += tPredict - tCrc;
        stats_->nsPredict += tStereo - tPredict;
        stats_->nsStereo += tEnd - tStereo;
    }
    head_ = h;
    used = end + 2;
    return Status::Ok;
}

void Decoder::write(Out fmt, int shift, void *const *dst) const
{
    const uint32_t n = head_.blockSize;
    const int ch = head_.channels;
    if (ch == 2 && (fmt == Out::S16Interleaved || fmt == Out::S32Interleaved) && shift >= 0 && shift <= 15) {
        interleave2(cur_[0], cur_[1], n, shift, fmt == Out::S32Interleaved, dst[0]);
        return;
    }
    switch (fmt) {
    case Out::S16Interleaved: {
        int16_t *o = static_cast<int16_t *>(dst[0]);
        for (int c = 0; c < ch; ++c) {
            const int32_t *s = cur_[c];
            for (uint32_t i = 0; i < n; ++i) o[size_t(i) * size_t(ch) + size_t(c)] = int16_t(uint16_t(uint32_t(s[i]) << shift));
        }
        break;
    }
    case Out::S16Planar:
        for (int c = 0; c < ch; ++c) {
            const int32_t *s = cur_[c];
            int16_t *o = static_cast<int16_t *>(dst[c]);
            for (uint32_t i = 0; i < n; ++i) o[i] = int16_t(uint16_t(uint32_t(s[i]) << shift));
        }
        break;
    case Out::S32Interleaved: {
        int32_t *o = static_cast<int32_t *>(dst[0]);
        for (int c = 0; c < ch; ++c) {
            const int32_t *s = cur_[c];
            for (uint32_t i = 0; i < n; ++i) o[size_t(i) * size_t(ch) + size_t(c)] = int32_t(uint32_t(s[i]) << shift);
        }
        break;
    }
    case Out::S32Planar:
        for (int c = 0; c < ch; ++c) {
            const int32_t *s = cur_[c];
            int32_t *o = static_cast<int32_t *>(dst[c]);
            for (uint32_t i = 0; i < n; ++i) o[i] = int32_t(uint32_t(s[i]) << shift);
        }
        break;
    }
}

}  // namespace orpheic::flac
