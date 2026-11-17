// Ο πυρήνας Rice του flacdecoder.cpp. Συμπεριλαμβάνεται μία φορά για κάθε μορφή του:
// η ίδια συνάρτηση μεταγλωττίζεται με άλλες εντολές του επεξεργαστή (βασικές, ή SHLX/SHRX/LZCNT όπου υπάρχουν) και για
// άλλο πλάτος δείγματος. Πριν από κάθε #include ορίζονται:
//   FLAC_RICE_NAME               το όνομα της συνάρτησης
//   FLAC_RICE_ATTR               ιδιότητες μεταγλώττισης (π.χ. target("bmi2,lzcnt")), ή κενό
//   FLAC_RICE_T                  ο τύπος του δείγματος εξόδου (int32_t ή int64_t)
//   FLAC_RICE_LEAD(v, z, onZero) βάζει στο z πόσα μηδενικά έχει μπροστά το v· αν το v είναι 0, είτε εκτελεί το onZero είτε
//                                δίνει z = 64 (και τα δύο οδηγούν στον ίδιο δρόμο)
//   FLAC_RICE_TOP(v, top, onZero) (αντί για το LEAD) βάζει στο top τη θέση του πιο σημαντικού 1 του v (0..63)· αν το v
//                                είναι 0 εκτελεί το onZero
//   FLAC_RICE_SHL(v, n), FLAC_RICE_SHR(v, n)   ολίσθηση 64 μπιτ, n στο 0..63
//   FLAC_RICE_UNFOLD(p, count)   το ξεδίπλωμα του προσήμου για `count` τιμές, επιτόπου
//   FLAC_RICE_KCONST             (προαιρετικό) η συνάρτηση γίνεται πρότυπο με το k ως σταθερά: μία συνάρτηση ανά k
//
// `count` σύμβολα Rice με παράμετρο k (RFC 9639, 9.2.7.2): z μηδενικά, ένα 1, και k μπιτ r. Η τιμή είναι (z << k) | r.
//
// ΠΩΣ ΚΕΡΔΙΖΕΙ ΧΡΟΝΟ. Κάθε σύμβολο πρέπει να ξέρει πού τελείωσε το προηγούμενο: αυτή η αλυσίδα ορίζει την ταχύτητα. Μέσα
// στην αλυσίδα μένουν ΔΥΟ πράξεις: το μέτρημα των μηδενικών και ΜΙΑ ολίσθηση κατά όλο το μήκος του συμβόλου
// (s = z + 1 + k). Η τιμή βγαίνει δίπλα της, από τα πάνω s μπιτ του αποθέματος: εκεί είναι z μηδενικά, το 1 και το r,
// δηλαδή ο αριθμός (1 << k) | r, άρα τιμή = αυτός + ((z - 1) << k), modulo 2^32. Το ξεδίπλωμα του προσήμου δεν γίνεται εδώ
// αλλά μετά, για όλο το τμήμα μαζί, 4 τιμές τη φορά. Το γέμισμα του αποθέματος (ένα διάβασμα 8 bytes) γίνεται μία φορά ανά
// ομάδα 1 ως 4 συμβόλων, όσα χωράνε συνήθως σε 56 μπιτ για αυτό το k· σύμβολο που δεν χωράει κλείνει την ομάδα νωρίτερα,
// και σύμβολο που δεν χωράει ούτε σε γεμάτο απόθεμα (μακριά σειρά μηδενικών) ή βρίσκεται στην άκρη των δεδομένων πηγαίνει
// από τον προσεκτικό δρόμο (riceSlow), που ελέγχει το τέλος.
#ifdef FLAC_RICE_KCONST   // το k είναι σταθερά της μεταγλώττισης: οι ολισθήσεις με το k γίνονται με σταθερό μετρητή
template <unsigned K>
FLAC_NOINLINE FLAC_RICE_ATTR bool FLAC_RICE_NAME(BitReader &br, FLAC_RICE_T *out, uint32_t count, uint64_t *slow)
{
    constexpr unsigned k = K;
#else
FLAC_NOINLINE FLAC_RICE_ATTR bool FLAC_RICE_NAME(BitReader &br, FLAC_RICE_T *out, uint32_t count, unsigned k, uint64_t *slow)
{
#endif
    const uint8_t *const base = br.base;
    const size_t n = br.n;
    size_t pos = br.pos;
    uint64_t buf = br.buf;
    unsigned cnt = unsigned(br.cnt);
    const unsigned perGroup = k <= 8 ? 4 : k <= 12 ? 3 : k <= 20 ? 2 : 1;
    FLAC_RICE_T *o = out;
    FLAC_RICE_T *const end = out + count;

// Ένα σύμβολο από το απόθεμα, χωρίς γέμισμα. Αν δεν χωράει στα `cnt` μπιτ που έχουμε (ή το απόθεμα είναι όλο μηδενικά)
// εκτελεί το `onShort`. Όταν χωράει, 1 ≤ s ≤ 63, άρα και οι δύο ολισθήσεις είναι μέσα στα όρια.
#ifdef FLAC_RICE_TOP
// Η ίδια αριθμητική γραμμένη με τη ΘΕΣΗ του πρώτου 1 (top = 63 - z), που είναι αυτό που δίνει η εντολή BSR των x86-64
// χωρίς LZCNT: s = 64 + k - top, τα πάνω s μπιτ είναι buf >> (top - k), και (z - 1) = 62 - top. Έτσι από το BSR ως την
// ολίσθηση του αποθέματος μεσολαβεί μία αφαίρεση.
#define FLAC_RICE_SYMBOL(onShort)                                                          \
    {                                                                                      \
        unsigned top;                                                                      \
        FLAC_RICE_TOP(buf, top, onShort)                                                   \
        const unsigned s = (64u + k) - top;                                                \
        if (FLAC_UNLIKELY(s > cnt)) { onShort; }                                           \
        const uint32_t u = uint32_t(FLAC_RICE_SHR(buf, top - k)) + ((62u - top) << k);     \
        buf = FLAC_RICE_SHL(buf, s);                                                       \
        cnt -= s;                                                                          \
        *o++ = FLAC_RICE_T(int32_t(u));                                                    \
    }
#else
#define FLAC_RICE_SYMBOL(onShort)                                                          \
    {                                                                                      \
        unsigned z;                                                                        \
        FLAC_RICE_LEAD(buf, z, onShort)                                                    \
        const unsigned s = z + 1 + k;                                                      \
        if (FLAC_UNLIKELY(s > cnt)) { onShort; }                                           \
        const uint32_t u = uint32_t(FLAC_RICE_SHR(buf, 64 - s)) + ((z - 1u) << k);         \
        buf = FLAC_RICE_SHL(buf, s);                                                       \
        cnt -= s;                                                                          \
        *o++ = FLAC_RICE_T(int32_t(u));                                                    \
    }
#endif

    // Το γέμισμα είναι το ίδιο με του BitReader::refill (flacdecoder.cpp· εκεί και η αναφορά στη δημοσίευση της τεχνικής).
    while (o < end) {
        if (FLAC_LIKELY(pos <= n)) {
            buf |= FLAC_RICE_SHR(be64(base + pos), cnt);
            pos += size_t(63 - cnt) >> 3;
            cnt |= 56;
            const size_t left = size_t(end - o);
            const unsigned group = left < perGroup ? unsigned(left) : perGroup;
            FLAC_RICE_SYMBOL(goto careful)   // το πρώτο της ομάδας: αν δεν χωράει ούτε τώρα, θέλει τον προσεκτικό δρόμο
            switch (group) {
            case 4: FLAC_RICE_SYMBOL(continue)
                [[fallthrough]];
            case 3: FLAC_RICE_SYMBOL(continue)
                [[fallthrough]];
            case 2: FLAC_RICE_SYMBOL(continue)
                [[fallthrough]];
            default: break;
            }
            continue;
        }
    careful:
        br.pos = pos;
        br.buf = buf;
        br.cnt = int(cnt);
        uint32_t u;
        if (!riceSlow(br, k, u)) return false;
        pos = br.pos;
        buf = br.buf;
        cnt = unsigned(br.cnt);
        *o++ = FLAC_RICE_T(int32_t(u));
        if (slow) ++*slow;
    }
#undef FLAC_RICE_SYMBOL
    br.pos = pos;
    br.buf = buf;
    br.cnt = int(cnt);
    FLAC_RICE_UNFOLD(out, count);
    return true;
}
