// Ο γρήγορος δρόμος του mp3huff.cpp. Συμπεριλαμβάνεται μία φορά για κάθε μορφή του: η ίδια συνάρτηση μεταγλωττίζεται με
// άλλες εντολές ολίσθησης. Πριν από κάθε #include ορίζονται:
//   HUFF_KERNEL_NAME             το όνομα της συνάρτησης
//   HUFF_KERNEL_ATTR             ιδιότητες μεταγλώττισης (π.χ. target("bmi2")), ή κενό
//   HUFF_SHL(v, n), HUFF_SHR(v, n)   ολίσθηση 64 μπιτ με μεταβλητό n στο 0 … 63
//
// ΤΟ ΠΑΡΑΘΥΡΟ. `buf`: τα επόμενα μπιτ της ροής, το πρώτο στο πιο σημαντικό μπιτ. Τα πάνω `cnt` είναι «δικά μας»· ό,τι
// βρίσκεται πιο κάτω είναι τα αμέσως επόμενα μπιτ της ροής (θα ξαναμπούν ίδια στο επόμενο γέμισμα) ή μηδέν. `p`: το επόμενο
// byte που δεν έχει μετρηθεί στο cnt. Ισχύει πάντα 8·p = θέση + cnt, και θέση = limit − left.
//
// ΠΟΤΕ ΕΙΝΑΙ ΑΣΦΑΛΕΣ ΤΟ ΓΕΜΙΣΜΑ. Το γέμισμα διαβάζει 8 bytes από το p. Γίνεται μόνο σε σημείο όπου η θέση είναι ≤ limit (στην
// αρχή, και μετά από βήμα που έγινε δεκτό), και cnt ≤ 63. Άρα 8·p ≤ limit + 63 ≤ 8·(bytes δεδομένων), δηλαδή το p είναι το
// πολύ το πρώτο byte του περιθωρίου, και τα 8 bytes από εκεί είναι τα 8 bytes που εγγυάται ο καλών.
//
// ΠΟΣΑ ΒΗΜΑΤΑ ΑΝΑ ΓΕΜΙΣΜΑ. Μετά το γέμισμα cnt ≥ 56. Κάθε βρόχος κάνει τόσα βήματα ανά γέμισμα ώστε (α) κανένα βήμα να μη
// θέλει περισσότερα μπιτ από όσα έμειναν, και (β) να μένουν αρκετά για τον δείκτη της ΠΡΩΤΗΣ αναζήτησης της επόμενης
// ομάδας, που υπολογίζεται ΠΡΙΝ από το γέμισμα. Το (β) κόβει την αλυσίδα που ορίζει τον χρόνο: η αναζήτηση δεν περιμένει το
// διάβασμα των 8 bytes, γιατί το γέμισμα δεν αλλάζει τα πάνω cnt μπιτ του buf. Ό,τι άλλο χρειάζεται το πρώτο βήμα
// (υποπίνακας, πρόσημα) διαβάζεται μετά το γέμισμα. Οι αριθμοί κάθε βρόχου είναι γραμμένοι πάνω του.
//
// ΤΙ ΕΙΝΑΙ ΕΝΑ ΒΗΜΑ. Αναζήτηση, γράψιμο των τιμών στο is, κατανάλωση των μπιτ, και ΤΟΤΕ ο έλεγχος: αν η θέση πέρασε το
// limit, το βήμα ξεγίνεται (γραμμή και θέση γυρίζουν πίσω) και ο δρόμος τελειώνει εκεί. Οι τιμές που πρόλαβαν να γραφτούν
// είναι μέσα στο is[0 … 575] και ο προσεκτικός δρόμος θα τις ξαναγράψει ή θα τις μηδενίσει.
HUFF_NOINLINE HUFF_KERNEL_ATTR void HUFF_KERNEL_NAME(const Tables &T, const uint8_t *data, int64_t limit, const Plan &plan, int16_t *is,
                                                    State &st)
{
    ptrdiff_t l = st.l;   // στο πλάτος των διευθύνσεων, ώστε το is[l] να μη θέλει επέκταση προσήμου σε κάθε βήμα
    int region = 0;       // οι περιοχές που έχουν ήδη τελειώσει (l ≥ τέλος τους) προσπερνιούνται μόνες τους
    unsigned nLin = st.nLin;
    int64_t left = limit - int64_t(st.pos);   // ≥ 0: το εγγυάται ο καλών
    // Το πρώτο διάβασμα: 8 bytes από το byte της θέσης (είναι byte δεδομένων ή το πρώτο του περιθωρίου, γιατί θέση ≤ limit).
    // Με θέση σε όριο byte δηλώνουμε 56 μπιτ δικά μας και όχι 64, ώστε το cnt να μένει ≤ 63 (η ολίσθηση του γεμίσματος δεν
    // ορίζεται για 64)· αλλιώς όσα έμειναν από το πρώτο byte και τα επόμενα 7.
    const unsigned bit = unsigned(st.pos & 7u);
    uint64_t buf = be64(data + (st.pos >> 3)) << bit;
    unsigned cnt = 63u - ((bit + 7u) & 7u);
    size_t p = (st.pos + cnt) >> 3;

// Το γέμισμα χωρίς διακλάδωση είναι δημοσιευμένη τεχνική: Fabian Giesen, «Reading bits in far too many ways (part 2)», ιστολόγιο
// «The ryg blog», Φεβρουάριος 2018.
#define HUFF_REFILL()                             \
    {                                             \
        buf |= HUFF_SHR(be64(data + p), cnt);     \
        p += (63u - cnt) >> 3;                    \
        cnt |= 56u;                               \
    }

// Κατανάλωση `total` μπιτ και `lines` γραμμών, και ο έλεγχος του βήματος. `undo`: ό,τι άλλο πρέπει να ξεγίνει.
// ΒΛΑΒΗ 7: ο γρήγορος δρόμος δέχεται βήμα που περνά το όριο κατά ένα μπιτ.
#define HUFF_TAKE(total, lines, undo)                     \
    buf = HUFF_SHL(buf, (total));                         \
    cnt -= (total);                                       \
    left -= int64_t(total);                               \
    l += ptrdiff_t(lines);                                \
    if (HUFF_UNLIKELY(left < (HUFF_MUT(7) ? -1 : 0))) {   \
        left += int64_t(total);                           \
        l -= ptrdiff_t(lines);                            \
        undo;                                             \
        goto handover;                                    \
    }

// Ένα ζεύγος χωρίς linbits από τον G, στο is[l], is[l + 1]. `i0`: ο δείκτης στη ρίζα· `rootBits`: το n της ρίζας. Αν η θέση
// είναι παραπομπή, ο υποπίνακας δείχνεται από τα s μπιτ μετά τα πρώτα n. Το πρώτο byte της θέσης είναι τα μπιτ του ζεύγους
// (κώδικας και πρόσημα), το δεύτερο το μήκος του κώδικα: η διαφορά τους είναι πόσες τιμές δεν είναι 0, και τα πρόσημα είναι
// τα μπιτ αμέσως μετά τον κώδικα (storePair).
#define HUFF_PLAIN_PAIR(i0, rootBits, total)                                                                      \
    {                                                                                                             \
        uint32_t e = G[i0];                                                                                       \
        HUFF_COV(gCov.gFast[set.gBase + (i0)]);                                                                   \
        if (HUFF_UNLIKELY(e & kLink)) {                                                                           \
            const uint32_t i1 = (e >> 16) + uint32_t(HUFF_SHR(HUFF_SHL(buf, (rootBits)), 64u - (e & 0xFFu)));     \
            e = G[i1];                                                                                            \
            HUFF_COV(gCov.gFast[set.gBase + i1]);                                                                 \
        }                                                                                                         \
        total = e & 0xFFu;                                                                                        \
        const unsigned len = (e >> 8) & 0xFFu;                                                                    \
        storePair(is + l, e, HUFF_SHL(buf, len), total - len);                                                    \
    }

// Ένα βήμα του D. Θέση με ζεύγη: τέσσερις τιμές με μία αποθήκευση (αν η θέση έχει μόνο ένα ζεύγος οι δύο τελευταίες είναι
// μηδέν και θα ξαναγραφτούν από το επόμενο βήμα), και προχωράμε 2 ή 4 γραμμές (2 + 2·μπιτ 31). Θέση 0 (το πρώτο ζεύγος
// δεν χωρά στον δείκτη): ένα ζεύγος από τον G, αφού πρώτα γεμίσει το παράθυρο (βλ. τον βρόχο).
#define HUFF_D_STEP()                                                  \
    {                                                                  \
        const uint32_t m = D[di];                                      \
        unsigned total = m & 0xFFu, lines = 2u + 2u * (m >> 31);       \
        HUFF_COV(gCov.d[set.d][di]);                                   \
        if (HUFF_LIKELY(m != 0u)) {                                    \
            expand4(is + l, m);                                        \
        } else {                                                       \
            HUFF_REFILL()                                              \
            const uint32_t i0 = uint32_t(HUFF_SHR(buf, 64u - gb));     \
            HUFF_PLAIN_PAIR(i0, gb, total)                             \
        }                                                              \
        HUFF_TAKE(total, lines, (void)0)                               \
    }

    for (; region < 3; ++region) {
        const ptrdiff_t end = plan.end[region];
        if (l >= end) continue;
        const unsigned sel = plan.sel[region];
        const unsigned kind = T.kind[sel];
        if (kind == kZero) {   // πίνακας 0: μηδενικά χωρίς κανένα μπιτ
            std::memset(is + l, 0, size_t(end - l) * sizeof(int16_t));
            l = end;
            continue;
        }
        if (kind == kBad) goto handover;   // ο προσεκτικός δρόμος θα φτάσει στο ίδιο ζεύγος και θα επιστρέψει −1
        const Set &set = T.set[T.setOf[sel]];
        const uint32_t *const G = T.g + set.gBase;
        HUFF_REFILL()   // κάθε περιοχή αρχίζει με cnt ≥ 56
        if (kind == kDirect) {
            // Τρία βήματα ανά γέμισμα. Ένα βήμα από τον D παίρνει ως kDBits μπιτ: 56 − 3·kDBits ≥ kDBits. Ένα βήμα από τον G
            // (θέση 0 του D) μπορεί να θέλει ως 21 (ο μακρύτερος κώδικας είναι 19, και δύο πρόσημα): γεμίζει πρώτα το ίδιο
            // το παράθυρο (είμαστε μετά από βήμα που έγινε δεκτό, άρα το γέμισμα είναι ασφαλές), και μετά από αυτό μένουν
            // 56 − 21 = 35, που φτάνουν για τα ως δύο βήματα που απομένουν στην ομάδα και τον επόμενο δείκτη: 35 ≥ 3·kDBits.
            const uint32_t *const D = T.d[set.d];
            const unsigned gb = set.gBits;
            while (l + 4 <= end) {   // χωράνε δύο ζεύγη στην περιοχή: η θέση του D χρησιμοποιείται ολόκληρη
                uint32_t di = uint32_t(buf >> (64u - kDBits));
                HUFF_REFILL()
                HUFF_D_STEP()
                if (l + 4 > end) break;
                di = uint32_t(buf >> (64u - kDBits));
                HUFF_D_STEP()
                if (l + 4 > end) break;
                di = uint32_t(buf >> (64u - kDBits));
                HUFF_D_STEP()
            }
            if (l < end) {
                // Μένει ακριβώς ένα ζεύγος στην περιοχή. Η θέση του D δεν λέει πόσα μπιτ είναι μόνο το πρώτο της ζεύγος,
                // οπότε αυτό το ένα ζεύγος διαβάζεται από τον G (το πολύ μία φορά ανά περιοχή).
                HUFF_REFILL()
                const uint32_t i0 = uint32_t(HUFF_SHR(buf, 64u - gb));
                unsigned total;
                HUFF_PLAIN_PAIR(i0, gb, total)
                HUFF_TAKE(total, 2, (void)0)
            }
        } else if (kind == kPlain) {
            // Δύο ζεύγη ανά γέμισμα: το καθένα ως 21 μπιτ, άρα το δεύτερο έχει τουλάχιστον 35 και μετά μένουν 14 ≥ kGBits.
            while (l < end) {
                uint32_t i0 = uint32_t(buf >> (64u - kGBits));
                unsigned total;
                HUFF_REFILL()
                HUFF_PLAIN_PAIR(i0, kGBits, total)
                HUFF_TAKE(total, 2, (void)0)
                if (l >= end) break;
                i0 = uint32_t(buf >> (64u - kGBits));
                HUFF_PLAIN_PAIR(i0, kGBits, total)
                HUFF_TAKE(total, 2, (void)0)
            }
        } else {
            // Με linbits: ένα ζεύγος ανά γέμισμα. Παίρνει το πολύ 56 − kGBits μπιτ (το ελέγχει το χτίσιμο των πινάκων για
            // κάθε αριθμό πίνακα· στους πίνακες του προτύπου το μακρύτερο είναι 36), άρα μένουν αρκετά για τον επόμενο δείκτη.
            // Το συνολικό μήκος βγαίνει από δύο πεδία της θέσης του πίνακα: το μήκος χωρίς linbits (κάτω byte), και το πόσο
            // ολισθαίνει το 2·linbits (πάνω μπιτ: δίνει 0, linbits ή 2·linbits). Έτσι το επόμενο ζεύγος δεν περιμένει ούτε
            // τις τιμές αυτού ούτε την απόφαση «έχει 15 ή όχι».
            // Ζεύγος χωρίς καμία τιμή 15: ίδιο με ένα ζεύγος χωρίς linbits (storePair). Είναι τα περισσότερα: ο πίνακας με
            // linbits διαλέγεται επειδή ΚΑΠΟΙΑ τιμή της περιοχής είναι πάνω από 15, και αυτές μαζεύονται στις χαμηλές
            // συχνότητες, οπότε η διακλάδωση προβλέπεται.
            // Ζεύγος με τιμή 15: όλα από το ίδιο παράθυρο, χωρίς άλλη διακλάδωση. Η σειρά των μπιτ (2.4.1.7): κώδικας,
            // linbits του x, πρόσημο του x, linbits του y, πρόσημο του y.
            // «(w >> 1) >> (63 − n)» = τα πάνω n μπιτ του w, και 0 για n = 0 (η ολίσθηση κατά 64 δεν ορίζεται).
            // «lin & (0 − σημαία)» = linbits αν η τιμή είναι 15, αλλιώς 0.
            const unsigned lin = T.lin[sel];
            const uint64_t lin2 = 2u * lin;
            while (l < end) {
                const uint32_t i0 = uint32_t(buf >> (64u - kGBits));
                HUFF_REFILL()
                uint32_t e = G[i0];
                unsigned total = (e & 0xFFu) + unsigned(HUFF_SHR(lin2, e >> 28));
                HUFF_COV(gCov.gFast[set.gBase + i0]);
                if (HUFF_UNLIKELY(e & kLink)) {
                    const uint32_t i1 = (e >> 16) + uint32_t(HUFF_SHR(HUFF_SHL(buf, kGBits), 64u - (e & 0xFFu)));
                    e = G[i1];
                    total = (e & 0xFFu) + unsigned(HUFF_SHR(lin2, e >> 28));
                    HUFF_COV(gCov.gFast[set.gBase + i1]);
                }
                const unsigned len = (e >> 8) & 31u;
                uint64_t w = HUFF_SHL(buf, len);
                if (HUFF_LIKELY((e & 0x6000u) == 0u)) {
                    storePair(is + l, e & 0x0FFFFFFFu, w, (e & 0xFFu) - len);
                    HUFF_TAKE(total, 2, (void)0)
                } else {
                    const uint32_t fx = (e >> 13) & 1u, fy = (e >> 14) & 1u;
                    const unsigned lx = lin & (0u - fx), ly = lin & (0u - fy);
                    uint32_t x = (e >> 16) & 15u, y = (e >> 24) & 15u;
                    x += uint32_t(HUFF_SHR(w >> 1, 63u - lx));
                    if (!HUFF_MUT(22)) w = HUFF_SHL(w, lx);   // ΒΛΑΒΗ 22: το πρόσημο του x διαβάζεται πριν από τα linbits του
                    const uint32_t sx = 0u - uint32_t(w >> 63);
                    if (HUFF_MUT(22)) w = HUFF_SHL(w, lx);
                    w = HUFF_SHL(w, unsigned(x != 0u));
                    y += uint32_t(HUFF_SHR(w >> 1, 63u - ly));
                    w = HUFF_SHL(w, ly);
                    const uint32_t sy = 0u - uint32_t(w >> 63);
                    is[l] = int16_t((x ^ sx) - sx);
                    is[l + 1] = int16_t((y ^ sy) - sy);
                    nLin += fx + fy;
                    HUFF_TAKE(total, 2, nLin -= fx + fy)
                }
            }
        }
    }
    // Οι τετράδες: ως τέσσερις ανά γέμισμα. Μία παίρνει το πολύ 6 + 4 = 10 μπιτ: 56 − 40 = 16 ≥ 6 για την επόμενη αναζήτηση.
    // Ο βρόχος του exact είναι «όσο l + 4 ≤ 576 και θέση < gend»· εδώ «θέση < limit» (left > 0), που είναι το ίδιο όσο
    // limit = gend, και αλλιώς σταματά νωρίτερα και συνεχίζει άλλος (η δεύτερη φάση ή ο προσεκτικός δρόμος).
    // Οι τέσσερις τιμές με τα πρόσημά τους: από τον κοινό πίνακα, με δείκτη ποιες είναι 1 και τα 4 μπιτ μετά τον κώδικα.
    {
        const uint16_t *const Q = T.quad[plan.count1];
#define HUFF_QUAD()                                                                                     \
    {                                                                                                   \
        const unsigned q = Q[qi];                                                                       \
        const unsigned total = q & 0xFFu;                                                               \
        const unsigned vi = ((q >> 8) & 0xF0u) | unsigned(HUFF_SHL(buf, (q >> 8) & 15u) >> 60);         \
        std::memcpy(is + l, T.quadVal[vi], 4 * sizeof(int16_t));                                        \
        HUFF_COV(gCov.quadFast[plan.count1][qi]);                                                       \
        HUFF_COV(gCov.quadValFast[vi]);                                                                 \
        HUFF_TAKE(total, 4, (void)0)                                                                    \
    }
        HUFF_REFILL()
        while (l + 4 <= 576 && left > 0) {
            unsigned qi = unsigned(buf >> 58);
            HUFF_REFILL()
            HUFF_QUAD()
            if (!(l + 4 <= 576 && left > 0)) break;
            qi = unsigned(buf >> 58);
            HUFF_QUAD()
            if (!(l + 4 <= 576 && left > 0)) break;
            qi = unsigned(buf >> 58);
            HUFF_QUAD()
            if (!(l + 4 <= 576 && left > 0)) break;
            qi = unsigned(buf >> 58);
            HUFF_QUAD()
        }
#undef HUFF_QUAD
    }
handover:
    st.l = int(l);
    st.region = region;
    st.pos = size_t(limit - left);
    st.nLin = nLin;
#undef HUFF_REFILL
#undef HUFF_TAKE
#undef HUFF_PLAIN_PAIR
#undef HUFF_D_STEP
}
