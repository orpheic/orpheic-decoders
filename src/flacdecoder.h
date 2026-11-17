// Ο δικός μας αποκωδικοποιητής FLAC, γραμμένος από το πρότυπο (RFC 9639).
//
// Μετρήθηκε απέναντι στον αποκωδικοποιητή του FFmpeg στα ίδια δεδομένα και βγάζει τα ίδια δείγματα, ένα προς ένα. Το FFmpeg
// μετριέται μόνο απ' έξω, τρέχοντάς το (χρόνος, έξοδος).
//
// ΤΙ ΚΑΝΕΙ. Ένα καρέ FLAC τη φορά: κεφαλίδα (με τον έλεγχο CRC-8 της), ένα υποκαρέ ανά κανάλι (σταθερό, αυτούσιο, σταθερός
// προβλέπτης τάξης 0-4, γραμμικός προβλέπτης τάξης 1-32), το υπόλοιπο σε κώδικα Rice με τμήματα, αχρησιμοποίητα μπιτ,
// αποσυσχέτιση στέρεο (αριστερό/διαφορά, διαφορά/δεξί, μέσο/διαφορά), και ο έλεγχος CRC-16 όλου του καρέ. Βάθη 4 ως 32
// μπιτ, 1 ως 8 κανάλια, κάθε μέγεθος μπλοκ ως 65535. Τα δείγματα βγαίνουν ως ακέραιοι 32 μπιτ με την αληθινή τους τιμή.
//
// ΑΚΡΙΒΕΙΑ. Όλη η αριθμητική είναι ακέραια. Ο γραμμικός προβλέπτης αθροίζει σε 32 μπιτ μόνο όταν το πρότυπο εγγυάται ότι
// χωράει (βάθος υποκαρέ + ακρίβεια συντελεστών + log2 της τάξης ≤ 32, παράρτημα A.3)· αλλιώς σε 64. Οι σταθεροί
// προβλέπτες δουλεύουν πάντα modulo 2^32: δεν έχουν ολίσθηση, άρα τα κάτω 32 μπιτ βγαίνουν σωστά όπως κι αν ξεχειλίσει το
// ενδιάμεσο, και το δείγμα ενός έγκυρου αρχείου χωράει εκεί. Το κανάλι διαφοράς ήχου 32 μπιτ θέλει 33: αυτά τα καρέ
// περνούν από δρόμο 64 μπιτ.
//
// ΑΣΦΑΛΕΙΑ. Τίποτα από όσα γράφει το αρχείο δεν γίνεται πιστευτό χωρίς έλεγχο: το μέγεθος του μπλοκ και τα κανάλια
// έρχονται από την κεφαλίδα ΤΟΥ ΚΑΡΕ (όχι από το STREAMINFO), κάθε τάξη και κάθε τμήμα ελέγχεται απέναντι στο μπλοκ, και ο
// αναγνώστης μπιτ δεν διαβάζει ποτέ πέρα από 8 bytes μετά το τέλος που του δόθηκε. Χαλασμένα δεδομένα δίνουν κατάσταση
// λάθους, όχι κατάρρευση.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace orpheic::flac {

struct StreamInfo {
    bool present = false;
    uint32_t minBlock = 0, maxBlock = 0;
    uint32_t minFrame = 0, maxFrame = 0;   // bytes· 0 = άγνωστο
    uint32_t sampleRate = 0;
    int channels = 0;
    int bits = 0;
    uint64_t totalSamples = 0;   // ανά κανάλι· 0 = άγνωστο
    uint8_t md5[16] = {};        // των αποκωδικοποιημένων δειγμάτων· όλο μηδέν = άγνωστο
    bool hasMd5() const
    {
        for (uint8_t b : md5)
            if (b) return true;
        return false;
    }
};

// Βρίσκει πού αρχίζει ο ήχος: ετικέτα ID3v2 αν υπάρχει μπροστά, «fLaC», τα μπλοκ μεταδεδομένων. Αρχείο που αρχίζει
// κατευθείαν με καρέ γίνεται δεκτό (info.present = false). false = δεν είναι FLAC ή τα μεταδεδομένα του δεν βγάζουν νόημα.
bool readHead(const uint8_t *data, size_t size, StreamInfo &info, size_t &audioStart);

// Μόνο το σώμα ενός μπλοκ STREAMINFO (τα 34 bytes του 8.2), όπως το δίνει π.χ. ένας αναγνώστης αρχείων ως «extradata».
// false = λιγότερα από 34 bytes ή τιμές που το πρότυπο δεν επιτρέπει (ρυθμός 0, μπλοκ κάτω από 16).
bool parseStreamInfo(const uint8_t *body, size_t size, StreamInfo &info);

enum class Status {
    Ok,
    Truncated,     // τα δεδομένα τελείωσαν πριν από το καρέ
    NotAFrame,     // εδώ δεν αρχίζει καρέ (συγχρονισμός, δεσμευμένες τιμές, CRC-8 της κεφαλίδας)
    Unsupported,   // η κεφαλίδα παραπέμπει στο STREAMINFO και δεν υπάρχει
    BadData,       // το καρέ παραβαίνει το πρότυπο
    BadCrc         // το CRC-16 του καρέ δεν ταιριάζει
};

struct FrameHead {
    uint32_t blockSize = 0;   // δείγματα ανά κανάλι
    uint32_t sampleRate = 0;
    int channels = 0;
    int bits = 0;
    int assignment = 0;       // 0 ανεξάρτητα κανάλια · 8 αριστερό/διαφορά · 9 διαφορά/δεξί · 10 μέσο/διαφορά
    bool variable = false;    // μεταβλητό μέγεθος μπλοκ: το `number` μετρά δείγματα, αλλιώς καρέ
    uint64_t number = 0;
    int headBytes = 0;
};

// Μόνο η κεφαλίδα (για αναζήτηση συγχρονισμού). Χρειάζονται ως 16 bytes.
Status parseFrameHead(const uint8_t *p, size_t n, const StreamInfo &si, FrameHead &h);

// Η πρώτη θέση ≥ from όπου διαβάζεται έγκυρη κεφαλίδα καρέ· n αν δεν υπάρχει.
size_t findFrame(const uint8_t *p, size_t n, size_t from, const StreamInfo &si);

// Ποιες μορφές των πυρήνων διάλεξε αυτός ο επεξεργαστής (για να γράφεται δίπλα στους αριθμούς).
const char *kernelText();

// ΔΟΚΙΜΗ: true = μόνο ό,τι έχει κάθε x86-64 (SSE2), σαν να ήταν παλιός ο επεξεργαστής· false = ό,τι έχει αυτός.
// Έτσι ο δρόμος των παλιών μηχανημάτων αποδεικνύεται και σε καινούργιο. Όχι από δύο νήματα μαζί.
void useBasicOnly(bool on);

// Τι περιέχει το αρχείο (για τις μετρήσεις): μετράει μόνο όταν του δοθεί.
struct Stats {
    uint64_t frames = 0, subframes = 0;
    uint64_t constant = 0, verbatim = 0, fixed[5] = {}, lpc[33] = {};
    uint64_t lpcWide = 0;                 // γραμμικοί προβλέπτες που θέλουν άθροισμα 64 μπιτ
    uint64_t precision[17] = {}, shift[32] = {};
    uint64_t wasted = 0;
    uint64_t partitionOrder[16] = {};
    uint64_t riceSamples[32] = {};        // δείγματα ανά παράμετρο Rice
    uint64_t escapedSamples = 0;
    uint64_t assignment[11] = {};
    uint64_t slowRice = 0;                // σύμβολα Rice που δεν χώρεσαν στον γρήγορο δρόμο
    uint64_t ringPairs = 0, ringOne = 0, ringRedone = 0;   // πρόβλεψη «σε δαχτυλίδι»: ζεύγη καναλιών, μόνα κανάλια, επαναλήψεις λόγω δείγματος 17 μπιτ
    double nsFrame = 0, nsResidual = 0, nsPredict = 0, nsCrc = 0, nsStereo = 0;   // πού πήγε ο χρόνος (μόνο όταν μετράει)
};

class Decoder {
public:
    explicit Decoder(const StreamInfo &si = StreamInfo());

    void setCheckCrc(bool on) { checkCrc_ = on; }
    void setStats(Stats *s) { stats_ = s; }

    // Αποκωδικοποιεί ΕΝΑ καρέ που αρχίζει στο p. Τα p[n .. n+8) πρέπει να ΔΙΑΒΑΖΟΝΤΑΙ (ό,τι κι αν έχουν). `used` = bytes του
    // καρέ. Μετά από Ok: head(), channel(c), write().
    Status decode(const uint8_t *p, size_t n, size_t &used);
    // Πόσα bytes κοίταξε η τελευταία κλήση του decode, ό,τι κι αν απάντησε (σε αποτυχία το `used` είναι 0): για όποιον μετρά τη
    // δουλειά που ξοδεύει ψάχνοντας καρέ μέσα σε χαλασμένα δεδομένα.
    size_t lastSpan() const { return span_; }
    // Και πόσα ΔΕΙΓΜΑΤΑ ζήτησαν τα υποκαρέ που άρχισε να διαβάζει (μέγεθος μπλοκ επί υποκαρέ). Τα bytes μόνα τους δεν μετρούν τη δουλειά:
    // ένα σταθερό υποκαρέ είναι 3 bytes και γεμίζει ως 65.535 δείγματα.
    size_t lastFilled() const { return filled_; }

    const FrameHead &head() const { return head_; }
    const int32_t *channel(int c) const { return cur_[c]; }   // head().blockSize δείγματα, αληθινές τιμές

    enum class Out { S16Interleaved, S16Planar, S32Interleaved, S32Planar };
    // Το τελευταίο καρέ στη μορφή που ζητείται, κάθε δείγμα ολισθημένο αριστερά κατά `shift`. Πλεγμένα: dst[0]. Χωριστά: dst[c].
    void write(Out fmt, int shift, void *const *dst) const;

private:
    StreamInfo si_;
    FrameHead head_;
    bool checkCrc_ = true;
    size_t span_ = 0;
    size_t filled_ = 0;
    Stats *stats_ = nullptr;
    void predict(const void *subs, const FrameHead &h);

    std::vector<int32_t> plane_[8];   // ό,τι διαβάστηκε από το καρέ: δείγματα προθέρμανσης και υπόλοιπο· ο απλός δρόμος προβλέπει επιτόπου
    std::vector<int32_t> alt_[8];     // η έξοδος των γρήγορων δρόμων πρόβλεψης (το υπόλοιπο μένει άθικτο, για επανάληψη)
    int32_t *cur_[8] = {};            // πού είναι τα δείγματα κάθε καναλιού μετά το τελευταίο καρέ
    std::vector<int64_t> wide_[2];
};

}  // namespace orpheic::flac
