// Δικός μας αποκωδικοποιητής MP3 (MPEG-1 Layer III), γραμμένος από το πρότυπο ISO/IEC 11172-3. Βλ. mp3decoder.h.
//
// Οι παραπομπές στα σχόλια («2.4.3.4.7», «πίνακας B.8», «σχήμα A.2») είναι σε εκείνο το έγγραφο.
#include "mp3decoder.h"

#include "mp3exact.h"
#include "mp3huff.h"
#include "mp3hybrid.h"
#include "mp3synth.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

namespace orpheic::mp3 {
namespace {

namespace tab = exact::tab;
using exact::Bits;
using exact::kPi;

// 2.4.2.3: ρυθμός μπιτ του Layer III ανά bitrate_index (kbit/s)· 0 = ελεύθερη μορφή, το 15 απαγορεύεται.
const int kBitrateL3[15] = {0, 32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256, 320};
// ISO/IEC 13818-3, 2.4.2.3: το ίδιο για τις χαμηλές συχνότητες δειγματοληψίας (ID = 0).
const int kBitrateL3Low[15] = {0, 8, 16, 24, 32, 40, 48, 56, 64, 80, 96, 112, 128, 144, 160};
// 2.4.2.3: συχνότητα δειγματοληψίας ανά sampling_frequency· το 3 είναι δεσμευμένο. Στο 13818-3 οι μισές: 22,05 · 24 · 16 kHz.
const int kSampleRate[3] = {44100, 48000, 32000};
// ISO/IEC 13818-3, 2.4.3.2: πόσοι συντελεστές κλίμακας έχει καθένα από τα τέσσερα τμήματα, ανά περιοχή τιμών του scalefac_compress
// (τρεις για τα συνηθισμένα κανάλια, τρεις για το δεξί κανάλι όταν υπάρχει στέρεο έντασης) και ανά είδος μπλοκ
// (μακρύ · κοντό · μικτό).
const uint8_t kPartition[6][3][4] = {
    {{6, 5, 5, 5}, {9, 9, 9, 9}, {6, 9, 9, 9}},       {{6, 5, 7, 3}, {9, 9, 12, 6}, {6, 9, 12, 6}},
    {{11, 10, 0, 0}, {18, 18, 0, 0}, {15, 18, 0, 0}}, {{7, 7, 7, 0}, {12, 12, 12, 0}, {6, 15, 12, 0}},
    {{6, 6, 6, 3}, {12, 9, 9, 6}, {6, 12, 9, 6}},     {{8, 8, 5, 0}, {15, 12, 9, 0}, {6, 18, 9, 0}}};
// 2.4.2.7: scalefac_compress → (slen1, slen2).
const uint8_t kSlen[16][2] = {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {3, 0}, {1, 1}, {1, 2}, {1, 3},
                              {2, 1}, {2, 2}, {2, 3}, {3, 1}, {3, 2}, {3, 3}, {4, 2}, {4, 3}};

// 2.4.3.1: CRC-16 με γεννήτορα X^16 + X^15 + X^2 + 1, αρχική κατάσταση όλο άσους, το πιο σημαντικό μπιτ πρώτο.
constexpr int kCrcTrust = 8;   // τόσα καρέ με σωστή λέξη ελέγχου πριν πιστέψουμε μια λάθος (βλ. Η ΛΕΞΗ ΕΛΕΓΧΟΥ ΠΙΣΤΕΥΕΤΑΙ…)

uint16_t crc16(uint16_t crc, const uint8_t *p, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        crc ^= uint16_t(p[i] << 8);
        for (int b = 0; b < 8; ++b) crc = uint16_t((crc & 0x8000) ? (crc << 1) ^ 0x8005 : (crc << 1));
    }
    return crc;
}

// Η πλευρική πληροφορία ενός granule ενός καναλιού (2.4.1.7).
struct Granule {
    int part23 = 0, bigValues = 0, globalGain = 0, scalefacCompress = 0;
    bool windowSwitching = false;
    int blockType = 0;
    bool mixed = false;
    int tableSelect[3] = {};
    int subblockGain[3] = {};
    int region0 = 0, region1 = 0;
    bool preflag = false, scaleOne = false;
    int count1Table = 0;
};

// Μια ζώνη συντελεστή κλίμακας με τη σειρά που έρχονται οι γραμμές από τους κώδικες Huffman (2.4.2.7, τέλος): στα μακριά
// μπλοκ κατά συχνότητα· στα κοντά κατά ζώνη, μέσα στη ζώνη κατά παράθυρο, μέσα στο παράθυρο κατά συχνότητα.
struct Band {
    uint16_t start, width;
    uint8_t sfb;
    int8_t window;   // −1 = μακρύ
};

int bandList(const Granule &g, int rate, Band *out)
{
    int n = 0;
    if (g.blockType == 2) {
        int i = 0, first = 0;
        if (g.mixed) {   // οι δύο κάτω υποζώνες (36 γραμμές) είναι μακριές: ζώνες 0 … 7 του μακριού πίνακα (0 … 5 στις χαμηλές
                         // συχνότητες, όπου οι ζώνες είναι πιο φαρδιές)· μετά οι κοντές από την 3
            for (int sfb = 0; tab::kBandLong[rate][sfb] < 36; ++sfb)
                out[n++] = {tab::kBandLong[rate][sfb], uint16_t(tab::kBandLong[rate][sfb + 1] - tab::kBandLong[rate][sfb]), uint8_t(sfb), -1};
            i = 36;
            first = 3;
        }
        for (int sfb = first; sfb < 13; ++sfb) {
            const int w = tab::kBandShort[rate][sfb + 1] - tab::kBandShort[rate][sfb];
            for (int win = 0; win < 3; ++win) {
                out[n++] = {uint16_t(i), uint16_t(w), uint8_t(sfb), int8_t(win)};
                i += w;
            }
        }
    } else {
        for (int sfb = 0; sfb < 22; ++sfb)
            out[n++] = {tab::kBandLong[rate][sfb], uint16_t(tab::kBandLong[rate][sfb + 1] - tab::kBandLong[rate][sfb]), uint8_t(sfb), -1};
    }
    return n;
}

// Αναγνώστης μπιτ του γρήγορου δρόμου για την πλευρική πληροφορία και τους συντελεστές κλίμακας: ως 25 μπιτ με μία ανάγνωση 4 bytes.
// Ίδια συμπεριφορά με τον exact::Bits: πέρα από το `end` δίνει μηδενικά και προχωρά τη θέση. Τα 4 bytes από τη θέση ανάγνωσης πρέπει
// να ΔΙΑΒΑΖΟΝΤΑΙ όσο pos < end (η δεξαμενή έχει 8 bytes περιθώριο· η πλευρική πληροφορία αντιγράφεται σε δικό της χώρο).
struct FastBits {
    const uint8_t *p = nullptr;
    size_t pos = 0, end = 0;   // σε μπιτ
    uint32_t get(int n)
    {
        if (n <= 0) return 0;
        if (pos + size_t(n) > end) {   // στο τέλος: μπιτ προς μπιτ, με μηδενικά πέρα από αυτό
            uint32_t v = 0;
            for (int k = 0; k < n; ++k, ++pos) v = (v << 1) | (pos < end ? (p[pos >> 3] >> (7 - (pos & 7))) & 1u : 0u);
            return v;
        }
        const uint8_t *q = p + (pos >> 3);
        const uint32_t w = uint32_t(q[0]) << 24 | uint32_t(q[1]) << 16 | uint32_t(q[2]) << 8 | uint32_t(q[3]);
        const uint32_t v = (w << (pos & 7)) >> (32 - n);
        pos += size_t(n);
        return v;
    }
    uint32_t get1() { return get(1); }
};

// Χρονόμετρο σταδίων: μετρά μόνο όταν του το ζητήσουν (Stats::timing).
struct Lap {
    bool on;
    std::chrono::steady_clock::time_point t;
    explicit Lap(bool want) : on(want) { if (on) t = std::chrono::steady_clock::now(); }
    void to(double &sum)
    {
        if (!on) return;
        const auto n = std::chrono::steady_clock::now();
        sum += std::chrono::duration<double, std::nano>(n - t).count();
        t = n;
    }
};

// Οι πίνακες του γρήγορου δρόμου.
struct FastTables {
    float pow43[8207];              // |is|^(4/3) για κάθε τιμή που μπορεί να βγει από τους κώδικες: ως 15 + 8191
    float pow2q[448];               // 2^((i − 402)/4): κάθε κέρδος του 2.4.3.4.7.1 είναι δύναμη του 2 με εκθέτη σε τέταρτα
    uint16_t permLong[576];         // πού πάει η γραμμή i (με τη σειρά που ήρθε) στη διάταξη του μετασχηματισμού, μακρύ μπλοκ
    uint16_t permShort[6][2][576];  // το ίδιο για κοντό μπλοκ, ανά πίνακα ζωνών και [1] για μικτό (η «αναδιάταξη» του 2.4.3.4.8)
    uint8_t bandOfLong[6][576];     // σε ποια ζώνη της λίστας (bandList) ανήκει η γραμμή i, μακρύ μπλοκ
    uint8_t bandOfShort[6][2][576]; // το ίδιο για κοντό και μικτό
    FastTables()
    {
        for (int i = 0; i < 8207; ++i) pow43[i] = float(consts::pow43(i));
        for (int i = 0; i < 448; ++i) pow2q[i] = float(consts::pow2q(i - 402));
        for (int i = 0; i < 576; ++i) permLong[i] = uint16_t(32 * (i % 18) + i / 18);
        for (int rate = 0; rate < 6; ++rate) {
            Band bands[40];
            Granule g;
            for (int kind = 0; kind < 3; ++kind) {   // μακρύ · κοντό · μικτό: η ίδια λίστα ζωνών που χρησιμοποιεί ο αποκωδικοποιητής
                g.blockType = kind == 0 ? 0 : 2;
                g.mixed = kind == 2;
                uint8_t *of = kind == 0 ? bandOfLong[rate] : bandOfShort[rate][kind - 1];
                const int n = bandList(g, rate, bands);
                for (int k = 0; k < n; ++k)
                    for (int i = bands[k].start; i < bands[k].start + bands[k].width; ++i) of[i] = uint8_t(k);
            }
        }
        for (int rate = 0; rate < 6; ++rate)
            for (int mixed = 0; mixed < 2; ++mixed) {
                uint16_t *perm = permShort[rate][mixed];
                int i = 0, first = 0;
                if (mixed) {
                    for (; i < 36; ++i) perm[i] = uint16_t(32 * (i % 18) + i / 18);
                    first = 3;
                }
                for (int sfb = first; sfb < 13; ++sfb) {
                    const int start = tab::kBandShort[rate][sfb], w = tab::kBandShort[rate][sfb + 1] - start;
                    for (int win = 0; win < 3; ++win)
                        for (int f = 0; f < w; ++f) {
                            const int line = start + f;
                            perm[i++] = uint16_t(32 * (6 * win + line % 6) + line / 6);
                        }
                }
            }
    }
};

const FastTables &fastTables()
{
    static const FastTables t;
    return t;
}

}  // namespace

Status parseHead(const uint8_t *p, size_t n, FrameHead &h)
{
    if (n < 4) return Status::Truncated;
    if (p[0] != 0xFF || (p[1] & 0xF0) != 0xF0) return Status::NotAFrame;   // syncword: δώδεκα άσοι
    const int id = (p[1] >> 3) & 1;
    const int layer = (p[1] >> 1) & 3;   // '11' Layer I · '10' Layer II · '01' Layer III · '00' δεσμευμένο
    const int bitrateIndex = p[2] >> 4;
    const int rateIndex = (p[2] >> 2) & 3;
    if (layer == 0 || bitrateIndex == 15 || rateIndex == 3) return Status::NotAFrame;
    h = FrameHead();
    h.version = id ? 1 : 2;
    h.layer = 4 - layer;
    h.crc = (p[1] & 1) == 0;   // protection_bit: 0 = έχει προστεθεί λέξη ελέγχου
    h.padding = (p[2] >> 1) & 1;
    h.mode = p[3] >> 6;
    h.channels = h.mode == 3 ? 1 : 2;
    const int ext = (p[3] >> 4) & 3;
    if (h.layer == 3 && h.mode == 1) {   // 2.4.2.3: '01' στέρεο έντασης · '10' μέσο/διαφορά · '11' και τα δύο
        h.intensity = ext & 1;
        h.midSide = ext & 2;
    }
    const bool low = h.version == 2;
    h.sampleRate = low ? kSampleRate[rateIndex] / 2 : kSampleRate[rateIndex];
    if (h.layer != 3) return Status::Ok;   // τα υπόλοιπα πεδία τα γεμίζει μόνο ό,τι αναλαμβάνουμε
    h.bitrate = (low ? kBitrateL3Low : kBitrateL3)[bitrateIndex] * 1000;
    // ISO/IEC 13818-3, 2.4.3.2: στις χαμηλές συχνότητες το καρέ έχει ένα granule (576 δείγματα) και η σταθερά του μήκους είναι 72.
    h.samples = low ? 576 : 1152;
    h.sideBytes = low ? (h.channels == 1 ? 9 : 17) : (h.channels == 1 ? 17 : 32);
    // 2.4.3.1: N = 144 · ρυθμός / συχνότητα, κομμένο προς τα κάτω, συν ένα byte όταν υπάρχει γέμισμα.
    h.frameBytes = h.bitrate ? (low ? 72 : 144) * h.bitrate / h.sampleRate + (h.padding ? 1 : 0) : 0;
    return Status::Ok;
}

struct Decoder::Impl {
    bool exact = true;
    Stats *stats = nullptr;
    FrameHead head;

    std::vector<uint8_t> reservoir;   // τα κύρια δεδομένα των τελευταίων καρέ (ως 511 bytes πίσω) και του τρέχοντος
    size_t reservoirLen = 0;
    // Κομμάτι της δεξαμενής που ΔΕΝ το έχουμε (η ουρά ενός κομμένου καρέ): [holeStart, holeEnd) σε bytes. Ό,τι granule πατά εκεί
    // δεν αποκωδικοποιείται. holeEnd == 0: κανένα.
    size_t holeStart = 0, holeEnd = 0;
    // Για το «το προηγούμενο καρέ πατά στο επόμενο»: πόσα bytes κύριων δεδομένων έχουν μπει ποτέ στη δεξαμενή, και πού (σε μπιτ, στην ίδια
    // αρίθμηση) τελειώνουν τα δεδομένα που δήλωσε το προηγούμενο καρέ· −1 = άγνωστο.
    uint64_t mainTotal = 0;
    long long prevEndAbs = -1;
    int crcGood = 0;                  // καρέ με σωστή λέξη ελέγχου που έχει δει αυτό το ρεύμα (ως kCrcTrust)· ΔΕΝ μηδενίζεται στο reset

    int scfLong[2][22] = {};          // οι συντελεστές κλίμακας του τελευταίου granule κάθε καναλιού (μένουν για το scfsi)
    // Το πρώτο granule αυτού του καρέ, ανά κανάλι: ήταν κοντό μπλοκ; μικτό; (βλ. readScalefactors, «απαγορευμένος συνδυασμός»).
    bool short0[2] = {false, false}, mixed0[2] = {false, false};
    int scfShort[2][13][3] = {};
    int illegalLong[22] = {};          // δεξί κανάλι: η τιμή του συντελεστή που σημαίνει «παράνομη θέση έντασης» (−1: καμία)
    int illegalShort[13][3] = {};
    exact::Hybrid hybridX[2];          // ο ακριβής δρόμος: επικάλυψη (2.4.3.4.10.4)
    exact::Synth synthX[2];            // και το ιστορικό της τράπεζας σύνθεσης (σχήμα A.2)
    HybridState hybridF[2];            // τα ίδια για τον γρήγορο δρόμο
    SynthState synthF[2];
    float out[2][1152] = {};

    void clearState()
    {
        reservoirLen = 0;
        holeStart = holeEnd = 0;
        mainTotal = 0;
        prevEndAbs = -1;
        std::memset(scfLong, 0, sizeof scfLong);
        std::memset(scfShort, 0, sizeof scfShort);
        short0[0] = short0[1] = mixed0[0] = mixed0[1] = false;
        for (int &v : illegalLong) v = -1;
        for (auto &band : illegalShort)
            for (int &v : band) v = -1;
        for (int ch = 0; ch < 2; ++ch) clearChannel(ch);
        lastRate = lastVersion = lastChannels = 0;
    }
    // Η επικάλυψη και το ιστορικό της τράπεζας σύνθεσης ενός καναλιού.
    void clearChannel(int ch)
    {
        hybridX[ch] = exact::Hybrid();
        synthX[ch] = exact::Synth();
        hybridReset(hybridF[ch]);
        synthReset(synthF[ch]);
    }
    int lastRate = 0, lastVersion = 0, lastChannels = 0;

    template <class B> bool readSideInfo(B &b, int &mainDataBegin, uint8_t scfsi[2][4], Granule gr[2][2]) const;
    template <class B> void readScalefactors(B &b, const Granule &g, int ch, int granule, const uint8_t *scfsi);
    template <class B> void readScalefactorsLow(B &b, Granule &g, int ch);
    int readHuffman(Bits &b, size_t gend, const Granule &g, int rate, int *is, long long frameStart);
    long debugFrame = -1;
    void requantize(const Granule &g, int ch, int rate, const int *is, double *xr) const;
    template <class T, class At0, class At1> void stereo(const Granule *g, int rate, T *x0, T *x1, At0 at0, At1 at1);
    exact::HuffmanIn huffmanIn(const Granule &g, int rate, long long frameStart) const;
    void noteHuffman(const exact::HuffmanNotes &notes);
    int readHuffmanFast(const uint8_t *data, size_t endBits, size_t &pos, size_t gend, const Granule &g, int rate, int16_t *is, long long frameStart);
    int requantizeFast(const Granule &g, int ch, int rate, const int16_t *is, int count, float *x) const;
    template <class T> int notSound(T *x, Stats *st) const;
    static void arrange(const Granule &g, int rate, const double *xr, double *x);
};

// 2.4.1.7. false = τιμές που το πρότυπο δεν επιτρέπει.
template <class B>
bool Decoder::Impl::readSideInfo(B &b, int &mainDataBegin, uint8_t scfsi[2][4], Granule gr[2][2]) const
{
    const int nch = head.channels;
    const bool low = head.version == 2;   // ISO/IEC 13818-3, 2.4.1.7: ένα granule, χωρίς scfsi και preflag, scalefac_compress 9 μπιτ
    if (low) {
        mainDataBegin = int(b.get(8));
        b.get(nch == 1 ? 1 : 2);   // private_bits
    } else {
        mainDataBegin = int(b.get(9));
        b.get(nch == 1 ? 5 : 3);   // private_bits
        for (int ch = 0; ch < nch; ++ch)
            for (int band = 0; band < 4; ++band) scfsi[ch][band] = uint8_t(b.get1());
    }
    bool ok = true;
    for (int g = 0; g < (low ? 1 : 2); ++g)
        for (int ch = 0; ch < nch; ++ch) {
            Granule &s = gr[g][ch];
            s = Granule();
            s.part23 = int(b.get(12));
            s.bigValues = int(b.get(9));
            s.globalGain = int(b.get(8));
            s.scalefacCompress = int(b.get(low ? 9 : 4));
            s.windowSwitching = b.get1();
            if (s.windowSwitching) {
                s.blockType = int(b.get(2));
                s.mixed = b.get1();
                for (int r = 0; r < 2; ++r) s.tableSelect[r] = int(b.get(5));
                for (int w = 0; w < 3; ++w) s.subblockGain[w] = int(b.get(3));
                // 2.4.2.7: με αλλαγή παραθύρου οι περιοχές δεν γράφονται: η περιοχή 0 έχει 8 ζώνες (9 όταν όλο το μπλοκ είναι
                // κοντό, όπου κάθε ζώνη μετρά τρεις φορές) και η περιοχή 1 όλα τα υπόλοιπα.
                s.region0 = (s.blockType == 2 && !s.mixed) ? 8 : 7;
                s.region1 = 36;
                if (s.blockType == 0) ok = false;   // «αν η σημαία δεν είναι ψηλά, ο τύπος είναι μηδέν»: ψηλά με τύπο μηδέν δεν ορίζεται
                if (s.blockType != 2) s.mixed = false;
            } else {
                for (int r = 0; r < 3; ++r) s.tableSelect[r] = int(b.get(5));
                s.region0 = int(b.get(4));
                s.region1 = int(b.get(3));
            }
            if (!low) s.preflag = b.get1();
            s.scaleOne = b.get1();
            s.count1Table = int(b.get1());
            if (s.bigValues > 288) ok = false;   // 576 γραμμές το πολύ
        }
    return ok;
}

// ISO/IEC 13818-3, 2.4.3.2 «Scalefactors»: στις χαμηλές συχνότητες οι συντελεστές έρχονται σε τέσσερα τμήματα· το πλήθος και το
// μήκος κάθε τμήματος, και το preflag, βγαίνουν από το scalefac_compress.
template <class B>
void Decoder::Impl::readScalefactorsLow(B &b, Granule &g, int ch)
{
    const int sc = g.scalefacCompress;
    int slen[4] = {0, 0, 0, 0}, row;
    g.preflag = false;
    const bool rightOfIntensity = head.intensity && ch == 1;
    if (!rightOfIntensity) {
        if (sc < 400) {
            slen[0] = (sc >> 4) / 5;
            slen[1] = (sc >> 4) % 5;
            slen[2] = (sc % 16) >> 2;
            slen[3] = sc % 4;
            row = 0;
        } else if (sc < 500) {
            slen[0] = ((sc - 400) >> 2) / 5;
            slen[1] = ((sc - 400) >> 2) % 5;
            slen[2] = (sc - 400) % 4;
            row = 1;
        } else {
            slen[0] = (sc - 500) / 3;
            slen[1] = (sc - 500) % 3;
            g.preflag = true;
            row = 2;
        }
    } else {
        const int isc = sc >> 1;   // το κάτω μπιτ είναι το intensity_scale
        if (isc < 180) {
            slen[0] = isc / 36;
            slen[1] = (isc % 36) / 6;
            slen[2] = (isc % 36) % 6;
            row = 3;
        } else if (isc < 244) {
            slen[0] = ((isc - 180) % 64) >> 4;
            slen[1] = ((isc - 180) % 16) >> 2;
            slen[2] = (isc - 180) % 4;
            row = 4;
        } else {
            slen[0] = (isc - 244) / 3;
            slen[1] = (isc - 244) % 3;
            row = 5;
        }
    }
    const uint8_t *nr = kPartition[row][g.blockType == 2 ? (g.mixed ? 2 : 1) : 0];
    // Με τη σειρά: όλοι οι συντελεστές του τμήματος 1, μετά του 2 … Όπου το μήκος είναι μηδέν ο συντελεστής είναι μηδέν («θέση
    // έντασης μηδέν»)· αλλιώς η μεγαλύτερη τιμή που χωρά στο μήκος σημαίνει, στο δεξί κανάλι του στέρεο έντασης, παράνομη θέση.
    int value[40], limit[40], n = 0;
    for (int part = 0; part < 4; ++part)
        for (int k = 0; k < nr[part]; ++k) {
            value[n] = slen[part] ? int(b.get(slen[part])) : 0;
            limit[n] = slen[part] ? (1 << slen[part]) - 1 : -1;
            ++n;
        }
    int i = 0;
    if (g.blockType == 2) {
        int first = 0;
        if (g.mixed) {
            for (int sfb = 0; sfb < 6; ++sfb, ++i) {
                scfLong[ch][sfb] = value[i];
                if (ch == 1) illegalLong[sfb] = limit[i];
            }
            first = 3;
        }
        for (int sfb = first; sfb < 12; ++sfb)
            for (int w = 0; w < 3; ++w, ++i) {
                scfShort[ch][sfb][w] = value[i];
                if (ch == 1) illegalShort[sfb][w] = limit[i];
            }
        for (int w = 0; w < 3; ++w) scfShort[ch][12][w] = 0;
    } else {
        for (int sfb = 0; sfb < 21; ++sfb, ++i) {
            scfLong[ch][sfb] = value[i];
            if (ch == 1) illegalLong[sfb] = limit[i];
        }
        scfLong[ch][21] = 0;
    }
}

// 2.4.1.7 main_data και 2.4.2.7 (slen1, slen2, scfsi).
template <class B>
void Decoder::Impl::readScalefactors(B &b, const Granule &g, int ch, int granule, const uint8_t *scfsi)
{
    const int slen1 = kSlen[g.scalefacCompress][0], slen2 = kSlen[g.scalefacCompress][1];
    if (g.blockType == 2) {
        int first = 0;
        if (g.mixed) {
            for (int sfb = 0; sfb < 8; ++sfb) scfLong[ch][sfb] = int(b.get(slen1));
            first = 3;
        }
        for (int sfb = first; sfb < 12; ++sfb)
            for (int w = 0; w < 3; ++w) scfShort[ch][sfb][w] = int(b.get(sfb < 6 ? slen1 : slen2));
        for (int w = 0; w < 3; ++w) scfShort[ch][12][w] = 0;
        if (granule == 0) {
            short0[ch] = true;
            mixed0[ch] = g.mixed;
        }
    } else {
        static const int kBound[5] = {0, 6, 11, 16, 21};
        for (int band = 0; band < 4; ++band) {
            if (granule == 1 && scfsi[band]) {   // ισχύουν και εδώ οι συντελεστές του πρώτου granule: δεν μεταδίδονται
                if (stats) ++stats->scfsiUsed;
                if (short0[ch]) {
                    // ΑΠΑΓΟΡΕΥΜΕΝΟΣ ΣΥΝΔΥΑΣΜΟΣ: το πρώτο granule ήταν κοντό μπλοκ. Το πρότυπο (2.4.2.7) ορίζει ότι τότε το scfsi είναι 0·
                    // μερικοί κωδικοποιητές το γράφουν 1 (βρέθηκε σε τρία από 16.000 αληθινά αρχεία). Συντελεστές μακριών
                    // ζωνών για το πρώτο granule ΔΕΝ ΥΠΑΡΧΟΥΝ (ένα κοντό μπλοκ δεν τους μεταδίδει), άρα «οι ίδιοι με του πρώτου» = μηδέν·
                    // σε μικτό μπλοκ οι οκτώ πρώτες μακριές ζώνες υπάρχουν και ισχύουν. Πριν από αυτό έμενε ό,τι είχε αφήσει ένα
                    // παλιότερο, άσχετο granule. Ποιο ακούγεται σωστό: σε 20 τέτοια σημεία το βήμα έντασης απέναντι στα γειτονικά ήταν
                    // τυπικά 0,6 dB με το μηδέν, 0,9 με τους συντελεστές του κοντού μπλοκ στη σειρά που μεταδόθηκαν, 1,5 με τα παλιά.
                    if (stats) ++stats->scfsiAfterShort;
                    for (int sfb = kBound[band]; sfb < kBound[band + 1]; ++sfb) {
                        if (!(mixed0[ch] && sfb < 8)) scfLong[ch][sfb] = 0;
                    }
                }
                continue;
            }
            for (int sfb = kBound[band]; sfb < kBound[band + 1]; ++sfb) scfLong[ch][sfb] = int(b.get(band < 2 ? slen1 : slen2));
        }
        scfLong[ch][21] = 0;
        if (granule == 0) short0[ch] = mixed0[ch] = false;
    }
}

// Οι κώδικες Huffman ενός granule: εδώ μόνο τα όρια των περιοχών (2.4.2.7 region0_count, region1_count)· η ανάγνωση στο mp3exact.h.
exact::HuffmanIn Decoder::Impl::huffmanIn(const Granule &g, int rate, long long frameStart) const
{
    exact::HuffmanIn in;
    in.big = g.bigValues * 2;
    if (g.blockType == 2 && !g.mixed) {   // κάθε κοντή ζώνη μετρά τρεις φορές: 9 «ζώνες» = οι τρεις πρώτες κοντές, επί τρία παράθυρα
        in.r1 = 3 * tab::kBandShort[rate][(g.region0 + 1) / 3];
        in.r2 = 576;
    } else {
        const int a = g.region0 + 1, c = a + g.region1 + 1;
        in.r1 = tab::kBandLong[rate][a < 22 ? a : 22];
        // Μικτό μπλοκ: η περιοχή 0 είναι το μακρύ του μέρος, οι 36 πρώτες γραμμές. Στο ISO/IEC 11172-3 αυτό βγαίνει και από τον πίνακα
        // (οι 8 πρώτες μακριές ζώνες είναι ακριβώς 36 γραμμές). Στις χαμηλές συχνότητες οι ζώνες είναι φαρδύτερες και το πρότυπο δεν
        // λέει πού πέφτει το όριο· παίρνουμε πάλι τις 36 γραμμές. Το ίδιο κάνει ο αποκωδικοποιητής του FFmpeg (μετρημένο
        // απ' έξω, τρέχοντάς τον σε ροές όπου γυρίσαμε εμείς το μπιτ του μικτού μπλοκ: κανένας κωδικοποιητής που έχουμε δεν παράγει τέτοια).
        if (g.mixed) in.r1 = 36;
        in.r2 = tab::kBandLong[rate][c < 22 ? c : 22];
    }
    if (in.r1 > in.big) in.r1 = in.big;
    if (in.r2 > in.big) in.r2 = in.big;
    for (int r = 0; r < 3; ++r) in.tableSelect[r] = g.tableSelect[r];
    in.count1Table = g.count1Table;
    (void)frameStart;
    return in;
}

void Decoder::Impl::noteHuffman(const exact::HuffmanNotes &notes)
{
    if (stats) {
        stats->linbitsValues += notes.linbitsValues;
        stats->stuffingBits += notes.stuffingBits;
        if (notes.bigOverrun) {
            ++stats->bigOverrun;
            stats->flags |= 2;
        }
        if (notes.quadOverrun) {
            ++stats->count1Overrun;
            stats->flags |= 1;
            stats->note[1] = notes.quadOverrunLine;
        }
        if (notes.stuffingBits) stats->flags |= 4;
    }
}

int Decoder::Impl::readHuffman(Bits &b, size_t gend, const Granule &g, int rate, int *is, long long frameStart)
{
    const exact::HuffmanIn in = huffmanIn(g, rate, frameStart);
    exact::HuffmanNotes notes;
    const int n = exact::huffman(b, gend, in, is, stats ? &notes : nullptr);
    noteHuffman(notes);
    return n;
}

int Decoder::Impl::readHuffmanFast(const uint8_t *data, size_t endBits, size_t &pos, size_t gend, const Granule &g, int rate, int16_t *is, long long frameStart)
{
    const exact::HuffmanIn in = huffmanIn(g, rate, frameStart);
    exact::HuffmanNotes notes;
    const int n = huffman(data, endBits, pos, gend, in, is, stats ? &notes : nullptr);
    noteHuffman(notes);
    return n;
}

// ΤΙ ΔΕΝ ΕΙΝΑΙ ΗΧΟΣ. Το πρότυπο επιτρέπει σε μια γραμμή του φάσματος τιμές ως 8206 με όποια ενίσχυση γράφει το καρέ· ένα
// χαλασμένο καρέ (ή ένας κωδικοποιητής με λάθος) δίνει έτσι γραμμές εκατοντάδες ή χιλιάδες φορές πάνω από τη μέγιστη στάθμη. Σε
// αληθινά αρχεία βρέθηκαν και τα δύο: ένα καρέ με κανονική μουσική και ΔΥΟ γραμμές 130 και 36 φορές πάνω από τη μέγιστη στάθμη,
// και καρέ όλο σκουπίδια που έβγαιναν 56.800 φορές πάνω. Όποιος τα παίζει τα κόβει στο όριο: ένα χτύπημα στη μέγιστη στάθμη.
//
// Το μέτρο: μια γραμμή με τιμή 1 δίνει ημίτονο με κορυφή 0,97 της μέγιστης στάθμης (μετρημένο με φτιαχτά καρέ).
// «Ενέργεια» ενός granule = το άθροισμα των τετραγώνων των 576 γραμμών ενός καναλιού. ΜΕΤΡΗΜΕΝΟ 14/11 στα 220 πιο δυνατά από 6.692 αληθινά
// MP3 (μουσική κομμένη ως 3,8 φορές πάνω από τη μέγιστη στάθμη): η μεγαλύτερη γραμμή 3,5, η μεγαλύτερη ενέργεια 26. Στα χαλασμένα
// σημεία έξι αρχείων: ενέργειες από 4.260 ως 2·10^13, γραμμές ως 4.220.000. Τα όρια μπαίνουν ανάμεσα, με περιθώριο από τη μουσική: γραμμή
// 16 (4,6 φορές η μεγαλύτερη αληθινή), ενέργεια 256 (9,8 φορές η μεγαλύτερη αληθινή, δηλαδή 10 dB).
//
// Ο κανόνας μας, ανά granule και κανάλι, μετά την αποκβάντιση: αν η ενέργεια περνά το kEnergyLimit, πρώτα μηδενίζονται οι γραμμές πάνω
// από το kLineLimit (η περίπτωση «μουσική με δύο παράλογες γραμμές»: η μουσική συνεχίζει)· αν και μετά η ενέργεια περνά το όριο, σιγεί
// ολόκληρο (η περίπτωση «σκουπίδια»). Και στις δύο περιπτώσεις το καρέ μετριέται ως χαλασμένο (BadData).
constexpr double kLineLimit = 16.0, kEnergyLimit = 256.0;
// ΤΙΜΕΣ ΠΑΝΩ Ή ΔΙΠΛΑ ΣΤΟ ΟΡΙΟ. Ο ακριβής δρόμος υπολογίζει σε διπλή ακρίβεια και ο γρήγορος σε απλή, άρα η ίδια γραμμή βγαίνει με
// διαφορά στα τελευταία ψηφία. Με σκέτο «μεγαλύτερο από το όριο» οι δύο πήραν ΑΛΛΗ απόφαση για το ίδιο καρέ σε μια δοκιμή με σκόπιμα
// χαλασμένα αρχεία (ο ένας μηδένισε, ο άλλος έπαιξε: διαφορά ίση με όλο το δείγμα). Το 16 και το 256 είναι επιπλέον «φυσικές» τιμές
// (8^(4/3) = 16, και μία γραμμή 16 έχει ενέργεια 256). Γι' αυτό η σύγκριση γίνεται ένα εκατομμυριοστό ΚΑΤΩ από το όριο: τιμή ίση με το
// όριο μετρά ως παράλογη και στους δύο δρόμους, και η στρογγύλευση (δέκα φορές μικρότερη) δεν φτάνει για να τους χωρίσει.
constexpr double kTie = 1.0 - 1e-6;
constexpr long long kOverlapFlawBits = 64;   // βλ. ΤΟ ΠΡΟΗΓΟΥΜΕΝΟ ΚΑΡΕ ΠΑΤΑ ΣΤΟ ΕΠΟΜΕΝΟ

// 0 = ήχος · 1 = μηδενίστηκαν γραμμές · 2 = σίγησε ολόκληρο. Η σειρά των γραμμών δεν παίζει ρόλο (ίδιο για τους δύο δρόμους).
template <class T> int Decoder::Impl::notSound(T *x, Stats *st) const
{
    // Οκτώ ανεξάρτητα αθροίσματα: ένα μόνο θα έκανε αλυσίδα 576 προσθέσεων που περιμένουν η μία την άλλη.
    T acc[8] = {};
    for (int i = 0; i < 576; i += 8)
        for (int k = 0; k < 8; ++k) acc[k] += x[i + k] * x[i + k];
    double e = 0;
    for (int k = 0; k < 8; ++k) e += double(acc[k]);
    if (e < kEnergyLimit * kTie) {
        if (st && e > st->loudestGranule) st->loudestGranule = e;
        return 0;
    }
    int zeroed = 0;
    e = 0;
    for (int i = 0; i < 576; ++i) {
        const double v = double(x[i]);
        if (v >= kLineLimit * kTie || v <= -kLineLimit * kTie || v != v) {   // και το «όχι αριθμός»
            x[i] = T(0);
            ++zeroed;
        } else {
            e += v * v;
        }
    }
    if (e >= kEnergyLimit * kTie) {
        for (int i = 0; i < 576; ++i) x[i] = T(0);
        if (st) ++st->absurdGranules;
        return 2;
    }
    if (st) st->absurdLines += uint64_t(zeroed);
    return 1;
}

// 2.4.3.4.7.1 για τον γρήγορο δρόμο: το ίδιο κέρδος από πίνακα (ο εκθέτης του 2 είναι πάντα ακέραιος αριθμός τετάρτων), το |is|^(4/3)
// από πίνακα, και η τιμή γράφεται κατευθείαν στη θέση της στη διάταξη του μετασχηματισμού. Επιστρέφει πόσες υποζώνες από την αρχή
// μπορεί να έχουν μη μηδενική τιμή (όλες οι άλλες είναι σίγουρα μηδέν).
//
// ΕΝΑΣ ΙΣΙΟΣ ΒΡΟΧΟΣ πάνω σε όλες τις γραμμές, όχι ένας μικρός ανά ζώνη: πρώτα το κέρδος κάθε ζώνης σε μικρό πίνακα, και μετά κάθε γραμμή
// βρίσκει τη ζώνη της από έτοιμο πίνακα. Με έναν βρόχο ανά ζώνη (4 ως 100 γραμμές ο καθένας) ο μεταγλωττιστής στο -O3 έφτιαχνε διανυσματικό
// κώδικα με αρχή και ουρά για κάθε ζώνη, και το στάδιο έβγαινε τέσσερις φορές πιο αργό απ' ό,τι στο -O2 (μετρημένο 13/11).
int Decoder::Impl::requantizeFast(const Granule &g, int ch, int rate, const int16_t *is, int count, float *x) const
{
    const FastTables &ft = fastTables();
    std::memset(x, 0, 576 * sizeof(float));
    if (count <= 0) return 0;
    Band bands[40];
    const int n = bandList(g, rate, bands);
    const bool isShort = g.blockType == 2;
    const uint16_t *perm = isShort ? ft.permShort[rate][g.mixed ? 1 : 0] : ft.permLong;
    const uint8_t *bandOf = isShort ? ft.bandOfShort[rate][g.mixed ? 1 : 0] : ft.bandOfLong[rate];
    const int step = g.scaleOne ? 4 : 2;   // 4 · scalefac_multiplier
    float gain[40];
    int subbands = 0;
    for (int k = 0; k < n; ++k) {
        const Band &bd = bands[k];
        int e = g.globalGain - 210 + 402;
        int top;   // η πρώτη υποζώνη που ΔΕΝ αγγίζει αυτή η ζώνη
        if (bd.window < 0) {
            e -= step * (scfLong[ch][bd.sfb] + (g.preflag ? tab::kPretab[bd.sfb] : 0));
            top = (bd.start + bd.width + 17) / 18;
        } else {
            e -= 8 * g.subblockGain[bd.window] + step * scfShort[ch][bd.sfb][bd.window];
            top = (tab::kBandShort[rate][bd.sfb + 1] + 5) / 6;
        }
        gain[k] = ft.pow2q[e < 0 ? 0 : e > 447 ? 447 : e];
        if (bd.start < count && top > subbands) subbands = top;
    }
    // Χωρίς διακλάδωση για το πρόσημο: τα πρόσημα των τιμών είναι τυχαία, και μια διακλάδωση εδώ μαντευόταν λάθος τις μισές φορές
    // (μετρημένο: το στάδιο αυτό ήταν το μισό του χρόνου όλου του αποκωδικοποιητή). Το μέτρο είναι ποτέ αρνητικό, άρα το μπιτ
    // προσήμου της ακέραιας τιμής μπαίνει κατευθείαν πάνω του.
    for (int i = 0; i < count; ++i) {
        const int32_t v = is[i];
        const int32_t neg = v >> 31;                     // 0 ή −1
        const float m = ft.pow43[(v ^ neg) - neg] * gain[bandOf[i]];
        // Τα μπιτ περνούν από τοπικές μεταβλητές και η τιμή γράφεται ως float: ένα memcpy κατευθείαν μέσα στον πίνακα γινόταν, με
        // τους ελέγχους ορίων του Ubuntu (_FORTIFY_SOURCE), κλήση συνάρτησης ανά τιμή μόλις η συνάρτηση ενσωματωνόταν στον καλούντα.
        uint32_t bits;
        std::memcpy(&bits, &m, sizeof bits);
        bits |= uint32_t(neg) & 0x80000000u;
        float signedValue;
        std::memcpy(&signedValue, &bits, sizeof signedValue);
        x[perm[i]] = signedValue;
    }
    if (!isShort) subbands = (count + 17) / 18;   // στα μακριά μπλοκ οι γραμμές είναι στη σειρά: φτάνει το πλήθος τους
    return subbands > 32 ? 32 : subbands;
}

// 2.4.3.4.7.1, κατά γράμμα. Οι τιμές μένουν στη σειρά που ήρθαν.
void Decoder::Impl::requantize(const Granule &g, int ch, int rate, const int *is, double *xr) const
{
    Band bands[40];
    const int n = bandList(g, rate, bands);
    // Σε τέταρτα της δύναμης του 2 (2.4.3.4.7.1): 2^(¼(global_gain − 210 − 8·subblock_gain)) · 2^(−scalefac_multiplier·…), με
    // scalefac_multiplier ½ ή 1 = 2 ή 4 τέταρτα. Ένα ακριβές 2^(q/4) (consts::pow2q) αντί για γινόμενο δύο std::pow: ίδιο σε κάθε σύστημα.
    const int quarters = g.scaleOne ? 4 : 2;
    for (int k = 0; k < n; ++k) {
        const Band &bd = bands[k];
        double gain;
        if (bd.window < 0)
            gain = consts::pow2q(g.globalGain - 210 - quarters * (scfLong[ch][bd.sfb] + (g.preflag ? tab::kPretab[bd.sfb] : 0)));
        else
            gain = consts::pow2q(g.globalGain - 210 - 8 * g.subblockGain[bd.window] - quarters * scfShort[ch][bd.sfb][bd.window]);
        for (int i = bd.start; i < bd.start + bd.width; ++i) {
            const int v = is[i];
            const double m = consts::pow43(v < 0 ? -v : v) * gain;   // ο πίνακας του προτύπου ως 8206 (15 + 2^13 − 1)
            xr[i] = v < 0 ? -m : m;
        }
    }
}

// 2.4.3.4.9: στέρεο μέσου/διαφοράς και στέρεο έντασης. Οι γραμμές μετριούνται με τη σειρά που ήρθαν από τους κώδικες Huffman: η
// γραμμή i του αριστερού ζευγαρώνει με τη γραμμή i του δεξιού. Το at0(i) και το at1(i) λένε πού βρίσκεται η γραμμή i μέσα στον
// πίνακα κάθε καναλιού. Έτσι ο ίδιος κώδικας, με τις ίδιες αποφάσεις, υπηρετεί τον ακριβή δρόμο (double, η θέση είναι το ίδιο το i)
// και τον γρήγορο (float, η θέση στη διάταξη του μετασχηματισμού, που εξαρτάται από το είδος μπλοκ ΤΟΥ ΚΑΘΕ καναλιού: σε χαλασμένα
// δεδομένα τα δύο κανάλια μπορεί να έχουν διαφορετικό).
template <class T, class At0, class At1>
void Decoder::Impl::stereo(const Granule *g, int rate, T *x0, T *x1, At0 at0, At1 at1)
{
    bool intensityLine[576] = {};
    const bool low = head.version == 2;
    if (head.intensity) {
        // Το «μηδενικό μέρος» του δεξιού καναλιού ορίζει από πού και πάνω ισχύει το στέρεο έντασης: η ζώνη με την τελευταία
        // μη μηδενική γραμμή του δεξιού είναι η τελευταία που ΔΕΝ είναι έντασης. Στα κοντά μπλοκ, χωριστά για κάθε παράθυρο.
        const Granule &gr = g[1];
        Band bands[40];
        const int n = bandList(gr, rate, bands);
        int lastBand[4] = {-1, -1, -1, -1};   // ανά παράθυρο (0 … 2) και [3] για τις μακριές ζώνες: η τελευταία ζώνη με μη μηδενική γραμμή
        // Η θέση έντασης κάθε ζώνης είναι ο συντελεστής κλίμακας του δεξιού καναλιού· η τελευταία ζώνη (χωρίς δικό της συντελεστή)
        // παίρνει της προηγούμενης. «Παράνομη» θέση: 7 στο ISO/IEC 11172-3 (2.4.3.4.9.3, βήμα 2)· στο ISO/IEC 13818-3, 2.4.3.2, η
        // μεγαλύτερη τιμή που χωρά στο μήκος του συντελεστή.
        int posOf[40];
        bool illegalOf[40];
        for (int k = 0; k < n; ++k) {
            const Band &bd = bands[k];
            int pos, illegal;
            if (bd.window < 0) {
                pos = scfLong[1][bd.sfb < 21 ? bd.sfb : 20];
                illegal = illegalLong[bd.sfb < 21 ? bd.sfb : 20];
            } else {
                pos = scfShort[1][bd.sfb < 12 ? bd.sfb : 11][bd.window];
                illegal = illegalShort[bd.sfb < 12 ? bd.sfb : 11][bd.window];
            }
            posOf[k] = pos;
            illegalOf[k] = low ? pos == illegal : pos >= 7;
        }
        for (int k = 0; k < n; ++k) {
            const Band &bd = bands[k];
            bool any = false;
            for (int i = bd.start; i < bd.start + bd.width; ++i)
                if (x1[at1(i)] != T(0)) {
                    any = true;
                    break;
                }
            if (any) lastBand[bd.window < 0 ? 3 : bd.window] = k;
        }
        if (debugFrame >= 0 && stats && long(stats->frames) - 1 == debugFrame) {
            std::fprintf(stderr, "frame %ld: right channel block type %d mixed %d scalefac_compress %d (intensity_scale %d) ms %d · last bands with sound: long %d, windows %d %d %d\n",
                         debugFrame, gr.blockType, int(gr.mixed), gr.scalefacCompress, gr.scalefacCompress & 1, int(head.midSide), lastBand[3], lastBand[0], lastBand[1], lastBand[2]);
            for (int k = 0; k < n; ++k) {
                double el = 0, er = 0;
                for (int i = bands[k].start; i < bands[k].start + bands[k].width; ++i) {
                    el += double(x0[at0(i)]) * double(x0[at0(i)]);
                    er += double(x1[at1(i)]) * double(x1[at1(i)]);
                }
                std::fprintf(stderr, "  band %2d (sfb %2d win %2d, lines %3d..%3d): pos %2d%s · energy L %.3g R %.3g\n", k, bands[k].sfb, bands[k].window, bands[k].start,
                             bands[k].start + bands[k].width - 1, posOf[k], illegalOf[k] ? " ILLEGAL" : "", el, er);
            }
        }
        const bool shortHasSound = lastBand[0] >= 0 || lastBand[1] >= 0 || lastBand[2] >= 0;
        bool used = false;
        for (int k = 0; k < n; ++k) {
            const Band &bd = bands[k];
            bool is;
            if (bd.window < 0)
                is = k > lastBand[3] && !(gr.blockType == 2 && shortHasSound);   // μικτό μπλοκ: οι μακριές ζώνες είναι κάτω από τις κοντές
            else
                is = k > lastBand[bd.window];
            if (!is) continue;
            // Η θέση είναι ο συντελεστής κλίμακας του δεξιού καναλιού· η τελευταία ζώνη (χωρίς συντελεστή) παίρνει της προηγούμενης.
            const int pos = posOf[k];
            if (stats && low) {
                const int lim = bd.window < 0 ? illegalLong[bd.sfb < 21 ? bd.sfb : 20] : illegalShort[bd.sfb < 12 ? bd.sfb : 11][bd.window];
                int bits = 0;
                while ((1 << bits) - 1 < lim) ++bits;
                if (lim >= 0 && bits < 6 && pos < 32 && !((bd.window < 0 && bd.sfb == 21) || (bd.window >= 0 && bd.sfb == 12))) ++stats->positionSeen[bits][pos];
            }
            if (illegalOf[k]) {   // η ζώνη δεν είναι στέρεο έντασης: μέσο/διαφορά αν ισχύει, αλλιώς τα κανάλια όπως είναι
                if (stats) stats->flags |= 512;
                continue;
            }
            used = true;
            if (stats) {
                ++stats->intensityBands;
                stats->flags |= 32;
                if (bd.window >= 0) stats->flags |= 64;                                   // στέρεο έντασης σε κοντό μπλοκ
                if ((bd.window < 0 && bd.sfb == 21) || (bd.window >= 0 && bd.sfb == 12)) stats->flags |= 128;   // στην τελευταία ζώνη, χωρίς δικό της συντελεστή
                if (pos == 0) stats->flags |= 256;
                if (pos >= 16) stats->flags |= 1024;
                else if (pos >= 8) stats->flags |= 2048;
            }
            if (low) {
                // ISO/IEC 13818-3, 2.4.3.2: R = L·kr, L = L·kl, με βάση i0 = 1/√2 όταν intensity_scale = 1, αλλιώς 1/⁴√2.
                // i0^n = 2^(−n/2) ή 2^(−n/4): ακριβές 2^(q/4) (consts::pow2q), όχι δύναμη ενός στρογγυλεμένου 1/√2.
                const int step = (gr.scalefacCompress & 1) ? 2 : 1;
                double kl = 1.0, kr = 1.0;
                if (pos > 0 && (pos & 1)) kl = consts::pow2q(-step * ((pos + 1) / 2));
                else if (pos > 0) kr = consts::pow2q(-step * (pos / 2));
                for (int i = bd.start; i < bd.start + bd.width; ++i) {
                    const T l = x0[at0(i)];
                    x0[at0(i)] = T(l * kl);
                    x1[at1(i)] = T(l * kr);
                    intensityLine[i] = true;
                }
                continue;
            }
            for (int i = bd.start; i < bd.start + bd.width; ++i) {
                const T l = x0[at0(i)];
                if (pos == 6) {   // is_ratio = tan(π/2): όλο στο αριστερό
                    x1[at1(i)] = T(0);
                } else {
                    const double ratio = consts::tanPi576(48 * pos);   // 2.4.3.4.9.3: tan(pos·π/12)
                    x0[at0(i)] = T(l * ratio / (1.0 + ratio));
                    x1[at1(i)] = T(l * 1.0 / (1.0 + ratio));
                }
                intensityLine[i] = true;
            }
        }
        if (stats && used) ++stats->intensityGranules;
    }
    if (head.midSide) {
        if (stats) ++stats->msGranules;
        const double r = 1.0 / std::sqrt(2.0);
        for (int i = 0; i < 576; ++i) {
            if (intensityLine[i]) continue;
            const T m = x0[at0(i)], s = x1[at1(i)];
            x0[at0(i)] = T((m + s) * r);
            x1[at1(i)] = T((m - s) * r);
        }
    }
    // Για τη διάγνωση: το τελικό φάσμα των δύο καναλιών αυτού του καρέ, γραμμή προς γραμμή (μαζί με τον ήχο που έβγαλε ένας άλλος
    // αποκωδικοποιητής λύνεται πίσω ποιες γραμμές έβγαλε αλλιώς).
    if (debugFrame >= 0 && stats && long(stats->frames) - 1 == debugFrame)
        for (int i = 0; i < 576; ++i) std::fprintf(stderr, "LINE %d %.9g %.9g\n", i, double(x0[at0(i)]), double(x1[at1(i)]));
}

// Από τη σειρά που ήρθαν οι τιμές στη διάταξη του μετασχηματισμού (θέση = 32·σειρά + υποζώνη, βλ. mp3exact.h). Στα κοντά μπλοκ
// αυτή είναι η «αναδιάταξη» του 2.4.3.4.8.
void Decoder::Impl::arrange(const Granule &g, int rate, const double *xr, double *x)
{
    if (g.blockType == 2) {
        int i = 0, first = 0;
        for (int k = 0; k < 576; ++k) x[k] = 0.0;
        if (g.mixed) {
            for (; i < 36; ++i) x[32 * (i % 18) + i / 18] = xr[i];
            first = 3;
        }
        for (int sfb = first; sfb < 13; ++sfb) {
            const int start = tab::kBandShort[rate][sfb], w = tab::kBandShort[rate][sfb + 1] - start;
            for (int win = 0; win < 3; ++win)
                for (int f = 0; f < w; ++f) {
                    const int line = start + f;   // γραμμή μέσα στο παράθυρο, 0 … 191
                    x[32 * (6 * win + line % 6) + line / 6] = xr[i++];
                }
        }
    } else {
        for (int i = 0; i < 576; ++i) x[32 * (i % 18) + i / 18] = xr[i];
    }
}

Decoder::Decoder() : d_(new Impl) { d_->clearState(); }
Decoder::~Decoder() = default;
void Decoder::setExact(bool on) { d_->exact = on; }
void Decoder::setDebugFrame(long frame) { d_->debugFrame = frame; }
void Decoder::setStats(Stats *s) { d_->stats = s; }
void Decoder::reset() { d_->clearState(); }
void Decoder::lostData()
{
    d_->reservoirLen = 0;
    d_->holeStart = d_->holeEnd = 0;
    d_->prevEndAbs = -1;
}
const FrameHead &Decoder::head() const { return d_->head; }
const float *Decoder::channel(int c) const { return d_->out[c]; }

Status Decoder::decode(const uint8_t *p, size_t n)
{
    Impl &d = *d_;
    double nsNone = 0;
    Lap lap(d.stats && d.stats->timing);
    FrameHead h;
    const Status hs = parseHead(p, n, h);
    if (hs != Status::Ok || h.layer != 3) {
        d.reservoirLen = 0;   // εδώ χάθηκε ένα καρέ άγνωστου μήκους: τα επόμενα δεν μπορούν να δείξουν πίσω από αυτό
        d.holeStart = d.holeEnd = 0;
        d.prevEndAbs = -1;
        return hs != Status::Ok ? hs : Status::Unsupported;
    }
    const bool low = h.version == 2;
    // Ελεύθερη μορφή: το καρέ είναι όσο μας έδωσαν. Καρέ κομμένο πριν από το μήκος που γράφει η κεφαλίδα του (τέλος κομμένου
    // αρχείου): αποκωδικοποιείται ό,τι υπάρχει, αρκεί να έχει ολόκληρη την πλευρική πληροφορία.
    const size_t frameBytes = h.frameBytes && size_t(h.frameBytes) < n ? size_t(h.frameBytes) : n;
    const size_t headBytes = 4 + (h.crc ? 2 : 0);
    // Καρέ κομμένο μέσα στην πλευρική του πληροφορία: δεν έχει τίποτα να αποκωδικοποιηθεί, αλλά ο χρόνος του υπάρχει: βγαίνει σιωπή
    // (με την ουρά του προηγούμενου καρέ, όπως σε κάθε σιωπηλό granule) και όλα τα κύρια δεδομένα του μετρούν ως τρύπα.
    const bool noSide = frameBytes < headBytes + size_t(h.sideBytes);
    // ΑΛΛΑΓΗ ΜΟΡΦΗΣ ΣΤΗ ΜΕΣΗ ΤΗΣ ΡΟΗΣ: άλλη συχνότητα, άλλο πρότυπο ή άλλος αριθμός καναλιών = άλλη κωδικοποίηση
    // κολλημένη στην προηγούμενη (δύο αρχεία ενωμένα). Η επικάλυψη και το ιστορικό της σύνθεσης είναι ήχος της προηγούμενης: καθαρίζουν,
    // και η συνέχεια βγαίνει ΑΚΡΙΒΩΣ όπως αν το νέο κομμάτι αποκωδικοποιούνταν μόνο του (μετρημένο: διαφορά 0).
    if (d.lastRate
        && (d.lastRate != h.sampleRate || d.lastVersion != h.version || d.lastChannels != h.channels)) {
        d.clearChannel(0);
        d.clearChannel(1);
        // Και η δεξαμενή: τα bytes της προηγούμενης μορφής δεν είναι «προηγούμενα δεδομένα» αυτής.
        d.reservoirLen = 0;
        d.holeStart = d.holeEnd = 0;
        d.prevEndAbs = -1;
    }
    d.lastRate = h.sampleRate;
    d.lastVersion = h.version;
    d.lastChannels = h.channels;
    d.head = h;
    Stats *st = d.stats;
    if (st) {
        ++st->frames;
        if (!h.bitrate) ++st->freeFormat;
        st->flags = 0;
        st->note[0] = st->note[1] = st->note[2] = st->note[3] = 0;
    }
    const int nch = h.channels;
    const int rateField = (p[2] >> 2) & 3;
    const int rate = low ? 3 + rateField : rateField;   // η σειρά των πινάκων ζωνών

    // 2.4.3.1: όταν η λέξη ελέγχου δεν ταιριάζει, η πλευρική πληροφορία είναι χαλασμένη και το πρότυπο συνιστά σίγαση του καρέ.
    bool crcBad = false;
    if (h.crc && !noSide) {   // πίνακας B.5: προστατεύονται τα μπιτ 16 … 31 της κεφαλίδας και όλη η πλευρική πληροφορία
        uint16_t crc = crc16(0xFFFF, p + 2, 2);
        crc = crc16(crc, p + 6, size_t(h.sideBytes));
        crcBad = crc != uint16_t(p[4] << 8 | p[5]);
        if (st) {
            ++st->crcFrames;
            if (crcBad) ++st->crcBad;
        }
        // Η ΛΕΞΗ ΕΛΕΓΧΟΥ ΠΙΣΤΕΥΕΤΑΙ ΜΟΝΟ ΑΦΟΥ ΤΟ ΡΕΥΜΑ ΕΔΕΙΞΕ ΟΤΙ ΤΗ ΓΡΑΦΕΙ ΣΩΣΤΑ. Ένας κωδικοποιητής που γράφει
        // λάθος λέξεις ελέγχου θα έβγαζε με τον κανόνα «λάθος λέξη = σιωπή» ολόκληρο το τραγούδι βουβό, ενώ ο ήχος του
        // αποκωδικοποιείται κανονικά. Γι' αυτό μια λάθος λέξη σιγεί το καρέ μόνο όταν έχουν προηγηθεί kCrcTrust καρέ με σωστή σε
        // αυτό το ρεύμα (μετριούνται σε όλη τη ζωή του αποκωδικοποιητή, όχι από το τελευταίο άλμα: είναι ιδιότητα του κωδικοποιητή).
        // Μετρημένο σε δείγμα 16.000 αληθινών αρχείων: τα 1.071 έχουν λέξη ελέγχου· σε 1.068 που αποκωδικοποιήθηκαν (12,2 εκατομμύρια
        // καρέ) βρέθηκαν δύο λάθος λέξεις, μία σε καθένα από δύο αρχεία, και κανένα αρχείο με συστηματικά λάθος λέξεις.
        if (!crcBad) {
            if (d.crcGood < kCrcTrust) ++d.crcGood;
        } else if (d.crcGood < kCrcTrust) {
            crcBad = false;   // δεν ξέρουμε ακόμη αν αυτό το ρεύμα γράφει σωστές λέξεις: το καρέ αποκωδικοποιείται όπως είναι
        }
    }

    int mainDataBegin = 0;
    uint8_t scfsi[2][4] = {};
    Granule gr[2][2];
    bool sideOk = false;
    if (!noSide) {
        if (d.exact) {
            Bits sb;
            sb.p = p + headBytes;
            sb.end = size_t(h.sideBytes) * 8;
            sideOk = d.readSideInfo(sb, mainDataBegin, scfsi, gr);
        } else {
            uint8_t side[40] = {};   // ως 32 bytes, με περιθώριο για τον γρήγορο αναγνώστη
            std::memcpy(side, p + headBytes, size_t(h.sideBytes));
            FastBits sb;
            sb.p = side;
            sb.end = size_t(h.sideBytes) * 8;
            sideOk = d.readSideInfo(sb, mainDataBegin, scfsi, gr);
        }
        if (crcBad) sideOk = false;
    }

    // Η δεξαμενή (2.4.2.7 main_data_begin, 2.4.3.4.3): τα κύρια δεδομένα αυτού του καρέ αρχίζουν main_data_begin bytes ΠΡΙΝ από
    // τα δικά του, μετρώντας μόνο κύρια δεδομένα (όχι κεφαλίδες και πλευρική πληροφορία των ενδιάμεσων καρέ).
    const size_t mainBytes = noSide ? 0 : frameBytes - headBytes - size_t(h.sideBytes);
    size_t lostBytes = 0;   // η ουρά που κόπηκε
    if (h.frameBytes && size_t(h.frameBytes) > n) lostBytes = size_t(h.frameBytes) - headBytes - size_t(h.sideBytes) - mainBytes;
    if (noSide && !h.frameBytes) {   // ελεύθερη μορφή: ούτε το μήκος του χαμένου ξέρουμε
        d.reservoirLen = 0;
        d.holeStart = d.holeEnd = 0;
        d.prevEndAbs = -1;
    }
    if (d.reservoirLen > 511) {
        const size_t drop = d.reservoirLen - 511;
        std::memmove(d.reservoir.data(), d.reservoir.data() + drop, 511);
        d.reservoirLen = 511;
        d.holeStart = d.holeStart > drop ? d.holeStart - drop : 0;
        d.holeEnd = d.holeEnd > drop ? d.holeEnd - drop : 0;
    }
    if (d.reservoir.size() < d.reservoirLen + mainBytes + lostBytes + 8) d.reservoir.resize(d.reservoirLen + mainBytes + lostBytes + 8);
    // Όταν λείπουν προηγούμενα καρέ (αρχή ροής, μετά από μετακίνηση) η αρχή των κύριων δεδομένων πέφτει πριν από ό,τι έχουμε:
    // `mainStart` αρνητικό. Ένα granule που αρχίζει μέσα σε αυτό που λείπει δεν αποκωδικοποιείται: σιωπή και στα δύο κανάλια
    // του (το ένα χωρίς το άλλο δεν είναι ο ήχος του καρέ, ιδίως στο στέρεο μέσου/διαφοράς)· το επόμενο granule, κανονικά.
    const bool haveEarlier = size_t(mainDataBegin) <= d.reservoirLen;
    const long long mainStart = (long long)(d.reservoirLen) * 8 - (long long)(mainDataBegin) * 8;
    if (mainBytes) std::memcpy(d.reservoir.data() + d.reservoirLen, p + headBytes + size_t(h.sideBytes), mainBytes);
    d.reservoirLen += mainBytes;
    if (lostBytes) {   // τα bytes που λείπουν μετρούν στη δεξαμενή (τα επόμενα καρέ δείχνουν πίσω με βάση το πλήρες μήκος), ως τρύπα
        std::memset(d.reservoir.data() + d.reservoirLen, 0, lostBytes);
        if (d.holeEnd == 0) d.holeStart = d.reservoirLen;
        d.reservoirLen += lostBytes;
        d.holeEnd = d.reservoirLen;
    }
    if (st && uint64_t(mainDataBegin) > st->reservoirMax) st->reservoirMax = uint64_t(mainDataBegin);

    // ΤΟ ΠΡΟΗΓΟΥΜΕΝΟ ΚΑΡΕ ΠΑΤΑ ΣΤΟ ΕΠΟΜΕΝΟ. Τα κύρια δεδομένα ενός καρέ πρέπει να τελειώνουν πριν αρχίσουν του επόμενου. Σε
    // χαλασμένα καρέ (και σε μερικά υγιή, από κωδικοποιητή που έγραψε 8 ή 19 μπιτ παραπάνω) το προηγούμενο δηλώνει περισσότερα απ' όσα
    // χωρούν. Εμείς ακολουθούμε το πρότυπο: κάθε καρέ διαβάζεται από εκεί που λέει το δικό του main_data_begin. Τα δεδομένα
    // του επόμενου καρέ διαβάζονται κανονικά και τελειώνουν ακριβώς εκεί που πρέπει, και η μουσική του παίζεται.
    //
    // ΚΑΙ ΕΙΝΑΙ ΖΗΜΙΑ ΟΤΑΝ ΕΙΝΑΙ ΜΕΓΑΛΟ. Μετρημένο 14/11 σε 329 αρχεία (τα 220 πιο δυνατά του δείγματος και 115 χαλασμένα ή περίεργα): πάτημα
    // υπάρχει σε 13, όλα από τα χαλασμένα. Τα δύο μικρά (8 και 19 μπιτ) είναι σε καρέ που κατά τα άλλα είναι υγιή· όλα τα άλλα είναι 64
    // μπιτ και πάνω (ως 2.014) και δίπλα τους υπάρχει πάντα καρέ με πλευρική πληροφορία-σκουπίδια (ενίσχυση 5, 22, 49) ή με λάθος
    // main_data_begin. Δύο καρέ που διεκδικούν τα ίδια 8 bytes δεν τα έγραψε έτσι κανένας κωδικοποιητής: το καρέ μετριέται ως χαλασμένο
    // (BadData), και παίζεται ό,τι αποκωδικοποιείται, γιατί δεν ξέρουμε ποιο από τα δύο λέει ψέματα.
    d.mainTotal += uint64_t(mainBytes + lostBytes);
    bool overlapFlaw = false, absurdFlaw = false, frameGarbage = false;
    {
        const long long startAbs = (long long)(d.mainTotal - d.reservoirLen) * 8 + mainStart;
        const long long overlap = d.prevEndAbs >= 0 && haveEarlier ? d.prevEndAbs - startAbs : 0;
        if (overlap > 0) {
            if (st) ++(overlap >= 8 ? st->overlapBig : st->overlapSmall);
            if (st && overlap >= 64) ++st->overlapHuge;
            overlapFlaw = overlap >= kOverlapFlawBits;
        }
        long long declared = 0;
        for (int g = 0; g < (low ? 1 : 2); ++g)
            for (int ch = 0; ch < nch; ++ch) declared += gr[g][ch].part23;
        d.prevEndAbs = sideOk ? startAbs + declared : -1;
    }

    Status result = Status::Ok;
    if (!sideOk) result = Status::BadData;
    if (st) {
        if (!haveEarlier) ++st->noEarlierData;
        else if (!sideOk) ++st->badFrames;
    }

    lap.to(st ? st->nsFront : nsNone);
    Bits mb;
    mb.p = d.reservoir.data();
    long long granuleStart = mainStart;   // σε μπιτ μέσα στη δεξαμενή· αρνητικό = μέσα σε ό,τι λείπει
    const size_t mainEnd = d.reservoirLen * 8;
    for (int g = 0; g < (low ? 1 : 2); ++g) {
        double xr[2][576];               // ο ακριβής δρόμος: οι τιμές με τη σειρά που ήρθαν
        alignas(64) float fx[2][576];    // ο γρήγορος: κατευθείαν στη διάταξη του μετασχηματισμού
        int16_t is16[2][592];
        int subbands[2] = {0, 0};
        int count[2] = {0, 0};
        bool quiet[2] = {false, false};
        // Κάποιο κανάλι του granule ΑΡΧΙΖΕΙ σε δεδομένα που δεν έχουμε (πριν από την αρχή, ή μέσα στην ουρά κομμένου καρέ): ούτε οι
        // συντελεστές κλίμακάς του δεν υπάρχουν, άρα όλο το granule είναι σιωπή. Κανάλι που αρχίζει σε δεδομένα που έχουμε και
        // συνεχίζει μέσα στην τρύπα αποκωδικοποιείται ως εκεί (οι χαμηλές συχνότητες έρχονται πρώτες) και τα υπόλοιπα είναι μηδέν.
        const long long holeFrom = (long long)(d.holeStart) * 8, holeTo = (long long)(d.holeEnd) * 8;
        bool missing = false;
        {
            long long at = granuleStart;
            for (int ch = 0; ch < nch; ++ch) {
                if (at < 0 || (d.holeEnd > 0 && at >= holeFrom && at < holeTo)) missing = true;
                at += gr[g][ch].part23;
            }
        }
        // ΤΟ ΔΕΥΤΕΡΟ GRANULE ΠΟΥ ΚΛΗΡΟΝΟΜΕΙ ΑΠΟ ΠΡΩΤΟ ΠΟΥ ΔΕΝ ΔΙΑΒΑΣΤΗΚΕ ΑΠΟΚΩΔΙΚΟΠΟΙΕΙΤΑΙ ΟΠΩΣ ΕΙΝΑΙ. Με το scfsi το δεύτερο granule δεν
        // μεταδίδει κάποιους συντελεστές κλίμακας (ισχύουν του πρώτου)· αν το πρώτο δεν αποκωδικοποιήθηκε, μένουν ό,τι υπήρχε (μηδέν μετά
        // από reset) και οι ζώνες αυτές βγαίνουν λίγο πιο δυνατά. Η εναλλακτική είναι η σιωπή· μετρήθηκε σε 398 αρχές από
        // τη μέση απέναντι στον ήχο της συνεχόμενης αποκωδικοποίησης: η σιωπή ήταν πιο μακριά από την αλήθεια σε 84 και πιο κοντά σε
        // καμία, και κανένα δείγμα δεν πέρασε τη στάθμη της μουσικής. Μένει έτσι.
        for (int ch = 0; ch < nch; ++ch) {
            Granule &s = gr[g][ch];
            int is[580];
            bool silent = result == Status::BadData || frameGarbage;
            if (missing) {
                silent = true;
                if (result == Status::Ok) result = Status::NoEarlierData;
            }
            bool empty = false;
            const long long thisStart = granuleStart;
            granuleStart += s.part23;   // το επόμενο granule αρχίζει εκεί που λέει το part2_3_length, ό,τι κι αν διαβαστεί
            if (!silent) {
                const size_t start = size_t(thisStart);
                mb.pos = start;
                size_t end = start + size_t(s.part23);
                if (end > mainEnd) {   // το granule λέει ότι έχει περισσότερα μπιτ από όσα υπάρχουν
                    end = mainEnd;
                    if (result != Status::BadData && st) ++st->badFrames;
                    result = Status::BadData;
                }
                if (d.holeEnd > 0 && (long long)(start) < holeFrom && (long long)(end) > holeFrom) {   // η ουρά του έχει χαθεί
                    end = size_t(holeFrom);
                    if (result == Status::Ok) result = Status::NoEarlierData;
                }
                mb.end = mainEnd;
                if (d.exact) {
                    if (low) d.readScalefactorsLow(mb, s, ch);
                    else d.readScalefactors(mb, s, ch, g, scfsi[ch]);
                } else {
                    FastBits fb;
                    fb.p = d.reservoir.data();
                    fb.pos = mb.pos;
                    fb.end = mainEnd;
                    if (low) d.readScalefactorsLow(fb, s, ch);
                    else d.readScalefactors(fb, s, ch, g, scfsi[ch]);
                    mb.pos = fb.pos;
                }
                lap.to(st ? st->nsScale : nsNone);
                if (mb.pos > end && s.part23 == 0) {
                    // ΑΔΕΙΟ GRANULE ΠΟΥ «ΖΗΤΑ» ΣΥΝΤΕΛΕΣΤΕΣ ΚΛΙΜΑΚΑΣ. Μηδέν μπιτ δεδομένων, αλλά scalefac_compress που θέλει μπιτ για
                    // συντελεστές: κατά γράμμα παραβαίνει το πρότυπο, όμως έτσι γράφει τη σιωπή ένας παλιός κωδικοποιητής (μετρημένο σε
                    // αληθινό αρχείο: πέντε τέτοια καρέ στην αρχή του τραγουδιού, που αλλιώς έβγαινε «χαλασμένο στο 0:00»). Δεν
                    // υπάρχει τίποτα να χαλάσει: το κανάλι είναι σιωπή σε αυτό το granule, και το καρέ ΔΕΝ μετριέται ως χαλασμένο.
                    // Οι «συντελεστές» που μόλις διαβάστηκαν είναι μπιτ του επόμενου granule: μηδενίζονται, για να μην τους κληρονομήσει το
                    // δεύτερο granule μέσω του scfsi.
                    empty = true;
                    // Στο στέρεο έντασης του MPEG-1 οι συντελεστές του δεξιού καναλιού είναι οι «θέσεις έντασης»: μηδέν θα έστελνε όλη τη
                    // μουσική στη μία πλευρά. Εκεί μπαίνει η θέση 7, που σημαίνει «αυτή η ζώνη δεν είναι στέρεο έντασης».
                    {
                        const int fill = (!low && ch == 1 && h.mode == 1 && h.intensity) ? 7 : 0;
                        for (int &v : d.scfLong[ch]) v = fill;
                        for (auto &band : d.scfShort[ch])
                            for (int &v : band) v = fill;
                    }
                } else if (mb.pos > end) {
                    silent = true;
                } else if (d.exact) {
                    count[ch] = d.readHuffman(mb, end, s, rate, is, mainStart);
                    if (count[ch] < 0) silent = true;
                } else {
                    size_t pos = mb.pos;
                    count[ch] = d.readHuffmanFast(d.reservoir.data(), mainEnd, pos, end, s, rate, is16[ch], mainStart);
                    if (count[ch] < 0) silent = true;
                }
                if (silent && result != Status::BadData) {
                    result = Status::BadData;
                    if (st) ++st->badFrames;
                }
            }
            if (silent || empty) count[ch] = 0;
            quiet[ch] = silent;
            lap.to(st ? st->nsHuffman : nsNone);
            if (!d.exact) {
                if (silent || empty) std::memset(is16[ch], 0, sizeof is16[ch]);
                subbands[ch] = d.requantizeFast(s, ch, rate, is16[ch], count[ch], fx[ch]);
                lap.to(st ? st->nsRequant : nsNone);
            } else if (silent || empty) {
                for (int i = 0; i < 576; ++i) xr[ch][i] = 0.0;
            } else {
                d.requantize(s, ch, rate, is, xr[ch]);
            }
            if (!silent) {
                // Το καρέ μετριέται ως χαλασμένο ΣΤΟ ΤΕΛΟΣ (absurdFlaw), όχι εδώ: ένα BadData εδώ θα σιγούσε και όλα τα επόμενα κανάλια και
                // granule του καρέ (το `silent` πιο πάνω), ενώ ο κανόνας λέει «η μουσική συνεχίζει».
                const int ns = d.exact ? d.notSound(xr[ch], st) : d.notSound(fx[ch], st);
                if (ns) absurdFlaw = true;
                // Granule που σίγησε ολόκληρο (σκουπίδια): το υπόλοιπο καρέ δεν είναι πιστευτό και σιγεί κι αυτό. Granule με λίγες παράλογες
                // γραμμές (μουσική): το υπόλοιπο καρέ παίζει.
                if (ns == 2) frameGarbage = true;
                if (ns == 2) {
                    count[ch] = 0;
                    subbands[ch] = 0;
                    quiet[ch] = true;
                }
            }
            if (st) {
                ++st->granules;
                ++st->blockType[s.blockType];
                if (s.blockType == 2) st->flags |= 16;
                if (s.mixed) ++st->mixed;
                if (s.preflag) ++st->preflag;
                if (s.scaleOne) ++st->scaleOne;
                ++st->count1Table[s.count1Table];
                for (int r = 0; r < (s.windowSwitching ? 2 : 3); ++r) ++st->table[s.tableSelect[r]];
            }
        }
        // Κοινό στέρεο με το ένα κανάλι σιωπηλό (χαλασμένο καρέ): το άλλο μόνο του δεν είναι ο ήχος του καρέ, και το στέρεο θα δούλευε
        // πάνω σε συντελεστές που δεν διαβάστηκαν. Και τα δύο σιωπηλά.
        if (nch == 2 && h.mode == 1 && (h.intensity || h.midSide) && (quiet[0] != quiet[1])) {
            for (int ch = 0; ch < 2; ++ch) {
                count[ch] = 0;
                subbands[ch] = 0;
                if (d.exact) for (int i = 0; i < 576; ++i) xr[ch][i] = 0.0;
                else std::memset(fx[ch], 0, sizeof fx[ch]);
            }
        }
        if (!d.exact) {
            if (nch == 2 && h.mode == 1 && (h.intensity || h.midSide)) {
                const int top = subbands[0] > subbands[1] ? subbands[0] : subbands[1];
                const FastTables &ft = fastTables();
                const Granule &left = gr[g][0], &right = gr[g][1];
                const uint16_t *perm0 = left.blockType == 2 ? ft.permShort[rate][left.mixed ? 1 : 0] : ft.permLong;
                const uint16_t *perm1 = right.blockType == 2 ? ft.permShort[rate][right.mixed ? 1 : 0] : ft.permLong;
                if (h.intensity || perm0 != perm1) {
                    // Σπάνιος δρόμος (στέρεο έντασης, ή κανάλια με διαφορετικό είδος μπλοκ): ο ίδιος κώδικας με τον ακριβή, γραμμή προς
                    // γραμμή, με τη θέση κάθε γραμμής στη διάταξη του δικού της καναλιού.
                    d.stereo(gr[g], rate, fx[0], fx[1], [perm0](int i) { return int(perm0[i]); }, [perm1](int i) { return int(perm1[i]); });
                } else {
                    // Μέσο/διαφορά σε όλο το φάσμα (2.4.3.4.9.2): μόνο στις υποζώνες που μπορεί να έχουν τιμές.
                    if (st) ++st->msGranules;
                    // Ολόκληρες σειρές, και στις υποζώνες που είναι σίγουρα μηδέν: ένας βρόχος σταθερού μήκους πάνω σε δύο χωριστούς
                    // πίνακες γίνεται διανυσματικός από τον μεταγλωττιστή, και το μηδέν μένει μηδέν.
                    const float r = 0.70710678118654752440f;
                    float *a = fx[0], *b = fx[1];
                    for (int i = 0; i < 576; ++i) {
                        const float m = a[i], sd = b[i];
                        a[i] = (m + sd) * r;
                        b[i] = (m - sd) * r;
                    }
                }
                // Μετά το στέρεο κάθε κανάλι μπορεί να έχει τιμές όπου είχε το άλλο. Όταν τα δύο κανάλια έχουν διαφορετικό είδος μπλοκ
                // (μόνο σε χαλασμένα δεδομένα) οι ίδιες γραμμές πέφτουν σε άλλες υποζώνες στο καθένα: τότε δεν υπόσχομαι τίποτα.
                subbands[0] = subbands[1] = perm0 != perm1 ? 32 : top;
            }
            lap.to(st ? st->nsStereo : nsNone);
            for (int ch = 0; ch < nch; ++ch) {
                alignas(64) float sub[576];
                hybrid(d.hybridF[ch], fx[ch], gr[g][ch].blockType, gr[g][ch].mixed, subbands[ch], sub);
                lap.to(st ? st->nsHybrid : nsNone);
                synth(d.synthF[ch], sub, d.out[ch] + 576 * g);
                lap.to(st ? st->nsSynth : nsNone);
            }
            continue;
        }
        if (nch == 2 && h.mode == 1 && (h.intensity || h.midSide))
            d.stereo(gr[g], rate, xr[0], xr[1], [](int i) { return i; }, [](int i) { return i; });
        for (int ch = 0; ch < nch; ++ch) {
            double x[576], sub[576], pcm[576];
            Impl::arrange(gr[g][ch], rate, xr[ch], x);
            exact::hybrid(d.hybridX[ch], x, gr[g][ch].blockType, gr[g][ch].mixed, sub);
            exact::synth(d.synthX[ch], sub, pcm);
            float *o = d.out[ch] + 576 * g;
            for (int i = 0; i < 576; ++i) o[i] = float(pcm[i]);
        }
    }
    if ((overlapFlaw || absurdFlaw) && result == Status::Ok) {
        result = Status::BadData;
        if (st) ++st->badFrames;
    }
    return result;
}

}  // namespace orpheic::mp3
