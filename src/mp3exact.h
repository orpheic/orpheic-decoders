// Οι πυρήνες του αποκωδικοποιητή MP3 ΚΑΤΑ ΓΡΑΜΜΑ από το ISO/IEC 11172-3, σε διπλή ακρίβεια, αργά: ο μάρτυρας.
//
// Κάθε γρήγορος δρόμος (mp3decoder.cpp και οι πυρήνες του) πρέπει να δίνει το ίδιο με αυτούς: τους ίδιους ακεραίους στους
// κώδικες Huffman, και την ίδια έξοδο μέσα σε ανοχή στρογγύλευσης στα υπόλοιπα. Εδώ δεν υπάρχει καμία βελτιστοποίηση, επίτηδες:
// κάθε συνάρτηση είναι ο τύπος ή το διάγραμμα του προτύπου, με την παραπομπή του.
//
// ΔΙΑΤΑΞΗ. Οι 576 τιμές ενός granule ζουν ως πίνακας 18 σειρών επί 32 υποζώνες: θέση = 32·σειρά + υποζώνη.
//   · πριν από τον μετασχηματισμό: σειρά = γραμμή συχνότητας μέσα στην υποζώνη (0 … 17) στα μακριά μπλοκ· στα κοντά
//     σειρά = 6·παράθυρο + γραμμή (η «αναδιάταξη» του 2.4.3.4.8: xr[υποζώνη][παράθυρο][γραμμή])
//   · μετά: σειρά = χρονική θέση (0 … 17), δηλαδή 32 δείγματα υποζωνών ανά χρονική θέση, όπως τα θέλει η τράπεζα σύνθεσης.
#pragma once

#include "mp3consts.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

namespace orpheic::mp3::exact {

namespace tab {
#include "mp3tables.inc"
}

inline constexpr double kPi = 3.14159265358979323846;

// Αναγνώστης μπιτ, ένα μπιτ τη φορά. Πέρα από το `end` (το τέλος των δεδομένων που υπάρχουν) δίνει μηδενικά και προχωρά τη
// θέση, ώστε όποιος διαβάζει να βλέπει μετά ότι ξεπέρασε.
struct Bits {
    const uint8_t *p = nullptr;
    size_t pos = 0, end = 0;   // σε μπιτ
    uint32_t get1()
    {
        if (pos >= end) {
            ++pos;
            return 0;
        }
        const uint32_t b = (p[pos >> 3] >> (7 - (pos & 7))) & 1u;
        ++pos;
        return b;
    }
    uint32_t get(int n)
    {
        uint32_t v = 0;
        while (n-- > 0) v = (v << 1) | get1();
        return v;
    }
};

// Δέντρο κωδίκων: node[2·i + μπιτ] → > 0 ο επόμενος κόμβος · < 0 φύλλο με τιμή −(v + 1) · 0 κανένας κώδικας.
struct Tree {
    std::vector<int16_t> node;
    void build(const tab::HuffCode *codes, int count, bool quad)
    {
        node.assign(2, 0);
        for (int c = 0; c < count; ++c) {
            const tab::HuffCode &h = codes[c];
            if (h.len == 0) continue;
            int n = 0;
            for (int b = h.len - 1; b >= 0; --b) {
                const size_t at = size_t(2 * n) + ((h.code >> b) & 1u);
                if (b == 0) {
                    node[at] = int16_t(-((quad ? h.x : h.x * 16 + h.y) + 1));
                } else {
                    if (node[at] == 0) {
                        node[at] = int16_t(node.size() / 2);
                        node.push_back(0);
                        node.push_back(0);
                    }
                    n = node[at];
                }
            }
        }
    }
    // Η τιμή του κώδικα που αρχίζει εδώ· −1 αν τα μπιτ δεν είναι κώδικας του πίνακα.
    int read(Bits &b) const
    {
        int n = 0;
        for (int depth = 0; depth < 32; ++depth) {
            const int v = node[size_t(2 * n) + b.get1()];
            if (v < 0) return -v - 1;
            if (v == 0) return -1;
            n = v;
        }
        return -1;
    }
};

struct Tables {
    Tree pair[32], quad[2];
    double synthN[64][32];      // 2.4.3.2.2: N[i][k] = cos((16 + i)(2k + 1)π/64)
    double synthD[512];         // πίνακας B.3
    double cs[8], ca[8];        // πίνακας B.9
    double imdctLong[36][18];   // 2.4.3.4.10.2 με n = 36
    double imdctShort[12][6];   // και με n = 12
    double window[4][36];       // 2.4.3.4.10.3, ανά block_type· το [2] είναι το παράθυρο 12 σημείων του κοντού μπλοκ
    Tables()
    {
        for (int t = 0; t < 32; ++t)
            if (tab::kHuffTable[t].codes) pair[t].build(tab::kHuffTable[t].codes, tab::kHuffTable[t].count, false);
        quad[0].build(tab::kHuffQuadA, 16, true);
        quad[1].build(tab::kHuffQuadB, 16, true);
        for (int i = 0; i < 64; ++i)
            for (int k = 0; k < 32; ++k) synthN[i][k] = consts::cosPi576(9L * (16 + i) * (2 * k + 1));   // cos((16+i)(2k+1)π/64)
        for (int i = 0; i < 512; ++i) synthD[i] = tab::kWindow[i] / 65536.0;
        for (int i = 0; i < 8; ++i) {
            cs[i] = 1.0 / std::sqrt(1.0 + tab::kAliasC[i] * tab::kAliasC[i]);
            ca[i] = tab::kAliasC[i] / std::sqrt(1.0 + tab::kAliasC[i] * tab::kAliasC[i]);
        }
        for (int i = 0; i < 36; ++i)
            for (int k = 0; k < 18; ++k) imdctLong[i][k] = consts::cosPi576(8L * (2 * i + 1 + 18) * (2 * k + 1));   // π/72·…
        for (int i = 0; i < 12; ++i)
            for (int k = 0; k < 6; ++k) imdctShort[i][k] = consts::cosPi576(24L * (2 * i + 1 + 6) * (2 * k + 1));   // π/24·…
        for (int i = 0; i < 36; ++i) {
            // sin(π/36·(i+½)) = sin(π(2i+1)/72) → m = 8(2i+1)· sin(π/12·(j+½)) = sin(π(2j+1)/24) → m = 24(2j+1)
            const double longSine = consts::sinPi576(8 * (2 * i + 1));
            const auto shortSine = [](int j) { return consts::sinPi576(24 * (2 * j + 1)); };
            window[0][i] = longSine;
            window[1][i] = i < 18 ? longSine : i < 24 ? 1.0 : i < 30 ? shortSine(i - 18) : 0.0;
            window[2][i] = i < 12 ? shortSine(i) : 0.0;
            window[3][i] = i < 6 ? 0.0 : i < 12 ? shortSine(i - 6) : i < 18 ? 1.0 : longSine;
        }
    }
};

inline const Tables &tables()
{
    static const Tables t;
    return t;
}

// ---------------------------------------------------------------------------------------------------------------------
// Οι κώδικες Huffman ενός granule ενός καναλιού (2.4.1.7 Huffmancodebits, 2.4.3.4.6).

struct HuffmanIn {
    int big = 0;                // big_values · 2: πόσες γραμμές είναι κωδικοποιημένες ως ζεύγη (≤ 576)
    int r1 = 0, r2 = 0;         // η γραμμή όπου αρχίζει η περιοχή 1 και η περιοχή 2 (ήδη κομμένες στο big)
    int tableSelect[3] = {};    // ο πίνακας κάθε περιοχής, 0 … 31
    int count1Table = 0;        // 0 = πίνακας A των τετράδων, 1 = πίνακας B
};

struct HuffmanNotes {
    uint64_t linbitsValues = 0;    // τιμές ≥ 15 που διάβασαν επιπλέον μπιτ
    bool bigOverrun = false;       // ένα ζεύγος πέρασε το τέλος του granule (και όλα τα επόμενα λείπουν)
    bool quadOverrun = false;      // μια τετράδα πέρασε το τέλος του granule
    int quadOverrunLine = 0;       // … σε αυτή τη γραμμή
    uint64_t stuffingBits = 0;     // μπιτ που περίσσεψαν πριν από το τέλος του granule
};

// Γεμίζει is[0 … 575]. Επιστρέφει πόσες γραμμές διάβασε (οι υπόλοιπες είναι μηδέν), ή −1 αν μια περιοχή διαλέγει πίνακα που
// «δεν χρησιμοποιείται» (4, 14) ή τα μπιτ δεν είναι κώδικας (τότε is = όλο μηδέν). `gend`: το τέλος των μπιτ αυτού του granule
// (≤ b.end). Το b.pos μένει εκεί που σταμάτησε η ανάγνωση (μπορεί να είναι μετά το gend).
inline int huffman(Bits &b, size_t gend, const HuffmanIn &in, int *is, HuffmanNotes *notes = nullptr)
{
    const Tables &t = tables();
    int l = 0;
    bool bad = false;
    for (; l < in.big; l += 2) {
        const int sel = in.tableSelect[l < in.r1 ? 0 : l < in.r2 ? 1 : 2];
        const tab::HuffTable &ht = tab::kHuffTable[sel];
        int x = 0, y = 0;
        if (sel != 0) {   // ο πίνακας 0 δεν έχει κώδικες: όλες οι τιμές της περιοχής είναι μηδέν, χωρίς κανένα μπιτ
            if (!ht.codes) {
                bad = true;
                break;
            }
            const int v = t.pair[sel].read(b);
            if (v < 0) {
                bad = true;
                break;
            }
            x = v >> 4;
            y = v & 15;
            if (x == 15 && ht.linbits) {
                x += int(b.get(ht.linbits));
                if (notes) ++notes->linbitsValues;
            }
            if (x && b.get1()) x = -x;
            if (y == 15 && ht.linbits) {
                y += int(b.get(ht.linbits));
                if (notes) ++notes->linbitsValues;
            }
            if (y && b.get1()) y = -y;
            if (b.pos > gend) {   // το ζεύγος πέρασε το τέλος των μπιτ του granule: δεν υπάρχει, ούτε όσα θα ακολουθούσαν
                if (notes) notes->bigOverrun = true;
                break;
            }
        }
        is[l] = x;
        is[l + 1] = y;
    }
    if (bad) {
        for (int i = 0; i < 576; ++i) is[i] = 0;
        return -1;
    }
    if (l >= in.big) {
        // Οι τετράδες: ώσπου να τελειώσουν τα μπιτ του granule ή οι 576 γραμμές, ό,τι έρθει πρώτο.
        const Tree &qt = t.quad[in.count1Table];
        while (l + 4 <= 576 && b.pos < gend) {
            const int q = qt.read(b);
            if (q < 0) break;
            int v[4] = {(q >> 3) & 1, (q >> 2) & 1, (q >> 1) & 1, q & 1};
            for (int k = 0; k < 4; ++k)
                if (v[k] && b.get1()) v[k] = -1;
            if (b.pos > gend) {
                if (notes) {
                    notes->quadOverrun = true;
                    notes->quadOverrunLine = l;
                }
                break;
            }
            for (int k = 0; k < 4; ++k) is[l + k] = v[k];
            l += 4;
        }
        if (notes && b.pos < gend) notes->stuffingBits += gend - b.pos;
    }
    for (int i = l; i < 576; ++i) is[i] = 0;
    return l;
}

// ---------------------------------------------------------------------------------------------------------------------
// Μείωση αναδίπλωσης, αντίστροφος MDCT, παράθυρο, επικάλυψη, αντιστροφή συχνότητας (2.4.3.4.10.1 ως 2.4.3.4.10.5).

struct Hybrid {
    double overlap[18 * 32] = {};   // το δεύτερο μισό του προηγούμενου μπλοκ: θέση = 32·i + υποζώνη
};

// x: 576 τιμές πριν από τον μετασχηματισμό (βλ. ΔΙΑΤΑΞΗ στην κορυφή)· δεν αλλάζει. out: 576 δείγματα υποζωνών, 32·χρονική θέση + υποζώνη.
// blockType 0 … 3 όπως στο πρότυπο. mixed (μόνο με blockType 2): οι υποζώνες 0 και 1 είναι μακριές με το κανονικό παράθυρο.
inline void hybrid(Hybrid &st, const double *x, int blockType, bool mixed, double *out)
{
    const Tables &t = tables();
    double v[576];
    std::memcpy(v, x, sizeof v);
    // Μείωση αναδίπλωσης στα σύνορα των υποζωνών: σε όλα στα μακριά μπλοκ· σε κανένα στα κοντά· στο μικτό μόνο ανάμεσα στις δύο
    // μακριές υποζώνες. Σειρά 17 − i της υποζώνης sb − 1 με σειρά i της υποζώνης sb.
    const int aliasEnd = blockType != 2 ? 32 : mixed ? 2 : 1;
    for (int sb = 1; sb < aliasEnd; ++sb)
        for (int i = 0; i < 8; ++i) {
            const double lo = v[32 * (17 - i) + sb - 1], hi = v[32 * i + sb];
            v[32 * (17 - i) + sb - 1] = lo * t.cs[i] - hi * t.ca[i];
            v[32 * i + sb] = hi * t.cs[i] + lo * t.ca[i];
        }
    for (int sb = 0; sb < 32; ++sb) {
        double z[36];
        if (blockType != 2 || (mixed && sb < 2)) {
            const double *win = t.window[blockType == 2 ? 0 : blockType];
            for (int i = 0; i < 36; ++i) {
                double s = 0.0;
                for (int k = 0; k < 18; ++k) s += v[32 * k + sb] * t.imdctLong[i][k];
                z[i] = s * win[i];
            }
        } else {
            for (int i = 0; i < 36; ++i) z[i] = 0.0;
            for (int w = 0; w < 3; ++w)
                for (int i = 0; i < 12; ++i) {
                    double s = 0.0;
                    for (int k = 0; k < 6; ++k) s += v[32 * (6 * w + k) + sb] * t.imdctShort[i][k];
                    z[6 + 6 * w + i] += s * t.window[2][i];
                }
        }
        for (int i = 0; i < 18; ++i) {
            double r = z[i] + st.overlap[32 * i + sb];
            st.overlap[32 * i + sb] = z[i + 18];
            if ((sb & 1) && (i & 1)) r = -r;   // 2.4.3.4.10.5: κάθε δεύτερο δείγμα κάθε δεύτερης υποζώνης αλλάζει πρόσημο
            out[32 * i + sb] = r;
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Η τράπεζα σύνθεσης (σχήμα A.2), για τις 18 χρονικές θέσεις ενός granule.

struct Synth {
    double v[1024] = {};
};

// in: 576 δείγματα υποζωνών, 32·χρονική θέση + υποζώνη. out: 576 δείγματα ήχου.
inline void synth(Synth &st, const double *in, double *out)
{
    const Tables &t = tables();
    for (int slot = 0; slot < 18; ++slot) {
        std::memmove(st.v + 64, st.v, sizeof(double) * 960);
        for (int i = 0; i < 64; ++i) {
            double s = 0.0;
            for (int k = 0; k < 32; ++k) s += t.synthN[i][k] * in[32 * slot + k];
            st.v[i] = s;
        }
        double w[512];
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 32; ++j) {
                w[i * 64 + j] = st.v[i * 128 + j] * t.synthD[i * 64 + j];
                w[i * 64 + 32 + j] = st.v[i * 128 + 96 + j] * t.synthD[i * 64 + 32 + j];
            }
        for (int j = 0; j < 32; ++j) {
            double s = 0.0;
            for (int i = 0; i < 16; ++i) s += w[j + 32 * i];
            out[32 * slot + j] = s;
        }
    }
}

}  // namespace orpheic::mp3::exact
