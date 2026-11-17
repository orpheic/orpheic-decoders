// dsthash: αποσυμπιέζει πλαίσια DST με τον αποσυμπιεστή DST του Orpheic και τυπώνει συνόψεις των bytes DSD.
//
//   dsthash αρχείο.dff | πλαίσια.dstf ...
//   dsthash --dsd-out έξοδος.dff αρχείο.dff      (γράφει και τον αποσυμπιεσμένο ήχο, ως DSDIFF χωρίς συμπίεση)
//
// Είσοδος: ένα αρχείο DSDIFF (.dff) με συμπίεση DST, η δημόσια μορφή αρχείου του DSD, ή η απλή απόθεση πλαισίων πιο κάτω.
// Το --dsd-out υπάρχει για τη σύγκριση με άλλον αποσυμπιεστή: το FFmpeg ανοίγει και τα δύο .dff, και αν τα bytes DSD είναι
// ίδια, η μετατροπή του σε PCM δίνει ίδια δείγματα (README, «Check it yourself»).
//
// Παράδειγμα χρήσης του dstdecoder.h, και έλεγχος: το ίδιο αρχείο πρέπει να δίνει την ίδια γραμμή σε κάθε μηχάνημα και
// μεταγλωττιστή. Χωρίς καμία εξάρτηση.
// Αυτό το αποθετήριο δεν έχει αναγνώστη δίσκων SACD: ο αποσυμπιεστής παίρνει ένα πλαίσιο τη φορά, και το πού βρίσκονται τα πλαίσια
// μέσα σε έναν δίσκο είναι δουλειά του καλούντος. Γι' αυτό η είσοδος εδώ είναι μια απλή απόθεση πλαισίων:
//   bytes 0-5    τα πέντε γράμματα "DSTF1" και αλλαγή γραμμής (0x0A)
//   byte  6      κανάλια (1 ως 6)
//   bytes 7-10   ρυθμός δειγματοληψίας του DSD σε Hz, μικρό byte πρώτο (2822400 για το DSD64 των SACD)
//   bytes 11-14  πλήθος πλαισίων, μικρό byte πρώτο
//   και για κάθε πλαίσιο: 4 bytes μήκος (μικρό byte πρώτο) και τα bytes του πλαισίου DST όπως είναι στον δίσκο.
// Ένα πλαίσιο είναι 1/75 του δευτερολέπτου: ρυθμός / 8 / 75 bytes DSD ανά κανάλι (4704 στο DSD64), πλεγμένα ανά κανάλι.
// Όπου ο αποσυμπιεστής αρνείται ένα πλαίσιο, στη θέση του μπαίνει σιωπή DSD (bytes 0x69), ώστε η έξοδος να κρατά τον χρόνο της.
// Τυπώνει SHA-256 όλων των bytes DSD, και FNV-1a 64 μπιτ των ίδιων bytes (μια μικρή σύνοψη, βολική για σύγκριση με το μάτι).
#include "dstdecoder.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <new>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Sha256 {
    uint32_t h[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    uint8_t block[64];
    size_t fill = 0;
    uint64_t total = 0;

    static uint32_t rotr(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }
    void compress(const uint8_t *p)
    {
        static const uint32_t k[64] = {
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be,
            0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa,
            0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85,
            0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
            0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f,
            0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
        uint32_t w[64];
        for (int i = 0; i < 16; ++i) w[i] = uint32_t(p[4 * i]) << 24 | uint32_t(p[4 * i + 1]) << 16 | uint32_t(p[4 * i + 2]) << 8 | p[4 * i + 3];
        for (int i = 16; i < 64; ++i) {
            const uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
        for (int i = 0; i < 64; ++i) {
            const uint32_t t1 = hh + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g)) + k[i] + w[i];
            const uint32_t t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
            hh = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
    }
    void add(const void *data, size_t n)
    {
        const uint8_t *p = static_cast<const uint8_t *>(data);
        total += n;
        while (n) {
            const size_t take = n < 64 - fill ? n : 64 - fill;
            std::memcpy(block + fill, p, take);
            fill += take; p += take; n -= take;
            if (fill == 64) { compress(block); fill = 0; }
        }
    }
    std::string hex()
    {
        const uint64_t bits = total * 8;
        const uint8_t one = 0x80, zero = 0;
        add(&one, 1);
        while (fill != 56) add(&zero, 1);
        uint8_t len[8];
        for (int i = 0; i < 8; ++i) len[i] = uint8_t(bits >> (56 - 8 * i));
        add(len, 8);
        char out[65];
        for (int i = 0; i < 8; ++i) std::snprintf(out + 8 * i, 9, "%08x", h[i]);
        return out;
    }
};

// Το αρχείο διαβάζεται κομμάτι-κομμάτι, ποτέ ολόκληρο: ένα DST .dff ενός ολόκληρου δίσκου είναι μερικά GB, και ο ήχος του
// αποσυμπιεσμένος ακόμη περισσότερα. Θέσεις 64 bit σε κάθε σύστημα (στα Windows το long του ftell είναι 32 bit).
struct File {
    std::FILE *f = nullptr;
    uint64_t size = 0;
    explicit File(const char *path)
    {
        f = std::fopen(path, "rb");
        if (f && seek(0, SEEK_END)) size = tell();
    }
    ~File() { if (f) std::fclose(f); }
    File(const File &) = delete;
    File &operator=(const File &) = delete;
    bool seek(uint64_t at, int whence = SEEK_SET)
    {
#if defined(_WIN32)
        return _fseeki64(f, int64_t(at), whence) == 0;
#else
        return fseeko(f, off_t(at), whence) == 0;
#endif
    }
    uint64_t tell()
    {
#if defined(_WIN32)
        const int64_t t = _ftelli64(f);
#else
        const int64_t t = int64_t(ftello(f));
#endif
        return t < 0 ? 0 : uint64_t(t);
    }
    // n bytes από τη θέση at, μόνο αν υπάρχουν όλα.
    bool read(uint64_t at, void *out, size_t n)
    {
        if (at > size || n > size - at) return false;
        return seek(at) && std::fread(out, 1, n, f) == n;
    }
};

uint32_t get32(const uint8_t *p) { return uint32_t(p[0]) | uint32_t(p[1]) << 8 | uint32_t(p[2]) << 16 | uint32_t(p[3]) << 24; }
uint64_t be(const uint8_t *p, int n)
{
    uint64_t v = 0;
    for (int i = 0; i < n; ++i) v = v << 8 | p[i];
    return v;
}

// Ό,τι χρειάζεται ο αποσυμπιεστής: κανάλια, ρυθμός, και πού είναι κάθε πλαίσιο μέσα στο αρχείο.
struct Input {
    int channels = 0;
    uint32_t rate = 0;
    std::vector<std::pair<uint64_t, uint32_t>> frames;   // θέση, μήκος
    bool cut = false;                                     // το αρχείο τελειώνει πριν από όσα δηλώνει
};

// Η απόθεση πλαισίων (η μορφή στην κορυφή του αρχείου).
const char *fromDump(File &file, Input &in)
{
    uint8_t h[15];
    if (!file.read(0, h, 15)) return "NOT-DFF-OR-FRAME-DUMP";
    in.channels = h[6];
    in.rate = get32(h + 7);
    const uint32_t count = get32(h + 11);
    uint64_t p = 15;
    for (uint32_t i = 0; i < count; ++i) {
        uint8_t l[4];
        if (!file.read(p, l, 4)) { in.cut = true; break; }
        const uint32_t len = get32(l);
        p += 4;
        if (len > file.size - p || len > uint32_t(INT32_MAX)) { in.cut = true; break; }
        in.frames.push_back({p, len});
        p += len;
    }
    return nullptr;
}

// DSDIFF: «FRM8», μέγεθος 8 bytes (μεγάλο byte πρώτο), «DSD », και μετά κομμάτια: όνομα 4 bytes, μέγεθος 8 bytes, περιεχόμενο,
// και ένα byte γέμισμα όταν το μέγεθος είναι μονό (το μέγεθος δεν το μετρά). Από όλα μας χρειάζονται τρία:
//   PROP («SND »), με μέσα του FS (ρυθμός, 4 bytes), CHNL (κανάλια, 2 bytes) και CMPR (η συμπίεση: «DST »)·
//   DST, με μέσα του FRTE (πλήθος πλαισίων) και ένα DSTF ανά πλαίσιο, με τα bytes του πλαισίου όπως είναι.
// Κάθε μέγεθος ελέγχεται απέναντι σε ό,τι μένει πριν χρησιμοποιηθεί· ένα κομμένο αρχείο δίνει όσα πλαίσια είναι ολόκληρα.
const char *fromDff(File &file, Input &in)
{
    uint8_t h[16];
    if (!file.read(0, h, 16) || std::memcmp(h + 12, "DSD ", 4) != 0) return "NOT-A-DSDIFF-FILE";
    uint64_t end = file.size;
    const uint64_t declared = be(h + 4, 8);
    if (declared < end - 12) end = declared + 12;
    else if (declared > end - 12) in.cut = true;
    if (end < 16) return "NOT-A-DSDIFF-FILE";   // ένα FRM8 που δηλώνει λιγότερα από το δικό του «DSD »
    bool dst = false;
    uint32_t declaredFrames = 0;
    uint8_t c[12], v[4];
    for (uint64_t p = 16; p <= end && end - p >= 12;) {
        if (!file.read(p, c, 12)) { in.cut = true; break; }
        const uint64_t size = be(c + 4, 8);
        const uint64_t body = p + 12;
        // Ένα κομμάτι μεγαλύτερο από όσα μένουν: το αρχείο κόπηκε. Διαβάζεται ως το τέλος του, ώστε ένα κομμένο DST να δίνει
        // όσα πλαίσια είναι ολόκληρα, και μετά σταματάμε.
        const bool last = size > end - body;
        if (last) in.cut = true;
        const uint64_t stop = last ? end : body + size;
        if (std::memcmp(c, "PROP", 4) == 0 && size >= 4 && file.read(body, v, 4) && std::memcmp(v, "SND ", 4) == 0) {
            for (uint64_t q = body + 4; q <= stop && stop - q >= 12;) {
                uint8_t k[12];
                if (!file.read(q, k, 12)) break;
                const uint64_t ks = be(k + 4, 8);
                const uint64_t kb = q + 12;
                if (ks > stop - kb) break;
                if (std::memcmp(k, "FS  ", 4) == 0 && ks >= 4 && file.read(kb, v, 4)) in.rate = uint32_t(be(v, 4));
                else if (std::memcmp(k, "CHNL", 4) == 0 && ks >= 2 && file.read(kb, v, 2)) in.channels = int(be(v, 2));
                else if (std::memcmp(k, "CMPR", 4) == 0 && ks >= 4 && file.read(kb, v, 4)) dst = std::memcmp(v, "DST ", 4) == 0;
                q = kb + ks + (ks & 1);
            }
        } else if (std::memcmp(c, "DST ", 4) == 0) {
            for (uint64_t q = body; q <= stop && stop - q >= 12;) {
                uint8_t k[12];
                if (!file.read(q, k, 12)) { in.cut = true; break; }
                const uint64_t ks = be(k + 4, 8);
                const uint64_t kb = q + 12;
                if (ks > stop - kb) { in.cut = true; break; }
                if (std::memcmp(k, "FRTE", 4) == 0 && ks >= 4 && file.read(kb, v, 4)) declaredFrames = uint32_t(be(v, 4));
                else if (std::memcmp(k, "DSTF", 4) == 0) {
                    if (ks > uint64_t(INT32_MAX)) { in.cut = true; break; }   // ο αποσυμπιεστής παίρνει μήκος int
                    in.frames.push_back({kb, uint32_t(ks)});
                }
                q = kb + ks + (ks & 1);
            }
        }
        if (last || stop + (size & 1) > end) break;
        p = stop + (size & 1);
    }
    if (!dst) return "NOT-DST-COMPRESSED";
    if (declaredFrames && in.frames.size() < declaredFrames) in.cut = true;
    return nullptr;
}

// DSDIFF χωρίς συμπίεση, με τα bytes DSD όπως βγαίνουν (πλεγμένα ανά κανάλι, MSB πρώτο, όπως το θέλει η μορφή). Η κεφαλίδα
// γράφεται πρώτη, με το μέγεθος του ήχου γνωστό από πριν (πλαίσια x bytes ανά πλαίσιο), και ο ήχος ακολουθεί πλαίσιο-πλαίσιο.
void put(std::vector<uint8_t> &o, const char *id, uint64_t size)
{
    o.insert(o.end(), id, id + 4);
    for (int i = 7; i >= 0; --i) o.push_back(uint8_t(size >> (8 * i)));
}
std::vector<uint8_t> dffHeader(int channels, uint32_t rate, uint64_t dsdBytes)
{
    static const char *const kIds[7][6] = {{}, {"C   "}, {"SLFT", "SRGT"}, {"MLFT", "MRGT", "C   "}, {"MLFT", "MRGT", "LS  ", "RS  "},
                                           {"MLFT", "MRGT", "C   ", "LS  ", "RS  "}, {"MLFT", "MRGT", "C   ", "LFE ", "LS  ", "RS  "}};
    static const char kNone[] = "not compressed";
    std::vector<uint8_t> prop;
    prop.insert(prop.end(), {'S', 'N', 'D', ' '});
    put(prop, "FS  ", 4);
    for (int i = 3; i >= 0; --i) prop.push_back(uint8_t(rate >> (8 * i)));
    put(prop, "CHNL", 2 + 4 * uint64_t(channels));
    prop.push_back(0);
    prop.push_back(uint8_t(channels));
    for (int c = 0; c < channels; ++c) prop.insert(prop.end(), kIds[channels][c], kIds[channels][c] + 4);
    const uint64_t cmpr = 4 + 1 + (sizeof kNone - 1);
    put(prop, "CMPR", cmpr);
    prop.insert(prop.end(), {'D', 'S', 'D', ' '});
    prop.push_back(uint8_t(sizeof kNone - 1));
    prop.insert(prop.end(), kNone, kNone + sizeof kNone - 1);
    if (cmpr & 1) prop.push_back(0);
    std::vector<uint8_t> o;
    const uint64_t body = 4 + (12 + 4) + (12 + prop.size()) + (12 + dsdBytes + (dsdBytes & 1));
    put(o, "FRM8", body);
    o.insert(o.end(), {'D', 'S', 'D', ' '});
    put(o, "FVER", 4);
    o.insert(o.end(), {1, 5, 0, 0});
    put(o, "PROP", prop.size());
    o.insert(o.end(), prop.begin(), prop.end());
    put(o, "DSD ", dsdBytes);
    return o;
}

}  // namespace

int usage()
{
    std::fprintf(stderr, "usage: dsthash file.dff|frames.dstf ...\n"
                         "       dsthash --dsd-out decoded.dff file.dff\n");
    return 2;
}

int main(int argc, char **argv)
{
    int failed = 0;
    const char *dsdOut = nullptr;
    int first = 1;
    if (argc < 2) return usage();
    if (std::strcmp(argv[1], "--dsd-out") == 0) {
        if (argc != 4) return usage();
        dsdOut = argv[2];
        first = 3;
        // Ποτέ πάνω στην ίδια την είσοδο, με όποιο όνομα κι αν γραφτεί (./album.dff, Album.dff, σύνδεσμος).
        std::error_code e1, e2;
        if (std::filesystem::exists(argv[2], e1) && std::filesystem::equivalent(argv[2], argv[3], e2)) {
            std::fprintf(stderr, "--dsd-out would write over its own input\n");
            return 2;
        }
    }
    for (int a = first; a < argc; ++a)
        if (argv[a][0] == '-' && argv[a][1] == '-') return usage();   // --help ή διακόπτης σε λάθος θέση
    for (int a = first; a < argc; ++a) try {
        File file(argv[a]);
        if (!file.f) {
            std::printf("CANNOT-READ %s\n", argv[a]);
            ++failed;
            continue;
        }
        uint8_t magic[6] = {};
        Input in;
        const char *why = nullptr;
        if (file.read(0, magic, 6) && std::memcmp(magic, "DSTF1\n", 6) == 0) why = fromDump(file, in);
        else if (std::memcmp(magic, "FRM8", 4) == 0) why = fromDff(file, in);
        else why = "NOT-DFF-OR-FRAME-DUMP";
        if (why) {
            std::printf("%s %s\n", why, argv[a]);
            ++failed;
            continue;
        }
        // Μόνο αληθινοί ρυθμοί DSD (πολλαπλάσια του 44100 x 64, ως το DSD512): ένας παράλογος ρυθμός σε φτιαγμένο αρχείο θα
        // ζητούσε τεράστιες ποσότητες ήχου για κάθε πλαίσιο των λίγων bytes.
        const bool realRate = in.rate >= 2822400 && in.rate <= 22579200 && in.rate % 2822400 == 0;
        const int perChannel = realRate ? int(in.rate / 8 / 75) : 0;
        DstDecoder dec;
        if (!realRate || in.channels < 1 || in.channels > 6 || !dec.init(in.channels) || perChannel <= 0 || perChannel > kDstMaxBytesPerChannel) {
            std::printf("BAD-SHAPE %s\n", argv[a]);
            ++failed;
            continue;
        }
        const size_t slot = size_t(in.channels) * size_t(perChannel);
        std::FILE *out = nullptr;
        bool writeOk = true;
        if (dsdOut) {
            out = std::fopen(dsdOut, "wb");
            if (!out) {
                std::printf("CANNOT-WRITE %s\n", dsdOut);
                ++failed;
                continue;
            }
            const std::vector<uint8_t> head = dffHeader(in.channels, in.rate, uint64_t(in.frames.size()) * slot);
            writeOk = std::fwrite(head.data(), 1, head.size(), out) == head.size();
        }
        std::vector<uint8_t> dsd(slot), frame;
        Sha256 sha;
        uint64_t fnv = 14695981039346656037ULL, decoded = 0, refused = 0;   // FNV-1a 64 (0xcbf29ce484222325)
        bool readError = false;
        for (const auto &fr : in.frames) {
            // Ένα πλαίσιο τη φορά στη μνήμη. Το μήκος του έχει ήδη ελεγχθεί απέναντι στο μέγεθος του αρχείου.
            frame.resize(fr.second);
            if (!file.read(fr.first, frame.data(), fr.second)) {   // ο δίσκος δεν έδωσε όσα έδειχνε το μέγεθος: σταματάμε
                readError = true;
                break;
            }
            if (dec.decodeFrame(frame.data(), int(fr.second), dsd.data(), perChannel)) {
                ++decoded;
            } else {
                std::memset(dsd.data(), 0x69, slot);
                ++refused;
            }
            sha.add(dsd.data(), slot);
            for (uint8_t b : dsd) { fnv ^= b; fnv *= 1099511628211ULL; }
            if (out && writeOk) writeOk = std::fwrite(dsd.data(), 1, slot, out) == slot;
        }
        if (out) {
            if (writeOk && (uint64_t(in.frames.size()) * slot) & 1) writeOk = std::fputc(0, out) == 0;
            if (std::fclose(out) != 0) writeOk = false;
        }
        if (in.cut || readError) ++failed;
        // «frames» = όσα αποσυμπιέστηκαν· «refused» = όσα αρνήθηκε ο αποσυμπιεστής (στη θέση τους μπήκε σιωπή DSD).
        std::printf("%s fnv %016llx %d ch %u Hz frames %llu refused %llu%s%s %s\n", sha.hex().c_str(), (unsigned long long)fnv, in.channels, in.rate,
                    (unsigned long long)decoded, (unsigned long long)refused, in.cut ? " FILE-CUT-SHORT" : "", readError ? " READ-ERROR" : "", argv[a]);
        if (dsdOut && (!writeOk || readError)) {
            std::remove(dsdOut);   // ένα μισό αρχείο θα έμοιαζε σωστό DSDIFF· καλύτερα κανένα
            std::printf("CANNOT-WRITE %s\n", dsdOut);
            ++failed;
        }
    } catch (const std::bad_alloc &) {
        std::printf("OUT-OF-MEMORY %s\n", argv[a]);
        ++failed;
    }
    return failed ? 1 : 0;
}
