// flachash: αποκωδικοποιεί ένα αρχείο .flac με τον αποκωδικοποιητή FLAC του Orpheic και τυπώνει μια σύνοψη SHA-256 των δειγμάτων.
//
//   flachash αρχείο.flac ...
//
// Παράδειγμα χρήσης του flacdecoder.h, και έλεγχος: το ίδιο αρχείο πρέπει να δίνει την ίδια γραμμή σε κάθε μηχάνημα και
// μεταγλωττιστή. Χωρίς καμία εξάρτηση.
// Τα δείγματα συνοψίζονται όπως τα ορίζει το ίδιο το FLAC για το MD5 του STREAMINFO (RFC 9639, 8.2): πλεγμένα ανά κανάλι, κάθε
// δείγμα σε όσα ολόκληρα bytes χρειάζεται το βάθος του, με το μικρό byte πρώτο και το πρόσημο επεκταμένο. Άρα:
//   · η σύνοψη SHA-256 είναι η ίδια με του `sha256sum` πάνω σε ακατέργαστο PCM του ίδιου βάθους (π.χ. s16le, s24le) από
//     οποιονδήποτε άλλον αποκωδικοποιητή·
//   · το MD5 των ίδιων bytes συγκρίνεται με αυτό που έγραψε ο κωδικοποιητής στο αρχείο, όταν υπάρχει («md5 same»).
// Ένα καρέ που δεν αποκωδικοποιείται μετριέται («bad») και η ανάγνωση συνεχίζει από την επόμενη έγκυρη κεφαλίδα καρέ.
#include "flacdecoder.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
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

// MD5 (RFC 1321), μόνο για τη σύγκριση με το STREAMINFO.
struct Md5 {
    uint32_t h[4] = {0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476};
    uint8_t block[64];
    size_t fill = 0;
    uint64_t total = 0;

    static uint32_t rotl(uint32_t x, int n) { return (x << n) | (x >> (32 - n)); }
    void compress(const uint8_t *p)
    {
        static const uint32_t k[64] = {
            0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501, 0x698098d8, 0x8b44f7af, 0xffff5bb1,
            0x895cd7be, 0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821, 0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa, 0xd62f105d, 0x02441453,
            0xd8a1e681, 0xe7d3fbc8, 0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed, 0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a, 0xfffa3942,
            0x8771f681, 0x6d9d6122, 0xfde5380c, 0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70, 0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05,
            0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665, 0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92, 0xffeff47d,
            0x85845dd1, 0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1, 0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391};
        static const int r[64] = {7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20,
                                  4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21};
        uint32_t w[16];
        for (int i = 0; i < 16; ++i) w[i] = uint32_t(p[4 * i]) | uint32_t(p[4 * i + 1]) << 8 | uint32_t(p[4 * i + 2]) << 16 | uint32_t(p[4 * i + 3]) << 24;
        uint32_t a = h[0], b = h[1], c = h[2], d = h[3];
        for (int i = 0; i < 64; ++i) {
            uint32_t f;
            int g;
            if (i < 16) { f = (b & c) | (~b & d); g = i; }
            else if (i < 32) { f = (d & b) | (~d & c); g = (5 * i + 1) & 15; }
            else if (i < 48) { f = b ^ c ^ d; g = (3 * i + 5) & 15; }
            else { f = c ^ (b | ~d); g = (7 * i) & 15; }
            const uint32_t t = d;
            d = c; c = b;
            b = b + rotl(a + f + k[i] + w[g], r[i]);
            a = t;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d;
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
    void finish(uint8_t out[16])
    {
        const uint64_t bits = total * 8;
        const uint8_t one = 0x80, zero = 0;
        add(&one, 1);
        while (fill != 56) add(&zero, 1);
        uint8_t len[8];
        for (int i = 0; i < 8; ++i) len[i] = uint8_t(bits >> (8 * i));
        add(len, 8);
        for (int i = 0; i < 16; ++i) out[i] = uint8_t(h[i / 4] >> (8 * (i % 4)));
    }
};

bool readFile(const char *path, std::vector<uint8_t> &data)
{
    std::FILE *f = std::fopen(path, "rb");
    if (!f) return false;
    uint8_t buf[1 << 16];
    size_t n;
    while ((n = std::fread(buf, 1, sizeof buf, f)) > 0) data.insert(data.end(), buf, buf + n);
    std::fclose(f);
    return true;
}

}  // namespace

int main(int argc, char **argv)
{
    using namespace orpheic::flac;
    int failed = 0;
    for (int a = 1; a < argc; ++a) {
        std::vector<uint8_t> data;
        if (!readFile(argv[a], data)) {
            std::printf("CANNOT-READ %s\n", argv[a]);
            ++failed;
            continue;
        }
        const size_t n = data.size();
        data.resize(n + 16, 0);   // ο αποκωδικοποιητής θέλει να ΔΙΑΒΑΖΟΝΤΑΙ 8 bytes μετά το τέλος (flacdecoder.h, decode)
        const uint8_t *p = data.data();
        StreamInfo si;
        size_t pos = 0;
        if (!readHead(p, n, si, pos)) {
            std::printf("NOT-FLAC %s\n", argv[a]);
            ++failed;
            continue;
        }
        Decoder dec(si);
        Sha256 sha;
        Md5 md5;
        uint64_t frames = 0, bad = 0, samples = 0;
        int rate = 0, channels = 0, bits = 0;
        std::vector<uint8_t> pcm;
        while (pos < n) {
            size_t used = 0;
            const Status s = dec.decode(p + pos, n - pos, used);
            if (s != Status::Ok || !used) {   // όχι καρέ εδώ: στο επόμενο σημείο όπου διαβάζεται έγκυρη κεφαλίδα
                ++bad;
                pos = findFrame(p, n, pos + 1, si);
                continue;
            }
            pos += used;
            ++frames;
            const FrameHead &h = dec.head();
            rate = int(h.sampleRate); channels = h.channels; bits = h.bits;
            const int width = (h.bits + 7) / 8;
            pcm.resize(size_t(h.blockSize) * size_t(h.channels) * size_t(width));
            uint8_t *o = pcm.data();
            for (uint32_t i = 0; i < h.blockSize; ++i)
                for (int c = 0; c < h.channels; ++c) {
                    const uint32_t v = uint32_t(dec.channel(c)[i]);   // συμπλήρωμα ως προς δύο: τα πάνω bytes είναι το πρόσημο
                    for (int b = 0; b < width; ++b) *o++ = uint8_t(v >> (8 * b));
                }
            sha.add(pcm.data(), pcm.size());
            md5.add(pcm.data(), pcm.size());
            samples += h.blockSize;
        }
        if (!frames) {
            std::printf("NO-FRAMES %s\n", argv[a]);
            ++failed;
            continue;
        }
        const char *verdict = "md5 absent";
        if (si.present && si.hasMd5()) {
            uint8_t got[16];
            md5.finish(got);
            const bool same = !std::memcmp(got, si.md5, 16);
            verdict = same ? "md5 same" : "MD5 DIFFERS";
            if (!same) ++failed;
        }
        std::printf("%s %d Hz %d ch %d bit frames %llu samples %llu", sha.hex().c_str(), rate, channels, bits, (unsigned long long)frames,
                    (unsigned long long)samples);
        if (si.present && si.totalSamples) std::printf(" of %llu", (unsigned long long)si.totalSamples);
        std::printf(" bad %llu %s %s\n", (unsigned long long)bad, verdict, argv[a]);
    }
    return failed ? 1 : 0;
}
