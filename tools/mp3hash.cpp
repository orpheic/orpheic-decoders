// mp3hash: decode an .mp3 file with the Orpheic MP3 decoder and print a SHA-256 of everything it produced.
//
//   mp3hash [--exact] file.mp3 ...
//
// For comparing two builds of the decoder: the same file must give the same line. No Qt, no FFmpeg.
// The frames are found here, in the simplest way that is the same for every build: an ID3v2 tag at the start is
// skipped, a frame is accepted when its header parses as Layer III and the header that should follow it parses too
// with the same version and sample rate (or the file ends there). Bytes between frames are skipped, and the decoder
// is told that data was lost. Free-format frames (no frame length in the header) are not framed: a file made only
// of them gets a line that says so.
// The hash covers, for every frame in order: the status the decoder returned, the number of channels and samples,
// and the samples of each channel as 32-bit floats.
#include "mp3decoder.h"

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
    using namespace orpheic::mp3;
    bool exact = false;
    int failed = 0;
    for (int a = 1; a < argc; ++a) {
        if (!std::strcmp(argv[a], "--exact")) { exact = true; continue; }
        std::vector<uint8_t> data;
        if (!readFile(argv[a], data)) {
            std::printf("CANNOT-READ %s\n", argv[a]);
            ++failed;
            continue;
        }
        const uint8_t *p = data.data();
        const size_t n = data.size();
        size_t pos = 0;
        if (n >= 10 && !std::memcmp(p, "ID3", 3))   // ID3v2: ten bytes of header, then a size of four 7-bit bytes
            pos = 10 + (size_t(p[6] & 0x7F) << 21 | size_t(p[7] & 0x7F) << 14 | size_t(p[8] & 0x7F) << 7 | size_t(p[9] & 0x7F));
        Decoder dec;
        dec.setExact(exact);
        Sha256 sha;
        uint64_t frames = 0, counts[6] = {}, samples = 0, skipped = 0;
        bool freeFormat = false, gap = false;
        while (pos + 4 <= n) {
            FrameHead h;
            bool ok = parseHead(p + pos, n - pos, h) == Status::Ok && h.layer == 3;
            if (ok && !h.frameBytes) {   // free format, or four bytes of something else that look like it: not a frame for this tool
                freeFormat = true;
                ok = false;
            }
            if (ok) {
                const size_t next = pos + size_t(h.frameBytes);
                if (next + 4 <= n) {
                    FrameHead h2;
                    ok = parseHead(p + next, n - next, h2) == Status::Ok && h2.layer == 3 && h2.version == h.version && h2.sampleRate == h.sampleRate;
                }
            }
            if (!ok) {
                ++pos;
                ++skipped;
                gap = true;
                continue;
            }
            if (gap && frames) dec.lostData();
            gap = false;
            const size_t len = pos + size_t(h.frameBytes) <= n ? size_t(h.frameBytes) : n - pos;
            const Status s = dec.decode(p + pos, len);
            pos += len;
            ++frames;
            const int code = int(s);
            if (code >= 0 && code < 6) ++counts[code];
            sha.add(&code, sizeof code);
            if (s == Status::Ok || s == Status::BadData || s == Status::NoEarlierData) {
                const FrameHead &hd = dec.head();
                const int shape[2] = {hd.channels, hd.samples};
                sha.add(shape, sizeof shape);
                for (int c = 0; c < hd.channels; ++c) sha.add(dec.channel(c), size_t(hd.samples) * sizeof(float));
                samples += uint64_t(hd.samples);
            }
        }
        if (!frames) {
            std::printf("%s %s\n", freeFormat ? "FREE-FORMAT (not framed by this tool)" : "NO-FRAMES", argv[a]);
            continue;
        }
        std::printf("%s %s frames %llu samples %llu ok %llu truncated %llu notaframe %llu unsupported %llu bad %llu noearlier %llu skippedbytes %llu %s\n",
                    sha.hex().c_str(), exact ? "exact" : "fast", (unsigned long long)frames, (unsigned long long)samples, (unsigned long long)counts[0],
                    (unsigned long long)counts[1], (unsigned long long)counts[2], (unsigned long long)counts[3], (unsigned long long)counts[4],
                    (unsigned long long)counts[5], (unsigned long long)skipped, argv[a]);
    }
    return failed ? 1 : 0;
}
