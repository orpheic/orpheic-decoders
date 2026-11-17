// mp3check: ο γρήγορος δρόμος του αποκωδικοποιητή MP3 απέναντι στον ΑΚΡΙΒΗ (τους τύπους του προτύπου σε διπλή ακρίβεια), δείγμα προς
// δείγμα, και πόσο γρήγορα τρέχει ο γρήγορος.
//
//   mp3check [--rounds N] [--time-exact] αρχείο.mp3 …      (--time-exact: χρονομετρά ΚΑΙ τον ακριβή δρόμο)
//
// Για κάθε αρχείο: η μεγαλύτερη απόσταση ενός δείγματος του γρήγορου από τον ακριβή, η μέση τετραγωνική απόσταση (rms), η στάθμη
// του ίδιου του σήματος (rms), και ο χρόνος του γρήγορου (ο καλύτερος από N περάσματα) ως πολλαπλάσιο της διάρκειας της μουσικής.
// Στο τέλος τα ίδια για όλα τα αρχεία μαζί. Πλήρης κλίμακα = 1.
#include "mp3decoder.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

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



struct Frame { size_t at, len; bool afterGap; };

// Τα καρέ ενός αρχείου, όπως τα βρίσκει και το mp3hash: κεφαλίδα Layer III με γνωστό μήκος, και η επόμενη κεφαλίδα να συμφωνεί.
std::vector<Frame> frameList(const std::vector<uint8_t> &data)
{
    using namespace orpheic::mp3;
    std::vector<Frame> out;
    const uint8_t *p = data.data();
    const size_t n = data.size();
    size_t pos = 0;
    if (n >= 10 && !std::memcmp(p, "ID3", 3))
        pos = 10 + (size_t(p[6] & 0x7F) << 21 | size_t(p[7] & 0x7F) << 14 | size_t(p[8] & 0x7F) << 7 | size_t(p[9] & 0x7F));
    bool gap = false;
    while (pos + 4 <= n) {
        FrameHead h;
        bool ok = parseHead(p + pos, n - pos, h) == Status::Ok && h.layer == 3 && h.frameBytes;
        if (ok) {
            const size_t next = pos + size_t(h.frameBytes);
            if (next + 4 <= n) {
                FrameHead h2;
                ok = parseHead(p + next, n - next, h2) == Status::Ok && h2.layer == 3 && h2.version == h.version && h2.sampleRate == h.sampleRate;
            }
        }
        if (!ok) {
            ++pos;
            gap = true;
            continue;
        }
        const size_t len = pos + size_t(h.frameBytes) <= n ? size_t(h.frameBytes) : n - pos;
        out.push_back({pos, len, gap && !out.empty()});
        gap = false;
        pos += len;
    }
    return out;
}

bool hasSamples(orpheic::mp3::Status s)
{
    using orpheic::mp3::Status;
    return s == Status::Ok || s == Status::BadData || s == Status::NoEarlierData;
}

}  // namespace

int main(int argc, char **argv)
{
    using namespace orpheic::mp3;
    int rounds = 5, files = 0;
    bool timeExact = false;
    double bestExactAll = 0;
    double worstAll = 0, errAll = 0, sigAll = 0, secondsAll = 0, bestAll = 0;
    unsigned long long nAll = 0;
    for (int a = 1; a < argc; ++a) {
        if (!std::strcmp(argv[a], "--rounds") && a + 1 < argc) { rounds = std::atoi(argv[++a]); continue; }
        if (!std::strcmp(argv[a], "--time-exact")) { timeExact = true; continue; }
        std::vector<uint8_t> data;
        if (!readFile(argv[a], data)) { std::printf("CANNOT-READ %s\n", argv[a]); continue; }
        const std::vector<Frame> frames = frameList(data);
        if (frames.empty()) { std::printf("NO-FRAMES %s\n", argv[a]); continue; }
        Decoder fast, exact;   // ένας καινούριος Decoder είναι ο ΑΚΡΙΒΗΣ· ο γρήγορος ζητιέται
        fast.setExact(false);
        exact.setExact(true);
        double worst = 0, err = 0, sig = 0, seconds = 0;
        unsigned long long n = 0;
        for (const Frame &f : frames) {
            if (f.afterGap) { fast.lostData(); exact.lostData(); }
            const Status sf = fast.decode(data.data() + f.at, f.len), se = exact.decode(data.data() + f.at, f.len);
            if (!hasSamples(sf) || !hasSamples(se)) continue;
            const FrameHead &h = exact.head();
            seconds += double(h.samples) / h.sampleRate;
            for (int c = 0; c < h.channels; ++c) {
                const float *x = fast.channel(c), *y = exact.channel(c);
                for (int i = 0; i < h.samples; ++i) {
                    const double d = double(x[i]) - double(y[i]);
                    if (std::fabs(d) > worst) worst = std::fabs(d);
                    err += d * d;
                    sig += double(y[i]) * y[i];
                    ++n;
                }
            }
        }
        double best = 1e30;
        for (int r = 0; r < rounds; ++r) {
            Decoder d;
            d.setExact(false);
            const auto t0 = std::chrono::steady_clock::now();
            for (const Frame &f : frames) {
                if (f.afterGap) d.lostData();
                d.decode(data.data() + f.at, f.len);
            }
            const double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (t < best) best = t;
        }
        double bestExact = 1e30;
        for (int r = 0; timeExact && r < (rounds < 2 ? rounds : 2); ++r) {
            Decoder d;
            d.setExact(true);
            const auto t0 = std::chrono::steady_clock::now();
            for (const Frame &f : frames) {
                if (f.afterGap) d.lostData();
                d.decode(data.data() + f.at, f.len);
            }
            const double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (t < bestExact) bestExact = t;
        }
        if (timeExact) bestExactAll += bestExact;
        if (!n) { std::printf("NO-SAMPLES %s\n", argv[a]); continue; }
        std::printf("largest %.3g · rms of difference %.3g · signal rms %.3g · %.1f s of sound in %.2f ms (%.0fx) · %s\n", worst, std::sqrt(err / n),
                    std::sqrt(sig / n), seconds, best * 1e3, seconds / best, argv[a]);
        if (worst > worstAll) worstAll = worst;
        errAll += err; sigAll += sig; nAll += n; secondsAll += seconds; bestAll += best; ++files;
    }
    if (files)
        std::printf("ALL %d files: largest %.3g · rms of difference %.3g · signal rms %.3g · %.1f s of sound in %.1f ms (%.0fx)\n", files, worstAll,
                    std::sqrt(errAll / nAll), std::sqrt(sigAll / nAll), secondsAll, bestAll * 1e3, secondsAll / bestAll);
    if (files && timeExact)
        std::printf("ALL %d files, exact path: %.1f s of sound in %.1f ms (%.0fx)\n", files, secondsAll, bestExactAll * 1e3, secondsAll / bestExactAll);
    return 0;
}
