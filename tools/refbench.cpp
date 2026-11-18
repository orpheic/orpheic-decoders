// refbench: ours vs libFLAC (FLAC) and vs libmpg123 (MP3). Both libraries used only through their public API, as black boxes.
// Method as speedtest: decoder alone, one thread, same input bytes, check round first, N timed rounds in alternating order, median.
//   FLAC: whole file in memory. libFLAC: init_stream + metadata untimed; timed = process_until_end_of_stream (its read callback
//         memcpy's the bytes in: unavoidable in its API, so the memcpy of the same bytes is timed separately and reported).
//         Ours: readHead untimed; timed = decode() frame after frame over the same contiguous bytes. Output on both sides: int32
//         planar in the decoder's own buffers, no conversion. CRC-16 checked on both sides. Check round: every frame memcmp.
//   MP3:  packets from libavformat (pure audio frames, held in memory). mpg123: feed mode, per packet mpg123_feed + decode_frame
//         until NEED_MORE, float32 output left in its buffer, gapless off. Ours: decode() per packet, float planar in its buffer.
#include "flacdecoder.h"
#include "mp3decoder.h"
#include <FLAC/stream_decoder.h>
#include <mpg123.h>
extern "C" {
#include <libavformat/avformat.h>
}
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using Clock = std::chrono::steady_clock;
static inline void escape(const void *p) { asm volatile("" : : "g"(p) : "memory"); }
static double now() { return std::chrono::duration<double>(Clock::now().time_since_epoch()).count(); }
static double median(std::vector<double> v) { std::sort(v.begin(), v.end()); size_t n = v.size(); return n % 2 ? v[n/2] : (v[n/2-1]+v[n/2])/2; }
static double spread(const std::vector<double> &v) { auto mm = std::minmax_element(v.begin(), v.end()); return (*mm.second-*mm.first)/median(v); }

// ---------------- FLAC ----------------
struct Mem { const uint8_t *p; size_t n, pos = 0; };
struct FlacCtx {
    Mem in; long frames = 0, errors = 0; uint64_t samples = 0;
    orpheic::flac::Decoder *ours = nullptr; const uint8_t *op = nullptr; size_t on = 0, opos = 0; bool same = true; long oursBad = 0;
};
static FLAC__StreamDecoderReadStatus rd(const FLAC__StreamDecoder *, FLAC__byte b[], size_t *bytes, void *c)
{
    Mem &m = static_cast<FlacCtx *>(c)->in;
    size_t k = std::min(*bytes, m.n - m.pos);
    if (!k) { *bytes = 0; return FLAC__STREAM_DECODER_READ_STATUS_END_OF_STREAM; }
    std::memcpy(b, m.p + m.pos, k); m.pos += k; *bytes = k;
    return FLAC__STREAM_DECODER_READ_STATUS_CONTINUE;
}
static FLAC__StreamDecoderWriteStatus wrTimed(const FLAC__StreamDecoder *, const FLAC__Frame *f, const FLAC__int32 *const buf[], void *c)
{
    for (unsigned i = 0; i < f->header.channels; ++i) escape(buf[i]);
    static_cast<FlacCtx *>(c)->samples += f->header.blocksize;
    return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE;
}
static FLAC__StreamDecoderWriteStatus wrCheck(const FLAC__StreamDecoder *, const FLAC__Frame *f, const FLAC__int32 *const buf[], void *c)
{
    FlacCtx &x = *static_cast<FlacCtx *>(c);
    ++x.frames; x.samples += f->header.blocksize;
    size_t used = 0;
    if (x.ours->decode(x.op + x.opos, x.on - x.opos, used) != orpheic::flac::Status::Ok) { ++x.oursBad; x.same = false; return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE; }
    x.opos += used;
    const auto &h = x.ours->head();
    if (h.blockSize != f->header.blocksize || unsigned(h.channels) != f->header.channels) { x.same = false; return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE; }
    for (unsigned ch = 0; ch < f->header.channels; ++ch)
        if (std::memcmp(x.ours->channel(int(ch)), buf[ch], size_t(h.blockSize) * 4)) x.same = false;
    return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE;
}
static void err(const FLAC__StreamDecoder *, FLAC__StreamDecoderErrorStatus, void *c) { ++static_cast<FlacCtx *>(c)->errors; }

static FLAC__StreamDecoder *libflacOpen(FlacCtx &x, bool check)
{
    FLAC__StreamDecoder *d = FLAC__stream_decoder_new();
    FLAC__stream_decoder_set_md5_checking(d, check);   // only in the check round: it confirms libFLAC's own output
    FLAC__stream_decoder_set_metadata_ignore_all(d);
    FLAC__stream_decoder_init_stream(d, rd, nullptr, nullptr, nullptr, nullptr, check ? wrCheck : wrTimed, nullptr, err, &x);
    FLAC__stream_decoder_process_until_end_of_metadata(d);
    return d;
}
static double libflacRound(const std::vector<uint8_t> &file, uint64_t *samples)
{
    FlacCtx x; x.in = {file.data(), file.size() - 8};
    FLAC__StreamDecoder *d = libflacOpen(x, false);
    double t0 = now();
    bool ok = FLAC__stream_decoder_process_until_end_of_stream(d);
    double t = now() - t0;
    FLAC__stream_decoder_finish(d); FLAC__stream_decoder_delete(d);
    if (samples) *samples = x.samples;
    return ok && !x.errors ? t : -1;
}
static double oursFlacRound(const std::vector<uint8_t> &file, const orpheic::flac::StreamInfo &si, size_t start, uint64_t *samples)
{
    orpheic::flac::Decoder dec(si); dec.setCheckCrc(true);
    const size_t n = file.size() - 8; size_t pos = start; uint64_t s = 0;
    double t0 = now();
    while (pos < n) {
        size_t used = 0;
        if (dec.decode(file.data() + pos, n - pos, used) != orpheic::flac::Status::Ok) return -1;
        for (int c = 0; c < dec.head().channels; ++c) escape(dec.channel(c));
        s += dec.head().blockSize; pos += used;
    }
    double t = now() - t0;
    if (samples) *samples = s;
    return t;
}
static double memcpyRound(const std::vector<uint8_t> &file)
{   // what libFLAC's read callback does, alone: copy the bytes in chunks into a fixed buffer
    static std::vector<uint8_t> b(65536);
    double t0 = now();
    for (size_t p = 0; p < file.size(); p += b.size()) { std::memcpy(b.data(), file.data() + p, std::min(b.size(), file.size() - p)); escape(b.data()); }
    return now() - t0;
}

static bool flacFile(const char *path, int rounds, std::vector<double> &ratios, std::vector<double> &ratiosNoCopy)
{
    FILE *f = std::fopen(path, "rb"); if (!f) return false;
    std::vector<uint8_t> file; uint8_t tmp[1 << 16]; size_t k;
    while ((k = std::fread(tmp, 1, sizeof tmp, f))) file.insert(file.end(), tmp, tmp + k);
    std::fclose(f);
    file.resize(file.size() + 8, 0);   // our decoder may read 8 bytes past the end
    orpheic::flac::StreamInfo si; size_t start = 0;
    if (!orpheic::flac::readHead(file.data(), file.size() - 8, si, start)) { std::printf("SKIPPED %s\n", path); return false; }
    // check round
    FlacCtx x; x.in = {file.data(), file.size() - 8};
    orpheic::flac::Decoder ours(si); ours.setCheckCrc(true);
    x.ours = &ours; x.op = file.data(); x.on = file.size() - 8; x.opos = start;
    FLAC__StreamDecoder *d = libflacOpen(x, true);
    bool ok = FLAC__stream_decoder_process_until_end_of_stream(d);
    bool md5 = FLAC__stream_decoder_finish(d); FLAC__stream_decoder_delete(d);
    const char *name = std::strrchr(path, '/') ? std::strrchr(path, '/') + 1 : path;
    bool same = ok && x.same && !x.errors && x.opos == file.size() - 8 && md5 && (!si.totalSamples || x.samples == si.totalSamples);
    std::printf("FLAC check %-16s frames %ld samples %llu (STREAMINFO %llu) libFLAC errors %ld, libFLAC MD5 %s, every frame %s, ours consumed all bytes %s\n",
                name, x.frames, (unsigned long long)x.samples, (unsigned long long)si.totalSamples, x.errors, md5 ? "ok" : "MISMATCH",
                x.same ? "IDENTICAL" : "DIFFERENT", x.opos == file.size() - 8 ? "yes" : "NO");
    if (!same) return false;
    std::vector<double> tl, to, tm;
    for (int r = 0; r < rounds; ++r) {
        if (r % 2 == 0) tl.push_back(libflacRound(file, nullptr));
        to.push_back(oursFlacRound(file, si, start, nullptr));
        if (r % 2 == 1) tl.push_back(libflacRound(file, nullptr));
        tm.push_back(memcpyRound(file));
    }
    if (*std::min_element(tl.begin(), tl.end()) < 0 || *std::min_element(to.begin(), to.end()) < 0) { std::printf("NOT COUNTED %s: a round failed\n", name); return false; }
    double ml = median(tl), mo = median(to), mm = median(tm);
    std::printf("FLAC %2d bit %5.1f kHz %6.1f s  libFLAC %.4f s ±%2.0f%%  ours %.4f s ±%2.0f%%  %5.2f times  (memcpy of input %.4f s; without it %5.2f)  %s\n",
                si.bits, si.sampleRate / 1000.0, double(si.totalSamples) / si.sampleRate, ml, 100 * spread(tl), mo, 100 * spread(to), ml / mo, mm, (ml - mm) / mo, name);
    ratios.push_back(ml / mo); ratiosNoCopy.push_back((ml - mm) / mo);
    return true;
}

// ---------------- MP3 ----------------
struct Pk { std::vector<uint8_t> b; };
static bool readMp3(const char *path, std::vector<Pk> &pk, int &rate, int &ch, long long &br)
{
    AVFormatContext *fmt = nullptr;
    if (avformat_open_input(&fmt, path, nullptr, nullptr) < 0) return false;
    avformat_find_stream_info(fmt, nullptr);
    int idx = av_find_best_stream(fmt, AVMEDIA_TYPE_AUDIO, -1, -1, nullptr, 0);
    rate = fmt->streams[idx]->codecpar->sample_rate; ch = fmt->streams[idx]->codecpar->ch_layout.nb_channels; br = fmt->streams[idx]->codecpar->bit_rate;
    AVPacket *p = av_packet_alloc();
    while (av_read_frame(fmt, p) >= 0) { if (p->stream_index == idx && p->size > 0) pk.push_back({std::vector<uint8_t>(p->data, p->data + p->size)}); av_packet_unref(p); }
    av_packet_free(&p); avformat_close_input(&fmt);
    return !pk.empty();
}
static mpg123_handle *mpgOpen()
{
    int e = 0; mpg123_handle *m = mpg123_new(nullptr, &e);
    mpg123_param(m, MPG123_ADD_FLAGS, MPG123_QUIET, 0);
    mpg123_param(m, MPG123_REMOVE_FLAGS, MPG123_GAPLESS, 0);
    mpg123_format_none(m);
    for (long r : {8000L, 11025L, 12000L, 16000L, 22050L, 24000L, 32000L, 44100L, 48000L}) mpg123_format(m, r, MPG123_MONO | MPG123_STEREO, MPG123_ENC_FLOAT_32);
    mpg123_open_feed(m);
    return m;
}
// one packet in; every frame that comes out goes to `each` (nullptr: just touch it)
template <class F> static long mpgFeed(mpg123_handle *m, const Pk &p, F each)
{
    long errs = 0;
    if (mpg123_feed(m, p.b.data(), p.b.size()) != MPG123_OK) ++errs;
    for (;;) {
        off_t num; unsigned char *audio; size_t bytes;
        int r = mpg123_decode_frame(m, &num, &audio, &bytes);
        if (r == MPG123_NEED_MORE) break;
        if (r == MPG123_NEW_FORMAT) continue;
        if (r != MPG123_OK) { ++errs; break; }
        each(audio, bytes);
    }
    return errs;
}
static double mpgRound(const std::vector<Pk> &pk)
{
    mpg123_handle *m = mpgOpen(); long errs = 0;
    double t0 = now();
    for (const Pk &p : pk) errs += mpgFeed(m, p, [](unsigned char *a, size_t) { escape(a); });
    double t = now() - t0;
    mpg123_delete(m);
    return errs ? -1 : t;
}
static double oursMp3Round(const std::vector<Pk> &pk, bool exact)
{
    orpheic::mp3::Decoder dec; dec.setExact(exact);
    double t0 = now();
    for (const Pk &p : pk) {
        auto st = dec.decode(p.b.data(), p.b.size());
        if (st == orpheic::mp3::Status::Ok || st == orpheic::mp3::Status::BadData || st == orpheic::mp3::Status::NoEarlierData)
            for (int c = 0; c < dec.head().channels; ++c) escape(dec.channel(c));
    }
    return now() - t0;
}
struct Gap { double mx = 0, sq = 0; uint64_t n = 0; void add(double a, double b) { double d = std::fabs(a - b); mx = std::max(mx, d); sq += d * d; ++n; } double rms() const { return n ? std::sqrt(sq / n) : 0; } };

static bool mp3File(const char *path, int rounds, std::vector<double> &ratios, std::vector<double> &gF, std::vector<double> &gX, std::vector<double> &gM)
{
    const char *name = std::strrchr(path, '/') ? std::strrchr(path, '/') + 1 : path;
    std::vector<Pk> pk; int rate = 0, chs = 0; long long br = 0;
    if (!readMp3(path, pk, rate, chs, br)) { std::printf("SKIPPED %s\n", name); return false; }
    // check round: full streams (float), then aligned comparison
    std::vector<float> sm, sf, sx;
    mpg123_handle *m = mpgOpen(); long merr = 0;
    orpheic::mp3::Decoder fast, exact; fast.setExact(false); exact.setExact(true);
    long bad = 0;
    for (const Pk &p : pk) {
        merr += mpgFeed(m, p, [&](unsigned char *a, size_t bytes) { const float *f = reinterpret_cast<const float *>(a); sm.insert(sm.end(), f, f + bytes / 4); });
        for (auto *d : {&fast, &exact}) {
            auto st = d->decode(p.b.data(), p.b.size());
            if (st == orpheic::mp3::Status::BadData && d == &fast) ++bad;
            if (st != orpheic::mp3::Status::Ok && st != orpheic::mp3::Status::BadData && st != orpheic::mp3::Status::NoEarlierData) { if (d == &fast) ++bad; continue; }
            auto &q = d == &fast ? sf : sx;
            for (int i = 0; i < d->head().samples; ++i) for (int c = 0; c < d->head().channels; ++c) q.push_back(d->channel(c)[i]);
        }
    }
    long mframes = 0; { const char *dn = mpg123_current_decoder(m); static bool said = false; if (!said) { std::printf("mpg123 decoder in use: %s\n", dn); said = true; } }
    mpg123_delete(m);
    (void)mframes;
    // alignment: mpg123 may leave out frames at the start (main data in a frame it never saw); find the offset (whole frames) that fits
    const size_t spf = size_t(chs) * (rate < 32000 ? 576 : 1152);
    size_t off = 0; double best = 1e9;
    for (size_t o = 0; o <= 4 * spf && o < sf.size(); o += spf) {
        double mx = 0; size_t n = std::min(sm.size(), sf.size() - o); n = std::min<size_t>(n, 200000);
        for (size_t i = 0; i < n; ++i) mx = std::max(mx, double(std::fabs(sf[o + i] - sm[i])));
        if (mx < best) { best = mx; off = o; }
    }
    Gap fm, xm, fx;
    size_t n = std::min(sm.size(), sf.size() - off);
    for (size_t i = 0; i < n; ++i) { fm.add(sf[off + i], sm[i]); xm.add(sx[off + i], sm[i]); }
    for (size_t i = 0; i < sf.size(); ++i) fx.add(sf[i], sx[i]);
    std::printf("MP3 check %-28s ours %zu samples, mpg123 %zu (mpg123 starts %zu frame(s) later, ends %ld samples apart), errors mpg123 %ld ours-bad %ld | "
                "largest gap: fast-vs-mpg123 %.3g (rms %.2g)  exact-vs-mpg123 %.3g (rms %.2g)  fast-vs-exact %.3g\n",
                name, sf.size() / chs, sm.size() / chs, off / spf, long(sf.size()) - long(off) - long(sm.size()), merr, bad, fm.mx, fm.rms(), xm.mx, xm.rms(), fx.mx);
    if (merr || bad) return false;
    gF.push_back(fm.mx); gX.push_back(xm.mx); gM.push_back(fx.mx);
    std::vector<double> tm, to;
    for (int r = 0; r < rounds; ++r) {
        if (r % 2 == 0) tm.push_back(mpgRound(pk));
        to.push_back(oursMp3Round(pk, false));
        if (r % 2 == 1) tm.push_back(mpgRound(pk));
    }
    if (*std::min_element(tm.begin(), tm.end()) < 0) { std::printf("NOT COUNTED %s\n", name); return false; }
    double a = median(tm), b = median(to);
    std::printf("MP3 %4lld kbps %5.1f kHz  mpg123 %.4f s ±%2.0f%%  ours fast %.4f s ±%2.0f%%  %5.2f times  %s\n", br / 1000, rate / 1000.0, a, 100 * spread(tm), b, 100 * spread(to), a / b, name);
    ratios.push_back(a / b);
    return true;
}

static void summary(const char *what, std::vector<double> v)
{
    if (v.empty()) return;
    auto mm = std::minmax_element(v.begin(), v.end());
    std::printf("  %-48s median %.3g · lowest %.3g · highest %.3g · files %zu\n", what, median(v), *mm.first, *mm.second, v.size());
}

int main(int argc, char **argv)
{
    int rounds = 5, first = 1;
    if (argc > 2 && !std::strcmp(argv[1], "--rounds")) { rounds = std::atoi(argv[2]); first = 3; }
    av_log_set_level(AV_LOG_ERROR);
    unsigned ma, mi, pa; mpg123_distversion(&ma, &mi, &pa);
    std::printf("libFLAC %s · libmpg123 %u.%u.%u · %d rounds, median, alternating order · one thread · flac kernels: %s\n", FLAC__VERSION_STRING, ma, mi, pa, rounds, orpheic::flac::kernelText());
    std::vector<double> fr, frn, mr, gF, gX, gM; int failed = 0;
    for (int a = first; a < argc; ++a) {
        std::string p = argv[a];
        bool ok = p.size() > 5 && p.substr(p.size() - 5) == ".flac" ? flacFile(argv[a], rounds, fr, frn) : mp3File(argv[a], rounds, mr, gF, gX, gM);
        failed += !ok;
        std::fflush(stdout);
    }
    std::printf("\nSUMMARY (times faster than the other library; gaps as fraction of full scale)\n");
    summary("FLAC vs libFLAC", fr);
    summary("FLAC vs libFLAC, its input memcpy subtracted", frn);
    summary("MP3 fast path vs mpg123", mr);
    summary("MP3 largest gap, fast path vs mpg123", gF);
    summary("MP3 largest gap, exact path vs mpg123", gX);
    summary("MP3 largest gap, fast path vs exact path", gM);
    if (failed) std::printf("%d file(s) not counted\n", failed);
    return failed != 0;
}
