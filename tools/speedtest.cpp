// speedtest: πόσο γρήγορα αποκωδικοποιούν οι αποκωδικοποιητές του Orpheic απέναντι στο FFmpeg του συστήματος, στα ίδια δεδομένα.
//
//   speedtest [--rounds N] αρχείο.mp3|αρχείο.flac|αρχείο.dff ...
//
// Το μόνο εργαλείο εδώ που χρειάζεται το FFmpeg (libavformat, libavcodec, libavutil, 5.1 ή νεότερο): χτίζεται όταν το βρίσκει
// το CMake.
//
// ΠΩΣ ΜΕΤΡΑ. Τα πακέτα κάθε αρχείου τα διαβάζει ΜΙΑ φορά το libavformat (τα καρέ MP3 και FLAC, τα πλαίσια DST ενός .dff) και
// μένουν στη μνήμη (άρα χρειάζεται μνήμη όσο το αρχείο). Τα ΙΔΙΑ πακέτα δίνονται και στους δύο αποκωδικοποιητές· στο FFmpeg ως
// πακέτα με δικό τους buffer, ώστε να μην τα αντιγράφει σε κάθε γύρο. Μετριέται μόνο η αποκωδικοποίηση ως τα δείγματα στη μορφή
// εξόδου: ο FLAC μας γράφει τα δείγματά του στην ΙΔΙΑ μορφή που βγάζει το FFmpeg (s16/s32, πλεγμένα ή χωριστά, με την ίδια
// ολίσθηση), και ο MP3 μας αφήνει float ανά κανάλι όπως το mp3float. Ένα νήμα ο καθένας. Το άνοιγμα του αποκωδικοποιητή δεν
// μετριέται. Πρώτα ένας γύρος ελέγχου (και ζέσταμα), μετά οι γύροι με σειρά που εναλλάσσεται (FFmpeg πρώτο στους μονούς, δικός
// μας πρώτος στους ζυγούς), και κρατιέται ο ΔΙΑΜΕΣΟΣ χρόνος του καθενός· ο λόγος = διάμεσος FFmpeg / διάμεσος δικού μας.
//
// ΤΙ ΕΛΕΓΧΕΙ, ΚΑΙ ΠΟΤΕ Η ΓΡΑΜΜΗ ΔΕΝ ΜΠΑΙΝΕΙ ΣΤΗ ΣΥΝΟΨΗ:
//   · FLAC: κάθε καρέ byte προς byte ίδιο με του FFmpeg. Ο δικός μας με τον έλεγχο CRC-16 κάθε καρέ ανοιχτό (όπως παίζει)· το
//     FFmpeg από προεπιλογή δεν τον κάνει, που είναι δουλειά ΚΑΤΑ του δικού μας.
//   · MP3: όσα δείγματα βγάζει ο καθένας (πρέπει να είναι ίσα, και κανένα καρέ να μην το προσπερνά ο δικός μας) και η μεγαλύτερη
//     διαφορά από το FFmpeg, για τον γρήγορο και για τον ακριβή δρόμο (το MP3 δεν έχει μία σωστή απάντηση σε μπιτ· βλ. README).
//   · MP3: εκτός σύνοψης όποιο αρχείο έχει καρέ που ο δικός μας δεν αποκωδικοποίησε ή βρήκε να παραβαίνουν το πρότυπο, ή απόκλιση
//     πάνω από kMp3GapLimit.
//   · DST: κανένα πλαίσιο να μην το αρνηθεί ο δικός μας, κανένα λάθος του FFmpeg, και το μέγεθος του πλαισίου να συμφωνεί με την
//     έξοδο του FFmpeg. Το FFmpeg βγάζει PCM: η μετατροπή DSD σε PCM είναι μέσα στον αποκωδικοποιητή του, ενώ ο δικός μας σταματά
//     στο DSD. Γι' αυτό μετριέται ΚΑΙ ο αποκωδικοποιητής DSD του FFmpeg (dsd_msbf, μόνο μετατροπή) πάνω στο ίδιο DSD, αφαιρείται,
//     και ο λόγος ΧΩΡΙΣ τη μετατροπή μπαίνει πρώτος (εκτίμηση: άλλος αποκωδικοποιητής με την ίδια μετατροπή).
//   · Ένας γύρος που απέτυχε βγάζει το αρχείο από τη σύνοψη.
#include "dstdecoder.h"
#include "flacdecoder.h"
#include "mp3decoder.h"

extern "C" {
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/avutil.h>
}

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <new>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

// Φράγμα για τον μεταγλωττιστή: «η μνήμη στο p διαβάζεται από κάπου». Χωρίς αυτό, με βελτιστοποίηση σε όλο το πρόγραμμα (LTO),
// θα μπορούσε να σβήσει δουλειά του δικού μας που δεν διαβάζει κανείς· η πλευρά του FFmpeg είναι κλήσεις βιβλιοθήκης και δεν
// σβήνεται. Δεν προσθέτει καμία δουλειά ανά δείγμα. Είναι το γνωστό δημόσιο μοτίβο του DoNotOptimize της Google Benchmark
// (και της ομιλίας του Chandler Carruth στο CppCon 2015, «Tuning C++: Benchmarks, and CPUs, and Compilers! Oh My!»).
inline void escape(const void *p)
{
#if defined(_MSC_VER)
    static const void *volatile seen;
    seen = p;
    _ReadWriteBarrier();
#else
    asm volatile("" : : "g"(p) : "memory");
#endif
}

struct Stream {
    AVCodecID codec = AV_CODEC_ID_NONE;
    AVCodecParameters *par = nullptr;
    std::vector<AVPacket *> packets;   // με δικό τους buffer και το περιθώριο μηδενικών του FFmpeg στο τέλος
    double seconds = 0;
    ~Stream()
    {
        for (AVPacket *p : packets) av_packet_free(&p);
        avcodec_parameters_free(&par);
    }
};

bool readPackets(const char *path, Stream &s, std::string &why)
{
    AVFormatContext *fmt = nullptr;
    if (avformat_open_input(&fmt, path, nullptr, nullptr) < 0) { why = "cannot open"; return false; }
    if (avformat_find_stream_info(fmt, nullptr) < 0) { avformat_close_input(&fmt); why = "no stream info"; return false; }
    const int idx = av_find_best_stream(fmt, AVMEDIA_TYPE_AUDIO, -1, -1, nullptr, 0);
    if (idx < 0) { avformat_close_input(&fmt); why = "no audio"; return false; }
    s.codec = fmt->streams[idx]->codecpar->codec_id;
    s.par = avcodec_parameters_alloc();
    avcodec_parameters_copy(s.par, fmt->streams[idx]->codecpar);
    if (fmt->duration > 0) s.seconds = double(fmt->duration) / AV_TIME_BASE;
    AVPacket *pkt = av_packet_alloc();
    while (av_read_frame(fmt, pkt) >= 0) {
        if (pkt->stream_index == idx && pkt->size > 0) {
            // Μόνο τα bytes, χωρίς τα παράπλευρα δεδομένα του περιέκτη (π.χ. «πέτα τα πρώτα N δείγματα»): και οι δύο
            // αποκωδικοποιητές παίρνουν ακριβώς το ίδιο.
            AVPacket *c = av_packet_alloc();
            if (av_new_packet(c, pkt->size) < 0) { av_packet_free(&c); break; }
            std::memcpy(c->data, pkt->data, size_t(pkt->size));
            s.packets.push_back(c);
        }
        av_packet_unref(pkt);
    }
    av_packet_free(&pkt);
    avformat_close_input(&fmt);
    if (s.packets.empty()) { why = "no packets"; return false; }
    return true;
}

// Κάθε γύρος πέτυχε (ένας αποτυχημένος γυρίζει −1 και δεν πρέπει να μπει στη διάμεσο).
bool good(const std::vector<double> &v)
{
    return std::all_of(v.begin(), v.end(), [](double x) { return x > 0; });
}

double median(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return n == 0 ? 0 : (n % 2 ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2);
}

// Πόσο σκόρπισαν οι γύροι: (μεγαλύτερος − μικρότερος) / διάμεσος.
double spread(const std::vector<double> &v)
{
    const auto mm = std::minmax_element(v.begin(), v.end());
    const double m = median(v);
    return m > 0 ? (*mm.second - *mm.first) / m : 0;
}

// ---- το FFmpeg ----------------------------------------------------------------------------------------------------------------

struct Ffmpeg {
    AVCodecContext *ctx = nullptr;
    AVFrame *fr = av_frame_alloc();
    long errors = 0;
    uint64_t samples = 0;   // ανά κανάλι
    bool open(const AVCodecParameters *par, AVCodecID id, std::string &why)
    {
        const AVCodec *codec = avcodec_find_decoder(id);
        if (!codec) { why = "FFmpeg has no decoder for this"; return false; }
        ctx = avcodec_alloc_context3(codec);
        avcodec_parameters_to_context(ctx, par);
        ctx->codec_id = id;
        ctx->thread_count = 1;
        if (avcodec_open2(ctx, codec, nullptr) < 0) { why = "FFmpeg cannot open its decoder"; return false; }
        return true;
    }
    ~Ffmpeg()
    {
        av_frame_free(&fr);
        avcodec_free_context(&ctx);
    }
    // Ένα πακέτο (ή nullptr = άδειασμα): ό,τι καρέ βγουν περνούν από το `each`, αν υπάρχει.
    void feed(const AVPacket *p, const std::function<void(const AVFrame *)> *each)
    {
        int r;
        while ((r = avcodec_send_packet(ctx, p)) == AVERROR(EAGAIN)) drain(each);
        if (r < 0 && r != AVERROR_EOF) ++errors;
        drain(each);
    }
    void drain(const std::function<void(const AVFrame *)> *each)
    {
        while (avcodec_receive_frame(ctx, fr) == 0) {
            samples += uint64_t(fr->nb_samples);
            if (each) (*each)(fr);
            av_frame_unref(fr);
        }
    }
};

// Ένας μετρημένος γύρος του FFmpeg: νέο πλαίσιο αποκωδικοποιητή (το άνοιγμα δεν μετριέται), όλα τα πακέτα, το άδειασμα.
double ffmpegRound(const std::vector<AVPacket *> &packets, const AVCodecParameters *par, AVCodecID id, long *errors = nullptr,
                   uint64_t *samples = nullptr)
{
    Ffmpeg f;
    std::string why;
    if (!f.open(par, id, why)) return -1;
    const auto t0 = Clock::now();
    for (const AVPacket *p : packets) f.feed(p, nullptr);
    f.feed(nullptr, nullptr);
    const double sec = std::chrono::duration<double>(Clock::now() - t0).count();
    if (errors) *errors = f.errors;
    if (samples) *samples = f.samples;
    return sec;
}

// ---- FLAC -----------------------------------------------------------------------------------------------------------------------

struct FlacOut {
    orpheic::flac::Decoder::Out fmt = orpheic::flac::Decoder::Out::S32Interleaved;
    int shift = 0, bytes = 4;
    bool planar = false;
};

double flacRound(const Stream &s, const orpheic::flac::StreamInfo &si, const FlacOut &o, std::string &why)
{
    orpheic::flac::Decoder dec(si);
    dec.setCheckCrc(true);
    const int ch = std::max(1, si.channels);
    const size_t block = std::max<uint32_t>(si.maxBlock, 65535);
    std::vector<uint8_t> buf(block * size_t(ch) * size_t(o.bytes));
    std::vector<void *> dst(static_cast<size_t>(ch));
    for (int c = 0; c < ch; ++c) dst[size_t(c)] = buf.data() + (o.planar ? size_t(c) * block * size_t(o.bytes) : 0);
    const auto t0 = Clock::now();
    for (const AVPacket *p : s.packets) {
        size_t used = 0;
        if (dec.decode(p->data, size_t(p->size), used) != orpheic::flac::Status::Ok) { why = "a frame did not decode"; return -1; }
        dec.write(o.fmt, o.shift, dst.data());
        escape(buf.data());
    }
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

// Ο γύρος ελέγχου: FFmpeg και δικός μας πάνω στο ίδιο πακέτο, και κάθε καρέ συγκρίνεται byte προς byte στη μορφή του FFmpeg.
// Μαθαίνει και ποια είναι αυτή η μορφή, για τους μετρημένους γύρους.
bool flacCheck(const Stream &s, const orpheic::flac::StreamInfo &si, FlacOut &o, std::string &why)
{
    Ffmpeg f;
    if (!f.open(s.par, AV_CODEC_ID_FLAC, why)) return false;
    orpheic::flac::Decoder dec(si);
    dec.setCheckCrc(true);
    bool same = true, known = false;
    std::vector<uint8_t> buf;
    std::vector<void *> dst;
    long frames = 0, theirs = 0;
    const std::function<void(const AVFrame *)> each = [&](const AVFrame *fr) {
        ++theirs;
        if (!known) {
            const auto fmt = AVSampleFormat(fr->format);
            o.planar = av_sample_fmt_is_planar(fmt);
            o.bytes = av_get_bytes_per_sample(fmt);
            const bool s16 = av_get_packed_sample_fmt(fmt) == AV_SAMPLE_FMT_S16;
            o.fmt = s16 ? (o.planar ? orpheic::flac::Decoder::Out::S16Planar : orpheic::flac::Decoder::Out::S16Interleaved)
                        : (o.planar ? orpheic::flac::Decoder::Out::S32Planar : orpheic::flac::Decoder::Out::S32Interleaved);
            o.shift = (s16 ? 16 : 32) - si.bits;
            known = true;
        }
        const auto &h = dec.head();
        const int ch = h.channels;
        if (fr->nb_samples != int(h.blockSize) || fr->ch_layout.nb_channels != ch) { same = false; return; }
        const size_t plane = size_t(h.blockSize) * size_t(o.bytes);
        buf.assign(plane * size_t(ch), 0);
        dst.assign(size_t(ch), nullptr);
        for (int c = 0; c < ch; ++c) dst[size_t(c)] = buf.data() + (o.planar ? size_t(c) * plane : 0);
        dec.write(o.fmt, o.shift, dst.data());
        if (o.planar) {
            for (int c = 0; c < ch && same; ++c) same = std::memcmp(fr->extended_data[c], dst[size_t(c)], plane) == 0;
        } else {
            same = std::memcmp(fr->extended_data[0], buf.data(), plane * size_t(ch)) == 0;
        }
    };
    for (const AVPacket *p : s.packets) {
        size_t used = 0;
        if (dec.decode(p->data, size_t(p->size), used) != orpheic::flac::Status::Ok) { why = "a frame did not decode (ours)"; return false; }
        ++frames;
        f.feed(p, &each);   // FLAC: ένα καρέ έξω για κάθε πακέτο, οπότε το καρέ του FFmpeg συγκρίνεται με το δικό μας της ίδιας στιγμής
        if (!same) break;
    }
    f.feed(nullptr, &each);
    if (f.errors) { why = "FFmpeg reported errors"; return false; }
    if (theirs != frames) same = false;
    if (!same) why = "samples DIFFERENT from FFmpeg";
    return same;
}

// ---- MP3 ------------------------------------------------------------------------------------------------------------------------

double mp3Round(const Stream &s, bool exact, long *skipped)
{
    orpheic::mp3::Decoder dec;
    dec.setExact(exact);
    long skip = 0;
    const auto t0 = Clock::now();
    for (const AVPacket *p : s.packets) {
        const auto st = dec.decode(p->data, size_t(p->size));
        if (st == orpheic::mp3::Status::Ok || st == orpheic::mp3::Status::BadData || st == orpheic::mp3::Status::NoEarlierData) {
            for (int c = 0; c < dec.head().channels; ++c) escape(dec.channel(c));
        } else {
            ++skip;
        }
    }
    if (skipped) *skipped = skip;
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

// Ο γύρος ελέγχου του MP3: τρεις ουρές δειγμάτων (FFmpeg, γρήγορος, ακριβής) που συγκρίνονται όσο γεμίζουν, χωρίς να κρατιέται
// όλο το αρχείο. Μετρά και τα δείγματα του καθενός και τα καρέ που ο δικός μας δεν αποκωδικοποίησε.
struct Mp3Check {
    double gapFast = 0, gapExact = 0;
    uint64_t ffSamples = 0, ourSamples = 0;
    long skipped = 0;   // καρέ που ο δικός μας δεν αποκωδικοποίησε καθόλου
    long badData = 0;   // καρέ που ο δικός μας βρήκε να παραβαίνουν το πρότυπο (βγήκε ό,τι μπορούσε)
};
// Πάνω από αυτή την απόκλιση από το FFmpeg (−80 dB της πλήρους κλίμακας) το αρχείο δεν μπαίνει στη σύνοψη: οι δύο δεν έκαναν τότε
// την ίδια δουλειά. Σε κανονικά αρχεία η απόκλιση είναι γύρω στο 0,000001.
constexpr double kMp3GapLimit = 1e-4;
bool mp3Check(const Stream &s, Mp3Check &r, std::string &why)
{
    Ffmpeg f;
    if (!f.open(s.par, AV_CODEC_ID_MP3, why)) return false;
    orpheic::mp3::Decoder fast, exact;
    fast.setExact(false);
    exact.setExact(true);
    std::vector<double> qf, qq, qx;   // πλεγμένα δείγματα, στην κλίμακα ±1
    size_t head = 0;                  // ό,τι συγκρίθηκε ήδη από την αρχή των τριών ουρών
    const auto compare = [&]() {
        const size_t n = std::min({qf.size(), qq.size(), qx.size()});
        for (size_t i = head; i < n; ++i) {
            r.gapFast = std::max(r.gapFast, std::fabs(qq[i] - qf[i]));
            r.gapExact = std::max(r.gapExact, std::fabs(qx[i] - qf[i]));
        }
        head = n;
        if (head > (1u << 20)) {   // ξεφόρτωμα ό,τι συγκρίθηκε
            qf.erase(qf.begin(), qf.begin() + long(head));
            qq.erase(qq.begin(), qq.begin() + long(head));
            qx.erase(qx.begin(), qx.begin() + long(head));
            head = 0;
        }
    };
    const std::function<void(const AVFrame *)> each = [&](const AVFrame *fr) {
        const auto fmt = AVSampleFormat(fr->format);
        const bool planar = av_sample_fmt_is_planar(fmt);
        const int ch = fr->ch_layout.nb_channels;
        for (int i = 0; i < fr->nb_samples; ++i)
            for (int c = 0; c < ch; ++c) {
                const uint8_t *b = planar ? fr->extended_data[c] : fr->extended_data[0];
                const int at = planar ? i : i * ch + c;
                double v = 0;
                switch (av_get_packed_sample_fmt(fmt)) {
                case AV_SAMPLE_FMT_FLT: v = reinterpret_cast<const float *>(b)[at]; break;
                case AV_SAMPLE_FMT_S16: v = reinterpret_cast<const int16_t *>(b)[at] / 32768.0; break;   // ο «σταθερής υποδιαστολής» mp3
                case AV_SAMPLE_FMT_S32: v = reinterpret_cast<const int32_t *>(b)[at] / 2147483648.0; break;
                default: break;
                }
                qf.push_back(v);
            }
    };
    const auto ours = [&](orpheic::mp3::Decoder &d, const AVPacket *p, std::vector<double> &q, bool count) {
        const auto st = d.decode(p->data, size_t(p->size));
        if (st != orpheic::mp3::Status::Ok && st != orpheic::mp3::Status::BadData && st != orpheic::mp3::Status::NoEarlierData) {
            if (count) ++r.skipped;
            return;
        }
        if (count && st == orpheic::mp3::Status::BadData) ++r.badData;
        const auto &h = d.head();
        for (int i = 0; i < h.samples; ++i)
            for (int c = 0; c < h.channels; ++c) q.push_back(d.channel(c)[i]);
        if (count) r.ourSamples += uint64_t(h.samples);
    };
    for (const AVPacket *p : s.packets) {
        f.feed(p, &each);
        ours(fast, p, qq, true);
        ours(exact, p, qx, false);
        compare();
    }
    f.feed(nullptr, &each);
    compare();
    r.ffSamples = f.samples;
    if (f.errors) { why = "FFmpeg reported errors"; return false; }
    return true;
}

// ---- DST ------------------------------------------------------------------------------------------------------------------------

double dstRound(const Stream &s, int channels, int perChannel, long *refused, std::vector<AVPacket *> *keepDsd = nullptr,
                size_t keepMax = 0)
{
    DstDecoder dec;
    if (!dec.init(channels)) return -1;
    std::vector<uint8_t> out(size_t(channels) * size_t(perChannel));
    long bad = 0;
    const auto t0 = Clock::now();
    for (const AVPacket *p : s.packets) {
        if (!dec.decodeFrame(p->data, p->size, out.data(), perChannel)) {
            std::memset(out.data(), 0x69, out.size());
            ++bad;
        }
        escape(out.data());
        if (keepDsd && keepDsd->size() < keepMax) {
            AVPacket *c = av_packet_alloc();
            av_new_packet(c, int(out.size()));
            std::memcpy(c->data, out.data(), out.size());
            keepDsd->push_back(c);
        }
    }
    if (refused) *refused = bad;
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

struct Result {
    std::string kind;
    double ratio;
};

int file(const char *path, int rounds, std::vector<Result> &all)
{
    const char *name = std::strrchr(path, '/');
    name = name ? name + 1 : path;
    Stream s;
    std::string why;
    if (!readPackets(path, s, why)) { std::printf("SKIPPED %s: %s\n", name, why.c_str()); return 1; }

    if (s.codec == AV_CODEC_ID_FLAC) {
        orpheic::flac::StreamInfo si;
        const uint8_t *x = s.par->extradata;
        int n = s.par->extradata_size;
        if (n >= 8 && std::memcmp(x, "fLaC", 4) == 0) { x += 8; n -= 8; }   // «fLaC» και κεφαλίδα μπλοκ μπροστά
        if (n < 34 || !orpheic::flac::parseStreamInfo(x, size_t(n), si)) { std::printf("SKIPPED %s: no STREAMINFO\n", name); return 1; }
        FlacOut o;
        if (!flacCheck(s, si, o, why)) { std::printf("NOT COUNTED %s: %s\n", name, why.c_str()); return 1; }
        std::vector<double> tf, to;
        for (int r = 0; r < rounds; ++r) {
            const bool ffFirst = r % 2 == 0;
            if (ffFirst) tf.push_back(ffmpegRound(s.packets, s.par, AV_CODEC_ID_FLAC));
            to.push_back(flacRound(s, si, o, why));
            if (!ffFirst) tf.push_back(ffmpegRound(s.packets, s.par, AV_CODEC_ID_FLAC));
        }
        if (!good(tf) || !good(to)) { std::printf("NOT COUNTED %s: a timed round failed (%s)\n", name, why.c_str()); return 1; }
        const double ratio = median(tf) / median(to);
        std::printf("FLAC %2d bit %5.1f kHz %7.1f s  ffmpeg %.4f s ±%2.0f%%  ours %.4f s ±%2.0f%%  %5.2f times  samples same  %s\n",
                    si.bits, si.sampleRate / 1000.0, s.seconds, median(tf), 100 * spread(tf), median(to), 100 * spread(to), ratio, name);
        all.push_back({"FLAC", ratio});
        return 0;
    }
    if (s.codec == AV_CODEC_ID_MP3) {
        Mp3Check c;
        if (!mp3Check(s, c, why)) { std::printf("NOT COUNTED %s: %s\n", name, why.c_str()); return 1; }
        std::vector<double> tf, tq, tx;
        for (int r = 0; r < rounds; ++r) {
            const bool ffFirst = r % 2 == 0;
            if (ffFirst) tf.push_back(ffmpegRound(s.packets, s.par, AV_CODEC_ID_MP3));
            tq.push_back(mp3Round(s, false, nullptr));
            tx.push_back(mp3Round(s, true, nullptr));
            if (!ffFirst) tf.push_back(ffmpegRound(s.packets, s.par, AV_CODEC_ID_MP3));
        }
        const double rq = median(tf) / median(tq), rx = median(tf) / median(tx);
        const bool ok = c.skipped == 0 && c.badData == 0 && c.ffSamples == c.ourSamples && c.gapFast <= kMp3GapLimit
                        && c.gapExact <= kMp3GapLimit && good(tf) && good(tq) && good(tx);
        std::printf("MP3 %4lld kbps %5.1f kHz %7.1f s  ffmpeg %.4f s ±%2.0f%%  fast path %.4f s %5.2f times  exact path %.4f s %5.2f times  "
                    "largest gap from ffmpeg %.2g (fast) %.2g (exact)%s  %s\n",
                    (long long)(s.par->bit_rate / 1000), s.par->sample_rate / 1000.0, s.seconds, median(tf), 100 * spread(tf), median(tq),
                    rq, median(tx), rx, c.gapFast, c.gapExact, ok ? "" : "  NOT COUNTED: ", name);
        if (!ok) {
            std::printf("    samples per channel: ffmpeg %llu, ours %llu; frames ours did not decode: %ld; frames against the standard: "
                        "%ld; largest gap allowed %.0e\n",
                        (unsigned long long)c.ffSamples, (unsigned long long)c.ourSamples, c.skipped, c.badData, kMp3GapLimit);
            return 1;
        }
        all.push_back({"MP3 fast path", rq});
        all.push_back({"MP3 exact path", rx});
        return 0;
    }
    if (s.codec == AV_CODEC_ID_DST) {
        const int ch = s.par->ch_layout.nb_channels;
        // Το libavformat δίνει για ένα .dff τον ρυθμό σε BYTES DSD ανά δευτερόλεπτο και κανάλι (μετρημένο: 352.800 στο DSD64,
        // δηλαδή 2.822.400 bit / 8). Ένα πλαίσιο = 1/75 δευτ. Αν κάποτε άλλαζε, ο αποκωδικοποιητής μας θα αρνιόταν τα πλαίσια
        // και ο έλεγχος πιο κάτω θα έβγαζε το αρχείο από τη σύνοψη.
        const int rate = s.par->sample_rate;
        const int perChannel = rate / 75;
        if (ch < 1 || ch > 6 || rate % 352800 != 0 || perChannel <= 0 || perChannel > kDstMaxBytesPerChannel) {
            std::printf("SKIPPED %s: unexpected shape (%d channels, %d Hz)\n", name, ch, rate);
            return 1;
        }
        // Έλεγχος και ζέσταμα: κανένα πλαίσιο αρνημένο, κανένα λάθος του FFmpeg. Κρατά το DSD που βγάζει ο δικός μας (ως 1 GB)
        // για τη μέτρηση της μετατροπής του FFmpeg.
        const size_t frameBytes = size_t(ch) * size_t(perChannel);
        const size_t keep = std::max<size_t>(1, (size_t(1) << 30) / frameBytes);
        long refused = 0, errors = 0;
        std::vector<AVPacket *> dsd;
        dstRound(s, ch, perChannel, &refused, &dsd, keep);
        uint64_t dstSamples = 0;
        ffmpegRound(s.packets, s.par, AV_CODEC_ID_DST, &errors, &dstSamples);
        // Και το μέγεθος του πλαισίου δεμένο με το FFmpeg: ο αποκωδικοποιητής του βγάζει ένα δείγμα PCM ανά 8 bit DSD κάθε
        // καναλιού, άρα (πλαίσια × bytes ανά κανάλι) δείγματα. Ένα λάθος μέγεθος ο δικός μας ΔΕΝ το αρνείται πάντα (δοκιμάστηκε):
        // αυτός ο έλεγχος το πιάνει.
        const uint64_t expect = uint64_t(s.packets.size()) * uint64_t(perChannel);
        const auto release = [&]() { for (AVPacket *p : dsd) av_packet_free(&p); };
        if (refused || errors || dstSamples != expect) {
            std::printf("NOT COUNTED %s: %ld frames refused by ours, %ld errors from FFmpeg, %llu samples from FFmpeg where %llu "
                        "were expected\n", name, refused, errors, (unsigned long long)dstSamples, (unsigned long long)expect);
            release();
            return 1;
        }
        // Η μετατροπή DSD σε PCM που κάνει ο αποκωδικοποιητής DST του FFmpeg, μετρημένη χωριστά: ο αποκωδικοποιητής DSD του FFmpeg
        // (dsd_msbf, μόνο μετατροπή) πάνω στο ΙΔΙΟ DSD που έβγαλε ο δικός μας (όλο, ως 1 GB· πέρα από αυτό ανακυκλώνεται). Είναι
        // εκτίμηση: άλλος αποκωδικοποιητής με την ίδια μετατροπή. Ένα ζέσταμα πριν από τους γύρους.
        std::vector<AVPacket *> conv;
        for (size_t i = 0; i < s.packets.size(); ++i) conv.push_back(dsd[i % dsd.size()]);
        uint64_t pcmSamples = 0;
        ffmpegRound(conv, s.par, AV_CODEC_ID_DSD_MSBF, nullptr, &pcmSamples);
        std::vector<double> tf, to, tc;
        for (int r = 0; r < rounds; ++r) {
            const bool ffFirst = r % 2 == 0;
            if (ffFirst) tf.push_back(ffmpegRound(s.packets, s.par, AV_CODEC_ID_DST));
            to.push_back(dstRound(s, ch, perChannel, nullptr));
            tc.push_back(ffmpegRound(conv, s.par, AV_CODEC_ID_DSD_MSBF));
            if (!ffFirst) tf.push_back(ffmpegRound(s.packets, s.par, AV_CODEC_ID_DST));
        }
        release();
        if (!good(tf) || !good(to) || !good(tc)) { std::printf("NOT COUNTED %s: a timed round failed\n", name); return 1; }
        const double with = median(tf) / median(to);
        const bool convOk = pcmSamples == dstSamples && median(tc) < median(tf);
        const double sec = s.seconds > 0 ? s.seconds : double(s.packets.size()) / 75;
        const std::string kind = ch > 2 ? "DST multichannel" : "DST stereo";
        if (!convOk) {
            std::printf("DST %d ch %7.1f s  ffmpeg %.4f s  ours %.4f s  %5.2f times with FFmpeg's DSD to PCM (the conversion alone could "
                        "not be measured, so not in the summary)  %s\n", ch, sec, median(tf), median(to), with, name);
            return 1;
        }
        // ΠΡΩΤΟΣ ο λόγος ΧΩΡΙΣ τη μετατροπή του FFmpeg (ίδια δουλειά και στις δύο πλευρές)· δεύτερος ο λόγος με αυτήν.
        const double without = (median(tf) - median(tc)) / median(to);
        std::printf("DST %d ch %7.1f s  ffmpeg %.4f s ±%2.0f%%, of it DSD to PCM about %.4f s  ours %.4f s ±%2.0f%%  %5.2f times "
                    "(%.2f with FFmpeg's DSD to PCM)  %s\n",
                    ch, sec, median(tf), 100 * spread(tf), median(tc), median(to), 100 * spread(to), without, with, name);
        all.push_back({kind, without});
        all.push_back({kind + ", incl. FFmpeg's DSD to PCM", with});
        return 0;
    }
    std::printf("SKIPPED %s: not MP3, FLAC or DST\n", name);
    return 1;
}

}  // namespace

int main(int argc, char **argv)
{
    int rounds = 5, first = 1;
    if (argc > 2 && std::strcmp(argv[1], "--rounds") == 0) {
        rounds = std::max(1, std::atoi(argv[2]));
        first = 3;
    }
    if (first >= argc) {
        std::fprintf(stderr, "usage: speedtest [--rounds N] file.mp3|file.flac|file.dff ...\n");
        return 2;
    }
    av_log_set_level(AV_LOG_ERROR);
    std::printf("FFmpeg %s (libavcodec %d.%d.%d) · %d rounds, median, alternating order · decoder alone, one thread, the same packets for both%s\n",
                av_version_info(), LIBAVCODEC_VERSION_MAJOR, LIBAVCODEC_VERSION_MINOR, LIBAVCODEC_VERSION_MICRO, rounds,
#ifdef NDEBUG
                ""
#else
                " · WARNING: this speedtest was built without optimisation (not a Release build); its numbers are not meaningful"
#endif
    );
    std::vector<Result> all;
    int failed = 0;
    for (int a = first; a < argc; ++a) {
        try {
            failed += file(argv[a], rounds, all);
        } catch (const std::bad_alloc &) {
            std::printf("SKIPPED %s: not enough memory (speedtest holds the whole compressed file in memory)\n", argv[a]);
            ++failed;
        }
    }
    std::map<std::string, std::vector<double>> by;
    for (const Result &r : all) by[r.kind].push_back(r.ratio);
    std::printf("\ntimes faster than FFmpeg (median · lowest · highest · files counted)\n");
    for (const auto &k : by) {
        const auto mm = std::minmax_element(k.second.begin(), k.second.end());
        std::printf("  %-44s %5.2f · %5.2f · %5.2f · %zu\n", k.first.c_str(), median(k.second), *mm.first, *mm.second, k.second.size());
    }
    if (failed) std::printf("%d file(s) not counted (see above)\n", failed);
    return failed ? 1 : 0;
}
