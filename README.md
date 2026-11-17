# orpheic-decoders

Three audio decoders from [Orpheic](https://orpheic.com), a music player for
the music you own: MP3, FLAC and DST (the lossless compression inside SACDs).
C++17, no dependencies, MIT licence. By Georgios Kakavelakis.

They are the decoders the player itself uses. They were written to be exact
first and fast second, and each claim below was measured.

| Decoder | Standard | Scope |
|---|---|---|
| MP3 | ISO/IEC 11172-3 and 13818-3, Layer III | 32, 44.1 and 48 kHz and the low sampling rates (16, 22.05, 24 kHz). Not Layer I or II, not the unofficial 8 to 12 kHz extension. |
| FLAC | RFC 9639 | 4 to 32 bits, 1 to 8 channels, every block size. Samples come out as 32-bit integers. |
| DST | ISO/IEC 14496-3:2009, Subpart 10 | Stereo and up to six channels. One compressed frame in, DSD bytes out. |

## What was measured

How many times faster each decoder runs than FFmpeg's decoder for the same
format. Decoder alone, one thread, the same packets for both, the median over
the files. Each number is against the FFmpeg that most people on that system
have.

| | Desktop PC (i9-14900K), FFmpeg 9.0.2 | MacBook (Apple M4), FFmpeg 9.0.2 from Homebrew | Raspberry Pi 5, FFmpeg 7.1.5 from the system |
|---|---|---|---|
| FLAC | 2.37 times (2.12 to 2.47) | 3.45 times (2.93 to 3.69) | 1.58 times (1.27 to 1.71) |
| MP3, fast path | 2.25 times (2.18 to 2.57) | 2.76 times (2.56 to 2.92) | 2.27 times (2.21 to 2.38) |
| DST, stereo | 2.01 times | 2.06 times | 2.89 times |
| DST, five or six channels | 1.81 times | 1.60 times | 2.74 times |

What a careful reader should know about this table:

- **FLAC and MP3** were measured on eight freely licensed recordings (8 FLAC
  files, and 32 MP3 files made from them). `tools/bench.sh` downloads the same
  recordings and repeats the measurement on your machine (see *Check it
  yourself* below).
- **MP3 has two paths, and the number is for the fast one.** The exact path,
  which Orpheic uses to play music, is written for accuracy and is several
  times slower than FFmpeg (0.14 to 0.27 times its speed on these machines).
  It still runs hundreds of times faster than the music.
- **DST** was measured on five excerpts of our own (two stereo, three with five
  or six channels) that we cannot give out, so this part cannot be repeated
  with our files, only with yours. FFmpeg's DST decoder also turns DSD into
  PCM and ours does not, so the numbers above leave that step out (it is an
  estimate, see below). With it left in they are 2.08 and 1.88 on the PC, 2.16
  and 1.70 on the MacBook, 3.01 and 2.85 on the Raspberry Pi 5. An earlier
  measurement on twelve discs, with that step in, gave 2.05 and 1.84 on the PC.
- **On the Raspberry Pi 5, FLAC gains least**, and 24-bit files least of all
  (1.27 to 1.45 times).
- **Against the FFmpeg that Ubuntu 24.04 ships (6.1.1)** the PC numbers are
  2.25 (FLAC), 2.38 (MP3), 1.89 and 1.79 (DST).
- **A 2007 laptop** (Core 2 Duo, Windows 10) was measured earlier with a
  different tool and other files, against the FFmpeg 6.1.2 that ships with
  Orpheic: MP3 2.45, FLAC 2.26, DST 2.09 and 1.86 times. Those are not
  comparable one to one with the table.
- FFmpeg 9.0.2 on the PC was built from the signed release with
  `./configure --enable-shared`, x86 assembly on.

- FLAC and DST output is identical to FFmpeg's, sample for sample (FLAC) and
  byte for byte (DST), on every file and disc measured. You can repeat this
  comparison yourself (see *Check it yourself* below).
- MP3 has no bit-exact definition: the standard is arithmetic on real numbers.
  The fast path is checked against a second path in this repository that
  executes the formulas of the standard literally in double precision
  (`mp3exact.h`). On the 32 MP3 files that `tools/bench.sh` makes, the largest
  gap was 0.42 millionths of full scale: about three and a half steps of 24-bit
  output, and 73 times smaller than one step of 16-bit output (the same on the
  PC and on the Raspberry Pi 5). `mp3check` repeats this on those files, or on
  any files you give it.

All the numbers, with charts, and how they were measured:
**[orpheic.com/decoders](https://orpheic.com/decoders/)**

## Build

```
cmake -B build
cmake --build build
```

This builds three static libraries, `orpheic_mp3`, `orpheic_flac` and
`orpheic_dst` (take one, two or all three; `orpheic_decoders` links them all),
and the small tools in `tools/`. Built and run with GCC 13 and Clang 18 on x86-64 Linux, GCC 14 on
ARM64 Linux (Raspberry Pi 5), Apple Clang 21 on ARM64 macOS and MSVC 19.44 on Windows.

The tools land in `build/`. With MSVC they land in `build/Debug`, and that build
is not optimised: for speed, build with `cmake --build build --config Release`
(they then land in `build/Release`). The
commands below name them without that path, so run them as `build/flachash`,
`build/speedtest` and so on, or add that folder to your `PATH`.

## The same bits everywhere

All three decoders give the same output, bit for bit, on every system and
compiler above. Measured on 15 November 2026 by hashing every sample:

| | Files | Linux GCC (x86-64) | Linux Clang (x86-64) | macOS Apple Clang (ARM64) | Raspberry Pi 5 GCC (ARM64) | Windows MSVC (x86-64) |
|---|---|---|---|---|---|---|
| MP3, exact path | 33 | reference | same | same | same | same |
| MP3, fast path | 33 | reference | same | same | same | same |
| FLAC | 11 | reference | same | same | same | same |
| DST, stereo and multichannel | 5 | reference | same | same | same | same |

FLAC and DST are integer arithmetic, so this holds by construction. MP3 is not:
it needs cosines, sines and powers, and each system's maths library answers
them a little differently (89 of the 13,288 values the exact path uses differed
in the last bit between Linux and Windows). So the MP3 decoder never asks the
system. Every such value is written in `mp3consts_tables.inc`, computed once to
60 digits and rounded correctly to the nearest double: the same on every system,
and the most accurate a double can hold. And the MP3 files are compiled with
`-ffp-contract=off`, so that no compiler fuses a multiply and an add on its own.

You can check this on your own machine. Each tool prints one SHA-256 line per
file; the same file must give the same line everywhere:

```
mp3hash --exact file.mp3 ...    # the exact path (without --exact: the fast path)
flachash file.flac ...          # also says whether the samples match the MD5 inside the file
dsthash file.dff ...            # a DSDIFF file with DST compression
```

## Check it yourself

**FLAC against FFmpeg.** The line `flachash` prints starts with the SHA-256 of
the samples as raw little-endian PCM, so for a 16-bit or 24-bit file these two
give the same hash:

```
flachash file.flac
ffmpeg -v error -i file.flac -f s16le - | sha256sum     # s24le for a 24-bit file
```

**DST against FFmpeg.** FFmpeg turns DST into PCM rather than DSD, so the
comparison goes through the same conversion on both sides: `dsthash` writes
what it decoded as an uncompressed DSDIFF file, and FFmpeg converts both files.

```
dsthash --dsd-out ours.dff file.dff
ffmpeg -v error -i file.dff -f f32le - | sha256sum
ffmpeg -v error -i ours.dff -f f32le - | sha256sum
```

(`sha256sum` is the Linux command. On macOS use `shasum -a 256`. On Windows,
write each output to a file instead of a pipe (`-f s16le out.pcm` for FLAC;
`-f f32le a.pcm` and `-f f32le b.pcm` for the two DST files) and hash it with
`certutil -hashfile out.pcm SHA256`.)

The two hashes are the same when the DSD bytes are the same (in our test, one
bit changed in the middle of the DSD changed the hash). Compare only a file that
`dsthash` reports as `refused 0` and without `FILE-CUT-SHORT`: a frame the
decoder refuses is written as DSD silence, and the comparison would then measure
the damage, not the decoder. `dsthash` reads DSDIFF files (`.dff`)
with DST compression, the usual file form of DST outside a disc. It also reads a
plain dump of frames (`.dstf`): the bytes `DSTF1` and a newline, the channel
count (1 byte), the DSD rate in Hz (4 bytes) and the number of frames (4 bytes),
then for each frame its length (4 bytes) and its bytes; numbers little-endian.

**MP3.** There is no bit-exact answer to compare with (see above): `mp3check`
compares the fast path with the exact one on any files you give it.

**Speed, on your machine.** `speedtest` times both decoders on the same data:
libavformat reads each file once, and the same packets go to FFmpeg's decoder
and to ours. Only decoding is timed, up to samples in FFmpeg's own output
format (our FLAC decoder writes its samples in the same format FFmpeg returns),
one thread each, in rounds of alternating order; the median of five is kept.
Before timing it checks every FLAC frame byte for byte against FFmpeg's, checks
that the MP3 decoders return as many samples as FFmpeg, decode every frame
without finding it against the standard and stay within one ten-thousandth of
full scale of it, and for DST that no frame was refused and that the frame size
agrees with FFmpeg's output. A file that fails a check, or a timing round that
fails, is reported and left out of the summary. Our FLAC decoder checks the CRC
of every frame; FFmpeg's, with its default settings, does not (measured: a frame
whose CRC we damaged on purpose came out of FFmpeg 6.1.1 and 9.0.2 with no
error, and ours refused it), so that part of the work is ours alone. FFmpeg's
DST decoder also turns DSD into PCM, which ours does not do: `speedtest` times
that conversion separately (FFmpeg's DSD decoder on the same DSD our decoder
produced), subtracts it and prints that ratio first, as an estimate, and the
ratio including the conversion after it.
It needs FFmpeg 5.1 or newer with its development files, and holds the
compressed file in memory.

```
tools/bench.sh                                # the free set below, against this machine's FFmpeg
speedtest file.flac file.mp3 file.dff ...     # or your own files
```

`tools/bench.sh` downloads eight freely licensed recordings from Wikimedia
Commons (listed with their licences in `tools/bench-corpus.tsv`: public domain,
CC0, CC BY and CC BY-SA; orchestra, string quartet, piano, baroque, rock,
post-rock, jazz; three of them 24-bit), checks their SHA-1, makes 32 MP3 files
from them with ffmpeg (128, 192 and 320 kbps and VBR `-q:a 2`), and runs
`speedtest`. The MP3 figure is for the fast path; `speedtest` times the exact
path too, which is the witness, written for accuracy, and is several times
slower than FFmpeg. No freely licensed DST recording is known to us, so DST is
measured on your own `.dff` files.

`tools/bench.sh` needs `curl`, the `ffmpeg` command with the libmp3lame encoder,
`sha1sum` or `shasum`, and `speedtest` built in `build/` (another build folder:
`BUILD=... tools/bench.sh`). The downloads stay in `bench-corpus/` for the next
run: about 260 MB of music, about 420 MB with the MP3 files made from it.

## Use

Each decoder takes one whole frame at a time and does no file reading of its
own; finding frames in a file is the caller's job (the tools show how).

- MP3: `orpheic::mp3::Decoder` in `mp3decoder.h`. `decode(frame, bytes)` leaves
  1152 or 576 float samples per channel in `channel(c)`. A new decoder runs the
  exact path; call `setExact(false)` for the fast one.
- FLAC: `orpheic::flac::Decoder` in `flacdecoder.h`.
- DST: `DstDecoder` in `dstdecoder.h`. `decodeFrame` writes interleaved DSD.

Damaged input gives an error status or silence, never a crash: nothing a file
says is believed before it is checked.

The comments in the code are in Greek, my language. If there is interest, they
will be translated into English in the future.

## How it was made

The code was written by Claude, an AI model from Anthropic, under my direction.
I set what to build and the rules it had to follow, and nothing was accepted on
its word: every claim here was measured against the standards, against a second
exact implementation, or against other decoders run as black boxes. The tools in
`tools/` repeat the main checks: `mp3check` (fast path against exact path),
`flachash` (against the MD5 each FLAC file carries), and `mp3hash`, `flachash`
and `dsthash` (the same bits on every system).

## Where the code comes from

The MP3 and FLAC decoders were written from their standards and checked against
other decoders that were only run, never read: their output and their speed were
measured from outside.

The DST decoder is the second one written for the player. The first was written
from the standard, and its output was compared byte for byte with the output of
the standard's reference decoder. The present one keeps the first one's header
parsing and rewrites the rest for speed. For that work no other DST decoder was
read, and its output is compared with the first decoder's and with FFmpeg's.

"Written from the standard" cannot be proved for code a model wrote: nobody can
show what a model remembers from its training. So the code was also compared
with other decoders for the same formats. The reviewers were AI models too:
separate ones, which had not written the code and were asked to look for
copying. The overlap in
text was at the level of unrelated code, and what is shared is what the
standards themselves print (tables, pseudo-code). The one part of the MP3
decoder they found organised like another decoder was rewritten.

## Licence

MIT, copyright Georgios Kakavelakis. See `LICENSE`.

Orpheic itself is at [orpheic.com](https://orpheic.com).
