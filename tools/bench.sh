#!/bin/sh
# bench.sh: η μέτρηση ταχύτητας σε μουσική που έχει ο καθένας. Κατεβάζει τα οκτώ κομμάτια του tools/bench-corpus.tsv (Wikimedia
# Commons· δημόσιο κτήμα, CC0, CC BY, CC BY-SA), ελέγχει το SHA-1 του καθενός, φτιάχνει από αυτά MP3 με το ffmpeg (128, 192 και
# 320 kbps και VBR -q:a 2, στα 44,1 kHz) και τρέχει το speedtest απέναντι στο FFmpeg αυτού του μηχανήματος.
#
#   tools/bench.sh [φάκελος]          (προεπιλογή: bench-corpus· τα κατεβασμένα μένουν εκεί για το επόμενο τρέξιμο)
#
# Θέλει: curl, ffmpeg με libmp3lame, sha1sum ή shasum, και το speedtest χτισμένο (cmake -B build && cmake --build build, με τα
# αρχεία ανάπτυξης του FFmpeg). Άλλος φάκελος χτισίματος: BUILD=… tools/bench.sh
set -eu

here=$(cd "$(dirname "$0")" && pwd)
build=${BUILD:-$here/../build}
dir=${1:-bench-corpus}
speedtest=$build/speedtest
[ -x "$speedtest" ] || { echo "no $speedtest: build first (cmake -B build && cmake --build build), with FFmpeg's development files installed" >&2; exit 2; }
command -v curl >/dev/null || { echo "curl is needed" >&2; exit 2; }
command -v ffmpeg >/dev/null || { echo "ffmpeg is needed (to make the MP3 files)" >&2; exit 2; }
ffmpeg -hide_banner -encoders 2>/dev/null | grep -q libmp3lame || { echo "this ffmpeg has no libmp3lame encoder" >&2; exit 2; }
if command -v sha1sum >/dev/null; then sha1() { sha1sum "$1" | cut -d' ' -f1; }
elif command -v shasum >/dev/null; then sha1() { shasum -a 1 "$1" | cut -d' ' -f1; }
else echo "sha1sum or shasum is needed" >&2; exit 2; fi

list=$here/bench-corpus.tsv
[ -f "$list" ] || { echo "no $list" >&2; exit 2; }
mkdir -p "$dir/flac" "$dir/mp3"
echo "1/3 the music (8 pieces from Wikimedia Commons, about 260 MB the first time; about 420 MB with the MP3 files)"
tab=$(printf '\t')
n=0
# Ανάγνωση από αρχείο, όχι από σωλήνα: ένα «exit» μέσα στον βρόχο σταματά όλο το script.
while IFS="$tab" read -r file url want what by licence page; do
    case $file in ''|'#'*) continue ;; esac
    f="$dir/flac/$file"
    if [ ! -f "$f" ] || [ "$(sha1 "$f")" != "$want" ]; then
        echo "  downloading $file  ($what · $by · $licence)"
        # Η Wikimedia ζητά περιγραφικό αναγνωριστικό και απαντά 429 σε γρήγορες διαδοχικές λήψεις: παύση και επαναλήψεις.
        if ! curl -fsSL --retry 6 --retry-delay 10 --retry-all-errors \
                -A "orpheic-decoders-bench/1.0 (https://github.com/orpheic/orpheic-decoders)" -o "$f.part" "$url"; then
            rm -f "$f.part"
            echo "  $file: download failed ($url)" >&2
            exit 1
        fi
        mv "$f.part" "$f"
        sleep 2
    fi
    if [ "$(sha1 "$f")" != "$want" ]; then
        echo "  $file: SHA-1 does not match the list (the file at $page may have been replaced); not used" >&2
        rm -f "$f"
        exit 1
    fi
    n=$((n + 1))
done < "$list"
[ "$n" -gt 0 ] || { echo "the list $list names no files" >&2; exit 2; }

echo "2/3 the MP3 files (ffmpeg $(ffmpeg -version | head -1 | cut -d' ' -f3), libmp3lame)"
for f in "$dir"/flac/*.flac; do
    b=$(basename "$f" .flac)
    for v in cbr128 cbr192 cbr320 vbr2; do
        o="$dir/mp3/$b-$v.mp3"
        [ -s "$o" ] && continue
        case $v in
            cbr128) opt="-b:a 128k" ;;
            cbr192) opt="-b:a 192k" ;;
            cbr320) opt="-b:a 320k" ;;
            vbr2)   opt="-q:a 2" ;;
        esac
        # shellcheck disable=SC2086
        ffmpeg -nostdin -v error -y -i "$f" -map 0:a -map_metadata -1 -ar 44100 -c:a libmp3lame $opt "$o.part.mp3"
        mv "$o.part.mp3" "$o"
    done
done

echo "3/3 the measurement (decoder alone, one thread, the same packets for both, median of 5 rounds)"
"$speedtest" "$dir"/flac/*.flac "$dir"/mp3/*.mp3
