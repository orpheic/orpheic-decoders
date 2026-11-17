#pragma once
// Ο ΑΠΟΣΥΜΠΙΕΣΤΗΣ DST (Direct Stream Transfer): κάθε συμπιεσμένο πλαίσιο μιας περιοχής DST ενός SACD σε DSD, ένα πλαίσιο τη
// φορά. Είναι ο γρήγορος: γράφτηκε για να τρέχει γρηγορότερα σε ΕΝΑ νήμα (σε Core 2 Duo ο προηγούμενος δεν προλάβαινε τη
// μουσική).
//
// DST (Direct Stream Transfer) decoder.
//
// Written from the specification ISO/IEC 14496-3:2009 Subpart 10 ("Technical
// description of lossless coding of oversampled audio"), clauses 10.5–10.7, and
// from our own previous decoder (June 2026),
// whose header parsing it copies. The previous decoder's decoded OUTPUT was
// compared byte-for-byte with the output of the MPEG-4 reference decoder; this
// one's output is compared with the previous decoder's and with FFmpeg's `dst`
// output (see below).
//
// DST = per-channel FIR sign prediction whose error is binary-arithmetic coded,
// driven by probability tables selected by the prediction value. Stereo + MCH.
//
// ΑΠΟ ΤΙ ΓΡΑΦΤΗΚΕ Ο ΓΡΗΓΟΡΟΣ (12/11/2026): μόνο από τον προηγούμενο δικό μας αποσυμπιεστή (η ανάγνωση της
// επικεφαλίδας είναι αντίγραφο εκείνου), από όσα ορίζει
// το ίδιο το πρότυπο (όπως ήταν ήδη υλοποιημένο εκεί) και από γενικές τεχνικές βελτιστοποίησης. Στη δουλειά αυτή δεν
// ανοίχτηκε πηγαίος κώδικας άλλου αποσυμπιεστή DST· άλλες υλοποιήσεις (το `dst` του FFmpeg) χρησιμοποιήθηκαν μόνο ως
// ΕΞΟΔΟΣ, για σύγκριση.
//
// ΠΩΣ ΕΙΝΑΙ ΓΡΗΓΟΡΟΣ (λεπτομέρειες και το επιχείρημα ισοδυναμίας στο .cpp): προβλέπει και τα οκτώ μπιτ ενός byte μαζί, από
// πίνακες με δείκτη ένα τετράμπιτο του ιστορικού, σε οκτώ «λωρίδες» των 16 μπιτ ανά θέση (8 KB ανά φίλτρο, μία πρόσθεση
// 128 μπιτ ανά γραμμή: SSE2, NEON ή δύο λέξεις των 64), με αριθμητικό αποκωδικοποιητή πάνω σε παράθυρο 64 μπιτ της εισόδου.
// Πλαίσιο με φίλτρο που το Σ|c| του περνά το 32767 (υγιής δίσκος δεν έχει) πηγαίνει από αργό, απλό δρόμο.
//
// ΤΙ ΒΓΑΖΕΙ. Σε κάθε πλαίσιο που δέχεται: τα ΙΔΙΑ bytes με τον προηγούμενο δικό μας αποσυμπιεστή,
// μπιτ προς μπιτ. Διαφέρει από εκείνον μόνο σε τρία είδη ΧΑΛΑΣΜΕΝΟΥ πλαισίου, που εκείνος «αποκωδικοποιούσε» σε θόρυβο και
// αυτός αρνείται (DstRefusal: ShortPlain, CutHeader, CodeStart), και σε πλαίσιο φτιαγμένο με το χέρι όπου ο προηγούμενος
// υπερχειλίζει έναν int. Και στα δύο ο προηγούμενος κάνει κάτι που η γλώσσα δεν ορίζει. Το μέτρο σύγκρισης σε ΚΑΘΕ είσοδο
// είναι ο διαιτητής (ένας απλός αποσυμπιεστής αναφοράς γραμμένος για τον έλεγχο): ίδια απάντηση, ίδιος λόγος άρνησης, ίδια
// bytes. Η απόδειξη: η σύγκριση σε ολόκληρους δίσκους, σε χαλασμένα και σε φτιαχτά πλαίσια, σε κάθε μεταγλωττιστή και
// μηχάνημα.

#include <cstdint>
#include <vector>

// Γιατί αρνήθηκε ένα πλαίσιο. Τα τρία τελευταία είναι οι αυστηρότεροι κανόνες: ο προηγούμενος εκεί δέχεται.
enum class DstRefusal : uint8_t {
  None,         // δεν αρνήθηκε
  Shape,        // πλήθος καναλιών ή μήκος πλαισίου που ο αποσυμπιεστής δεν δέχεται (δεν έρχεται από τον δίσκο)
  Header,       // πεδίο της επικεφαλίδας που κανένας κωδικοποιητής δεν γράφει (οι έλεγχοι που είχε και ο προηγούμενος)
  ShortPlain,   // απλό (ασυμπίεστο) πλαίσιο με λιγότερα bytes από όσα έχει ένα πλαίσιο
  CutHeader,    // το πλαίσιο τελείωσε μέσα στην επικεφαλίδα του
  CodeStart,    // ο αριθμητικός κώδικας αρχίζει εκτός διαστήματος (τα πρώτα 13 μπιτ >= 4095)
};
constexpr int kDstRefusals = 6;
// Πάνω από αυτό το μήκος (bytes ανά κανάλι) άρνηση: το DSD512 θέλει 37.632, και από 2^28 η αριθμητική σε int δεν χωρά.
constexpr int kDstMaxBytesPerChannel = 1 << 20;

// Μία θέση των πινάκων πρόβλεψης: οκτώ «λωρίδες» των 16 μπιτ σε δύο λέξεις (δες το .cpp). Στοιχισμένη στα 16 bytes, γιατί
// διαβάζεται ολόκληρη με μία εντολή 128 μπιτ.
struct DstLanes {
  alignas(16) uint64_t w[2];
};

// ΜΟΝΟ ΓΙΑ ΤΗ ΜΕΤΡΗΣΗ: τι είδους πλαίσια πέρασαν. Τα «bytes» είναι bytes DSD ενός καναλιού (ή ενός διαστήματος όπου
// κανένα κανάλι δεν αλλάζει φίλτρο), ώστε τα ποσοστά να λένε πόση από τη ΔΟΥΛΕΙΑ είναι κάθε είδος, όχι πόσα πλαίσια.
// ΧΩΡΙΣ #ifdef, ΕΠΙΤΗΔΕΣ, μαζί με το countInto και το noteFrame του .cpp: ο κώδικας που ελέγχουν τα εργαλεία
// πρέπει να είναι αυτός που στέλνεται. Δεν έχουν κείμενα, και όποιος δεν δίνει χώρο για μετρητές
// πληρώνει έναν έλεγχο δείκτη ανά πλαίσιο.
struct DstStats {
  long long refused[kDstRefusals] = {};   // ανά DstRefusal· το [None] είναι όσα δέχτηκε
  long long plain = 0;                    // από όσα δέχτηκε, τα απλά (ασυμπίεστα)
  long long coded = 0;                    // από όσα δέχτηκε, τα συμπιεσμένα
  long long multiSeg = 0;                 // συμπιεσμένα όπου κάποιο κανάλι έχει πάνω από ένα τμήμα (φίλτρων ή πινάκων)
  long long halfProb = 0;                 // συμπιεσμένα με Half_Probability σε κάποιο κανάλι
  long long narrow = 0;                   // συμπιεσμένα όπου κανένα φίλτρο σε χρήση δεν έχει τάξη πάνω από 64
  long long slow = 0;                     // συμπιεσμένα που πήγαν από τον αργό δρόμο (φίλτρο με Σ|c| > 32767)
  long long headCut = 0;                  // συμπιεσμένα όπου κάποιο τμήμα αρχίζει μέσα στα πρώτα 16 bytes (bytes 1..15)
  long long maxFilters = 0;               // τα περισσότερα φίλτρα που δήλωσε ένα πλαίσιο
  long long stretches = 0;                // διαστήματα σε όλα τα συμπιεσμένα
  long long stretchBytes = 0;             // bytes διαστημάτων (ένα ανά byte του πλαισίου, όχι ανά κανάλι)
  long long orderBytes[8] = {};           // bytes καναλιού με φίλτρο τάξης 1-16, 17-32, ..., 113-128
  long long distinctBytes[7] = {};        // bytes διαστήματος με 0..6 ΔΙΑΦΟΡΕΤΙΚΑ φίλτρα σε χρήση
  long long longBytes[7] = {};            // bytes διαστήματος με 0..6 διαφορετικά φίλτρα τάξης πάνω από 64 σε χρήση
  void add(const DstStats& o) {
    for (int i = 0; i < kDstRefusals; ++i) refused[i] += o.refused[i];
    plain += o.plain; coded += o.coded; multiSeg += o.multiSeg; halfProb += o.halfProb; narrow += o.narrow; slow += o.slow;
    headCut += o.headCut;
    if (o.maxFilters > maxFilters) maxFilters = o.maxFilters;
    stretches += o.stretches; stretchBytes += o.stretchBytes;
    for (int i = 0; i < 8; ++i) orderBytes[i] += o.orderBytes[i];
    for (int i = 0; i < 7; ++i) { distinctBytes[i] += o.distinctBytes[i]; longBytes[i] += o.longBytes[i]; }
  }
};

class DstDecoder {
 public:
  bool init(int channels) { channels_ = channels; return channels >= 1 && channels <= 6; }
  // Ένα πλαίσιο DST σε channels * bytesPerChannel bytes πλεγμένου DSD (MSB πρώτο). false = άρνηση· τότε ο καλών δεν
  // χρησιμοποιεί την έξοδο και αποφασίζει ο ίδιος τι θα βάλει στη θέση της. (Προσοχή: bytes 0x00 στο DSD δεν είναι σιωπή
  // αλλά συνεχής τάση σε πλήρη κλίμακα· σιωπή είναι το 0x69.)
  bool decodeFrame(const uint8_t* in, int inBytes, uint8_t* dsdOut, int bytesPerChannel);
  // Γιατί αρνήθηκε το τελευταίο πλαίσιο (None αν το δέχτηκε).
  DstRefusal lastRefusal() const { return why_; }
  // Μόνο για τη μέτρηση: από εδώ και πέρα κάθε πλαίσιο μετριέται στο `s` (nullptr = όχι). Ο χώρος είναι του καλούντος.
  void countInto(DstStats* s) { stats_ = s; }

 private:
  int channels_ = 2;
  DstRefusal why_ = DstRefusal::None;
  DstStats* stats_ = nullptr;
  // Πρόχειρος χώρος που μεγαλώνει μόνο (ποτέ ανά πλαίσιο όταν έχει φτάσει το μέγεθος που θέλει το πλαίσιο):
  std::vector<int> coef_;       // kMaxElem * 128, οι συντελεστές κάθε φίλτρου
  std::vector<int> pone_;       // kMaxElem * 64, οι πίνακες πιθανοτήτων όπως διαβάζονται
  std::vector<uint8_t> p8_;     // οι ίδιοι σε bytes (τιμές 1..128), για την κρίσιμη διαδρομή
  std::vector<DstLanes> tab_;   // 32 γραμμές των 16 θέσεων ανά φίλτρο (8 KB): οι πίνακες πρόβλεψης
};
