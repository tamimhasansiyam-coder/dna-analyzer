import os
import re
from kivy.app import App
from kivy.core.clipboard import Clipboard
from kivy.core.window import Window
from kivy.metrics import dp
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.button import Button
from kivy.uix.label import Label
from kivy.uix.popup import Popup
from kivy.uix.scrollview import ScrollView
from kivy.uix.switch import Switch
from kivy.uix.textinput import TextInput
from kivy.utils import platform

Window.clearcolor = (0.95, 0.97, 0.99, 1)

CODON_TABLE = {
    'UUU': ('F', 'Phe'), 'UUC': ('F', 'Phe'), 'UUA': ('L', 'Leu'), 'UUG': ('L', 'Leu'),
    'CUU': ('L', 'Leu'), 'CUC': ('L', 'Leu'), 'CUA': ('L', 'Leu'), 'CUG': ('L', 'Leu'),
    'AUU': ('I', 'Ile'), 'AUC': ('I', 'Ile'), 'AUA': ('I', 'Ile'), 'AUG': ('M', 'Met'),
    'GUU': ('V', 'Val'), 'GUC': ('V', 'Val'), 'GUA': ('V', 'Val'), 'GUG': ('V', 'Val'),
    'UCU': ('S', 'Ser'), 'UCC': ('S', 'Ser'), 'UCA': ('S', 'Ser'), 'UCG': ('S', 'Ser'),
    'CCU': ('P', 'Pro'), 'CCC': ('P', 'Pro'), 'CCA': ('P', 'Pro'), 'CCG': ('P', 'Pro'),
    'ACU': ('T', 'Thr'), 'ACC': ('T', 'Thr'), 'ACA': ('T', 'Thr'), 'ACG': ('T', 'Thr'),
    'GCU': ('A', 'Ala'), 'GCC': ('A', 'Ala'), 'GCA': ('A', 'Ala'), 'GCG': ('A', 'Ala'),
    'UAU': ('Y', 'Tyr'), 'UAC': ('Y', 'Tyr'), 'UAA': ('*', 'STOP'), 'UAG': ('*', 'STOP'),
    'CAU': ('H', 'His'), 'CAC': ('H', 'His'), 'CAA': ('Q', 'Gln'), 'CAG': ('Q', 'Gln'),
    'AAU': ('N', 'Asn'), 'AAC': ('N', 'Asn'), 'AAA': ('K', 'Lys'), 'AAG': ('K', 'Lys'),
    'GAU': ('D', 'Asp'), 'GAC': ('D', 'Asp'), 'GAA': ('E', 'Glu'), 'GAG': ('E', 'Glu'),
    'UGU': ('C', 'Cys'), 'UGC': ('C', 'Cys'), 'UGA': ('*', 'STOP'), 'UGG': ('W', 'Trp'),
    'CGU': ('R', 'Arg'), 'CGC': ('R', 'Arg'), 'CGA': ('R', 'Arg'), 'CGG': ('R', 'Arg'),
    'AGU': ('S', 'Ser'), 'AGC': ('S', 'Ser'), 'AGA': ('R', 'Arg'), 'AGG': ('R', 'Arg'),
    'GGU': ('G', 'Gly'), 'GGC': ('G', 'Gly'), 'GGA': ('G', 'Gly'), 'GGG': ('G', 'Gly')
}

AA_WEIGHTS = {
    'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
    'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17,
    'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
}
WATER_WEIGHT = 18.015


class DNAAnalyzerRoot(BoxLayout):
    def __init__(self, **kwargs):
        super().__init__(orientation='vertical', padding=dp(12), spacing=dp(10), **kwargs)
        self.last_txt_report = ""
        self.last_fasta_report = ""

        header = Label(
            text="[b]DNA to Protein Analyzer[/b]",
            markup=True,
            font_size=dp(20),
            color=(0, 0.35, 0.75, 1),
            size_hint_y=None,
            height=dp(36)
        )
        self.add_widget(header)

        top_btns = BoxLayout(size_hint_y=None, height=dp(40), spacing=dp(8))
        btn_paste = Button(text="Paste", background_color=(0.2, 0.5, 0.85, 1), bold=True)
        btn_paste.bind(on_release=self.paste_dna)

        btn_sample = Button(text="Sample DNA", background_color=(0.15, 0.6, 0.45, 1), bold=True)
        btn_sample.bind(on_release=self.load_sample)

        btn_clear = Button(text="Clear", background_color=(0.75, 0.25, 0.25, 1), bold=True)
        btn_clear.bind(on_release=self.clear_all)

        top_btns.add_widget(btn_paste)
        top_btns.add_widget(btn_sample)
        top_btns.add_widget(btn_clear)
        self.add_widget(top_btns)

        self.dna_input = TextInput(
            hint_text="Enter DNA sequence (A, T, G, C) or FASTA...",
            multiline=True,
            size_hint_y=None,
            height=dp(110),
            font_size=dp(14),
            background_color=(1, 1, 1, 1),
            foreground_color=(0.1, 0.1, 0.1, 1)
        )
        self.add_widget(self.dna_input)

        orf_row = BoxLayout(size_hint_y=None, height=dp(38), spacing=dp(8))
        orf_lbl = Label(
            text="ORF Mode Only (First AUG to Stop):",
            color=(0.15, 0.15, 0.15, 1),
            font_size=dp(14),
            halign='left',
            valign='middle'
        )
        orf_lbl.bind(size=orf_lbl.setter('text_size'))
        self.orf_switch = Switch(active=False, size_hint_x=None, width=dp(65))
        orf_row.add_widget(orf_lbl)
        orf_row.add_widget(self.orf_switch)
        self.add_widget(orf_row)

        btn_analyze = Button(
            text="Analyze Sequence",
            size_hint_y=None,
            height=dp(46),
            background_color=(0, 0.45, 0.85, 1),
            bold=True,
            font_size=dp(16)
        )
        btn_analyze.bind(on_release=self.analyze_sequence)
        self.add_widget(btn_analyze)

        scroll = ScrollView(size_hint=(1, 1))
        self.output_label = Label(
            text="Results will appear here...",
            markup=True,
            font_size=dp(13.5),
            color=(0.1, 0.12, 0.18, 1),
            size_hint_y=None,
            halign='left',
            valign='top',
            padding=(dp(8), dp(8))
        )
        self.output_label.bind(
            width=lambda instance, value: setattr(instance, 'text_size', (value, None)),
            texture_size=lambda instance, value: setattr(instance, 'height', value[1])
        )
        scroll.add_widget(self.output_label)
        self.add_widget(scroll)

        bottom_btns = BoxLayout(size_hint_y=None, height=dp(42), spacing=dp(8))
        btn_copy = Button(text="Copy Report", background_color=(0.25, 0.35, 0.5, 1), bold=True)
        btn_copy.bind(on_release=self.copy_report)

        btn_txt = Button(text="Save .TXT", background_color=(0.1, 0.55, 0.35, 1), bold=True)
        btn_txt.bind(on_release=lambda x: self.save_file("txt"))

        btn_fasta = Button(text="Save .FASTA", background_color=(0.1, 0.55, 0.35, 1), bold=True)
        btn_fasta.bind(on_release=lambda x: self.save_file("fasta"))

        bottom_btns.add_widget(btn_copy)
        bottom_btns.add_widget(btn_txt)
        bottom_btns.add_widget(btn_fasta)
        self.add_widget(bottom_btns)

    def show_alert(self, title, message):
        content = BoxLayout(orientation='vertical', padding=dp(12), spacing=dp(10))
        lbl = Label(text=message, halign='center', valign='middle')
        lbl.bind(size=lbl.setter('text_size'))
        close_btn = Button(text="OK", size_hint_y=None, height=dp(40))
        content.add_widget(lbl)
        content.add_widget(close_btn)
        popup = Popup(title=title, content=content, size_hint=(0.85, 0.38))
        close_btn.bind(on_release=popup.dismiss)
        popup.open()

    def paste_dna(self, instance):
        clip = Clipboard.paste()
        if clip:
            self.dna_input.text = clip
        else:
            self.show_alert("Clipboard Empty", "No text found in clipboard.")

    def load_sample(self, instance):
        self.dna_input.text = (
            ">Sample_Human_Beta_Globin\n"
            "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC"
            "GTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAG"
        )
        self.analyze_sequence(None)

    def clear_all(self, instance):
        self.dna_input.text = ""
        self.output_label.text = "Results cleared."
        self.last_txt_report = ""
        self.last_fasta_report = ""

    def analyze_sequence(self, instance):
        raw = self.dna_input.text
        lines = [line for line in raw.splitlines() if not line.strip().startswith(">")]
        cleaned = re.sub(r"\s+", "", "".join(lines)).upper()

        if not cleaned:
            self.show_alert("Empty Input", "Please enter a DNA sequence first.")
            return

        invalid = sorted(set(re.findall(r"[^ATGC]", cleaned)))
        if invalid:
            self.show_alert(
                "Invalid DNA Characters",
                f"Only A, T, G, C allowed.\nInvalid found: {', '.join(invalid)}"
            )
            return

        comp_map = str.maketrans("ATGC", "TACG")
        comp_3_to_5 = cleaned.translate(comp_map)
        rev_comp_5_to_3 = comp_3_to_5[::-1]
        mrna = cleaned.replace("T", "U")

        dna_len = len(cleaned)
        c_a, c_t = cleaned.count("A"), cleaned.count("T")
        c_g, c_c = cleaned.count("G"), cleaned.count("C")
        gc_pct = ((c_g + c_c) / dna_len) * 100.0
        tm = (c_a + c_t) * 2.0 + (c_g + c_c) * 4.0 if dna_len < 14 else 64.9 + 41.0 * ((c_g + c_c - 16.4) / dna_len)

        orf_only = self.orf_switch.active
        start_idx = 0
        if orf_only:
            start_idx = mrna.find("AUG")
            if start_idx == -1:
                self.show_alert("No ORF Found", "No Start Codon (AUG) found. Turn off ORF Mode to translate.")
                return

        codons_data = []
        for i in range(start_idx, len(mrna) - 2, 3):
            codon = mrna[i:i + 3]
            aa1, aa3 = CODON_TABLE[codon]
            codons_data.append((codon, aa1, aa3))
            if orf_only and aa1 == "*":
                break

        valid_aas = [aa1 for _, aa1, _ in codons_data if aa1 != "*"]
        if valid_aas:
            mw_da = sum(AA_WEIGHTS[aa] for aa in valid_aas) - ((len(valid_aas) - 1) * WATER_WEIGHT)
        else:
            mw_da = 0.0
        mw_kda = mw_da / 1000.0

        def colorize(text, codon, aa1):
            if codon == "AUG":
                return f"[b][color=0b6e28]{text}[/color][/b]"
            if aa1 == "*":
                return f"[b][color=b71c1c]{text}[/color][/b]"
            return text

        codons_markup = " ".join(colorize(c, c, a1) for c, a1, _ in codons_data)
        aa3_markup = "-".join(colorize(a3, c, a1) for c, a1, a3 in codons_data)
        aa1_markup = "-".join(colorize(a1, c, a1) for c, a1, _ in codons_data)

        self.output_label.text = (
            f"[b][color=005ac1]=== SEQUENCE ANALYTICS ===[/color][/b]\n"
            f"Length: [b]{dna_len} bp[/b]  |  GC Content: [b]{gc_pct:.2f}%[/b]  |  Tm: [b]{tm:.1f} C[/b]\n"
            f"Protein MW: [b]{mw_da:.2f} Da ({mw_kda:.3f} kDa)[/b]  |  Residues: [b]{len(valid_aas)} aa[/b]\n"
            f"Counts: A={c_a}, T={c_t}, G={c_g}, C={c_c}\n\n"
            f"[b][color=005ac1]=== NUCLEIC ACIDS ===[/color][/b]\n"
            f"[b]5'->3' Coding DNA:[/b]\n{cleaned}\n\n"
            f"[b]3'->5' Complementary DNA:[/b]\n{comp_3_to_5}\n\n"
            f"[b]5'->3' Reverse Complement:[/b]\n{rev_comp_5_to_3}\n\n"
            f"[b]5'->3' Transcribed mRNA:[/b]\n{mrna}\n\n"
            f"[b][color=005ac1]=== PROTEIN TRANSLATION ===[/color][/b]\n"
            f"([b][color=0b6e28]Green = Start AUG[/color][/b], [b][color=b71c1c]Red = Stop[/color][/b])\n\n"
            f"[b]Codons:[/b]\n{codons_markup}\n\n"
            f"[b]3-Letter Code:[/b]\n{aa3_markup}\n\n"
            f"[b]1-Letter Code:[/b]\n{aa1_markup}"
        )

        self.last_txt_report = (
            f"=== DNA TO PROTEIN ANALYSIS REPORT ===\n"
            f"Length: {dna_len} bp | GC Content: {gc_pct:.2f}% | Tm: {tm:.1f} C\n"
            f"Protein MW: {mw_da:.2f} Da ({mw_kda:.3f} kDa)\n\n"
            f"5' -> 3' Coding DNA:\n{cleaned}\n\n"
            f"3' -> 5' Complementary DNA:\n{comp_3_to_5}\n\n"
            f"5' -> 3' Reverse Complement:\n{rev_comp_5_to_3}\n\n"
            f"5' -> 3' Transcribed mRNA:\n{mrna}\n\n"
            f"Protein (3-Letter): {'-'.join(a3 for _, _, a3 in codons_data)}\n"
            f"Protein (1-Letter): {'-'.join(a1 for _, a1, _ in codons_data)}\n"
        )

        self.last_fasta_report = (
            f">Translated_Protein | Length={len(valid_aas)}aa | MW={mw_da:.2f}Da\n"
            f"{''.join(valid_aas)}\n"
        )

    def copy_report(self, instance):
        if not self.last_txt_report:
            self.show_alert("Nothing to Copy", "Analyze a sequence first.")
            return
        Clipboard.copy(self.last_txt_report)
        self.show_alert("Copied", "Analysis report copied to clipboard!")

    def save_file(self, ext):
        if not self.last_txt_report:
            self.show_alert("Nothing to Save", "Analyze a sequence first.")
            return

        content = self.last_fasta_report if ext == "fasta" else self.last_txt_report
        filename = f"DNA_Protein_Report.{ext}"

        if platform == "android":
            from android.storage import primary_external_storage_path
            download_dir = os.path.join(primary_external_storage_path(), "Download")
            os.makedirs(download_dir, exist_ok=True)
            file_path = os.path.join(download_dir, filename)
        else:
            file_path = os.path.abspath(filename)

        try:
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(content)
            self.show_alert("File Saved", f"Saved to:\n{file_path}")
        except Exception as e:
            self.show_alert("Save Error", str(e))


class DNAProteinApp(App):
    def build(self):
        self.title = "DNA to Protein Analyzer"
        return DNAAnalyzerRoot()


if __name__ == "__main__":
    DNAProteinApp().run()
      
