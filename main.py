import re
import flet as ft

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


def main(page: ft.Page):
    page.title = "DNA to Protein Analyzer"
    page.theme_mode = ft.ThemeMode.LIGHT
    page.scroll = ft.ScrollMode.AUTO
    page.padding = 16

    last_report = {"txt": ""}

    dna_input = ft.TextField(
        label="DNA Sequence (A, T, G, C) or FASTA",
        multiline=True,
        min_lines=4,
        max_lines=6,
        border_radius=10
    )

    orf_switch = ft.Switch(label="ORF Mode Only (First AUG to Stop)", value=False)
    result_text = ft.Text("Results will appear here...", selectable=True, size=14)

    def show_msg(msg):
        page.snack_bar = ft.SnackBar(ft.Text(msg))
        page.snack_bar.open = True
        page.update()

    def load_sample(e):
        dna_input.value = (
            ">Sample_Human_Beta_Globin\n"
            "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC"
            "GTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAG"
        )
        analyze(None)

    def clear_all(e):
        dna_input.value = ""
        result_text.value = "Results cleared."
        last_report["txt"] = ""
        page.update()

    def copy_report(e):
        if not last_report["txt"]:
            show_msg("Analyze a sequence first.")
            return
        page.set_clipboard(last_report["txt"])
        show_msg("Report copied to clipboard!")

    def analyze(e):
        raw = dna_input.value or ""
        lines = [line for line in raw.splitlines() if not line.strip().startswith(">")]
        cleaned = re.sub(r"\s+", "", "".join(lines)).upper()

        if not cleaned:
            show_msg("Please enter a DNA sequence first.")
            return

        invalid = sorted(set(re.findall(r"[^ATGC]", cleaned)))
        if invalid:
            show_msg(f"Invalid characters found: {', '.join(invalid)}")
            return

        comp_3_to_5 = cleaned.translate(str.maketrans("ATGC", "TACG"))
        rev_comp = comp_3_to_5[::-1]
        mrna = cleaned.replace("T", "U")

        dna_len = len(cleaned)
        c_a, c_t, c_g, c_c = cleaned.count("A"), cleaned.count("T"), cleaned.count("G"), cleaned.count("C")
        gc_pct = ((c_g + c_c) / dna_len) * 100.0
        tm = (c_a + c_t) * 2.0 + (c_g + c_c) * 4.0 if dna_len < 14 else 64.9 + 41.0 * ((c_g + c_c - 16.4) / dna_len)

        start_idx = 0
        if orf_switch.value:
            start_idx = mrna.find("AUG")
            if start_idx == -1:
                show_msg("No Start Codon (AUG) found.")
                return

        codons_data = []
        for i in range(start_idx, len(mrna) - 2, 3):
            codon = mrna[i:i + 3]
            aa1, aa3 = CODON_TABLE[codon]
            codons_data.append((codon, aa1, aa3))
            if orf_switch.value and aa1 == "*":
                break

        valid_aas = [aa1 for _, aa1, _ in codons_data if aa1 != "*"]
        mw_da = sum(AA_WEIGHTS[aa] for aa in valid_aas) - ((len(valid_aas) - 1) * WATER_WEIGHT) if valid_aas else 0.0

        report = (
            f"=== SEQUENCE ANALYTICS ===\n"
            f"Length: {dna_len} bp | GC: {gc_pct:.2f}% | Tm: {tm:.1f} C\n"
            f"Protein MW: {mw_da:.2f} Da ({mw_da / 1000.0:.3f} kDa) | Residues: {len(valid_aas)} aa\n"
            f"Counts: A={c_a}, T={c_t}, G={c_g}, C={c_c}\n\n"
            f"5'->3' Coding DNA:\n{cleaned}\n\n"
            f"3'->5' Complementary DNA:\n{comp_3_to_5}\n\n"
            f"5'->3' Reverse Complement:\n{rev_comp}\n\n"
            f"5'->3' mRNA:\n{mrna}\n\n"
            f"=== PROTEIN TRANSLATION ===\n"
            f"Codons: {' '.join(c for c, _, _ in codons_data)}\n\n"
            f"3-Letter: {'-'.join(a3 for _, _, a3 in codons_data)}\n\n"
            f"1-Letter: {'-'.join(a1 for _, a1, _ in codons_data)}"
        )

        last_report["txt"] = report
        result_text.value = report
        page.update()

    page.add(
        ft.Text("DNA to Protein Analyzer", size=22, weight=ft.FontWeight.BOLD, color=ft.Colors.BLUE_800),
        ft.Row([
            ft.ElevatedButton("Sample DNA", on_click=load_sample),
            ft.ElevatedButton("Clear", on_click=clear_all),
            ft.ElevatedButton("Copy Report", on_click=copy_report),
        ], wrap=True),
        dna_input,
        orf_switch,
        ft.ElevatedButton(
            "Analyze Sequence",
            on_click=analyze,
            bgcolor=ft.Colors.BLUE_700,
            color=ft.Colors.WHITE,
            width=400,
            height=45
        ),
        ft.Divider(),
        result_text
    )


ft.app(target=main)
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
      
