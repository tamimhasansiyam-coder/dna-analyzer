
import re
import flet as ft

APP_NAME = "DNA to Protein Analyzer"

CODON_TABLE = {
    "UUU": ("F", "Phe"), "UUC": ("F", "Phe"),
    "UUA": ("L", "Leu"), "UUG": ("L", "Leu"),
    "CUU": ("L", "Leu"), "CUC": ("L", "Leu"),
    "CUA": ("L", "Leu"), "CUG": ("L", "Leu"),
    "AUU": ("I", "Ile"), "AUC": ("I", "Ile"),
    "AUA": ("I", "Ile"), "AUG": ("M", "Met"),
    "GUU": ("V", "Val"), "GUC": ("V", "Val"),
    "GUA": ("V", "Val"), "GUG": ("V", "Val"),

    "UCU": ("S", "Ser"), "UCC": ("S", "Ser"),
    "UCA": ("S", "Ser"), "UCG": ("S", "Ser"),
    "CCU": ("P", "Pro"), "CCC": ("P", "Pro"),
    "CCA": ("P", "Pro"), "CCG": ("P", "Pro"),
    "ACU": ("T", "Thr"), "ACC": ("T", "Thr"),
    "ACA": ("T", "Thr"), "ACG": ("T", "Thr"),
    "GCU": ("A", "Ala"), "GCC": ("A", "Ala"),
    "GCA": ("A", "Ala"), "GCG": ("A", "Ala"),

    "UAU": ("Y", "Tyr"), "UAC": ("Y", "Tyr"),
    "UAA": ("*", "STOP"), "UAG": ("*", "STOP"),
    "CAU": ("H", "His"), "CAC": ("H", "His"),
    "CAA": ("Q", "Gln"), "CAG": ("Q", "Gln"),
    "AAU": ("N", "Asn"), "AAC": ("N", "Asn"),
    "AAA": ("K", "Lys"), "AAG": ("K", "Lys"),
    "GAU": ("D", "Asp"), "GAC": ("D", "Asp"),
    "GAA": ("E", "Glu"), "GAG": ("E", "Glu"),

    "UGU": ("C", "Cys"), "UGC": ("C", "Cys"),
    "UGA": ("*", "STOP"), "UGG": ("W", "Trp"),
    "CGU": ("R", "Arg"), "CGC": ("R", "Arg"),
    "CGA": ("R", "Arg"), "CGG": ("R", "Arg"),
    "AGU": ("S", "Ser"), "AGC": ("S", "Ser"),
    "AGA": ("R", "Arg"), "AGG": ("R", "Arg"),
    "GGU": ("G", "Gly"), "GGC": ("G", "Gly"),
    "GGA": ("G", "Gly"), "GGG": ("G", "Gly")
}

AA_WEIGHTS = {
    "A": 89.09, "R": 174.20, "N": 132.12,
    "D": 133.10, "C": 121.16, "E": 147.13,
    "Q": 146.15, "G": 75.07, "H": 155.16,
    "I": 131.17, "L": 131.17, "K": 146.19,
    "M": 149.21, "F": 165.19, "P": 115.13,
    "S": 105.09, "T": 119.12, "W": 204.23,
    "Y": 181.19, "V": 117.15
}

WATER_WEIGHT = 18.015

SAMPLE_DNA = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC"
    "GTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAG"
)


def clean_sequence(raw):
    lines = []

    for line in raw.splitlines():
        line = line.strip()

        if not line or line.startswith(">"):
            continue

        lines.append(line)

    sequence = re.sub(r"\s+", "", "".join(lines)).upper()

    return sequence


def analyze_sequence(sequence, orf_mode=False):
    if not sequence:
        raise ValueError("Please enter a DNA sequence.")

    invalid = sorted(set(re.findall(r"[^ATGC]", sequence)))

    if invalid:
        raise ValueError(
            "Invalid DNA characters: " + ", ".join(invalid)
        )

    length = len(sequence)

    count_a = sequence.count("A")
    count_t = sequence.count("T")
    count_g = sequence.count("G")
    count_c = sequence.count("C")

    gc_count = count_g + count_c
    gc_percent = (gc_count / length) * 100

    if length < 14:
        tm = 2 * (count_a + count_t) + 4 * gc_count
        tm_method = "Wallace rule"
    else:
        tm = 64.9 + 41 * ((gc_count - 16.4) / length)
        tm_method = "Long-sequence approximation"

    complementary = sequence.translate(
        str.maketrans("ATGC", "TACG")
    )

    reverse_complement = complementary[::-1]

    mrna = sequence.replace("T", "U")

    start_index = 0

    if orf_mode:
        start_index = mrna.find("AUG")

        if start_index == -1:
            raise ValueError("No start codon (AUG) found.")

    codons = []
    amino_acids = []
    three_letter = []
    stop_found = False

    for i in range(start_index, len(mrna) - 2, 3):
        codon = mrna[i:i + 3]
        aa1, aa3 = CODON_TABLE[codon]

        codons.append(codon)
        amino_acids.append(aa1)
        three_letter.append(aa3)

        if aa1 == "*":
            stop_found = True

            if orf_mode:
                break

    valid_aas = [
        aa for aa in amino_acids if aa != "*"
    ]

    residue_count = len(valid_aas)

    molecular_weight = sum(
        AA_WEIGHTS[aa] for aa in valid_aas
    )

    if residue_count > 1:
        molecular_weight -= (
            (residue_count - 1) * WATER_WEIGHT
        )

    protein = "".join(valid_aas)

    report = f"""
DNA TO PROTEIN ANALYSIS REPORT
================================

SEQUENCE STATISTICS

DNA Length: {length} bp
GC Content: {gc_percent:.2f}%
Melting Temperature: {tm:.2f} °C
Tm Method: {tm_method}

Nucleotide Count:
A = {count_a}
T = {count_t}
G = {count_g}
C = {count_c}

DNA SEQUENCES

5' to 3' Coding DNA:
{sequence}

3' to 5' Complementary DNA:
{complementary}

5' to 3' Reverse Complement:
{reverse_complement}

5' to 3' mRNA:
{mrna}

PROTEIN TRANSLATION
================================

Codons:
{" ".join(codons) if codons else "None"}

Three-Letter Amino Acids:
{"-".join(three_letter) if three_letter else "None"}

One-Letter Protein Sequence:
{protein if protein else "None"}

Protein Residues: {residue_count}

Stop Codon Found: {"Yes" if stop_found else "No"}

MOLECULAR WEIGHT
================================

Molecular Weight:
{molecular_weight:.2f} Da

Molecular Weight:
{molecular_weight / 1000:.4f} kDa

NOTE:
Tm is an approximate estimate.
This tool analyzes the supplied sequence as a
coding strand in the 5' to 3' direction.
ORF mode starts at the first AUG and stops
at the first in-frame stop codon.
"""

    return report.strip()


def main(page: ft.Page):

    page.title = APP_NAME
    page.theme_mode = ft.ThemeMode.LIGHT
    page.bgcolor = "#F3F6FB"
    page.padding = 16
    page.scroll = ft.ScrollMode.AUTO

    report_data = {"text": ""}

    title = ft.Text(
        "DNA Analyzer",
        size=27,
        weight=ft.FontWeight.BOLD,
        color="#123B70"
    )

    subtitle = ft.Text(
        "DNA Sequence and Protein Analysis",
        size=13,
        color="#60738B"
    )

    dna_input = ft.TextField(
        label="Enter DNA Sequence or FASTA",
        hint_text="Example: ATGGTGCATCTGACT...",
        multiline=True,
        min_lines=5,
        max_lines=9,
        border_radius=12,
        bgcolor=ft.Colors.WHITE,
        border_color="#C9D5E5"
    )

    orf_switch = ft.Switch(
        label="ORF Mode (First AUG to Stop)",
        value=False
    )

    status_text = ft.Text(
        "Enter a DNA sequence to begin.",
        size=12,
        color="#60738B"
    )

    result_text = ft.Text(
        "Your analysis report will appear here.",
        selectable=True,
        size=13,
        color="#26364A"
    )

    result_card = ft.Container(
        content=result_text,
        padding=16,
        bgcolor=ft.Colors.WHITE,
        border_radius=12,
        border=ft.Border.all(1, "#DCE4EF")
    )

    def show_status(message, error=False):
        status_text.value = message
        status_text.color = "#C62828" if error else "#19734B"
        page.update()

    def analyze_click(e):
        try:
            sequence = clean_sequence(dna_input.value or "")

            report = analyze_sequence(
                sequence,
                orf_switch.value
            )

            report_data["text"] = report
            result_text.value = report

            show_status("Analysis completed successfully.")

        except ValueError as error:
            show_status(str(error), True)

    def load_sample(e):
        dna_input.value = SAMPLE_DNA
        analyze_click(None)

    def clear_all(e):
        dna_input.value = ""
        result_text.value = (
            "Your analysis report will appear here."
        )
        report_data["text"] = ""

        status_text.value = "All data cleared."
        status_text.color = "#60738B"

        page.update()

    def copy_report(e):
        if not report_data["text"]:
            show_status("Analyze a sequence first.", True)
            return

        page.set_clipboard(report_data["text"])
        show_status("Report copied to clipboard.")

    page.add(
        ft.Container(
            content=ft.Column(
                [
                    ft.Text("🧬", size=38),
                    title,
                    subtitle
                ],
                spacing=4
            ),
            padding=ft.Padding.only(bottom=15)
        ),

        ft.Container(
            content=ft.Column(
                [
                    ft.Text(
                        "Sequence Input",
                        size=18,
                        weight=ft.FontWeight.BOLD,
                        color="#123B70"
                    ),

                    dna_input,

                    orf_switch,

                    ft.Row(
                        [
                            ft.OutlinedButton(
                                "Sample DNA",
                                icon=ft.Icons.SCIENCE,
                                on_click=load_sample
                            ),

                            ft.OutlinedButton(
                                "Clear",
                                icon=ft.Icons.DELETE_OUTLINE,
                                on_click=clear_all
                            )
                        ],
                        wrap=True
                    ),

                    ft.Button(
                        "Analyze Sequence",
                        icon=ft.Icons.ANALYTICS,
                        on_click=analyze_click,
                        width=float("inf"),
                        height=48,
                        style=ft.ButtonStyle(
                            bgcolor="#1557A0",
                            color=ft.Colors.WHITE
                        )
                    ),

                    status_text
                ],
                spacing=14
            ),
            padding=16,
            bgcolor=ft.Colors.WHITE,
            border_radius=14,
            border=ft.Border.all(1, "#DCE4EF")
        ),

        ft.Container(height=12),

        ft.Row(
            [
                ft.Text(
                    "Analysis Results",
                    size=19,
                    weight=ft.FontWeight.BOLD,
                    color="#123B70"
                ),

                ft.Container(expand=True),

                ft.IconButton(
                    icon=ft.Icons.CONTENT_COPY,
                    tooltip="Copy Report",
                    on_click=copy_report
                )
            ]
        ),

        result_card,

        ft.Container(
            content=ft.Text(
                "DNA to Protein Analyzer | Version 1.0",
                size=11,
                color="#718096",
                text_align=ft.TextAlign.CENTER
            ),
            padding=12,
            alignment=ft.Alignment.CENTER
        )
    )


if __name__ == "__main__":
    ft.run(main)
