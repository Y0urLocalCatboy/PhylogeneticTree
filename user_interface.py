from tkinter import *
from tkinter import ttk, filedialog
from algorithm import *
"""This script implements a GUI for the Needleman-Wunsch algorithm using Tkinter.
It allows users to input two sequences and their scoring parameters 
and then calculates the optimal alignment using the Needleman-Wunsch algorithm from algorithm.py.
The GUI includes:
- Input fields for the first and second sequences
- Input fields for match score, mismatch penalty, and gap penalty
- A button to calculate the alignment
- A button to load sequences from FASTA files
"""

root = Tk()
root.title("Needleman Wunsch Algorithm")
mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

first_sentence = StringVar()
second_sentence = StringVar()
missmatch_penalty = StringVar()
gap_penalty = StringVar()
match_score = StringVar()

first_sentence_entry = ttk.Entry(mainframe, width=7, textvariable=first_sentence)
second_sentence_entry = ttk.Entry(mainframe, width=7, textvariable=second_sentence)
missmatch_penalty_entry = ttk.Entry(mainframe, width=7, textvariable=missmatch_penalty)
gap_penalty_entry = ttk.Entry(mainframe, width=7, textvariable=gap_penalty)
match_score_entry = ttk.Entry(mainframe, width=7, textvariable=match_score)

ttk.Label(mainframe, text="First sentence").grid(column=1, row=1, sticky=W)
ttk.Label(mainframe, text="Second sentence").grid(column=3, row=1, sticky=W)
ttk.Label(mainframe, text="Match score").grid(column=1, row=2, sticky=W)
ttk.Label(mainframe, text="Missmatch penalty").grid(column=3, row=2, sticky=W)
ttk.Label(mainframe, text="Gap penalty").grid(column=1, row=3, sticky=W)

first_sentence_entry.grid(column=2, row=1, sticky=(W, E))
second_sentence_entry.grid(column=4, row=1, sticky=(W, E))
missmatch_penalty_entry.grid(column=4, row=2, sticky=(W, E))
gap_penalty_entry.grid(column=2, row=3, sticky=(W, E))
match_score_entry.grid(column=2, row=2, sticky=(W, E))

ttk.Separator(mainframe, orient='horizontal').grid(column=1, row=4, columnspan=4, sticky=(W, E), pady=10)
ttk.Label(mainframe, text="From Fasta files:").grid(column=1, row=4, columnspan=4, sticky=W, pady=10)

(ttk.Button(mainframe,
           text="Calculate",
           command=lambda: main(first_sentence.get(),
                                second_sentence.get(),
                                int(match_score.get()),
                                int(missmatch_penalty.get()),
                                int(gap_penalty.get())))
                                .grid(column=1, row=5, columnspan=5, sticky=(W, E), pady=10))


def load_fasta_file(target_var):
    """Loads a FASTA format file and extract the sequence

    Args:
        target_var (StringVar): The variable to store the sequence
        """

    file_path = filedialog.askopenfilename(
        title="Select FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
    )
    if not file_path:
        return

    try:
        with open(file_path, "r") as fasta_file:
            lines = fasta_file.readlines()

            sequence = ""
            for line in lines:
                if not line.startswith(">"):
                    sequence += line.strip()

            target_var.set(sequence)
    except Exception as e:
        print(f"Error loading FASTA file: {e}")


ttk.Button(mainframe, text="Load sequence 1",
           command=lambda: load_fasta_file(first_sentence)).grid(column=2, row=4, sticky=W)

ttk.Button(mainframe, text="Load sequence 2",
           command=lambda: load_fasta_file(second_sentence)).grid(column=4, row=4, sticky=W)


for child in mainframe.winfo_children():
    child.grid_configure(padx=5, pady=5)

root.mainloop()