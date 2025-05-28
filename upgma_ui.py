from tkinter import *
from tkinter import ttk, filedialog, messagebox
import numpy as np
import os
from upgma import *

"""
This script implements a GUI for the UPGMA algorithm using Tkinter.
It allows users to:
1. Load multiple sequences from FASTA files
2. Load a distance matrix directly from a file
3. Set parameters for the alignment (match score, mismatch penalty, gap penalty)
4. Generate and display the UPGMA tree
5. Save the tree to a file
"""

class UPGMAApp:
    def __init__(self, root):
        self.root = root
        self.root.title("UPGMA Phylogenetic Tree Generator")
        
        # Variables
        self.sequences = {}  # Dictionary to store sequences
        self.distance_matrix = None
        self.sequence_names = None
        self.match_score = StringVar(value="1")
        self.mismatch_penalty = StringVar(value="-1")
        self.gap_penalty = StringVar(value="-1")
        self.output_file = StringVar(value="phylogenetic_tree.txt")
        self.output_image = StringVar(value="phylogenetic_tree.png")
        self.input_type = StringVar(value="sequences")
        
        # Create main frame
        mainframe = ttk.Frame(root, padding="10 10 10 10")
        mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        
        # Input type selection
        ttk.Label(mainframe, text="Input Type:").grid(column=1, row=1, sticky=W)
        ttk.Radiobutton(mainframe, text="Sequences", variable=self.input_type, 
                        value="sequences", command=self.toggle_input_type).grid(column=2, row=1, sticky=W)
        ttk.Radiobutton(mainframe, text="Distance Matrix", variable=self.input_type, 
                        value="matrix", command=self.toggle_input_type).grid(column=3, row=1, sticky=W)
        
        # Sequences frame
        self.sequences_frame = ttk.LabelFrame(mainframe, text="Sequences", padding="5 5 5 5")
        self.sequences_frame.grid(column=1, row=2, columnspan=3, sticky=(W, E), pady=5)
        
        # Sequence list
        self.sequence_listbox = Listbox(self.sequences_frame, width=50, height=5)
        self.sequence_listbox.grid(column=1, row=1, columnspan=2, sticky=(W, E))
        
        # Sequence buttons
        ttk.Button(self.sequences_frame, text="Add Sequence", 
                   command=self.add_sequence).grid(column=1, row=2, sticky=W, pady=5)
        ttk.Button(self.sequences_frame, text="Remove Selected", 
                   command=self.remove_sequence).grid(column=2, row=2, sticky=E, pady=5)
        
        # Distance matrix frame
        self.matrix_frame = ttk.LabelFrame(mainframe, text="Distance Matrix", padding="5 5 5 5")
        self.matrix_frame.grid(column=1, row=3, columnspan=3, sticky=(W, E), pady=5)
        self.matrix_frame.grid_remove()  # Hide initially
        
        ttk.Label(self.matrix_frame, text="No distance matrix loaded").grid(column=1, row=1, sticky=W)
        ttk.Button(self.matrix_frame, text="Load Distance Matrix", 
                   command=self.load_distance_matrix).grid(column=1, row=2, sticky=W, pady=5)
        
        # Parameters frame
        params_frame = ttk.LabelFrame(mainframe, text="Alignment Parameters", padding="5 5 5 5")
        params_frame.grid(column=1, row=4, columnspan=3, sticky=(W, E), pady=5)
        
        ttk.Label(params_frame, text="Match Score:").grid(column=1, row=1, sticky=W)
        ttk.Entry(params_frame, width=5, textvariable=self.match_score).grid(column=2, row=1, sticky=W)
        
        ttk.Label(params_frame, text="Mismatch Penalty:").grid(column=1, row=2, sticky=W)
        ttk.Entry(params_frame, width=5, textvariable=self.mismatch_penalty).grid(column=2, row=2, sticky=W)
        
        ttk.Label(params_frame, text="Gap Penalty:").grid(column=1, row=3, sticky=W)
        ttk.Entry(params_frame, width=5, textvariable=self.gap_penalty).grid(column=2, row=3, sticky=W)
        
        # Output frame
        output_frame = ttk.LabelFrame(mainframe, text="Output Options", padding="5 5 5 5")
        output_frame.grid(column=1, row=5, columnspan=3, sticky=(W, E), pady=5)
        
        ttk.Label(output_frame, text="Text Output File:").grid(column=1, row=1, sticky=W)
        ttk.Entry(output_frame, width=30, textvariable=self.output_file).grid(column=2, row=1, sticky=W)
        
        ttk.Label(output_frame, text="Image Output File:").grid(column=1, row=2, sticky=W)
        ttk.Entry(output_frame, width=30, textvariable=self.output_image).grid(column=2, row=2, sticky=W)
        
        # Generate button
        ttk.Button(mainframe, text="Generate Phylogenetic Tree", 
                   command=self.generate_tree).grid(column=1, row=6, columnspan=3, sticky=(W, E), pady=10)
        
        # Status label
        self.status_var = StringVar(value="Ready")
        ttk.Label(mainframe, textvariable=self.status_var).grid(column=1, row=7, columnspan=3, sticky=W)
        
        # Configure grid
        for child in mainframe.winfo_children():
            child.grid_configure(padx=5, pady=5)
    
    def toggle_input_type(self):
        """Toggle between sequence input and distance matrix input"""
        if self.input_type.get() == "sequences":
            self.sequences_frame.grid()
            self.matrix_frame.grid_remove()
        else:
            self.sequences_frame.grid_remove()
            self.matrix_frame.grid()
    
    def add_sequence(self):
        """Add a sequence from a FASTA file"""
        file_paths = filedialog.askopenfilenames(
            title="Select FASTA file(s)",
            filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
        )
        
        if not file_paths:
            return
        
        for file_path in file_paths:
            try:
                new_sequences = load_sequences_from_fasta(file_path)
                self.sequences.update(new_sequences)
                
                # Update the listbox
                self.update_sequence_listbox()
                
                self.status_var.set(f"Loaded {len(new_sequences)} sequences from {os.path.basename(file_path)}")
            except Exception as e:
                messagebox.showerror("Error", f"Error loading FASTA file: {e}")
    
    def update_sequence_listbox(self):
        """Update the sequence listbox with current sequences"""
        self.sequence_listbox.delete(0, END)
        for name in self.sequences.keys():
            seq = self.sequences[name]
            display_seq = seq[:20] + "..." if len(seq) > 20 else seq
            self.sequence_listbox.insert(END, f"{name}: {display_seq}")
    
    def remove_sequence(self):
        """Remove the selected sequence"""
        selection = self.sequence_listbox.curselection()
        if not selection:
            return
        
        index = selection[0]
        name = list(self.sequences.keys())[index]
        del self.sequences[name]
        
        self.update_sequence_listbox()
        self.status_var.set(f"Removed sequence: {name}")
    
    def load_distance_matrix(self):
        """Load a distance matrix from a file"""
        file_path = filedialog.askopenfilename(
            title="Select distance matrix file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        
        if not file_path:
            return
        
        try:
            self.distance_matrix, self.sequence_names = load_distance_matrix_from_file(file_path)
            
            # Update the label
            for widget in self.matrix_frame.winfo_children():
                widget.destroy()
            
            ttk.Label(self.matrix_frame, 
                      text=f"Loaded matrix with {len(self.sequence_names)} sequences").grid(column=1, row=1, sticky=W)
            ttk.Button(self.matrix_frame, text="Load Different Matrix", 
                       command=self.load_distance_matrix).grid(column=1, row=2, sticky=W, pady=5)
            
            self.status_var.set(f"Loaded distance matrix from {os.path.basename(file_path)}")
        except Exception as e:
            messagebox.showerror("Error", f"Error loading distance matrix: {e}")
    
    def generate_tree(self):
        """Generate the phylogenetic tree"""
        try:
            # Get parameters
            match_score = int(self.match_score.get())
            mismatch_penalty = int(self.mismatch_penalty.get())
            gap_penalty = int(self.gap_penalty.get())
            output_file = self.output_file.get()
            output_image = self.output_image.get()
            
            # Check input
            if self.input_type.get() == "sequences":
                if not self.sequences:
                    messagebox.showerror("Error", "No sequences loaded")
                    return
                
                # Generate tree from sequences
                self.status_var.set("Generating distance matrix...")
                self.root.update()
                
                root = main_upgma(
                    sequences=self.sequences,
                    match_score=match_score,
                    mismatch_penalty=mismatch_penalty,
                    gap_penalty=gap_penalty,
                    output_file=output_file,
                    output_image=output_image
                )
            else:
                if self.distance_matrix is None:
                    messagebox.showerror("Error", "No distance matrix loaded")
                    return
                
                # Generate tree from distance matrix
                self.status_var.set("Generating tree from distance matrix...")
                self.root.update()
                
                root = main_upgma(
                    distance_matrix=self.distance_matrix,
                    sequence_names=self.sequence_names,
                    output_file=output_file,
                    output_image=output_image
                )
            
            # Show success message
            self.status_var.set("Tree generated successfully!")
            messagebox.showinfo("Success", 
                               f"Phylogenetic tree generated successfully!\n\n"
                               f"Text output: {output_file}\n"
                               f"Image output: {output_image}")
            
            # Open the image
            os.startfile(output_image)
            
        except Exception as e:
            messagebox.showerror("Error", f"Error generating tree: {e}")
            self.status_var.set("Error generating tree")

if __name__ == "__main__":
    root = Tk()
    app = UPGMAApp(root)
    root.mainloop()