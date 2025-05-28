from tkinter import *
from tkinter import ttk
import sys
import os

"""
Main entry point for the Phylogenetic Tree application.
This script provides a launcher that allows users to choose between:
1. Needleman-Wunsch Pairwise Alignment
2. UPGMA Phylogenetic Tree Construction

Author: [Your Name]
Date: [Current Date]
"""

class LauncherApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Phylogenetic Analysis Tools")
        
        # Set window size and position
        window_width = 400
        window_height = 300
        screen_width = root.winfo_screenwidth()
        screen_height = root.winfo_screenheight()
        center_x = int(screen_width/2 - window_width/2)
        center_y = int(screen_height/2 - window_height/2)
        root.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')
        
        # Create main frame
        mainframe = ttk.Frame(root, padding="20 20 20 20")
        mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        
        # Title
        ttk.Label(mainframe, text="Phylogenetic Analysis Tools", 
                 font=("Helvetica", 16)).grid(column=1, row=1, columnspan=2, pady=20)
        
        # Description
        ttk.Label(mainframe, text="Select a tool to launch:", 
                 font=("Helvetica", 12)).grid(column=1, row=2, columnspan=2, pady=10)
        
        # Buttons
        ttk.Button(mainframe, text="Needleman-Wunsch Alignment", 
                  command=self.launch_needleman_wunsch,
                  width=30).grid(column=1, row=3, columnspan=2, pady=5)
        
        ttk.Button(mainframe, text="UPGMA Phylogenetic Tree", 
                  command=self.launch_upgma,
                  width=30).grid(column=1, row=4, columnspan=2, pady=5)
        
        # Exit button
        ttk.Button(mainframe, text="Exit", 
                  command=root.destroy,
                  width=10).grid(column=1, row=5, columnspan=2, pady=20)
        
        # Configure grid
        for child in mainframe.winfo_children():
            child.grid_configure(padx=5, pady=5)
    
    def launch_needleman_wunsch(self):
        """Launch the Needleman-Wunsch alignment interface"""
        self.root.destroy()  # Close the launcher
        os.system(f"{sys.executable} user_interface.py")
    
    def launch_upgma(self):
        """Launch the UPGMA phylogenetic tree interface"""
        self.root.destroy()  # Close the launcher
        os.system(f"{sys.executable} upgma_ui.py")

def main():
    root = Tk()
    app = LauncherApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()