from tkinter import *
from tkinter import ttk, filedialog, messagebox
import sys
import os
from upgma import *

"""
Main entry point for the Phylogenetic Tree application.
This script directly launches the UPGMA Phylogenetic Tree Construction.
Gabriel Pankowski
"""

# Import the UPGMAApp class from upgma_ui.py
from upgma_ui import UPGMAApp

def main():
    """
    Main function to launch the UPGMA Phylogenetic Tree application directly
    """
    root = Tk()
    app = UPGMAApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
