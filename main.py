from tkinter import *

"""
Main application.
Directly launches the UPGMA Phylogenetic Tree Construction.
Created by Gabriel Pankowski
"""

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
