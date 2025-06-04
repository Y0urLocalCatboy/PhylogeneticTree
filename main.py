from tkinter import *

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
