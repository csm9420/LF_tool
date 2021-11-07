# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:08:37 2020

@author: micha
"""

import tkinter as tk
import tkinter.messagebox


def check_overwrite(project_id):
    """ Creates a GUI window; asking the user if overwriting data is acceptable.

    :param project_id: Name of the project which could be overwritten.
    :return:
    """
    root = tk.Tk()
    root.withdraw()
    msg_box = tk.messagebox.askquestion('Warning',
                                        f'A project with the ID "{project_id}" already exists. Data will be '
                                        f'overwritten. This cannot be undone.\n\n'
                                        f'Do you want to delete existing data and continue?',
                                        icon="warning")
    root.destroy()

    if msg_box == "yes":
        # Overwriting data is okay
        return True
    else:
        # Overwriting data is not okay
        return False


def convergenceError(process):
    root = tk.Tk()
    root.withdraw()
    msgBox = tk.messagebox.showerror("Convergence Error", "Calculation for " + process + " did not converge!")
    root.destroy()


def nyquistFreqWarning():
    root = tk.Tk()
    root.withdraw()
    tk.messagebox.showwarning("Nyquist Diagram Frequency Warning", \
                              "Starting frequency is not 0! Nyquist Diagram requires a starting frequnency of 0 and is automaticilly extended.")
    root.destroy()


def chooseRunMode():
    root = tk.Tk()

    ws = root.winfo_screenwidth()
    hs = root.winfo_screenheight()
    w = root.winfo_width()
    h = root.winfo_height()
    x = (ws / 2) - (w / 2)
    y = (hs / 2) - (h / 2)
    root.geometry('+%d+%d' % (x, y))

    feedSys = tk.BooleanVar()
    feedSys.set(False)
    modFeedSys = tk.BooleanVar()
    modFeedSys.set(False)

    # root.withdraw()

    class ModeChooser:

        def __init__(self, master):
            frame = tk.Canvas(master)
            frame.pack()
            self.Button1 = tk.Button(frame, text="Nyquist and Bode Diagrams", command=self.nyqMode)
            self.Button1.pack(side=tk.LEFT)
            self.Button2 = tk.Button(frame, text="Pressure Drop Diagrams", command=self.presDropMode)
            self.Button2.pack(side=tk.LEFT)
            self.CheckBox1 = tk.Checkbutton(frame, text="FeedSystemOff", var=feedSys, command=self.feedSysOff)
            self.CheckBox1.pack(side=tk.BOTTOM)
            self.CheckBox2 = tk.Checkbutton(frame, text="Modular Feed System", var=modFeedSys, command=self.modularFeed)
            self.CheckBox2.pack(side=tk.BOTTOM)

        def nyqMode(self):
            print("Nyquist Mode")
            global nyquist
            nyquist = True
            root.quit()

        def presDropMode(self):
            print("Pressure Drop Diagrams")
            global nyquist
            nyquist = False
            root.quit()

        def feedSysOff(self):
            if (feedSys.get()):
                feedSys.set(True)
                print("Feed system Off")
            else:
                feedSys.set(False)
                print("Feed system On")

        def modularFeed(self):
            if (modFeedSys.get()):
                modFeedSys.set(True)
                print("Modular Feed System selected")
            else:
                modFeedSys.set(False)
                print("Modular Feed System deselected")

    app = ModeChooser(root)
    root.mainloop()
    root.destroy()

    ##Inverting FeedSys, nessesary in many program as bool feedOn while button is more useful as feedOff
    return nyquist, not feedSys.get(), modFeedSys.get()
